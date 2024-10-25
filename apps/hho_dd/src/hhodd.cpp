#include <iostream>
#include <regex>
#include <set>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

template<typename Mesh>
using sdmap_t = std::map<size_t, std::vector<typename Mesh::cell_type>>;

using flagmap_t = std::map<size_t, std::vector<bool>>;

template<typename Mesh>
static flagmap_t
make_cell_flagmap(const Mesh& msh) {
    flagmap_t ret;
    for (size_t clnum = 0; clnum < msh.cells_size(); clnum++) {
        auto cl = msh.cell_at(clnum);
        auto di = msh.domain_info(cl);
        auto& flags = ret[di.tag()];
        if (flags.size() == 0)
            flags.resize( msh.cells_size() );
        flags[clnum] = true; 
    }

    return ret;
}

template<typename Mesh>
static flagmap_t
make_face_flagmap(const Mesh& msh) {
    flagmap_t ret;
    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        auto& flags = ret[di.tag()];
        if (flags.size() == 0)
            flags.resize( msh.faces_size() );
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs) {
            auto bi = msh.boundary_info(fc);
            if (bi.is_boundary() /*and bi.is_internal()*/)
                continue;
            flags[offset(msh, fc)] = true;
        } 
    }

    return ret;
}

template<typename Mesh>
auto
make_overlapping_subdomains(const Mesh& msh, size_t overlap_layers)
{
    using cell_type = typename Mesh::cell_type;
    using cvi = typename std::vector<cell_type>::iterator;
    auto cvf = connectivity_via_faces(msh);

    auto subdomain_cells = make_cell_flagmap(msh);
    for (size_t ol = 0; ol < overlap_layers; ol++) {
        std::map<size_t, std::set<size_t>> layer_cells;
        /* Iterate on the currently identified subdomains */
        for (auto& [tag, present] : subdomain_cells) {
            /* for each cell */
            for (size_t clnum = 0; clnum < present.size(); clnum++) {
                if (not present[clnum])
                    continue;

                auto cl = msh.cell_at(clnum);
                auto fcs = faces(msh, cl);
                for (auto& fc : fcs) {
                    /* for each neighbour */
                    auto neigh = cvf.neighbour_via(msh, cl, fc);
                    if (not neigh)
                        continue;

                    auto neighnum = offset(msh, neigh.value());
                    if (not present[neighnum])
                        layer_cells[tag].insert(neighnum);
                }
            }
        }

        /* move the found neighbours to the subdomain */
        for (auto& [tag, cellnums] : layer_cells) {
            for (auto& cellnum : cellnums)
                subdomain_cells[tag][cellnum] = true;
        }
    }

    auto subdomain_faces = make_face_flagmap(msh);
    for (auto& [tag, cell_present] : subdomain_cells) {
        for(size_t clnum = 0; clnum < cell_present.size(); clnum++) {
            if (not cell_present[clnum])
                continue;
            auto cl = msh.cell_at(clnum);
            auto fcs = faces(msh, cl);
            for (auto& fc : fcs) {
                auto neigh = cvf.neighbour_via(msh, cl, fc);
                if (not neigh)
                    continue;

                auto neigh_id = offset(msh, neigh.value());
                if (cell_present[neigh_id] /*or include_dd_bnd*/) {
                    subdomain_faces[tag][offset(msh, fc)] = true;
                }
            }
        }
    }

    return std::pair(subdomain_cells, subdomain_faces);
}

template<typename Mesh>
void
hhodd(const Mesh& msh, size_t levels)
{
    auto [sd_cells, sd_faces] = make_overlapping_subdomains(msh, levels);
    
    disk::silo_database silo;
    silo.create("hhodd.silo");
    silo.add_mesh(msh, "mesh");

    for (auto& [tag, cell_present] : sd_cells)
    {
        std::vector<double> yesno(msh.cells_size());
        std::transform(cell_present.begin(), cell_present.end(),
            yesno.begin(), [](bool x) { return double(x); } );

        std::stringstream ss;
        ss << "domain" << tag;
        silo.add_variable("mesh", ss.str(), yesno, disk::zonal_variable_t);
    }

    silo.close();

    if constexpr (Mesh::dimension == 2) {
        for (auto& [tag, ifcs] : sd_faces ) {
            std::stringstream ss;
            ss << "faces_" << tag << ".m";
            std::ofstream ofs(ss.str());
            for (size_t i = 0; i < ifcs.size(); i++) {
                if (not ifcs[i])
                    continue;
                auto fc = msh.face_at(i);
                auto pts = points(msh, fc);
                ofs << "line([" << pts[0].x() << "," << pts[1].x() << "],";
                ofs << "[" << pts[0].y() << "," << pts[1].y() << "]);" << std::endl;
            }
        }
    }

}

int main(int argc, char **argv)
{
    int levels = 1;
    int opt;
    while ((opt = getopt(argc, argv, "l:")) != -1) {
        switch (opt) {
        case 'l': {
            int l = atoi(optarg);
            if (l < 0) {
                std::cerr << "Levels must be positive, resetting to 1." << std::endl;
                l = 1;
            }
            levels = l;
            } break;
        }
    }
    argc -= optind;
    argv += optind;

    if (argc < 1) {
        std::cout << "missing filename" << std::endl;
        return 1;
    }

    const char *mesh_filename = argv[0];

    using T = double;

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<T,2>;
        mesh_type msh;
        disk::gmsh_geometry_loader< mesh_type > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        hhodd(msh, levels);

        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<T,3>;
        mesh_type msh;
        disk::gmsh_geometry_loader< mesh_type > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        hhodd(msh, levels);

        return 0;
    }

}
