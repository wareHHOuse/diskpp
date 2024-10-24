#include <iostream>
#include <regex>
#include <set>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

template<typename Mesh>
using sdmap_t = std::map<size_t, std::vector<typename Mesh::cell_type>>;

template<typename Mesh>
auto
make_overlapping_subdomains(const Mesh& msh, size_t overlap_layers)
{
    /* What this should actually do is a BFS, but Disk++ is missing
     * some infrastructure to do a proper BFS. Therefore this code is
     * quite inefficient and sooner or later must disappear. */
    using cell_type = typename Mesh::cell_type;

    using cvi = typename std::vector<cell_type>::iterator;

    auto cvf = connectivity_via_faces(msh);

    /* Collect the cells in the various subdomains:
     * the cells are stored ordered in msh, therefore each
     * vector in the map will be sorted too. this is necessary
     * precondition for std::lower_bound() below. */
    sdmap_t<Mesh> subdomains;
    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        auto tag = di.tag();
        subdomains[tag].push_back(cl);
    }

    for (size_t ol = 0; ol < overlap_layers; ol++) {
        std::map<size_t, std::set<cell_type>> layer_cells;
        /* Iterate on the currently identified subdomains */
        for (auto& [tag, cells] : subdomains) {
            /* for each cell */
            for (auto& cl : cells) {
                auto fcs = faces(msh, cl);
                for (auto& fc : fcs) {
                    /* for each neighbour */
                    auto neigh = cvf.neighbour_via(msh, cl, fc);
                    if (not neigh)
                        continue;

                    /* and check if it is already in the current subdomain */
                    bool present = std::binary_search(cells.begin(),
                        cells.end(), neigh.value());

                    /* if not, save it */
                    if (not present)
                        layer_cells[tag].insert(neigh.value());
                }
            }
        }

        /* move the found neighbours to the subdomain */
        for (auto& [tag, cells] : layer_cells) {
            subdomains[tag].insert(subdomains[tag].end(), cells.begin(), cells.end());
            std::sort(subdomains[tag].begin(), subdomains[tag].end());
        }
    }

    return subdomains;
}

template<typename Mesh>
auto
detect_dd_faces(const Mesh& msh, const sdmap_t<Mesh>& subdomains, bool include_dd_bnd) {
    std::map<size_t, std::vector<bool>> sd_flags;
    for (auto& [tag, cells] : subdomains) {
        auto& flags = sd_flags[tag];
        if (flags.size() == 0)
            flags.resize(msh.cells_size());
        for (auto& cl : cells)
            flags[offset(msh, cl)] = true; 
    }

    auto cvf = connectivity_via_faces(msh);

    std::map<size_t, std::vector<bool>> f_flags;
    for (auto& [tag, cells] : subdomains) {
        if (f_flags[tag].size() == 0)
            f_flags[tag].resize( msh.faces_size() );

        for (auto& cl : cells) {
            auto fcs = faces(msh, cl);
            for (auto& fc : fcs) {
                auto neigh = cvf.neighbour_via(msh, cl, fc);
                if (not neigh)
                    continue;

                auto neigh_id = offset(msh, neigh.value());
                if (sd_flags[tag][neigh_id] or include_dd_bnd) {
                    f_flags[tag][offset(msh, fc)] = true;
                }
            }
        }
    }

    if constexpr (Mesh::dimension == 2) {
        for (auto& [tag, ifcs] : f_flags ) {
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

template<typename Mesh>
void
hhodd(const Mesh& msh, size_t levels)
{
    auto subdomains = make_overlapping_subdomains(msh, levels);
    detect_dd_faces(msh, subdomains, false);
    
    disk::silo_database silo;
    silo.create("hhodd.silo");
    silo.add_mesh(msh, "mesh");

    for (auto& [tag, cells] : subdomains)
    {
        std::vector<double> yesno(msh.cells_size());
        for (auto& cl : cells)
            yesno[offset(msh, cl)] += 1.0;

        std::stringstream ss;
        ss << "domain" << tag;
        silo.add_variable("mesh", ss.str(), yesno, disk::zonal_variable_t);
    }

    silo.close();
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
