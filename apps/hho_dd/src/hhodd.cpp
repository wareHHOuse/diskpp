#include <iostream>
#include <regex>
#include <set>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

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

    /* Collect the cells in the various subdomains */
    std::map<size_t, std::vector<cell_type>> subdomains;
    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        subdomains[di.tag()].push_back(cl);
    }
    /* the cells are stored ordered in msh, therefore each
     * vector in the map is sorted too. this is necessary
     * precondition for std::lower_bound() below. */

    for (size_t ol = 0; ol < overlap_layers; ol++) {
        std::map<size_t, std::set<cell_type>> layer_cells;
        /* Iterate on the currently identified subdomains */
        for (auto& [tag, cells] : subdomains) {
            /* for each cell */
            for (auto& cl : cells) {
                auto fcs = faces(msh, cl);
                for (auto& fc : fcs) {
                    /* for each neighbour */
                    auto neigh_opt = cvf.neighbour_via(msh, cl, fc);
                    if (neigh_opt) {
                        /* and check if it is already in the current subdomain */
                        bool present = std::binary_search(cells.begin(),
                            cells.end(), neigh_opt.value());

                        /* if not, save it */
                        if (not present)
                            layer_cells[tag].insert(neigh_opt.value());
                    }
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
void
hhodd(const Mesh& msh)
{
    auto subdomains = make_overlapping_subdomains(msh, 5);
    

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
    if (argc < 2) {
        std::cout << "missing filename" << std::endl;
        return 1;
    }

    const char *mesh_filename = argv[1];

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

        hhodd(msh);

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

        hhodd(msh);

        return 0;
    }

}
