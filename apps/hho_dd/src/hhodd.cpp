#include <iostream>
#include <regex>
#include <set>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

template<typename Mesh>
auto
subdomains(const Mesh& msh)
{
    std::set<size_t> sds;
    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        sds.insert(di.tag());
    }

    return sds; 
}

template<typename Mesh>
void
hhodd(const Mesh& msh)
{
    auto subdoms = subdomains(msh);
    auto cvf = connectivity_via_faces(msh);

    std::map<size_t, std::vector<double>> cellmaps;
    for (auto& sd : subdoms) {
        cellmaps[sd].resize(msh.cells_size());
    }

    size_t clnum = 0;
    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        cellmaps[di.tag()][clnum] = 1.0;

        auto fcs = faces(msh, cl);
        for (auto& fc : fcs) {

            auto bi = msh.boundary_info(fc);
            if (!bi.is_boundary())
                continue;
            
            auto [neigh, valid] = cvf.neighbour_via(msh, cl, fc);
            if (!valid)
                continue;

            auto nofs = offset(msh, neigh);
            cellmaps[di.tag()][nofs] = 1;
        }

        clnum++;
    }


    disk::silo_database silo;
    silo.create("hhodd.silo");
    silo.add_mesh(msh, "mesh");

    for (auto& [tag, cellmap] : cellmaps)
    {
        std::stringstream ss;
        ss << "domain" << tag;
        silo.add_variable("mesh", ss.str(), cellmap, disk::zonal_variable_t);
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
