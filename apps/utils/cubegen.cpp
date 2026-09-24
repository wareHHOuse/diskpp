#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/loaders/loader.hpp"

int main(int argc, char **argv)
{
    using T = double;

    size_t      levels = 0;
    size_t      degree = 0;
    char *      mesh_filename = nullptr;

    int ch;
    while ( (ch = getopt(argc, argv, "r:k:m:")) != -1 )
    {
        switch(ch)
        {
            case 'r':
                levels = std::stoull(optarg);
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    using mesh_type = disk::cartesian_mesh<T,3>;

    //using mesh_type = disk::generic_mesh<T,2>;

    mesh_type msh;
    auto mesher = disk::make_simple_mesher(msh);

    std::vector<double> dids;
    size_t i = 0;
    for (i = 0; i < levels; i++) {
        mesher.refine();
    }

    //std::cout << "Diameter: " << disk::average_diameter(msh) << std::endl;
    msh.statistics();
    std::string silo_fn = "cubegen.silo";
    disk::silo_database db;
    db.create(silo_fn);
    db.add_mesh(msh, "mesh");
    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        dids.push_back( di.tag() );

        
        std::cout << cl << std::endl;
        std::cout << "  - measure = " << measure(msh, cl) << std::endl;
        std::cout << "  - barycenter = " << barycenter(msh, cl) << std::endl;
    
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs) {
            std::cout << "  " << fc << std::endl;
            std::cout << "    - measure = " << measure(msh, fc) << std::endl;
            std::cout << "    - barycenter = " << barycenter(msh, fc) << std::endl;
        }
    }
    db.add_variable("mesh", "domain_ids", dids, disk::zonal_variable_t);
}