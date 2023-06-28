#include <iostream>

#include <cstdio>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/geometry/geometry.hpp"


int main(int argc, char **argv)
{
    using T = double;

    if (argc < 2)
    {
        std::cout << "Please specify the CSV of the polygon" << std::endl;
        return 1;
    }

    disk::generic_mesh<T,2> msh;
    disk::load_single_element_csv(msh, argv[1]);
    
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) {
        std::cout << "Can't open a pipe to gnuplot" << std::endl;
        return -1;
    }

    for (auto& cl : msh)
    {
        auto tris = disk::triangulate_nonconvex_polygon(msh, cl);
        fprintf(gp, "$DATA << EOD\n");
        for (auto& tri : tris)
        {
            auto p0 = tri.p0; auto p1 = tri.p1; auto p2 = tri.p2;
            fprintf(gp, "%f %f %f %f\n", p0.x(), p0.y(), (p1-p0).x(), (p1-p0).y());
            fprintf(gp, "%f %f %f %f\n", p1.x(), p1.y(), (p2-p1).x(), (p2-p1).y());
            fprintf(gp, "%f %f %f %f\n", p2.x(), p2.y(), (p0-p2).x(), (p0-p2).y());
        }
        fprintf(gp, "EOD\n");
        
        auto bar = barycenter(msh, cl);
        
        fprintf(gp, "plot $DATA with vectors nohead title 'Triangulation',");
        fprintf(gp, " '-' w p lc rgb 'red' pt 4 title 'Barycenter'\n%f %f\ne\n",
            bar.x(), bar.y());

        break;
    }
    fflush(gp);
    pclose(gp);
    return 0;
}