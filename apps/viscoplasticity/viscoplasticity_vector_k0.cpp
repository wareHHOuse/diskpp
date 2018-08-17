#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"

#include "viscoplasticity_vector_solver_k0.hpp"
template<typename T>
auto
run_viscoplasticity(size_t degree,
                    const T & alpha,
                    const problem_type& problem,
                    const std::string & other_info)
{
    /* Medit 2d*/
    std::cout << "Guessed mesh format: Medit format" << std::endl;
    typedef disk::generic_mesh<T, 2>  mesh_type;

    T tolerance = 1.e-10, Ninf = 10.e+5;
    size_t max_iters = 1; //50000;

    std::string name, filename;
    switch (problem)
    {
        case DRIVEN:
            name = "driven";
            filename = "../../../diskpp/meshes/2D_quads/medit/square_h0025.medit2d";
            break;
        case COUETTE:
            name = "couette";
            filename = "../../../diskpp/meshes/2D_triangles/medit/couronne_01.medit2d";
            break;
        case VANE:
            name = "vane";
            filename = "../../../diskpp/meshes/2D_triangles/medit/vane05.medit2d";
            break;
        default:
            std::cout << "wrong arguments" << std::endl;
            exit(1);
    }

    mesh_type                     msh;
    disk::medit_mesh_loader<T, 2> loader;
    loader.read_mesh(filename);
    loader.populate_mesh(msh);

    std::string info = name + "_k" + tostr(degree) + "_a" + tostr(alpha) + other_info;
    std::ofstream ofs("errors_" + info + ".data");

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    typename revolution::hho_degree_info hdi(degree, degree);
    augmented_lagrangian_viscoplasticity<mesh_type> als(msh, hdi, alpha);

    std::cout << "before define problem" << std::endl;
    auto assembler = als.define_problem(msh, problem);
    std::cout << "after define problem" << std::endl;
    als.initialize(msh, assembler);

    size_t i;
    for(i = 0; i < max_iters; i++)
    {
        als.run_stokes_like(msh, assembler, i);

        //als.update_multiplier(msh, assembler);
        //auto error = als.compute_errors(msh, assembler, false);

        T cvg_total (0.), cvg_stress(0.), cvg_gamma(0.);
        std::tie(cvg_total, cvg_stress, cvg_gamma) = als.convergence;

        if(i % 1000 == 0)
        {
            std::cout << "  i : "<< i<<"  - " << std::sqrt(cvg_total)<<std::endl;
            als.post_processing( msh, assembler, info +"_i" + tostr(i), problem);
            std::cout << "done" << std::endl;
        }

        assert(std::sqrt(cvg_total) < Ninf);
        if( std::sqrt(cvg_total)  < tolerance)
            break;
    }
    ofs.close();

    std::cout << "Finish" << std::endl;
    auto final_error = als.compute_errors(msh, assembler, true);
    als.post_processing( msh, assembler, info +"_i" + tostr(i), problem);

    return final_error;
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *word      = nullptr;
    int ch;
    size_t degree = 0;
    RealType alpha = 1.;

    problem_type problem = DRIVEN;

    while ( (ch = getopt(argc, argv, "k:a:dcv")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case 'a':
                alpha = atof(optarg);
                if (alpha <= 0)
                {
                    std::cout << "alpha must be >=0. Falling back to 1." << std::endl;
                    alpha = 1.;
                }
            case 'd':
                problem = DRIVEN;
                break;

            case 'c':
                problem = COUETTE;
                std::cout << "couette chosen" << std::endl;
                break;
            case 'v':
                problem = VANE;
                std::cout << "vane chosen" << std::endl;
                break;
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    word = argv[0];

    if (word == nullptr)
    {
        std::cout << "no word specified" << std::endl;
        return 1;
    }

    run_viscoplasticity(degree, alpha, problem, word);

    return 0;
}
//#endif
