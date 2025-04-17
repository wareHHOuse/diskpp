/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <cstddef>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <filesystem>

#include "diskpp/common/util.h"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho"
#include "rasdd.hpp"

#include "diskpp_git_revision.h"


namespace disk {


template<typename Mesh>
struct source;

template<mesh_2D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return 2.0*M_PI*M_PI*sx*sy;
    }
};


template<mesh_3D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return 3.0*M_PI*M_PI*sx*sy*sz;
    }
};


template<typename Mesh>
struct solution;

template<mesh_2D Mesh>
struct solution<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return sx*sy;
    }
};


template<mesh_3D Mesh>
struct solution<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return sx*sy*sz;
    }
};


} // namespace disk



enum class hho_variant {
    equal_order,
    mixed_order
};

enum class ras_mode {
    bicgstab,
    iterate
};

enum class mesh_type {
    triangles,
    cartesian,
    hexas,
    tetras,
    undefined,
};
struct solver_config {
    size_t          overlap;
    size_t          ras_maxiter;
    size_t          degree;
    bool            ras_debug;
    hho_variant     variant;
    ras_mode        mode;
    std::string     fn_bicg_hist;
    std::string     fn_silo;
    std::string     fn_err_hist;
    std::string     outdir;
    mesh_type       meshtype;
    size_t          imesh_reflevels;
    size_t          imesh_partitions;
};

struct iterdata {
    double error;
    double residual;
};

std::ostream&
operator<<(std::ostream& os, const iterdata& id)
{
    os << id.error << " " << id.residual;
    return os;
};

template<typename Mesh>
void
diffusion_solver(const Mesh& msh, const solver_config& scfg)
{
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using hho_space = hho_space<Mesh>;
    using T = typename hho_space::scalar_type;
    using cell_basis_type = typename hho_space::cell_basis_type;
    using face_basis_type = typename hho_space::face_basis_type;

    bool mixed_order = (scfg.variant == hho_variant::mixed_order);

    auto face_degree = scfg.degree;
    auto cell_degree = mixed_order ? scfg.degree + 1 : scfg.degree;

    degree_info di(cell_degree, face_degree);

    disk::source<Mesh> f;

    auto assm = make_assembler(msh, di);

    auto sizeF = face_basis_type::size_of_degree(di.face);

    using dm = disk::dynamic_matrix<T>;
    using dv = disk::dynamic_vector<T>;
    std::vector<std::pair<dm, dv>> local_contribs;
    local_contribs.reserve(msh.cells_size());

    /* MATRIX ASSEMBLY */
    std::cout << "Assembling global matrix..." << std::endl;
    for (auto& cl : msh)
    {
        dm lhs;
        if (mixed_order) {
            auto [R, A] = local_operator(msh, cl, di);
            auto S = local_stabilization_hdg(msh, cl, di);
            lhs = A+S;
        } else {
            auto [R, A] = local_operator(msh, cl, di);
            auto S = local_stabilization(msh, cl, di, R);
            lhs = A+S;
        }

        auto phiT = hho_space::cell_basis(msh, cl, di.cell);
        dv rhs = integrate(msh, cl, f, phiT);
        auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);
        assm.assemble(msh, cl, lhsc, rhsc);
        local_contribs.push_back({lhs, rhs});
    }
    assm.finalize();

    std::cout << "Unknowns: " << assm.LHS.rows() << " ";
    std::cout << "Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    /* PREPARE RAS */
    rasdd ras(msh, sizeF, scfg.overlap);
    ras.prepare( assm.LHS, assm.dirichlet_faces_flags() );
    if (scfg.ras_debug)
        ras.save_debug_data();

    /* HELPER FUNCTION FOR POSTPRO */
    auto postpro = [&](const std::string& filename,
        const dv& sol, const dv& r) {
        std::vector<T> u_data;
        std::vector<T> r_data;
        T error = 0.0;
        disk::solution<Mesh> u_sol;
        for(size_t cell_i = 0; cell_i < msh.cells_size(); cell_i++)
        {
            auto cl = msh.cell_at(cell_i);   
            auto& [lhs, rhs] = local_contribs[cell_i];
            auto phiT = hho_space::cell_basis(msh, cl, di.cell);
            dv sol_ana = local_reduction(msh, cl, di, u_sol);
            auto locsolF = assm.take_local_solution(msh, cl, sol);
            dv locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
            u_data.push_back(locsol(0));
            dv diff = locsol - sol_ana;
            error += diff.dot(lhs*diff);

            /* This does not make much sense, as the residual is not
             * in the solution space. It is however useful since it
             * allows to check that things are zero where they should
             * be zero. */
            auto locresF = assm.take_local_solution(msh, cl, r);
            dv zero = dv::Zero(rhs.size());
            dv locres = disk::hho::deschur(lhs, zero, locresF, phiT);
            r_data.push_back(locres(0));
        }
        if (filename != "") {
            disk::silo_database silo_db;
            silo_db.create(filename);
            silo_db.add_mesh(msh, "mesh");
            silo_db.add_variable("mesh", "u", u_data, disk::zonal_variable_t);
            silo_db.add_variable("mesh", "residual", r_data, disk::zonal_variable_t);
        }
        return std::sqrt(error);
    };
    
    /* SOLVE */
    if (scfg.mode == ras_mode::iterate) {
        std::vector<iterdata> ids;
        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
        for (size_t niter = 0; niter < scfg.ras_maxiter; niter++) {
            iterdata id;
            disk::dynamic_vector<T> r = assm.RHS - assm.LHS*sol;
            sol = sol + ras(r);
            id.residual = r.norm();
            std::cout << "Postprocessing..." << std::endl;
            std::stringstream ss;
            
            if (scfg.fn_silo != "") {
                ss << scfg.fn_silo << "_iter" << niter << ".silo";
                std::cout << "  Saving solution to " << ss.str() << std::endl;
            }
            
            id.error = postpro(ss.str(), sol, r);

            std::cout << "  A-norm error: " << id.error << std::endl;
            ids.push_back(id);
        }

        std::cout << "  Saving error history to " << scfg.fn_err_hist << std::endl;
        std::ofstream ofs(scfg.fn_err_hist);
        for (auto& id : ids)
            ofs << id << std::endl;
    }
    else {
        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
        bicgstab_io<T> bio;
        bio.verbose = true;
        sol = bicgstab(assm.LHS, assm.RHS, ras, bio);
        std::cout << "Postprocessing..." << std::endl;
        std::cout << "  Saving BiCGStab history to " << scfg.fn_bicg_hist << std::endl;
        std::ofstream ofs(scfg.fn_bicg_hist);
        for (auto& rr : bio.history)
            ofs << rr << std::endl;

        if (scfg.fn_silo != "") {
            std::cout << "  Saving solution to " << scfg.fn_silo << std::endl;
        }
        
        auto err = postpro(scfg.fn_silo, sol, bio.residual);
        std::cout << "  A-norm error: " << err << std::endl;
    }
}

template<typename Mesh>
void
partition_unit_square_mesh(Mesh& msh, size_t np)
{
    if (np < 2)
        return;

    auto storage = msh.backend_storage();
    for (size_t cell_i = 0; cell_i < msh.cells_size(); cell_i++) {
        auto cl = msh.cell_at(cell_i);
        auto bar = barycenter(msh, cl);
        auto domxy = bar * np;

        size_t subdom = 0;
        if constexpr (Mesh::dimension == 1)
            subdom = size_t(domxy.x());
        if constexpr (Mesh::dimension == 2)
            subdom = size_t(domxy.x()) + np*size_t(domxy.y());
        if constexpr (Mesh::dimension == 3)
            subdom = size_t(domxy.x()) + np*(size_t(domxy.y()) + np*size_t(domxy.z()));

        storage->subdomain_info[cell_i] = disk::subdomain_descriptor(subdom);
    }

    disk::make_interpartition_boundaries(msh);
    disk::renumber_hypercube_boundaries(msh);
}

void
run_on_internal_mesh(const solver_config& scfg)
{
    using T = double;

    if (scfg.meshtype == mesh_type::triangles) {
        disk::simplicial_mesh<T,2> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < scfg.imesh_reflevels; l++)
            mesher.refine();
        
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver(msh, scfg);
    }

    if (scfg.meshtype == mesh_type::cartesian) {
        disk::cartesian_mesh<T,2> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < scfg.imesh_reflevels; l++)
            mesher.refine();
        
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver(msh, scfg);
    }

    if (scfg.meshtype == mesh_type::hexas) {
        disk::generic_mesh<T,2> msh;
        auto mesher = make_fvca5_hex_mesher(msh);
        mesher.make_level(scfg.imesh_reflevels);

        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver(msh, scfg);
    }

    if (scfg.meshtype == mesh_type::tetras) {
        disk::simplicial_mesh<T,3> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < scfg.imesh_reflevels; l++)
            mesher.refine();
        
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver(msh, scfg);
    }
}

int main(int argc, char **argv)
{
    std::cout << "HHO-DD solver, revision " << GIT_REVISION << std::endl;
    
    rusage_monitor rm(false);

    solver_config scfg;
    scfg.overlap = 1;
    scfg.ras_maxiter = 20;
    scfg.degree = 0;
    scfg.ras_debug = false;
    scfg.variant = hho_variant::equal_order;
    scfg.mode = ras_mode::bicgstab;
    scfg.fn_bicg_hist = "bicgstab.txt";
    scfg.fn_err_hist = "error.txt";
    scfg.fn_silo = "";
    scfg.outdir = "";
    scfg.meshtype = mesh_type::undefined;
    scfg.imesh_partitions = 2;
    scfg.imesh_reflevels = 2;

    int opt;
    while ((opt = getopt(argc, argv, "b:d:e:i:k:m:o:p:r:s:DIMR")) != -1) {
        switch (opt) {

        case 'b': /* BiCGStab history filename (bicgstab mode) */
            scfg.fn_bicg_hist = optarg;
            break;

        case 'd': /* output data directory */
            scfg.outdir = optarg;
            break;

        case 'e': /* Error history filename (iter mode) */
            scfg.fn_err_hist = optarg;
            break;

        case 'i': {
            int optval = std::max(0, atoi(optarg));
            scfg.ras_maxiter = optval;
            } break;
        
        case 'k': { /* hho method degree */
            int optval = std::max(0, atoi(optarg));
            scfg.degree = optval;
            } break;

        case 'm': { /* internal mesh type */
            if ( std::string(optarg) == "tri" )
                scfg.meshtype = mesh_type::triangles;
            else if ( std::string(optarg) == "quad" )
                scfg.meshtype = mesh_type::cartesian;
            else if ( std::string(optarg) == "hex" )
                scfg.meshtype = mesh_type::hexas;
            else if ( std::string(optarg) == "tet" )
                scfg.meshtype = mesh_type::tetras;
            else
                std::cout << "Warning: Wrong mesh type '" << optarg << "'" << std::endl;
            } break;

        case 'o': { /* number of overlap layers */
            int optval = std::max(0, atoi(optarg));
            scfg.overlap = optval;
            } break;

        case 'p': { /* partitions for internal mesh */
            int optval = std::max(0, atoi(optarg));
            scfg.imesh_partitions = optval;
            } break;

        case 'r': { /* refinement levels for internal mesh */
            int optval = std::max(0, atoi(optarg));
            scfg.imesh_reflevels = optval;
            } break;
            
        case 's': /* silo filename or prefix */
            scfg.fn_silo = optarg;
            break;

        case 'D': /* Output subdomain debugging data */
            scfg.ras_debug = true;
            break;
        
        case 'I': /* Iterate instead of using BiCGStab */
            scfg.mode = ras_mode::iterate;
            break;

        case 'M': /* enable mixed order */
            scfg.variant = hho_variant::mixed_order;
            break;

        case 'R': /* Enable resource usage reporting */
            rm.enabled(true);
            break;
        }
    }
    argc -= optind;
    argv += optind;

    if ( (argc < 1) and (scfg.meshtype == mesh_type::undefined) ) {
        std::cout << "Usage: hhodd [options] [filename]\n\n";
        std::cout << " If a mesh is provided via [filename] is, the currently supported formats are:\n";
        std::cout << "   * GMSH .geo (not .msh): if a simplicial 2D mesh is wanted the GMSH script\n";
        std::cout << "     must have extension `.geo2s`, whereas if a 3D simplicial mesh is wanted\n";
        std::cout << "     the script must have `.geo3s` extension\n";
        std::cout << "   * FVCA5-style meshes in .typ1 format\n";
        std::cout << "   * FVCA6-style meshes in .msh format\n\n";
        std::cout << " Meshes can be also generated by the internal meshers, see -m below. In any\n";
        std::cout << " case the code assumes that the domain is the unit hypercube [0,1]^d.\n\n";
        std::cout << " Other options:\n";
        std::cout << "   -b <str> : BiCGStab history filename (default: bicgstab.txt)\n";
        std::cout << "   -d <str> : Directory for all the output files\n";
        std::cout << "   -e <str> : Error history filename (default: error.txt)\n";
        std::cout << "   -i <int> : RAS iterations (in iter mode)\n";
        std::cout << "   -k <int> : HHO degree\n";
        std::cout << "   -m <std> : Use internal meshers (tri, quad, hex or tet)\n";
        std::cout << "   -o <int> : Number of overlap layers\n";
        std::cout << "   -p <int> : Partitions for internal mesh\n";
        std::cout << "   -r <int> : Refinement levels for internal mesh\n";
        std::cout << "   -s <str> : Silo filename (bicgstab mode) or prefix (iter mode)\n";
        std::cout << "   -D       : Output subdomain debugging data\n";
        std::cout << "   -I       : Iterate instead of using BiCGStab\n";
        std::cout << "   -M       : Enable mixed-order HHO\n";
        std::cout << "   -R       : Report resource usage at exit\n";
        return 1;
    }

    if (scfg.outdir != "") {
        std::filesystem::create_directory(scfg.outdir);
        scfg.fn_bicg_hist = scfg.outdir + "/" + scfg.fn_bicg_hist;
        scfg.fn_err_hist = scfg.outdir + "/" + scfg.fn_err_hist;
        if (scfg.fn_silo != "")
            scfg.fn_silo = scfg.outdir + "/" + scfg.fn_silo;
    }

    const char *mesh_filename = argv[0];
    if ( mesh_filename and (scfg.meshtype != mesh_type::undefined) ) {
        std::cout << "Warning: both '-m' and a filename were specified. Ignoring -m.\n";
    }

    using T = double;

    if (scfg.meshtype != mesh_type::undefined)
    {
        run_on_internal_mesh(scfg);
        return 0;
    }

    assert(argc > 0);

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<T,2>;
        mesh_type msh;
        disk::gmsh_geometry_loader< mesh_type > loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        diffusion_solver(msh, scfg);
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
        diffusion_solver(msh, scfg);
        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        disk::generic_mesh<T,2> msh;
        disk::load_mesh_fvca5_2d<T>(mesh_filename, msh);
        //partition_unit_square_mesh(msh, scfg.imesh_partitions);
        disk::simplicial_mesh<T,2> fmsh;
        submesh_via_gmsh(msh, fmsh, 0.1);
        diffusion_solver(fmsh, scfg);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver(msh, scfg);
        return 0;
    }

    std::cout << "Unrecognized mesh format. Exiting." << std::endl;
    return 1;
}
