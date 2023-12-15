/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <vector>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <map>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "mumps.hpp"
#include "diskpp/output/silo.hpp"

#include "diskpp/methods/implementation_hho/curl.hpp"

#include "diskpp/common/timecounter.hpp"

#include "sol/sol.hpp"
#include "diskpp/solvers/feast.hpp"

#include "sgr.hpp"

#include <Eigen/Eigenvalues>


#include "lua/mesh.hpp"
#include "lua/hho.hpp"
#include "lua/feast.hpp"


template<typename T>
struct steklov_solver_configuration {
    disk::lua::mesh_parameters          mesh;
    disk::lua::hho_parameters           hho;
    disk::feast_eigensolver_params<T>   feast;
};




template<typename T>
void register_my_usertypes(sol::state& lua, steklov_solver_configuration<T>& config)
{
    disk::lua::register_mesh_usertypes(lua);
    disk::lua::register_hho_usertypes(lua);
    disk::lua::register_feast_usertypes<T>(lua);

    using ssc_t = steklov_solver_configuration<T>;
    sol::usertype<ssc_t> ssct = lua.new_usertype<ssc_t>(
        "steklov_solver_configuration",
        sol::constructors<ssc_t()>()
    );
    ssct["mesh"] = &ssc_t::mesh;
    ssct["hho"] = &ssc_t::hho;
    ssct["feast"] = &ssc_t::feast;
    lua["config"] = &config;
}

template<typename T>
bool
lua_init_context(sol::state& lua, steklov_solver_configuration<T>& config)
{
    lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::io, sol::lib::table, sol::lib::string);

    register_my_usertypes(lua, config);

    lua["mesh"] = lua.create_table();
    lua["hho"] = lua.create_table();
    lua["feast"] = lua.create_table();
    lua["boundary"] = lua.create_table();

    auto sf = lua.safe_script_file("steklov.lua");

    return sf.valid();
}


template<typename Mesh>
class hho_assembler_steklov
{
    using coordinate_type = typename Mesh::coordinate_type;
    using scalar_type = coordinate_type;
    using mesh_type = Mesh;
    using face_type = typename Mesh::face_type;
    using cell_type = typename Mesh::cell_type;

    //Mesh                                msh;
    disk::hho_degree_info               hdi;

    std::vector<Triplet<scalar_type>>   tripletsL;
    std::vector<Triplet<scalar_type>>   tripletsR;
    std::vector<bool>                   is_dirichlet;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    const size_t INVALID_OFFSET = (size_t) ~0;

    size_t face_basis_size(void) const {
        return disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    }

    /* Get the offset of a face in the linear system */
    size_t get_system_offset(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(compress_table.size() == msh.faces_size());
        auto cnum = compress_table[face_num];
        assert(cnum != INVALID_OFFSET);
        auto fbs = face_basis_size();
        return cnum*fbs;
    }

    /* Determine if a face should be assembled */
    bool is_in_system(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto bi = msh.boundary_info(fc);
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(is_dirichlet.size() == msh.faces_size());
        return not (bi.is_boundary() and is_dirichlet[face_num]);
    }

    void make_tables(mesh_type& msh)
    {
        compress_table.resize( msh.faces_size(), INVALID_OFFSET);
        expand_table.resize( sysfcs );

        size_t face_i = 0;
        size_t compressed_ofs = 0;
        for (auto& fc : faces(msh))
        {
            assert(compressed_ofs <= face_i);
            if ( is_in_system(msh, fc) )
            {
                assert(face_i < compress_table.size());
                compress_table[face_i] = compressed_ofs;
                assert(compressed_ofs < expand_table.size());
                expand_table[compressed_ofs] = face_i;
                compressed_ofs++;
            }

            face_i++;
        }

        assert(face_i == msh.faces_size());
    }

public:
    
    SparseMatrix<scalar_type>           LHS;
    SparseMatrix<scalar_type>           RHS;

    size_t                              syssz;
    size_t                              sysfcs;


    hho_assembler_steklov()
    {}

    hho_assembler_steklov(Mesh& msh, const disk::hho_degree_info& hdi,
        const std::vector<bool>& is_dirichlet)
    {
        initialize(msh, hdi, is_dirichlet);
    }

    void clear()
    {
        tripletsL.clear();
        tripletsR.clear();
        is_dirichlet.clear();
        compress_table.clear();
        expand_table.clear();
        syssz = 0;
        sysfcs = 0;
    }

    void initialize(Mesh& msh, const disk::hho_degree_info& p_hdi,
                    const std::vector<bool>& p_is_dirichlet)
    {
        clear();

        is_dirichlet = p_is_dirichlet;
        hdi = p_hdi;

        auto fbs = face_basis_size();

        auto in_system = [&](const face_type& fc) -> bool {
            auto ofs = offset(msh, fc);
            assert(ofs < is_dirichlet.size());
            return not (msh.is_boundary(fc) and is_dirichlet[ofs]);
        };

        sysfcs = std::count_if(msh.faces_begin(), msh.faces_end(), in_system);
        syssz = fbs*sysfcs;

        make_tables(msh);

        LHS = SparseMatrix<scalar_type>(syssz, syssz);
        RHS = SparseMatrix<scalar_type>(syssz, syssz);

        std::cout << "Assembler initialized: " << sysfcs << " faces in system, ";
        std::cout << syssz << " DoFs" << std::endl;
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<scalar_type, Dynamic, Dynamic>& lhsc,
             const Matrix<scalar_type, Dynamic, Dynamic>& rhs)
    {
        auto fbs = face_basis_size();

        auto fcs = faces(msh, cl);
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( not is_in_system(msh, fcs[fi]) )
                continue;

            auto cofsi = get_system_offset(msh, fcs[fi]);
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                auto lofsi = fi*fbs;
                auto lofsj = fj*fbs;

                if ( not is_in_system(msh, fcs[fj]) ) 
                    continue;

                auto cofsj = get_system_offset(msh, fcs[fj]);
                for (size_t i = 0; i < fbs; i++)
                    for(size_t j = 0; j < fbs; j++)
                        tripletsL.push_back( Triplet<scalar_type>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
            }

            auto lofs = fi*fbs;
            for (size_t i = 0; i < fbs; i++)
                for(size_t j = 0; j < fbs; j++)
                    tripletsR.push_back( Triplet<scalar_type>(cofsi+i, cofsi+j, rhs(lofs+i, lofs+j)) );
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(tripletsL.begin(), tripletsL.end());
        RHS.setFromTriplets(tripletsR.begin(), tripletsR.end());
        tripletsL.clear();
        tripletsR.clear();
        std::cout << "LHS has " << LHS.nonZeros() << " nonzeros." << std::endl; 
        std::cout << "RHS has " << RHS.nonZeros() << " nonzeros." << std::endl; 
    }

    disk::dynamic_vector<scalar_type>
    get_expanded_solution(const Mesh& msh, disk::dynamic_vector<scalar_type>& sol)
    {
        auto fbs = face_basis_size();

        disk::dynamic_vector<scalar_type> ret = 
            disk::dynamic_vector<scalar_type>::Zero( fbs*msh.faces_size() );

        for (size_t i = 0; i < sysfcs; i++)
        {
            auto in_offset = i*fbs;
            auto out_offset = expand_table.at(i)*fbs;
            ret.segment(out_offset, fbs) = sol.segment(in_offset, fbs);
        }

        return ret;
    }

    disk::dynamic_matrix<scalar_type>
    get_element_dofs(const Mesh& msh, const typename Mesh::cell& cl,
        disk::dynamic_matrix<scalar_type>& sol)
    {
        auto fbs = face_basis_size();
        auto fcs = faces(msh, cl);
        disk::dynamic_matrix<scalar_type> ret = 
            disk::dynamic_matrix<scalar_type>::Zero( fbs*fcs.size(), sol.cols() );

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            if ( not is_in_system(msh, fc) )
                continue;

            auto ofs = get_system_offset(msh, fc);
            ret.block(i*fbs, 0, fbs, sol.cols()) = sol.block(ofs, 0, fbs, sol.cols());
        }

        return ret;
    }
};


enum class boundary_type {
    internal,
    dirichlet,
    neumann,
    robin,
    steklov
};

template<typename Mesh>
struct hho_steklov_solver_state
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename Mesh::coordinate_type;
    using assembler_type = hho_assembler_steklov<Mesh>;
    using vector_type = disk::dynamic_vector<scalar_type>;

    mesh_type                           msh;
    assembler_type                      assm;
    disk::dynamic_vector<scalar_type>   eigvals;
    disk::dynamic_matrix<scalar_type>   eigvecs;
    disk::dynamic_matrix<scalar_type>   eigvecs_full;

    std::vector<boundary_type>          bndtype;
    std::vector<bool>                   is_dirichlet;

    disk::hho_degree_info               hdi;
};

boundary_type
lua_get_boundary_type(sol::state& lua, size_t bndid)
{
    sol::optional<std::string> bndtype_opt = lua["boundary"][bndid];
    if (not bndtype_opt)
        return boundary_type::dirichlet;

    std::string bndtype = bndtype_opt.value();
    if (bndtype == "dirichlet")
        return boundary_type::dirichlet;
    if (bndtype == "neumann")
        return boundary_type::neumann;
    if (bndtype == "robin")
        return boundary_type::robin;
    if (bndtype == "steklov")
        return boundary_type::steklov;

    return boundary_type::dirichlet;
}


template<typename Config, typename State>
static int
assemble(const Config& config, State& state)
{
    using mesh_type = typename State::mesh_type;
    using scalar_type = typename State::scalar_type;
    using dm = disk::dynamic_matrix<scalar_type>;
    using dv = disk::dynamic_vector<scalar_type>;

    auto cd = state.hdi.cell_degree();
    auto fd = state.hdi.face_degree();
    auto rd = state.hdi.reconstruction_degree();
    auto cbs = disk::scalar_basis_size(cd, mesh_type::dimension);
    auto fbs = disk::scalar_basis_size(fd, mesh_type::dimension-1);
    auto rbs = disk::scalar_basis_size(rd, mesh_type::dimension);

    state.assm.initialize(state.msh, state.hdi, state.is_dirichlet);

    for (auto& cl : state.msh)
    {
        /* Make standard HHO Laplacian */
        auto [GR, A] = make_scalar_hho_laplacian(state.msh, cl, state.hdi);
        dm S = make_scalar_hho_stabilization(state.msh, cl, GR, state.hdi);
        dm L = A+S;

        /* Make face mass matrices */
        auto fcs = faces(state.msh, cl);

        dm BFF = dm::Zero(fcs.size()*fbs, fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto& fc = fcs[face_i];
            auto gfnum = offset(state.msh, fc);
            assert(gfnum < state.bndtype.size());
            auto bndtype = state.bndtype[gfnum];

            switch (bndtype) {
                case boundary_type::dirichlet:
                case boundary_type::neumann:
                    continue;
                    break;

                case boundary_type::robin: {
                    auto fb = disk::make_scalar_monomial_basis(state.msh, fc, fd);
                    auto idx = cbs+fbs*face_i;
                    L.block(idx, idx, fbs, fbs) +=
                        disk::make_mass_matrix(state.msh, fc, fb);
                } break;

                case boundary_type::steklov: {
                    auto fb = disk::make_scalar_monomial_basis(state.msh, fc, fd);
                    auto idx = fbs*face_i;
                    L.block(cbs+idx, cbs+idx, fbs, fbs) +=
                        disk::make_mass_matrix(state.msh, fc, fb);
                    BFF.block(idx, idx, fbs, fbs) =
                        disk::make_mass_matrix(state.msh, fc, fb);
                } break;

                case boundary_type::internal:
                    break;

                default:
                    break;
            }
        }

        dm LC = disk::static_condensation(L, cbs);
        state.assm.assemble(state.msh, cl, LC, BFF);
    }

    state.assm.finalize();

    return 0;
}

template<typename Config, typename State>
static int
solve(const Config& config, State& state)
{
    using mesh_type = typename State::mesh_type;
    using scalar_type = typename State::scalar_type;
    using dm = disk::dynamic_matrix<scalar_type>;
    using dv = disk::dynamic_vector<scalar_type>;

    auto cd = state.hdi.cell_degree();
    auto cbs = disk::scalar_basis_size(cd, mesh_type::dimension);

    auto ret = disk::feast(config.feast, state.assm.LHS, state.assm.RHS,
        state.eigvecs, state.eigvals);

    if (disk::feast_status::success != ret) {
        std::cout << "FEAST algorithm did not converge: "<< ret << std::endl;
        return -1;
    }

    std::cout << "Eigenvalues found: " << state.eigvals.size() << std::endl;
    auto found_eigs = state.eigvals.size();

    if (found_eigs == 0)
        return 0;

    
    state.eigvecs_full = dm::Zero(cbs*state.msh.cells_size(), found_eigs);
    size_t cl_i = 0;
    for (auto& cl : state.msh)
    {
        auto [GR, A] = make_scalar_hho_laplacian(state.msh, cl, state.hdi);
        dm S = make_scalar_hho_stabilization(state.msh, cl, GR, state.hdi);
        dm L = A+S;
        dm leigs = state.assm.get_element_dofs(state.msh, cl, state.eigvecs);
        dm leigs_full = disk::static_decondensation(L, leigs);

        auto row = cl_i*cbs;
        state.eigvecs_full.block(row, 0, cbs, found_eigs) =
            leigs_full.block(0, 0, cbs, found_eigs);

        cl_i++;
    }

    return 0;
}

template<typename Config, typename State>
static int
export_to_visit(const Config& config, State& state)
{
    using mesh_type = typename State::mesh_type;
    using scalar_type = typename State::scalar_type;
    using dm = disk::dynamic_matrix<scalar_type>;
    using dv = disk::dynamic_vector<scalar_type>;

    auto found_eigs = state.eigvals.size();

    if (0 == found_eigs) {
        std::cout << "Warning: No eigenvalues found, not exporting to Silo DB.";
        std::cout << std::endl;
        return 1; 
    }

    disk::silo_database db;
    db.create("steklov.silo");
    db.add_mesh(state.msh, "mesh");

    auto cd = state.hdi.cell_degree();
    auto cbs = disk::scalar_basis_size(cd, mesh_type::dimension);

    dm plot_eigvecs = dm::Zero(state.msh.cells_size(), found_eigs);

    size_t cl_i = 0;
    for (auto& cl : state.msh) {
        plot_eigvecs.row(cl_i) = state.eigvecs_full.row(cl_i*cbs);
        cl_i++;
    }

    for (size_t i = 0; i < found_eigs; i++) {
        std::string vname = "eig_" + std::to_string(i);
        dv eigvec = plot_eigvecs.col(i);
        db.add_variable("mesh", vname, eigvec, disk::zonal_variable_t);
    }

    return 0;
}

template<typename Mesh>
static void
detect_boundaries(sol::state& lua, Mesh& msh, std::vector<boundary_type>& bndtypes)
{
    bndtypes.reserve( msh.faces_size() );
    for (auto& fc : faces(msh))
    {
        auto bi = msh.boundary_info(fc);
        if (not bi.is_boundary()) {
            bndtypes.push_back(boundary_type::internal);
            continue;
        }

        if (bi.is_internal()) {
            bndtypes.push_back(boundary_type::internal);
            continue;
        }

        auto type = lua_get_boundary_type(lua, bi.tag());
        bndtypes.push_back(type);
    }

    assert(bndtypes.size() == msh.faces_size());
}

template<typename Config, typename State>
static int
check_error(sol::state& lua, const Config& config, State& state)
{
    using mesh_type = typename State::mesh_type;
    using scalar_type = typename State::scalar_type;
    using dm = disk::dynamic_matrix<scalar_type>;
    using dv = disk::dynamic_vector<scalar_type>;

    auto found_eigs = state.eigvals.size();

    dv ones = dv::Ones(found_eigs);
    dv num_eigs = state.eigvals.segment(0,found_eigs)-ones;
    dv ana_eigs = dv::Zero(found_eigs);
    
    if constexpr (mesh_type::dimension == 2) {
        for (size_t i = 0; i < found_eigs; i++)
            ana_eigs(i) = (i+1) * M_PI * std::tanh((i+1)*M_PI);
    }
    else {
        const size_t tmp_eigs = 5;
        const size_t tmp_len = ((tmp_eigs+2)*(tmp_eigs+1))/2;
        dv ana_eigs_tmp = dv::Zero(tmp_len);
        size_t pos = 0;
        for (size_t im = 0; im <= tmp_eigs; im++) {
            for (size_t n = 0; n <= im; n++) {
                auto m = im - n;
                auto l = std::sqrt(m*m+n*n);
                ana_eigs_tmp(pos++) = l * M_PI * std::tanh(l*M_PI);
            }
        }
        std::sort(ana_eigs_tmp.begin(), ana_eigs_tmp.end());
        ana_eigs = ana_eigs_tmp.head(found_eigs);
    }
    
    std::ios coutfmt(NULL);
    coutfmt.copyfmt(std::cout);
    std::cout << "Mesh h = " << disk::average_diameter(state.msh) << std::endl;
    std::cout << "Num: " << std::setprecision(10) << num_eigs.transpose() << std::endl;
    std::cout << "Ana: " <<  ana_eigs.transpose() << std::endl;
    dv errs = (num_eigs - ana_eigs).transpose().cwiseAbs();
    std::cout << std::setprecision(4) << std::setw(7) << std::scientific;
    for (size_t i = 0; i < errs.size(); i++)
        std::cout << errs(i) << "\t";
    std::cout << std::endl;
    std::cout.copyfmt(coutfmt);

    return 0;
}

template<typename Config, typename State>
static int
steklov_solver(sol::state& lua, const Config& config, State& state)
{
    auto cd = config.hho.order;
    auto fd = config.hho.order;
    auto rd = config.hho.order+1;
    state.hdi.cell_degree(cd);
    state.hdi.face_degree(fd);
    state.hdi.reconstruction_degree(rd);

    state.is_dirichlet.resize( state.msh.faces_size(), false );
    detect_boundaries(lua, state.msh, state.bndtype);
    std::transform(state.bndtype.begin(), state.bndtype.end(), state.is_dirichlet.begin(),
        [](const boundary_type& b) {
            return boundary_type::dirichlet == b;
        }
    );

    assemble(config, state);
    solve(config, state);
    export_to_visit(config, state);
    check_error(lua, config, state);

    return 0;
}

template<typename T>
static int
run_mesh_internal(sol::state& lua, const steklov_solver_configuration<T>& config)
{
    using namespace disk::lua;

    switch (config.mesh.type) {

        case internal_mesh_type::triangles: {
            using mesh_type = disk::simplicial_mesh<T,2>;
            hho_steklov_solver_state<mesh_type> state;
            auto mesher = disk::make_simple_mesher(state.msh);
            for (size_t i = 0; i < config.mesh.reflevel; i++)
                mesher.refine();
            return steklov_solver(lua, config, state);
        } break;

        case internal_mesh_type::quadrangles: {
            using mesh_type = disk::cartesian_mesh<T,2>;
            hho_steklov_solver_state<mesh_type> state;
            auto mesher = disk::make_simple_mesher(state.msh);
            for (size_t i = 0; i < config.mesh.reflevel; i++)
                mesher.refine();
            return steklov_solver(lua, config, state);
        } break;

        case internal_mesh_type::hexagons: {
            using mesh_type = disk::generic_mesh<T,2>;
            hho_steklov_solver_state<mesh_type> state;
            auto mesher = disk::make_fvca5_hex_mesher(state.msh);
            mesher.make_level(config.mesh.reflevel);
            return steklov_solver(lua, config, state);
        } break;

        case internal_mesh_type::tetrahedra: {
            using mesh_type = disk::simplicial_mesh<T,3>;
            hho_steklov_solver_state<mesh_type> state;
            auto mesher = disk::make_simple_mesher(state.msh);
            for (size_t i = 0; i < config.mesh.reflevel; i++)
                mesher.refine();
            return steklov_solver(lua, config, state);
        } break;

        default:
            std::cout << "config.mesh.type: invalid value" << std::endl;
            return -1;
    }

    return 0;
}

template<typename T>
static int
run_mesh_from_file(sol::state& lua, const steklov_solver_configuration<T>& config)
{
    auto& fname = config.mesh.filename;

    if (std::regex_match(fname, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        using mesh_type = disk::generic_mesh<T, 2>;
        hho_steklov_solver_state<mesh_type> state;
        disk::load_mesh_fvca5_2d<T>(fname.c_str(), state.msh);
        return steklov_solver(lua, config, state);
    }
    
#ifdef HAVE_GMSH
    /* GMSH 2D simplicials */
    if (std::regex_match(fname, std::regex(".*\\.geo2s$") ))
    {
        using mesh_type = disk::simplicial_mesh<T,2>;
        hho_steklov_solver_state<mesh_type> state;
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::gmsh_geometry_loader<mesh_type> loader;
        loader.read_mesh(fname);
        loader.populate_mesh(state.msh);
        return steklov_solver(lua, config, state);
    }

    /* GMSH 3D simplicials */
    if (std::regex_match(fname, std::regex(".*\\.geo3s$") ))
    {
        using mesh_type = disk::simplicial_mesh<T,3>;
        hho_steklov_solver_state<mesh_type> state;
        std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
        disk::gmsh_geometry_loader<mesh_type> loader;
        loader.read_mesh(fname);
        loader.populate_mesh(state.msh);
        return steklov_solver(lua, config, state);
    }
#endif

    return -1;
}

template<typename T>
int
run(sol::state& lua, const steklov_solver_configuration<T>& config)
{
    using namespace disk::lua;

    int ret = 0;
    switch (config.mesh.source)
    {
        case mesh_source::internal:
            ret = run_mesh_internal(lua, config);
            break;

        case mesh_source::file:
            ret = run_mesh_from_file(lua, config);
            break;

        case mesh_source::invalid:
            std::cout << "Mesh was not specified." << std::endl;
            return -1;
    }

    return ret;
}

int main(int argc, char **argv)
{
    using T = double;

    sol::state lua;
    steklov_solver_configuration<T> config;

    if ( not lua_init_context(lua, config) )
    {
        std::cout << "Problem in Lua config" << std::endl;
        return 1;
    }

    lua["run"] = [&](){ run(lua, config); };

    lua["solution_process"]();

    return 0;
}