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

template<typename T>
void
lua_get_feast_params(sol::state& lua, disk::feast_eigensolver_params<T>& fep)
{
    sol::optional<int> tolerance_opt = lua["feast"]["tolerance"];
    if (tolerance_opt)
        fep.tolerance = tolerance_opt.value();

    sol::optional<int> subspace_size_opt = lua["feast"]["subspace_size"];
    if (subspace_size_opt)
        fep.subspace_size = subspace_size_opt.value();

    sol::optional<T> eigval_min_opt = lua["feast"]["eigval_min"];
    if (eigval_min_opt)
        fep.min_eigval = eigval_min_opt.value();

    sol::optional<T> eigval_max_opt = lua["feast"]["eigval_max"];
    if (eigval_max_opt)
        fep.max_eigval = eigval_max_opt.value();

    sol::optional<bool> verbose_opt = lua["feast"]["verbose"];
    if (verbose_opt)
        fep.verbose = verbose_opt.value();

    sol::optional<std::string> inner_solver_opt = lua["feast"]["inner_solver"];
    if (inner_solver_opt) {
        std::string inner_solver = inner_solver_opt.value();
        if (inner_solver == "eigen_sparselu")
            fep.fis = disk::feast_inner_solver::eigen_sparselu;
        if (inner_solver == "mumps")
            fep.fis = disk::feast_inner_solver::mumps;
    }
}

void
lua_get_hho_params(sol::state& lua, size_t& order)
{
    sol::optional<size_t> order_opt = lua["hho"]["order"];
    if (order_opt)
        order = order_opt.value();
}

template<typename Mesh>
void
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

template<typename Mesh>
void
steklov_solver(sol::state& lua, Mesh& msh)
{
    using T = typename Mesh::coordinate_type;

    hho_assembler_steklov<Mesh> assm;

    size_t degree = 0;
    lua_get_hho_params(lua, degree);
    auto cd = degree;
    auto fd = degree;
    auto rd = degree+1;
    disk::hho_degree_info hdi;
    hdi.cell_degree(cd);
    hdi.face_degree(fd);
    hdi.reconstruction_degree(rd);

    auto cbs = disk::scalar_basis_size(cd, Mesh::dimension);
    auto fbs = disk::scalar_basis_size(fd, Mesh::dimension-1);
    auto rbs = disk::scalar_basis_size(rd, Mesh::dimension);

    std::vector<bool> is_dirichlet( msh.faces_size(), false );
    std::vector<boundary_type> bndtypes;
    detect_boundaries(lua, msh, bndtypes);
    assert(bndtypes.size() == is_dirichlet.size());
    for (size_t i = 0; i < bndtypes.size(); i++)
        if (bndtypes[i] == boundary_type::dirichlet)
            is_dirichlet[i] = true;

    assm.initialize(msh, hdi, is_dirichlet);

    for (auto& cl : msh)
    {
        auto [GR, A] = make_scalar_hho_laplacian(msh, cl, hdi);
        disk::dynamic_matrix<T> S = make_scalar_hho_stabilization(msh, cl, GR, hdi);
        disk::dynamic_matrix<T> L = A+S;

        auto fcs = faces(msh, cl);

        disk::dynamic_matrix<T> BFF =
            disk::dynamic_matrix<T>::Zero(fcs.size()*fbs, fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto& fc = fcs[face_i];
            auto gfnum = offset(msh, fc);
            auto bndtype = bndtypes.at(gfnum);

            switch (bndtype) {
                case boundary_type::dirichlet:
                case boundary_type::neumann:
                    continue;
                    break;

                case boundary_type::robin: {
                    auto fb = disk::make_scalar_monomial_basis(msh, fc, fd);
                    auto idx = cbs+fbs*face_i;
                    L.block(idx, idx, fbs, fbs) +=
                        disk::make_mass_matrix(msh, fc, fb);
                } break;

                case boundary_type::steklov: {
                    auto fb = disk::make_scalar_monomial_basis(msh, fc, fd);
                    auto idx = fbs*face_i;
                    L.block(cbs+idx, cbs+idx, fbs, fbs) +=
                        disk::make_mass_matrix(msh, fc, fb);
                    BFF.block(idx, idx, fbs, fbs) =
                        disk::make_mass_matrix(msh, fc, fb);
                } break;

                case boundary_type::internal :
                    break;

                default:
                    break;
            }
        }

        disk::dynamic_matrix<T> LC = disk::static_condensation(L, cbs);
        assm.assemble(msh, cl, LC, BFF);
    }

    assm.finalize();


    disk::dynamic_matrix<T> eigvecs;
    disk::dynamic_vector<T> eigvals;


    disk::feast_eigensolver_params<T> fep;
    fep.verbose = true;
    fep.tolerance = 8;
    fep.min_eigval = 1;
    fep.max_eigval = 6;
    fep.subspace_size = 30;
    fep.max_iter = 10;
    fep.fis = disk::feast_inner_solver::eigen_sparselu;
    lua_get_feast_params(lua, fep);

    int retries = lua["feast"]["retries"].get_or(20);
    bool use_mkl_feast = lua["feast"]["use_mkl"].get_or(false);

    if (use_mkl_feast) {
        bool s = false;
        size_t retry = 0;
        do {
            int ret = disk::generalized_eigenvalue_solver(fep, assm.LHS, assm.RHS, eigvecs, eigvals);
            s = ret != -2;
        } while (!s && retry++ < retries);
    }
    else {
        auto ret = disk::feast(fep, assm.LHS, assm.RHS, eigvecs, eigvals);
        std::cout << ret << std::endl;
    }

    std::cout << "Found eigs: " << fep.eigvals_found << std::endl;
    auto found_eigs = fep.eigvals_found;

    disk::dynamic_vector<T> ones = disk::dynamic_vector<T>::Ones(found_eigs);
    disk::dynamic_vector<T> num_eigs = eigvals.segment(0,found_eigs)-ones;
    disk::dynamic_vector<T> ana_eigs = disk::dynamic_vector<T>::Zero(found_eigs);
    //for (size_t i = 0; i < found_eigs; i++)
    //    ana_eigs(i) = i * M_PI * std::tanh(i*M_PI);

    const size_t tmp_eigs = 5;
    const size_t tmp_len = ((tmp_eigs+2)*(tmp_eigs+1))/2;
    disk::dynamic_vector<T> ana_eigs_tmp = disk::dynamic_vector<T>::Zero(tmp_len);
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
    
    std::cout << "Mesh h = " << disk::average_diameter(msh) << std::endl;
    std::cout << "Num: " << num_eigs.transpose() << std::endl;
    std::cout << "Ana: " <<  ana_eigs.transpose() << std::endl;
    std::cout << "Err: " <<  (num_eigs - ana_eigs).transpose().cwiseAbs() << std::endl;

    if (found_eigs == 0)
        return;

    disk::silo_database db;
    db.create("steklov.silo");
    db.add_mesh(msh, "mesh");

    disk::dynamic_matrix<T> es = disk::dynamic_matrix<T>::Zero(msh.cells_size(), found_eigs);

    size_t cl_i = 0;
    for (auto& cl : msh)
    {
        auto [GR, A] = make_scalar_hho_laplacian(msh, cl, hdi);
        disk::dynamic_matrix<T> S = make_scalar_hho_stabilization(msh, cl, GR, hdi);
        disk::dynamic_matrix<T> L = A+S;
        disk::dynamic_matrix<T> leigs = assm.get_element_dofs(msh, cl, eigvecs);
        disk::dynamic_matrix<T> leigs_full = disk::static_decondensation(L, leigs);

        es.row(cl_i) = leigs_full.row(0);
        cl_i++;
    }

    for (size_t i = 0; i < found_eigs; i++)
    {
        std::string vname = "eig_" + std::to_string(i);
        disk::dynamic_vector<T> e = es.col(i);
        db.add_variable("mesh", vname, e, disk::zonal_variable_t);
    }

}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <mesh filename>" << std::endl;
        return 1;
    }

    sol::state lua;

    lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::io, sol::lib::table, sol::lib::string);
    lua["boundary"] = lua.create_table();
    lua["feast"] = lua.create_table();
    lua["hho"] = lua.create_table();
    auto sf = lua.safe_script_file("steklov.lua");
    if (not sf.valid())
    {
        std::cout << "Problem in Lua config" << std::endl;
        return 1;
    }

    using T = double;

    std::string mesh_filename = argv[1];
        
#ifdef HAVE_GMSH
    /* GMSH 2D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        using mesh_type = disk::simplicial_mesh<T,2>;
        mesh_type msh;
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::gmsh_geometry_loader<mesh_type> loader;
    
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        steklov_solver(lua, msh);
    }

    /* GMSH 3D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        using mesh_type = disk::simplicial_mesh<T,3>;
        mesh_type msh;
        std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
        disk::gmsh_geometry_loader<mesh_type> loader;
    
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        steklov_solver(lua, msh);
    }
#else
    std::cout << "GMSH support not compiled. Exiting." << std::endl;
    return 1;
#endif


    return 0;
}