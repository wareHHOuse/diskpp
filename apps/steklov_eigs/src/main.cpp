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

        if (sysfcs == msh.faces_size()) {
            std::cout << "Pure Neumann problem" << std::endl;
        }

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




template<typename Mesh>
void
steklov_solver(Mesh& msh)
{
    using T = typename Mesh::coordinate_type;

    hho_assembler_steklov<Mesh> assm;

    const size_t degree = 0;
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

    std::vector<bool> is_dirichlet;
    is_dirichlet.resize( msh.faces_size() );
    assm.initialize(msh, hdi, is_dirichlet);

    for (auto& cl : msh)
    {
        auto [GR, A] = make_scalar_hho_laplacian(msh, cl, hdi);
        disk::dynamic_matrix<T> S = make_scalar_hho_stabilization(msh, cl, GR, hdi);
        disk::dynamic_matrix<T> L = A+S;
        disk::dynamic_matrix<T> LC = disk::static_condensation(L, cbs);

        auto fcs = faces(msh, cl);

        disk::dynamic_matrix<T> BFF =
            disk::dynamic_matrix<T>::Zero(fcs.size()*fbs, fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto& fc = fcs[face_i];
            auto bi = msh.boundary_info(fc);
            if (not bi.is_boundary())
                continue;

            auto boundary_id = bi.tag();

            auto fb = disk::make_scalar_monomial_basis(msh, fc, fd);
            auto idx = fbs*face_i;
            BFF.block(idx, idx, fbs, fbs) =
                disk::make_mass_matrix(msh, fc, fb);
        }

        assm.assemble(msh, cl, LC, BFF);
    }

    assm.finalize();


    disk::dynamic_matrix<T> eigvecs;
    disk::dynamic_vector<T> eigvals;
    disk::feast_eigensolver_params<T> fep;
    fep.verbose = true;
    fep.tolerance = 8;
    fep.min_eigval = 10;
    fep.max_eigval = 60;
    fep.subspace_size = 30;
    disk::generalized_eigenvalue_solver(fep, assm.LHS, assm.RHS, eigvecs, eigvals);
    std::cout << "Found eigs: " << fep.eigvals_found << std::endl;
    std::cout << eigvals.transpose() << std::endl;


    disk::silo_database db;
    db.create("steklov.silo");
    db.add_mesh(msh, "mesh");

    disk::dynamic_matrix<T> es = disk::dynamic_matrix<T>::Zero(msh.cells_size(), fep.eigvals_found);

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

    for (size_t i = 0; i < fep.eigvals_found; i++)
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

        steklov_solver(msh);
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

        steklov_solver(msh);
    }
#else
    std::cout << "GMSH support not compiled. Exiting." << std::endl;
    return 1;
#endif


    return 0;
}