#include <iostream>
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/vem"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"

#include "mumps.hpp"

template<disk::mesh_2D Mesh>
struct rhs_functor {
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type; 

    rhs_functor(const Mesh& msh)
    {}

    T operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI * pt.x());
        auto sy = std::sin(M_PI * pt.y());
        return 2.0 * M_PI * M_PI * sx * sy;
    }

};

int main(void)
{
    using T = double;
    size_t degree = 2;

    //disk::generic_mesh<T,2> msh;
    //auto mesher = make_fvca5_hex_mesher(msh);
    //mesher.make_level(4);

    disk::cartesian_mesh<T,2> msh;
    auto mesher = make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();
    mesher.refine();

    rhs_functor f(msh);

    auto nvtx = msh.points_size();
    auto nfcs = msh.faces_size();
    auto totdofs = nvtx + (degree - 1)*nfcs;

    std::vector<bool> is_dirichlet(totdofs, false);

    disk::dof_mapper dofmap(msh, degree);

    size_t nbndfcs = 0;
    for (auto& fc : faces(msh)) {
        if ( msh.is_boundary(fc) ) {
            nbndfcs++;
            auto b2g = dofmap.bnd_to_global( offset(msh, fc) );
            for (auto& d : b2g)
                is_dirichlet.at(d) = true;
        }
    }

    auto dirdofs = std::count(is_dirichlet.begin(),
        is_dirichlet.end(), true);

    auto dirdofs_chk = degree*nbndfcs;
    assert(dirdofs == dirdofs_chk);

    std::cout << "nvtx = " << nvtx << ", nfcs = " << nfcs << ", ";
    std::cout << "totdofs = " << totdofs << ", nbndfcs = ";
    std::cout << nbndfcs << ", dirdofs = " << dirdofs << std::endl;

    std::vector<size_t> compress_map(totdofs);

    for (size_t i = 0, cnum = 0; i < totdofs; i++) {
        if ( is_dirichlet[i] )
            continue;

        compress_map[i] = cnum++;
    }
    size_t realdofs = totdofs - dirdofs;

    disk::sparse_matrix<T> LHS(realdofs, realdofs);
    disk::dynamic_vector<T> RHS = disk::dynamic_vector<T>::Zero(realdofs);
    disk::dynamic_vector<T> csol = disk::dynamic_vector<T>::Zero(realdofs);

    using triplet_t = Eigen::Triplet<T>;
    std::vector<triplet_t> trips;

    for (size_t cl_i = 0; cl_i < msh.cells_size(); cl_i++)
    {
        auto cl = msh[cl_i];
        auto pts = points(msh, cl);
        auto npts = pts.size();
        auto nbdofs = npts*degree;
        
        auto [L, R] = disk::vem_2d::compute_local(msh, cl, degree, f);
        auto [Lc, Rc] = disk::vem::schur(L, R, nbdofs);

        auto l2g = dofmap.cell_to_global(cl_i);
        assert(Lc.rows() == Lc.cols());
        assert(Lc.cols() == l2g.size());

        for (size_t i = 0; i < l2g.size(); i++) {
            auto gi = l2g[i];
            if (is_dirichlet[gi])
                continue;
            gi = compress_map[gi];
            for (size_t j = 0; j < l2g.size(); j++) {
                auto gj = l2g[j];
                if (is_dirichlet[gj])
                    continue;
                gj = compress_map[gj];
                trips.push_back( triplet_t(gi, gj, Lc(i,j)) );
            }

            RHS[gi] += Rc[i];
        }
    }

    LHS.setFromTriplets( trips.begin(), trips.end() );

    csol = mumps_lu(LHS, RHS);

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(msh.points_size());

    for (size_t i = 0; i < msh.points_size(); i++) {
        if (is_dirichlet[i])
            continue;
        sol[i] = csol[ compress_map[i] ];
    }

    disk::silo_database db;
    db.create("vem.silo");
    db.add_mesh(msh, "mesh");
    db.add_variable("mesh", "u", sol, disk::nodal_variable_t);

    return 0;
}