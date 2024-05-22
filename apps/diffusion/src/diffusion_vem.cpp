#include <iostream>
#include <iomanip>
#include <optional>

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

std::vector<std::optional<size_t>>
compress(const std::vector<size_t>& l2g,
    const std::vector<std::optional<size_t>>& compress_map)
{
    std::vector<std::optional<size_t>> ret;
    for (size_t i = 0; i < l2g.size(); i++) {
        assert(l2g[i] < compress_map.size());
        ret.push_back(compress_map[l2g[i]]);
    }

    return ret;
}

template<disk::mesh_2D Mesh>
void
test_vem_diffusion(const Mesh& msh, size_t degree)
{
    using T = typename Mesh::coordinate_type;

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

    std::vector<std::optional<size_t>> compress_map(totdofs);

    for (size_t i = 0, cnum = 0; i < totdofs; i++) {
        if ( is_dirichlet[i] )
            continue;

        compress_map[i] = cnum++;
    }
    size_t realdofs = totdofs - dirdofs;

    std::cout << "Linear system DoFs: " << realdofs << std::endl;

    disk::sparse_matrix<T> LHS(realdofs, realdofs);
    disk::dynamic_vector<T> RHS = disk::dynamic_vector<T>::Zero(realdofs);
    disk::dynamic_vector<T> csol = disk::dynamic_vector<T>::Zero(realdofs);

    using triplet_t = Eigen::Triplet<T>;
    std::vector<triplet_t> trips;

    for (size_t cl_i = 0; cl_i < msh.cells_size(); cl_i++)
    {
        std::cout << "\r Assembling cell " << cl_i+1 << "/";
        std::cout << msh.cells_size();
        std::cout.flush();
        auto cl = msh[cl_i];
        auto pts = points(msh, cl);
        auto npts = pts.size();
        auto nbdofs = npts*degree;
        
        auto [L, R] = disk::vem_2d::compute_local(msh, cl, degree, f);
        auto [Lc, Rc] = disk::vem::schur(L, R, nbdofs);

        auto l2g = compress(dofmap.cell_to_global(cl_i),
            compress_map);

        assert(Lc.rows() == Lc.cols());
        assert(Lc.cols() == l2g.size());

        for (size_t i = 0; i < l2g.size(); i++) {
            if (not l2g[i])
                continue;
            
            auto gi = l2g[i].value();
            for (size_t j = 0; j < l2g.size(); j++) {
                if (not l2g[j])
                    continue;
                auto gj = l2g[j].value();
                trips.push_back( triplet_t(gi, gj, Lc(i,j)) );
            }

            RHS[gi] += Rc[i];
        }
    }
    std::cout << std::endl;

    LHS.setFromTriplets( trips.begin(), trips.end() );
    trips.clear();

    std::cout << "Linear solver..." << std::flush;
    csol = mumps_lu(LHS, RHS);
    std::cout << "done" << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(msh.points_size());

    for (size_t i = 0; i < msh.points_size(); i++) {
        if (not compress_map[i])
            continue;
        sol[i] = csol[ compress_map[i].value() ];
    }

    disk::silo_database db;
    db.create("vem.silo");
    db.add_mesh(msh, "mesh");
    db.add_variable("mesh", "u", sol, disk::nodal_variable_t);
}

enum class mesh_type {
    triangles,
    cartesian,
    hexas,
};

int main(int argc, char **argv)
{
    using T = double;
    size_t degree = 2;
    size_t level = 2;
    mesh_type mt = mesh_type::hexas;

    int opt;
    while ((opt = getopt(argc, argv, "k:m:l:")) != -1) {
        switch (opt) {
        case 'k': {
            int k = atoi(optarg);
            if (k < 1) {
                std::cerr << "Degree must be greater than zero." << std::endl;
                return 1;
            }
            degree = k;
            } break;
        case 'm': {
            if ( std::string(optarg) == "tri" )
                mt = mesh_type::triangles;
            else if ( std::string(optarg) == "quad" )
                mt = mesh_type::cartesian;
            else if ( std::string(optarg) == "hex" )
                mt = mesh_type::hexas;
            else {
                std::cerr << "Mesh type: tri, quad or hex" << std::endl;
                return 1;
            }
            } break;
        case 'l': {
            int l = atoi(optarg);
            if (l < 0) {
                std::cerr << "Level must be positive." << std::endl;
                return 1;
            }
            level = l;
            } break;
        }
    }

    if (mt == mesh_type::triangles) {
        disk::simplicial_mesh<T,2> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < level; l++)
            mesher.refine();

        test_vem_diffusion(msh, degree);
    }

    if (mt == mesh_type::cartesian) {
        disk::cartesian_mesh<T,2> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < level; l++)
            mesher.refine();

        test_vem_diffusion(msh, degree);
    }

    if (mt == mesh_type::hexas) {
        disk::generic_mesh<T,2> msh;
        auto mesher = make_fvca5_hex_mesher(msh);
        mesher.make_level(level);

        test_vem_diffusion(msh, degree);
    }
    

    return 0;
}