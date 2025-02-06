/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <set>

/* ENABLE_CCW_FACES is a horrible hack to circumvent some limitations
 * of the original polygonal element implementation of DiSk++.
 * See diskpp/geometry/geometry_generic.hpp for the code enabled
 * by this define.
 *
 * in generic_mesh<> the cells deduce their faces by pointing
 * to the corresponding objects in the mesh storage. In turn, the
 * faces are stored as a pair (a,b) such that a < b where a and b
 * are the indices of the points of the segment. This does not allow
 * to compute normals only by looking at the segment itself, therefore
 * the implementation of normal() found in geometry_all.hpp.
 *
 * The code that this define enables essentially gives you a faces()
 * that returns the faces 1) sorted in CCW order and 2) with the proper
 * orientation, so normals can be computed trivially. This however breaks
 * the VEM solver and possibly other things.
 *
 * The only way to really fix this is to reimplement generic_mesh<T,2>
 * properly in something like polygonal_mesh<T> and gradually migrate
 * everything to polygonal_mesh<T>.
 *
 * The reason why we need ENABLE_CCW_FACES in this program is that the
 * agglomeration process generates polygons on which the original
 * implementation of the normal computation fails because the direction
 * of the normal cannot be deduced from the barycenter of the polygon
 * and the barycenter of the face.
 */

#define ENABLE_CCW_FACES

#include "diskpp/common/eigen.hpp"
#include "diskpp/mesh/mesh.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/common/timecounter.hpp"

#include "diskpp/methods/hho"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/dg"
#include "mumps.hpp"

#include "sgr.hpp"

template<typename Mesh>
struct source_functor;

template<disk::mesh_2D Mesh>
struct source_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return 2.0*M_PI*M_PI*sx*sy;
    }
};

template<disk::mesh_3D Mesh>
struct source_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return 3.0*M_PI*M_PI*sx*sy*sz;
    }
};


template<typename Mesh>
struct solution_functor;

template<disk::mesh_2D Mesh>
struct solution_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return sx*sy;
    }
};

template<disk::mesh_3D Mesh>
struct solution_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return sx*sy*sz;
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return source_functor<Mesh>();
}

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

template<typename Mesh>
void
hho_diffusion_solver(const Mesh& msh, size_t degree, disk::silo_database& silo)
{
    std::cout << "HHO solver" << std::endl;
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;

    degree_info di(degree);

    auto f = make_rhs_function(msh);

    auto assm = make_assembler(msh, di);

    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    std::vector<std::pair<MT, VT>> lcs;

    timecounter tc;
    tc.tic();
    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);
        disk::dynamic_matrix<T> lhs = A+S;
        auto phiT = hho_space<mesh_type>::cell_basis(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
        lcs.push_back({lhs, rhs});
        auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);
        assm.assemble(msh, cl, lhsc, rhsc);
    }
    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;
    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;
    tc.tic();
    disk::dynamic_vector<T> sol = mumps_lu(assm.LHS, assm.RHS);
    std::cout << " Solver time: " << tc.toc() << std::endl;
    std::vector<T> u_data;
    
    T error = 0.0;
    T L2error = 0.0;
    auto u_sol = make_solution_function(msh);
    tc.tic();
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto phiT = hho_space<mesh_type>::cell_basis(msh, cl, di.cell);
        auto MMe = integrate(msh, cl, phiT, phiT);
        const auto& [lhs, rhs] = lcs[cell_i++];
        disk::dynamic_vector<T> sol_ana = local_reduction(msh, cl, di, u_sol);
        auto locsolF = assm.take_local_solution(msh, cl, sol);
        disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
        u_data.push_back(locsol(0));
        disk::dynamic_vector<T> diff = locsol - sol_ana;
        error += diff.dot(lhs*diff);
        disk::dynamic_vector<T> diffT = diff.head(phiT.size());
        L2error += diffT.transpose() * (MMe*diffT);
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error: " << std::sqrt(L2error) << ", ";
    std::cout << "A-norm error: " << std::sqrt(error) << std::endl;

    silo.add_variable("dstmesh", "u_hho", u_data, disk::zonal_variable_t);
}

template<typename Mesh>
auto
dg_diffusion_solver(Mesh& msh, size_t degree,
    const typename Mesh::coordinate_type eta,
    disk::silo_database& silo)
{   
    std::cout << "DG solver" << std::endl;
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto basis_rescaling = disk::basis::rescaling_strategy::inertial;

    auto f = make_rhs_function(msh);

    auto cbs = disk::scalar_basis_size(degree, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_assembler(msh, cbs);
    
    timecounter tc;
    tc.tic();
    for (auto& tcl : msh)
    {
        auto tbasis = disk::basis::scaled_monomial_basis(msh, tcl, degree, basis_rescaling);
        
        matrix_type K = integrate(msh, tcl, grad(tbasis), grad(tbasis));
        vector_type loc_rhs = integrate(msh, tcl, f, tbasis);

        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {   
            auto n     = normal(msh, tcl, fc);
            auto eta_l = eta / diameter(msh, fc);
            
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            if (nv) {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
                
                auto ncl = nv.value();
                auto nbasis = disk::basis::scaled_monomial_basis(msh, ncl, degree, basis_rescaling);
                assert(tbasis.size() == nbasis.size());

                Att += + eta_l * integrate(msh, fc, tbasis, tbasis);
                Att += - 0.5 * integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - 0.5 * integrate(msh, fc, tbasis, grad(tbasis).dot(n));

                Atn += - eta_l * integrate(msh, fc, nbasis, tbasis);
                Atn += - 0.5 * integrate(msh, fc, grad(nbasis).dot(n), tbasis);
                Atn += + 0.5 * integrate(msh, fc, nbasis, grad(tbasis).dot(n));

                assm.assemble(msh, tcl, tcl, Att);
                assm.assemble(msh, tcl, ncl, Atn);
            }
            else {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                Att += + eta_l * integrate(msh, fc, tbasis, tbasis);
                Att += - integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - integrate(msh, fc, tbasis, grad(tbasis).dot(n));
                assm.assemble(msh, tcl, tcl, Att);
            }   
        }
        
        assm.assemble(msh, tcl, K, loc_rhs);
    }

    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;

    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    tc.tic();
    sol = mumps_lu(assm.LHS, assm.RHS);
    std::cout << " Solver time: " << tc.toc() << std::endl;

    auto sol_fun = make_solution_function(msh);

    std::vector<double> u;
    u.reserve(msh.cells_size());

    T err = 0.0; size_t cell_i = 0;
    tc.tic();
    for (auto& cl : msh)
    {
        auto cb = disk::basis::scaled_monomial_basis(msh, cl, degree, basis_rescaling);
        auto MMe = integrate(msh, cl, cb, cb);
        auto arhs = integrate(msh, cl, sol_fun, cb);

        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());
        Matrix<T, Dynamic, 1> diff = lsol - asol;
        err += diff.dot(MMe*diff);
        u.push_back( lsol(0) );
        cell_i++;
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error: " << std::sqrt(err) << std::endl;;
    silo.add_variable("dstmesh", "u_dg", u, disk::zonal_variable_t);

    return sol;
}


template<disk::mesh_2D FineMesh>
using coarse_mesh_t = disk::generic_mesh<typename FineMesh::coordinate_type, 2>;

/* This function agglomerates elements in a 2D mesh to polygons
 * of a new 2D mesh. Agglomeration is based on the subdomain
 * numbering: elements with the same subdomain id get agglomerated
 * into the same polygon. This works only for connected patches
 * of elements. If elements with the same subdomain id are in two
 * different, disconnected patches, the result is not defined.
 */
template<typename FineMesh>
void agglomerate_by_subdomain(const FineMesh& srcmsh,
    coarse_mesh_t<FineMesh>& dstmsh)
{
    struct adjlist {
        size_t nodes[2];
        size_t used;
        adjlist() : used(0) {}
        void insert(size_t n) {
            if (used > 1)
                throw std::logic_error("Adjacency list full.");
            nodes[used++] = n;
        }
    };

    using src_mesh_type = FineMesh;
    using coord_type = typename FineMesh::coordinate_type;
    using dst_mesh_type = coarse_mesh_t<FineMesh>;
    using src_face_type = typename src_mesh_type::face_type;
    using dst_node_type = typename dst_mesh_type::node_type;
    using dst_edge_type = typename dst_mesh_type::edge_type;
    using dst_surf_type = typename dst_mesh_type::surface_type;

    using polygraph_t = std::map<size_t, adjlist>;
    std::map<size_t, polygraph_t> pgs;
    std::vector<std::optional<size_t>> compress_map;
    compress_map.resize(srcmsh.points_size());

    /* Collect all the boundary edges of the subdomains we want to
     * agglomerate in a polygon. Make a graph out of them. */
    for (auto& cl : srcmsh) {
        auto di = srcmsh.domain_info(cl);
        auto& pg = pgs[di.tag()];
        auto fcs = faces(srcmsh, cl);
        for (auto& fc : fcs) {
            auto bi = srcmsh.boundary_info(fc);
            if (not bi.is_boundary())
                continue;
            auto [pi1, pi2] = fc.point_ids();
            pg[pi1].insert(pi2);
            pg[pi2].insert(pi1);
            compress_map[pi1] = 1;
            compress_map[pi2] = 1;
        }
    }
    
    /* Make the original mesh to new mesh node mapping. Collect the
     * necessary points and node numbers. */
    auto srcstor = srcmsh.backend_storage();
    auto dststor = dstmsh.backend_storage();
    for (size_t i = 0, ci = 0; i < compress_map.size(); i++) {
        if (compress_map[i]) {
            dststor->points.push_back(srcstor->points[i]);
            dststor->nodes.push_back(typename dst_mesh_type::node_type{ci});
            compress_map[i] = ci++;
        }
    }

    /* Do a DFS to build the polygon out of its edges. This is a special
     * case of DFS as each node is guaranteed to have only two neighbours.*/
    std::map<size_t, std::vector<size_t>> paths;
    for (auto& [tag, pg] : pgs) {
        auto nvtx = pg.size();
        assert(nvtx >= 3);
        std::vector<size_t> path;
        path.reserve(nvtx);
        auto [start, adj] = *pg.begin();
        assert(adj.used == 2);
        size_t visiting = adj.nodes[0];
        path.push_back(start);
        for (size_t i = 1; i < nvtx; i++) {
            path.push_back(visiting);
            adj = pg.at(visiting);
            visiting = (adj.nodes[0] == path[i-1]) ? adj.nodes[1] : adj.nodes[0];
        }
        assert(visiting == start);
        
        /* Reverse the path if vertices are in clockwise order 
         * (uses shoelace formula).*/
        double dir = 0.0;
        for (size_t i = 0; i < path.size(); i++) {
            auto p0 = srcstor->points[ path[i] ];
            auto p1 = srcstor->points[ path[(i+1)%path.size()] ];
            dir += (p1.x() - p0.x())*(p1.y() + p0.y());
        }

        if (dir > 0) std::reverse(path.begin(), path.end());

        /* and translate the node numbers to the new mesh numbering */
        for (auto & n : path)
            n = compress_map[ n ].value();

        /* Now put all the edges in the mesh */
        for (size_t i = 0; i < path.size(); i++) {
            auto p0 = path[i];
            auto p1 = path[(i+1)%path.size()];
            auto node1 = typename dst_node_type::id_type(p0);
            auto node2 = typename dst_node_type::id_type(p1);
            auto e = dst_edge_type(node1, node2);
            dststor->edges.push_back(e);
        }

        paths[tag] = std::move(path);
    }

    /* We have all the edges: sort them and make them unique */
    std::sort(dststor->edges.begin(), dststor->edges.end());
    auto last = std::unique(dststor->edges.begin(), dststor->edges.end());
    dststor->edges.erase(last, dststor->edges.end());

    using newsurf_pair_t = std::pair<size_t, dst_surf_type>;
    using newsurfs_vec_t = std::vector<newsurf_pair_t>;
    newsurfs_vec_t newsurfs;
    newsurfs.reserve(paths.size());

    /* and finally all the elements */
    for (auto& [tag, path] : paths) {
        std::vector<typename dst_edge_type::id_type> surfedges;
        for (size_t i = 0; i < path.size(); i++) {
            auto p0 = path[i];
            auto p1 = path[(i+1)%path.size()];
            auto node1 = typename dst_node_type::id_type(p0);
            auto node2 = typename dst_node_type::id_type(p1);
            auto e = dst_edge_type(node1, node2);
            auto eid = find_element_id(dststor->edges.begin(),
                dststor->edges.end(), e);
            surfedges.push_back(eid.second);
        }

        /* We store in [tag, surf] pairs to not lose the association
         * during sorting. Later we unzip the pair and store the
         * data in the mesh. */
        dst_surf_type newsurf(surfedges);
        newsurf.set_point_ids(path.begin(), path.end());
        newsurfs.push_back({tag, newsurf});
    }

    /* Sort all the tag-surface pairs according to the surface ordering */
    std::sort(newsurfs.begin(), newsurfs.end(),
        [](const newsurf_pair_t& a, const newsurf_pair_t& b) {
            return a.second < b.second;
        }
    );

    /* Unzip the newsurf array in the appropriate storage arrays */
    dststor->surfaces.clear();
    dststor->subdomain_info.clear();
    dststor->surfaces.reserve( newsurfs.size() );
    dststor->subdomain_info.reserve( newsurfs.size() );
    for (const auto& [tag, newsurf] : newsurfs) {
        dststor->surfaces.push_back(newsurf);
        auto di = disk::subdomain_descriptor(tag);
        dststor->subdomain_info.push_back(di);
    }
    assert( dststor->surfaces.size() == newsurfs.size() );
    assert( dststor->subdomain_info.size() == newsurfs.size() );

    dststor->boundary_info.resize( dststor->edges.size() );

    /* Count the elements adjacent to each face */
    std::vector<size_t> counts(dststor->edges.size());
    for (auto& cl : dstmsh) {
        auto fcs = faces(dstmsh, cl);
        for (auto& fc : fcs) {
            counts[offset(dstmsh, fc)]++;
        }
    }

    /* All the faces adjacent to only one element are
     * boundary faces. This should preserve the original
     * boundaries, but we leave this for later. Not
     * necessary for now. */
    for ( size_t i = 0; i < counts.size(); i++) {
        if (counts[i] == 1) {
            dststor->boundary_info[i].is_boundary(true);
            dststor->boundary_info[i].id(0);
        }
    }

    disk::mark_internal_faces(dstmsh);
}

template<typename FM>
using cc2ff_t = std::map<typename coarse_mesh_t<FM>::cell_type, std::set<typename FM::face_type>>;

template<typename FineMesh>
auto make_cc2ff(const FineMesh& fmsh,
    const coarse_mesh_t<FineMesh>& cmsh)
{
    using fine_mesh_type = FineMesh;
    using coord_type = typename fine_mesh_type::coordinate_type;
    using coarse_mesh_type = coarse_mesh_t<FineMesh>;
    using fine_face_type = typename fine_mesh_type::face_type;
    using coarse_cell_type = typename coarse_mesh_type::cell_type;

    /* This maps from Coarse Cells to Fine Faces */
    cc2ff_t<fine_mesh_type> cc2ff;

    for (const auto& fcl : fmsh) {
        auto di = fmsh.domain_info(fcl);
        auto ccl = cmsh.cell_at(di.tag());
        auto ffcs = faces(fmsh, fcl);
        cc2ff[ccl].insert(ffcs.begin(), ffcs.end());
    }

    return cc2ff;
}

template<typename FineMesh>
auto make_projectors(const FineMesh& fmsh, const coarse_mesh_t<FineMesh>& cmsh,
    cc2ff_t<FineMesh>& cc2ff, size_t coarse_degree, size_t fine_degree)
{
    using T = typename FineMesh::coordinate_type;
    using fine_mesh_type = FineMesh;
    using coarse_mesh_type = coarse_mesh_t<FineMesh>;
    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using SM = Eigen::SparseMatrix<T>;
    using triplet = Eigen::Triplet<T>;

    auto cbs = disk::scalar_basis_size(coarse_degree, coarse_mesh_type::dimension);
    auto fbs = disk::scalar_basis_size(fine_degree, fine_mesh_type::dimension-1);

    auto nrows = fbs * fmsh.faces_size();
    auto ncols = cbs * cmsh.cells_size();

    SM gproj = SM(nrows, ncols);

    std::vector<triplet> triplets;
    for (auto& ccl : cmsh) {
        /* Coarse cell basis */
        auto cphi = disk::basis::scaled_monomial_basis(cmsh, ccl, coarse_degree);

        auto col = cbs * offset(cmsh, ccl);

        const auto& ffcs = cc2ff.at(ccl);
        for (const auto& ffc : ffcs) {
            /* Fine face basis */
            auto fphi = disk::basis::scaled_monomial_basis(fmsh, ffc, fine_degree);

            auto row = fbs * offset(fmsh, ffc);

            auto C2F = integrate(fmsh, ffc, cphi, fphi);
            auto mass = integrate(fmsh, ffc, fphi, fphi);
            decltype(C2F) P = mass.llt().solve(C2F);

            for (size_t i = 0; i < fbs; i++) {
                for (size_t j = 0; j < cbs; j++) {
                    triplets.push_back( {row+i, col+j, P(i,j)} );
                }
            }
        }
    }
    gproj.setFromTriplets( triplets.begin(), triplets.end() );
    triplets.clear();

    return gproj;
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
}

template<disk::mesh_2D Mesh>
void dump_mesh(const Mesh& msh)
{
    for (size_t cli = 0; cli < msh.cells_size(); cli++) {
        auto cl = msh.cell_at(cli);

        std::stringstream ss;
        ss << "cell_" << cli << ".m";
        std::ofstream ofs(ss.str());
        auto bar = barycenter(msh, cl);
        ofs << "plot(" << bar.x() << "," << bar.y() << ",'*');\n";
        ofs << "hold on;";
        auto pts = points(msh, cl);
        for (size_t i = 0; i < pts.size(); i++) {
            auto p0 = pts[i];
            auto p1 = pts[(i+1)%pts.size()];
            if (i == 0)
            {
                ofs << "line([" << p0.x() << ", " << p1.x() << "], "
                            "[" << p0.y() << ", " << p1.y() << "], 'Color', 'red');\n";
            }
            else
            {
                ofs << "line([" << p0.x() << ", " << p1.x() << "], "
                            "[" << p0.y() << ", " << p1.y() << "], 'Color', 'blue');\n";
            }
        }
    }
}

/* If called with `base = 0`, simply make the subdomain
 * numbering zero-based. Otherwise, make the subdomain
 * numbering start from `base`. */
template<typename Mesh>
void rebase_subdomain_numbering(Mesh& msh, size_t base = 0)
{
    if (msh.cells_size() == 0)
        return;

    auto di0 = msh.domain_info(msh.cell_at(0));
    auto min_tag = di0.tag();

    for (size_t i = 1; i < msh.cells_size(); i++) {
        const auto& cl = msh.cell_at(i);
        auto di = msh.domain_info(cl);
        min_tag = std::min(min_tag, di.tag());
    }

    auto storage = msh.backend_storage();
    for (auto& di : storage->subdomain_info)
        di = disk::subdomain_descriptor( (di.tag() - min_tag) + base);
}

template<typename FineMesh>
void dg_to_hho(const FineMesh& msh, size_t degree,
    const Eigen::SparseMatrix<typename FineMesh::coordinate_type>& proj,
    const disk::dynamic_vector<typename FineMesh::coordinate_type>& dg_sol,
    disk::silo_database& silo)
{
    using namespace disk::hho::slapl;

    auto f = make_rhs_function(msh);

    degree_info di(degree);

    using mesh_type = FineMesh;
    using T = typename mesh_type::coordinate_type;

    disk::dynamic_vector<T> hho_sol = proj * dg_sol;

    auto fbs = disk::scalar_basis_size(degree, mesh_type::dimension - 1);

    std::vector<double> u;

    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);
        disk::dynamic_matrix<T> lhs = A+S;
        auto phiT = hho_space<mesh_type>::cell_basis(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);

        auto fcs = faces(msh, cl);
        disk::dynamic_vector<T> locsolF = disk::dynamic_vector<T>::Zero(fbs * fcs.size());
        for (size_t iF = 0; iF < fcs.size(); iF++) {
            const auto& fc = fcs[iF];
            auto gofs = offset(msh, fc);
            locsolF.segment(fbs*iF, fbs) = hho_sol.segment(fbs*gofs, fbs);
        }

        disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
        u.push_back(locsol(0));
    }

    silo.add_variable("srcmesh", "u_dg2hho", u, disk::zonal_variable_t);
}

int main(int argc, char **argv)
{
    using T = double;
    
    if (argc < 2) {
        std::cout << "Please specify a 2D GMSH mesh" << std::endl;
        return 1;
    }
    
    #if 0
    //disk::simplicial_mesh<T,2> srcmsh;
    //auto mesher = make_simple_mesher(srcmsh);
    //for (size_t l = 0; l < 1; l++)
    //    mesher.refine();
    

    disk::generic_mesh<T,2> srcmsh;
    auto mesher = make_fvca5_hex_mesher(srcmsh);
    mesher.make_level(4);

    partition_unit_square_mesh(srcmsh, 4);
    #endif

    timecounter tc;
    
    tc.tic();
    using mesh_type = disk::simplicial_mesh<T,2>;
    mesh_type finemsh;
    disk::gmsh_geometry_loader< mesh_type > loader;
    loader.read_mesh(argv[1]);
    loader.populate_mesh(finemsh);
    std::cout << "GMSH: " << tc.toc() << " seconds" << std::endl;

    rebase_subdomain_numbering(finemsh);

    disk::generic_mesh<T, 2> coarsemsh;
    
    tc.tic();
    agglomerate_by_subdomain(finemsh, coarsemsh);
    std::cout << "Agglomeration: " << tc.toc() << " seconds" << std::endl;
    

    tc.tic();
    auto cc2ff = make_cc2ff(finemsh, coarsemsh);
    std::cout << "cc2ff: " << tc.toc() << " seconds" << std::endl;

    std::vector<double> cell_partnum;
    for (auto& cl : finemsh) {
        auto di = finemsh.domain_info(cl);
        cell_partnum.push_back(di.tag());
    }

    std::vector<double> elemnum;
    for (const auto& cl : coarsemsh) {
        auto di = coarsemsh.domain_info(cl);
        elemnum.push_back(di.tag());
    }

    disk::silo_database silo;
    silo.create("agglo.silo");
    silo.add_mesh(finemsh, "srcmesh");
    silo.add_mesh(coarsemsh, "dstmesh");
    silo.add_variable("srcmesh", "partnum", cell_partnum, disk::zonal_variable_t);
    silo.add_variable("dstmesh", "cellnum", elemnum, disk::zonal_variable_t);
    

    //dump_mesh(coarsemsh);

    auto degree = 1;

    size_t min_faces = faces(coarsemsh, coarsemsh.cell_at(0)).size();
    size_t max_faces = faces(coarsemsh, coarsemsh.cell_at(0)).size();

    for (auto& cl : coarsemsh) {
        auto fcs = faces(coarsemsh, cl);
        min_faces = std::min(min_faces, fcs.size());
        max_faces = std::max(max_faces, fcs.size());
    }

    std::cout << "Min. number of faces in a cell: " << min_faces << std::endl;
    std::cout << "Max. number of faces in a cell: " << max_faces << std::endl;

    tc.tic();
    auto proj = make_projectors(finemsh, coarsemsh, cc2ff, degree+1, degree);
    std::cout << "projectors: " << tc.toc() << " seconds" << std::endl;

    /*
    for (degree = 0; degree < 5; degree++) {
        std::vector<double> etas = {0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0};
        for(auto& eta : etas) {
            std::cout << sgr::BRedfg << "Degree " << degree+1 << ", eta = " << eta << sgr::reset << std::endl; 
            dg_diffusion_solver(coarsemsh, degree+1, eta, silo);
        }
    }
    */
    
    auto sol = dg_diffusion_solver(coarsemsh, degree+1, 1.0, silo);
    
    std::cout << "DG to HHO" << std::endl;
    tc.tic();
    dg_to_hho(finemsh, degree, proj, sol, silo);
    std::cout << "  " << tc.toc() << " seconds" << std::endl;
    hho_diffusion_solver(coarsemsh, degree, silo);

    return 0;
}