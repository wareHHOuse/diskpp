/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

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
 */
#define ENABLE_CCW_FACES

#include "diskpp/common/eigen.hpp"
#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/dg"
#include "mumps.hpp"

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

    timecounter tc;
    tc.tic();
    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);
        disk::dynamic_matrix<T> lhs = A+S;
        auto phiT = typename hho_space<mesh_type>::cell_basis_type(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
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
    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);       
        disk::dynamic_matrix<T> lhs = A+S;
        auto phiT = typename hho_space<mesh_type>::cell_basis_type(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
        auto MMe = integrate(msh, cl, phiT, phiT);
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
void
dg_diffusion_solver(Mesh& msh, size_t degree,
    const typename Mesh::coordinate_type eta,
    disk::silo_database& silo)
{   
    std::cout << "DG solver" << std::endl;
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto f = make_rhs_function(msh);

    auto cbs = disk::scalar_basis_size(degree, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_assembler(msh, cbs);
    
    timecounter tc;
    tc.tic();
    for (auto& tcl : msh)
    {
        auto tbasis = disk::make_scalar_monomial_basis(msh, tcl, degree);
        auto qps = disk::integrate(msh, tcl, 2*degree);
        
        matrix_type K = matrix_type::Zero(tbasis.size(), tbasis.size());
        vector_type loc_rhs = vector_type::Zero(tbasis.size());
        for (auto& qp : qps)
        {
            auto ep = qp.point();
            auto phi = tbasis.eval_functions(ep);
            auto dphi = tbasis.eval_gradients(ep);
            
            K += qp.weight() * dphi * dphi.transpose();
            loc_rhs += qp.weight() * phi * f(qp.point());
        }
        
        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            
            auto n     = normal(msh, tcl, fc);
            auto eta_l = eta / diameter(msh, fc);
            auto f_qps = disk::integrate(msh, fc, 2*degree);
            
            if (nv) {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
                
                auto ncl = nv.value();
                auto nbasis = disk::make_scalar_monomial_basis(msh, ncl, degree);
                assert(tbasis.size() == nbasis.size());
                
                for (auto& fqp : f_qps) {
                    auto ep     = fqp.point();
                    auto tphi   = tbasis.eval_functions(ep);
                    auto tdphi  = tbasis.eval_gradients(ep);
                
                    /* NOT on a boundary */
                    Att += + fqp.weight() * eta_l * tphi * tphi.transpose();
                    Att += - fqp.weight() * 0.5 * tphi * (tdphi*n).transpose();
                    Att += - fqp.weight() * 0.5 * (tdphi*n) * tphi.transpose();
                
                    auto nphi   = nbasis.eval_functions(ep);
                    auto ndphi  = nbasis.eval_gradients(ep);
                
                    Atn += - fqp.weight() * eta_l * tphi * nphi.transpose();
                    Atn += - fqp.weight() * 0.5 * tphi * (ndphi*n).transpose();
                    Atn += + fqp.weight() * 0.5 * (tdphi*n) * nphi.transpose();
                }
                assm.assemble(msh, tcl, tcl, Att);
                assm.assemble(msh, tcl, ncl, Atn);
            }
            else {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                for (auto& fqp : f_qps) {
                    auto ep     = fqp.point();
                    auto tphi   = tbasis.eval_functions(ep);
                    auto tdphi  = tbasis.eval_gradients(ep);
                    
                    /* On a boundary*/
                    Att += + fqp.weight() * eta_l * tphi * tphi.transpose();
                    Att += - fqp.weight() * tphi * (tdphi*n).transpose();
                    Att += - fqp.weight() * (tdphi*n) * tphi.transpose();
                    
                    //loc_rhs -= fqp.weight() * (tdphi*n);
                    //loc_rhs += fqp.weight() * eta_l * tphi;
                }
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

    std::vector<double> data;

    T err = 0.0; size_t cell_i = 0;
    tc.tic();
    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, degree);
        Matrix<T, Dynamic, Dynamic> MMe = disk::make_mass_matrix(msh, cl, cb);
        Matrix<T, Dynamic, 1> arhs = disk::make_rhs(msh, cl, cb, sol_fun);
        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());
        Matrix<T, Dynamic, 1> diff = lsol - asol;
        err += diff.dot(MMe*diff);
        data.push_back( lsol(0) );
        cell_i++;
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error: " << std::sqrt(err) << std::endl;
    disk::silo_zonal_variable<T> u("u_dg", data);
    silo.add_variable("dstmesh", u);
}

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

template<typename SrcMesh>
void agglomerate_by_subdomain(const SrcMesh& srcmsh,
    disk::generic_mesh<typename SrcMesh::coordinate_type, 2>& dstmsh)
{
    using src_mesh_type = SrcMesh;
    using coord_type = typename SrcMesh::coordinate_type;
    using dst_mesh_type = disk::generic_mesh<coord_type, 2>;
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
        
        /* Reverse the path if vertices are in clockwise order */
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

    /* Sort the edges and make them unique */
    std::sort(dststor->edges.begin(), dststor->edges.end());
    auto last = std::unique(dststor->edges.begin(), dststor->edges.end());
    dststor->edges.erase(last, dststor->edges.end());

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

        dst_surf_type newsurf(surfedges);
        newsurf.set_point_ids(path.begin(), path.end());
        dststor->surfaces.push_back(newsurf);
    }

    std::sort(dststor->surfaces.begin(), dststor->surfaces.end());
    dststor->subdomain_info.resize( dststor->surfaces.size() );
    dststor->boundary_info.resize( dststor->edges.size() );


    std::vector<size_t> counts(dststor->edges.size());

    for (auto& cl : dstmsh) {
        auto fcs = faces(dstmsh, cl);
        for (auto& fc : fcs) {
            counts[offset(dstmsh, fc)]++;
        }
    }

    for ( size_t i = 0; i < counts.size(); i++) {
        if (counts[i] == 1) {
            dststor->boundary_info[i].is_boundary(true);
            dststor->boundary_info[i].id(0);
        }
    }

    disk::mark_internal_faces(dstmsh);
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

    
    using mesh_type = disk::simplicial_mesh<T,2>;
    mesh_type srcmsh;
    disk::gmsh_geometry_loader< mesh_type > loader;
    loader.read_mesh(argv[1]);
    loader.populate_mesh(srcmsh);
    

    

    std::vector<double> cp;
    for (auto& cl : srcmsh) {
        auto di = srcmsh.domain_info(cl);
        cp.push_back(di.tag());
    }


    disk::generic_mesh<T, 2> dstmsh;

    agglomerate_by_subdomain(srcmsh, dstmsh);


    disk::silo_database silo;
    silo.create("agglo.silo");
    silo.add_mesh(srcmsh, "srcmesh");
    silo.add_mesh(dstmsh, "dstmesh");
    silo.add_variable("srcmesh", "partnum", cp, disk::zonal_variable_t);
    
    std::map<size_t, double> areas;
    for (auto& scl : srcmsh) {
        auto tag = srcmsh.domain_info(scl).tag();
        areas[tag] += measure(srcmsh, scl);
    }


    //for (auto& [tag, area] : areas)
    //    std::cout << tag << " " << area << std::endl;

    double tot_area = 0.0;
    for (auto& dcl : dstmsh) {
        tot_area += measure(dstmsh, dcl);
    }

    std::cout << "Total mesh area: " << tot_area << std::endl;

    dump_mesh(dstmsh);

    auto degree = 2;

    size_t min_faces = faces(dstmsh, dstmsh.cell_at(0)).size();
    size_t max_faces = faces(dstmsh, dstmsh.cell_at(0)).size();

    for (auto& cl : dstmsh) {
        auto fcs = faces(dstmsh, cl);
        min_faces = std::min(min_faces, fcs.size());
        max_faces = std::max(max_faces, fcs.size());
    }

    std::cout << "Min. number of faces in a cell: " << min_faces << std::endl;
    std::cout << "Max. number of faces in a cell: " << max_faces << std::endl;

    dg_diffusion_solver(dstmsh, degree+1, 100.0, silo);
    hho_diffusion_solver(dstmsh, degree, silo);

    return 0;
}