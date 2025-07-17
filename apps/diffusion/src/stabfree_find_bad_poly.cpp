#include <iostream>
#include <complex>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"

//#if 0
namespace disk {

template<typename T>
std::vector<typename generic_mesh<T, 2>::face>
faces(const generic_mesh<T, 2>& msh,
      const typename generic_mesh<T, 2>::cell& cl)
{
    using face_type = typename generic_mesh<T, 2>::face;
    using nodeid_type = typename generic_mesh<T, 2>::node_type::id_type;
    std::vector<face_type> ret;
    auto ptids = cl.point_ids();
    for (size_t i = 0; i < ptids.size(); i++) {
        auto p0 = nodeid_type(ptids[i]);
        auto p1 = nodeid_type(ptids[(i+1)%ptids.size()]);
        ret.push_back(face_type(p0,p1));
    }
    return ret;
}

template<typename T>
static_vector<T, 2>
normal(const generic_mesh<T, 2>& msh,
       const typename generic_mesh<T, 2>::cell& cl,
       const typename generic_mesh<T, 2>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() == 2);

    auto v = pts[1] - pts[0];
    auto n = (point<T,2>({v.y(), -v.x()})).to_vector();

    return n/n.norm();
}

}

//#endif

#include "diskpp/methods/hho"
#include "diskpp/bases/bases_utils.hpp"
#include "diskpp/output/silo.hpp"











enum class hho_variant {
    mixed_order_low,
    equal_order,
    mixed_order_high
};

template<disk::mesh_2D Mesh>
void
adjust_stabfree_recdeg(const Mesh& msh, const typename Mesh::cell_type& cl,
    disk::hho_degree_info& hdi)
{
    size_t cd = hdi.cell_degree();
    size_t fd = hdi.face_degree();
    bool is_mixed_high = (hdi.cell_degree() > hdi.face_degree());
    size_t n = faces(msh, cl).size();   
    size_t rpd = cd+2;

    /* HHO space dofs */
    size_t from = ((cd+2)*(cd+1))/2 + n*(fd+1);
    /* Reconstruction dofs */
    size_t to = ((rpd+2)*(rpd+1))/2;

    if (from <= to) {
        hdi.reconstruction_degree(rpd);
    }
    else {
        /* Every harmonic degree provides 2 additional dofs, therefore
         * we need an increment that it is sufficient to accomodate
         * (from-to) dofs => ((from - to) + (2-1))/2 */
        size_t incr = (from - to + 1)/2;
        hdi.reconstruction_degree(rpd+incr);
    }
}

template<disk::mesh_2D Mesh>
auto test(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;

    hho_variant hv = hho_variant::mixed_order_high;

    Eigen::Matrix<T,Mesh::dimension,Mesh::dimension> Id =
        Eigen::Matrix<T,Mesh::dimension,Mesh::dimension>::Identity();
    
    auto cl = msh[0];
    disk::hho_degree_info hdi(1);

    //#define STABFREE
    #ifdef STABFREE

    adjust_stabfree_recdeg(msh, cl, hdi);

    if (hv == hho_variant::mixed_order_high) {
        /* Build the standard reconstruction + projection on the cells */
        auto oper = make_shl_face_proj_harmonic(msh, cl, hdi, Id);
        A = oper.second;
        GR = oper.first;
        //std::cout << "+GR rows: " << GR.rows() << std::endl;
    } else {
        /* Build the nested reconstruction */
        auto oper = make_sfl(msh, cl, hdi, Id);
        A = oper.second;
        GR = oper.first;
        //std::cout << "-GR rows: " << GR.rows() << std::endl;
    }

    #else

    auto oper = make_scalar_hho_laplacian(msh, cl, hdi);
    A = oper.second;
    GR = oper.first;
    A = A + make_scalar_hho_stabilization(msh, cl, GR, hdi);
    
    #endif

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigs = A.eigenvalues();
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigsa = eigs.cwiseAbs();

    std::sort(eigsa.begin(), eigsa.end());

    //for (size_t i = 0; i < eigsa.size(); i++)
    //    std::cout << eigsa[i] << " ";
    //std::cout << std::endl;

    auto min_nonzero = eigsa[1];
    for (size_t i = 1; i < eigsa.size(); i++)
        min_nonzero = std::min(min_nonzero, eigsa[i]);

    return min_nonzero;
}


struct constraints {
    double  Rmin;
    double  Rmax;
    double  dmin;

    constraints()
        : Rmin(0.1), Rmax(1.0), dmin(0.1)
    {}
};

enum move_status {
    moved,
    hit_Rmin,
    hit_Rmax,
    too_close
};

template<disk::mesh_2D Mesh>
move_status
move_node(Mesh& msh, size_t which, const typename Mesh::point_type& delta)
{
    using point_type = Mesh::point_type;
    auto storage = msh.backend_storage();
    auto& mpts = storage->points;

    assert(which < mpts.size());
    mpts[which] = mpts[which] + delta;

    return move_status::moved;
}

template<disk::mesh_2D Mesh>
void
set_node(Mesh& msh, size_t which, const typename Mesh::point_type& p)
{
    using point_type = Mesh::point_type;
    auto storage = msh.backend_storage();
    auto& mpts = storage->points;

    assert(which < mpts.size());
    mpts[which] = p;
}

template<disk::mesh_2D Mesh>
struct candidate {
    double      value;
    using point_type = Mesh::point_type;
    point_type  point;

    bool operator<(const candidate& other) {
        return value < other.value;
    }
};

template<disk::mesh_2D Mesh>
std::vector<candidate<Mesh>>
minimize_step(Mesh& msh)
{
    double eps = 0.05;

    using point_type = Mesh::point_type;
    auto storage = msh.backend_storage();
    auto& mpts = storage->points;

    double origeig = test(msh);

    std::vector<candidate<Mesh>> candidates( mpts.size() );

    for (size_t i = 0; i < mpts.size(); i++)
    {
        point_type orig = mpts[i];

        // Move north
        mpts[i] = point_type(orig.x(), orig.y() + eps);
        double mxp = test(msh);
        // Move south
        mpts[i] = point_type(orig.x(), orig.y() - eps);
        double mxm = test(msh);
        // Move east
        mpts[i] = point_type(orig.x() + eps, orig.y());
        double myp = test(msh);
        // Move west
        mpts[i] = point_type(orig.x() - eps, orig.y());
        double mym = test(msh);
        // Restore
        mpts[i] = orig;
        
        double dx = 0.5*(mxp - mxm)/eps;
        double dy = 0.5*(myp - mym)/eps;

        point_type newp = orig - eps*point_type(dx,dy);

        mpts[i] = newp;
        double eig = test(msh);
        mpts[i] = orig;
        
        std::cout << i << ": " << dx << " " << dy << ", eig: " << eig << ", orig eig: " << origeig << std::endl;

        candidates[i].value = eig;
        candidates[i].point = newp;
    }

    std::sort(candidates.begin(), candidates.end());

    return candidates;
}

template<disk::mesh_2D Mesh>
void explore(Mesh& msh, size_t pid)
{
    disk::silo_database silo;
    std::string fn = "explore.silo";
    silo.create(fn);
    silo.add_mesh(msh, "polygon");

    size_t N = 201;
    double eps = 0.02;

    using point_type = Mesh::point_type;
    auto storage = msh.backend_storage();
    auto& mpts = storage->points;
    auto p0 = mpts[pid];
    point_type base(p0.x() - (N/2)*eps, p0.y() - (N/2)*eps);

    std::vector<point_type> coords;
    std::vector<double> vals;

    for (size_t i = 0; i < N; i++) {
        std::cout << "**** " << i << std::endl;
        for (size_t j = 0; j < N; j++) {
            point_type ofs(eps*i, eps*j);
            coords.push_back( base+ofs );
            mpts[pid] = base+ofs;
            vals.push_back( test(msh) );
        }
    }
    mpts[pid] = p0;

    silo.add_mesh(coords, "pmsh");
    silo.add_variable("pmsh", "mineig", vals);
}

int main(void)
{
    using T = double;
    using mesh_type = disk::generic_mesh<T,2>;
    using point_type = mesh_type::point_type;

    for (size_t i = 6; i < 11; i++) {
        std::cout << " **** Faces = " << i << " ****" << std::endl;
        mesh_type msh_gen;
        double scale = 1.0;
        disk::make_single_element_mesh(msh_gen, scale, i);

        disk::silo_database silo;
        std::string fn = "badpolys_" + std::to_string(i) + ".silo";
        silo.create(fn);
        silo.add_mesh(msh_gen, "mesh");

        minimize_step(msh_gen);

        explore(msh_gen, 0);

        break;
    }

    

    return 0;
}