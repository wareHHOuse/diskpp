#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"

//#if 0
namespace disk {

/*
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
*/

template<typename T>
static_vector<T, 2>
normal(const generic_mesh<T, 2>& msh,
       const typename generic_mesh<T, 2>::cell& cl,
       const typename generic_mesh<T, 2>::face& fc)
{
    auto pts = points(msh, cl);
    Eigen::Matrix<T, 3, 3> OM;
    OM << 1.0, pts[0].x(), pts[0].y(),
          1.0, pts[1].x(), pts[1].y(),
          1.0, pts[2].x(), pts[2].y();

    double orientation = OM.determinant();
    orientation /= std::abs(orientation);

    auto clptids = cl.point_ids();
    auto fcpts = fc.point_ids();
    assert(fcpts.size() == 2);

    for (size_t i = 0; i < clptids.size(); i++) {
        auto p0 = clptids[i];
        auto p1 = clptids[(i+1)%clptids.size()];
        if ( (p0 == fcpts[0] and p1 == fcpts[1]) or
             (p0 == fcpts[1] and p1 == fcpts[0]) ) {

            auto v = pts[(i+1)%clptids.size()] - pts[i];
            auto n = (point<T,2>({v.y(), -v.x()})).to_vector();

            return orientation*n/n.norm();
        }
    }

    throw std::logic_error("face not found");
}

}

//#endif

#include "diskpp/methods/hho"
#include "diskpp/bases/bases_utils.hpp"
#include "diskpp/output/silo.hpp"



#define SOL_ALL_SAFETIES_ON 1
#include <sol/sol.hpp>

struct config {
    size_t          degree;
    bool            use_stabfree;
    size_t          num_vertices;
    size_t          vertex;
    std::string     silo_fn;
    size_t          nsamples;
    double          eps;
    int             variant;
    bool            multi_thread;

    config() :
        degree(0), use_stabfree(true), num_vertices(3),
        vertex(0), nsamples(101), eps(0.02),
        variant(0), multi_thread(false)
    {}
};

template<typename T>
std::vector<T> make_T_vector() {
    return {};
}

sol::table init_config(sol::this_state L) {
    sol::state_view lua(L);
    sol::table module = lua.create_table();
    module.new_usertype<config>("config",
        sol::constructors<config()>(),
        "degree", &config::degree,
        "use_stabfree", &config::use_stabfree,
        "num_vertices", &config::num_vertices,
        "vertex", &config::vertex,
        "silo_fn", &config::silo_fn,
        "nsamples", &config::nsamples,
        "eps", &config::eps,
        "variant", &config::variant,
        "multi_thread", &config::multi_thread
    );

    using T = double;
    using point_type = disk::point<T,2>;
    module.new_usertype<point_type>("point",
        sol::constructors<point_type(),
            point_type(const T&, const T&)
            >(),
        "x", sol::resolve<T(void) const>(&point_type::template x<T>),
        "y", sol::resolve<T(void) const>(&point_type::template y<T>)
    );

    lua["make_point_vector"] = make_T_vector<point_type>;

    return module;
}

inline std::ostream&
operator<<(std::ostream& os, const config& cfg)
{
    os << "degree:       " << cfg.degree << std::endl;
    os << "use_stabfree: " << std::boolalpha << cfg.use_stabfree << std::endl;
    os << "num_vertices: " << cfg.num_vertices << std::endl;
    os << "vertex:       " << cfg.vertex << std::endl;
    os << "silo_fn:      " << cfg.silo_fn << std::endl;
    os << "nsamples:     " << cfg.nsamples << std::endl;
    os << "eps:          " << cfg.eps << std::endl;
    if (cfg.variant < 0)
        os << "variant:      mixed-order low" << std::endl;
    if (cfg.variant == 0)
        os << "variant:      equal-order" << std::endl;
    if (cfg.variant > 0)
        os << "variant:      mixed-order high" << std::endl;

    os << "multi_thread: " << std::boolalpha << cfg.multi_thread << std::endl;    
    
    return os;
}



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
auto test(const Mesh& msh, const config& cfg)
{
    using T = typename Mesh::coordinate_type;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GR;

    hho_variant hv = hho_variant::equal_order;
    if (cfg.variant < 0)
        hv = hho_variant::mixed_order_low;
    if (cfg.variant == 0)
        hv = hho_variant::equal_order;
    if (cfg.variant > 0)
        hv = hho_variant::mixed_order_high;

    Eigen::Matrix<T,Mesh::dimension,Mesh::dimension> Id =
        Eigen::Matrix<T,Mesh::dimension,Mesh::dimension>::Identity();
    
    auto cl = msh[0];
    disk::hho_degree_info hdi(cfg.degree);

    if (cfg.use_stabfree) {
        adjust_stabfree_recdeg(msh, cl, hdi);

        if (hv == hho_variant::mixed_order_high) {
            auto oper = make_shl_face_proj_harmonic(msh, cl, hdi, Id);
            A = oper.second;
            GR = oper.first;
        } else {
            auto oper = make_sfl(msh, cl, hdi, Id);
            A = oper.second;
            GR = oper.first;
        }
    }
    else {
        auto oper = make_scalar_hho_laplacian(msh, cl, hdi);
        A = oper.second;
        GR = oper.first;
        A = A + make_scalar_hho_stabilization(msh, cl, GR, hdi);
    }

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
    using point_type = typename Mesh::point_type;
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
    using point_type = typename Mesh::point_type;
    auto storage = msh.backend_storage();
    auto& mpts = storage->points;

    assert(which < mpts.size());
    mpts[which] = p;
}

template<disk::mesh_2D Mesh>
struct candidate {
    double      value;
    using point_type = typename Mesh::point_type;
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

    using point_type = typename Mesh::point_type;
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


    #if 0
    auto nthreads = std::thread::hardware_concurrency();

    std::mutex mtx;
    auto compute = [&](int tnum) {
        Mesh mymsh;
        msh.copy_to(mymsh);
        auto& my_mpts = mymsh.backend_storage()->points;

        size_t lines_per_proc = std::ceil( double(N)/nthreads );
        size_t from = tnum * lines_per_proc;
        size_t to = std::min( (tnum+1)*lines_per_proc, N );

        std::vector<point_type> my_coords;
        my_coords.reserve( N*lines_per_proc );
        std::vector<double> my_vals;
        my_vals.reserve( N*lines_per_proc );

        for (size_t i = from; i < to; i++) {
            for (size_t j = 0; j < N; j++) {
                point_type ofs(eps*i, eps*j);
                my_coords.push_back( base+ofs );
                my_mpts[cfg.vertex] = base+ofs;
                my_vals.push_back( test(mymsh, cfg) );
            }
        }
        my_mpts[cfg.vertex] = p0;

        mtx.lock();
        coords.insert( coords.end(), my_coords.begin(), my_coords.end() );
        vals.insert( vals.end(), my_vals.begin(), my_vals.end() );
        mtx.unlock();
    };

    std::cout << "Started computation with " << nthreads << " threads";
    std::cout << std::endl;

    std::vector<std::thread> threads;
    for (size_t i = 0; i < nthreads; i++)
        threads.push_back( std::thread(compute, i) );

    for (auto& t : threads)
        t.join();
    #endif


template<typename T>
bool
is_convex(const std::vector<disk::point<T,2>>& pts)
{
    if (pts.size() < 4)
        return false;

    for (size_t i = 0; i < pts.size(); i++) {
        auto v0 = pts[(i+1) % pts.size()] - pts[i];
        auto v1 = pts[(i+2) % pts.size()] - pts[(i+1) % pts.size()];

        if ( det(v0,v1) <= 0 )
            return true;
    }

    return false;
}


template<typename T>
int orientation(const disk::point<T,2>& p, const disk::point<T,2>& q, const disk::point<T,2>& r) {
    double val = (q.y() - p.y()) * (r.x() - q.x()) -
                 (q.x() - p.x()) * (r.y() - q.y());
    if (val == 0.0) return 0;  // colinear
    return (val > 0.0) ? 1 : 2; // clock or counterclockwise
}

// Check if point q lies on segment pr
template<typename T>
bool onSegment(const disk::point<T,2>& p, const disk::point<T,2>& q, const disk::point<T,2>& r) {
    return q.x() <= std::max(p.x(), r.x()) && q.x() >= std::min(p.x(), r.x()) &&
           q.y() <= std::max(p.y(), r.y()) && q.y() >= std::min(p.y(), r.y());
}

// Check if segments p1q1 and p2q2 intersect
template<typename T>
bool segmentsIntersect(const disk::point<T,2>& p1, const disk::point<T,2>& q1,
                       const disk::point<T,2>& p2, const disk::point<T,2>& q2) {
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special cases
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false;
}

// Check if a polygon is self-intersecting
template<typename T>
bool is_self_intersecting(const std::vector<disk::point<T,2>>& polygon) {
    int n = polygon.size();
    if (n < 4) return false; // Triangle or less cannot self-intersect

    for (int i = 0; i < n; ++i) {
        disk::point<T,2> p1 = polygon[i];
        disk::point<T,2> q1 = polygon[(i + 1) % n];

        for (int j = i + 1; j < n; ++j) {
            // Skip adjacent segments and same segment
            if (j == i || (j + 1) % n == i || (i + 1) % n == j)
                continue;

            disk::point<T,2> p2 = polygon[j];
            disk::point<T,2> q2 = polygon[(j + 1) % n];

            if (segmentsIntersect(p1, q1, p2, q2))
                return true;
        }
    }
    return false;
}













template<disk::mesh_2D Mesh>
void explore(Mesh& msh, const config& cfg)
{
    disk::silo_database silo;
    if (cfg.silo_fn.length() > 0) {
        silo.create(cfg.silo_fn);
        silo.add_mesh(msh, "polygon");
    }

    size_t N = cfg.nsamples;
    double eps = cfg.eps;

    using point_type = typename Mesh::point_type;
    auto storage = msh.backend_storage();
    auto& mpts = storage->points;
    
    is_self_intersecting(mpts);
    for(size_t i = 0; i < mpts.size(); i++)
    {
        std::vector<point_type> coords;
        std::vector<double> vals;

        auto p0 = mpts[i];
        point_type base(p0.x() - (N/2)*eps, p0.y() - (N/2)*eps);

        if (cfg.multi_thread) {
            auto nthreads = std::thread::hardware_concurrency();

            std::mutex mtx;
            auto compute = [&](int tnum) {
                Mesh mymsh;
                msh.copy_to(mymsh);
                auto& my_mpts = mymsh.backend_storage()->points;

                size_t lines_per_proc = std::ceil( double(N)/nthreads );
                size_t from = tnum * lines_per_proc;
                size_t to = std::min( (tnum+1)*lines_per_proc, N );

                std::vector<point_type> my_coords;
                my_coords.reserve( N*lines_per_proc );
                std::vector<double> my_vals;
                my_vals.reserve( N*lines_per_proc );

                for (size_t j = from; j < to; j++) {
                    for (size_t k = 0; k < N; k++) {
                        point_type ofs(eps*j, eps*k);
                        my_mpts[i] = base+ofs;
                        if ( is_self_intersecting(my_mpts) )
                            continue;
                        my_coords.push_back( base+ofs );
                        my_vals.push_back( test(mymsh, cfg) );
                    }
                }
                my_mpts[i] = p0;

                mtx.lock();
                coords.insert( coords.end(), my_coords.begin(), my_coords.end() );
                vals.insert( vals.end(), my_vals.begin(), my_vals.end() );
                mtx.unlock();
            };

            std::cout << "Started computation with " << nthreads << " threads";
            std::cout << std::endl;

            std::vector<std::thread> threads;
            for (size_t i = 0; i < nthreads; i++)
                threads.push_back( std::thread(compute, i) );

            for (auto& t : threads)
                t.join();
        }
        else {
            for (size_t j = 0; j < N; j++) {
                std::cout << "\rVertex " << i << ": " << j+1 << "/" << N << std::flush;
                for (size_t k = 0; k < N; k++) {
                    point_type ofs(eps*j, eps*k);
                    
                    mpts[i] = base+ofs;
                    if ( is_self_intersecting(mpts) )
                        continue;
                    coords.push_back( base+ofs );
                    vals.push_back( test(msh, cfg) );
                }
            }
            mpts[i] = p0;
            std::cout << std::endl;
        }


        if (cfg.silo_fn.length() > 0) {
            std::string meshname = "pmsh" + std::to_string(i);
            silo.add_mesh(coords, meshname);
            silo.add_variable(meshname, "mineig" + std::to_string(i), vals);
        }
    }
}


