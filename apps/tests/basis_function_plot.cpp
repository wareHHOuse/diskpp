#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/output/silo.hpp"

template<typename T>
static auto
make_test_points(const disk::cartesian_mesh<T,2>& msh,
    const typename disk::cartesian_mesh<T,2>::cell_type& cl, size_t N)
{
    using mesh_type = disk::cartesian_mesh<T,2>;
    using cell_type = typename mesh_type::cell_type;
    using point_type = typename mesh_type::point_type;
    
    std::vector<point_type> ret;

    auto pts = points(msh, cl);
    auto hx = (pts[1] - pts[0]).to_vector().norm()/(N+1);
    auto hy = (pts[3] - pts[0]).to_vector().norm()/(N+1);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            auto px = pts[0].x() + j*hx;
            auto py = pts[0].y() + i*hy;
            ret.push_back( point_type(px, py) );
        }
    }

    return ret;
}

int main(void)
{
    using T = double;
    using mesh_type = disk::cartesian_mesh<T,2>;

    mesh_type msh;
    disk::make_single_element_mesh(msh);

    auto cl = *msh.cells_begin();

    auto cb = disk::make_scalar_harmonic_top_basis(msh, cl, 3);
    cb.maximum_polynomial_degree(1);
    
    std::cout << cb.size() << std::endl;

    std::vector<std::vector<T>> funs;
    std::vector<std::vector<T>> grads_x;
    std::vector<std::vector<T>> grads_y;

    funs.resize( cb.size() );
    grads_x.resize( cb.size() );
    grads_y.resize( cb.size() );

    auto tps = make_test_points(msh, cl, 100);
    for (auto& tp : tps)
    {
        auto phi = cb.eval_functions(tp);
        auto dphi = cb.eval_gradients(tp);

        for (size_t i = 0; i < cb.size(); i++)
        {
            funs[i].push_back( phi(i) );  
            grads_x[i].push_back( dphi(i,0) ); 
            grads_y[i].push_back( dphi(i,1) ); 
        }
    }

    disk::silo_database db;
    db.create("basisfunctions.silo");
    db.add_mesh(tps, "pointmesh");

    for (size_t i = 0; i < cb.size(); i++) {
        std::string fun_name = "fun_" + std::to_string(i);
        db.add_variable("pointmesh", fun_name, funs[i]);

        std::string grad_x_name = "grad_x_" + std::to_string(i);
        db.add_variable("pointmesh", grad_x_name, grads_x[i]);

        std::string grad_y_name = "grad_y_" + std::to_string(i);
        db.add_variable("pointmesh", grad_y_name, grads_y[i]);
        
        std::string grad_name = "grad" + std::to_string(i);
        db.add_expression(grad_name, "{" + grad_x_name + ", " + grad_y_name + "}", DB_VARTYPE_VECTOR);
    }

    return 0;
}