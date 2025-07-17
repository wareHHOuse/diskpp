#include <iostream>
#include <complex>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"
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
    disk::hho_degree_info hdi(0);
    adjust_stabfree_recdeg(msh, cl, hdi);

    if (hv == hho_variant::mixed_order_high) {
        /* Build the standard reconstruction + projection on the cells */
        auto oper = make_shl_face_proj_harmonic(msh, cl, hdi, Id);
        A = oper.second;
        GR = oper.first;
        std::cout << "+GR rows: " << GR.rows() << std::endl;
    } else {
        /* Build the nested reconstruction */
        auto oper = make_sfl(msh, cl, hdi, Id);
        A = oper.second;
        GR = oper.first;
        std::cout << "-GR rows: " << GR.rows() << std::endl;
    }

    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> eigs = A.eigenvalues();
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigsa = eigs.cwiseAbs();

    std::sort(eigsa.begin(), eigsa.end());

    for (size_t i = 0; i < eigsa.size(); i++)
        std::cout << eigsa[i] << " ";
    std::cout << std::endl;

    auto min_nonzero = eigsa[1];
    for (size_t i = 1; i < eigsa.size(); i++)
        min_nonzero = std::min(min_nonzero, eigsa[i]);

    return min_nonzero;
}




int main(void)
{
    using T = double;
    using mesh_type = disk::generic_mesh<T,2>;
    using point_type = mesh_type::point_type;

    for (size_t i = 3; i < 11; i++) {
        std::cout << " **** Faces = " << i << " ****" << std::endl;
        mesh_type msh_gen;
        double scale = 1.0;
        disk::make_single_element_mesh(msh_gen, scale, i);

        disk::silo_database silo;
        std::string fn = "badpolys_" + std::to_string(i) + ".silo";
        silo.create(fn);
        silo.add_mesh(msh_gen, "mesh");

        std::cout << test(msh_gen) << std::endl;

        auto storage = msh_gen.backend_storage();
        auto& mpts = storage->points;

        mpts[0] = mpts[0] + point_type(-2, -2);

        std::cout << test(msh_gen) << std::endl;

        break;
    }

    return 0;
}