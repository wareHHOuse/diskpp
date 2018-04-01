
#include "core/loaders/loader.hpp"

#include "core/revolution/bases"
#include "core/revolution/quadratures"

#include "timecounter.h"

template<typename Mesh>
void testperf_pre(const Mesh& msh, size_t degree)
{
	using T = typename Mesh::coordinate_type;

    std::vector< std::vector< std::pair<T, Matrix<T, Dynamic, Dynamic> > > > pre_bases;

	for (auto& cl : msh)
	{
		std::vector< std::pair<T, Matrix<T, Dynamic, Dynamic> > > p;

		auto cb = revolution::make_scalar_monomial_basis(msh, cl, degree);

		Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(cb.size(), cb.size());

		auto qps = revolution::integrate(msh, cl, 2*degree);
		for (auto& qp : qps)
		{
			auto dphi = cb.eval_functions(qp.point());
			p.push_back( std::make_pair(qp.weight(), dphi) );
		}

		pre_bases.push_back(p);
	}

	timecounter tc;
    tc.tic();

	size_t i = 0;
	for (auto& cl : msh)
	{
		auto cb = revolution::make_scalar_monomial_basis(msh, cl, degree);
		Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(cb.size(), cb.size());

		auto pb = pre_bases[i++];

		for (auto& p : pb)
		{
			stiff += p.first * p.second * p.second.transpose();
		}
	}

	tc.toc();

    std::cout << "Precomputed, degree " << degree << ": " << tc << std::endl;
}

template<typename Mesh>
void testperf_otf(const Mesh& msh, size_t degree)
{
	using T = typename Mesh::coordinate_type;

	timecounter tc;
    tc.tic();

	for (auto& cl : msh)
	{
		auto cb = revolution::make_scalar_monomial_basis(msh, cl, degree);

		Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(cb.size(), cb.size());

		auto qps = revolution::integrate(msh, cl, 2*degree);
		for (auto& qp : qps)
		{
			auto dphi = cb.eval_functions(qp.point());
			stiff += qp.weight() * dphi * dphi.transpose();
		}
	}

	tc.toc();

    std::cout << "On the fly, degree " << degree << ": " << tc << std::endl;
}

int main(void)
{
	using T = double;

    typedef disk::generic_mesh<T, 2>  mesh_type;

    std::vector< mesh_type > ret;

    mesh_type msh;
    disk::fvca5_mesh_loader<T, 2> loader;

    if ( !loader.read_mesh("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1") )
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);

    for (size_t k = 0; k < 8; k++)
    	testperf_otf(msh, k);

    for (size_t k = 0; k < 8; k++)
    	testperf_pre(msh, k);
    
    return 0;
}
