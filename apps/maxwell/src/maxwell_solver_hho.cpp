/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020, 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "diskpp/methods/hho"
#include "diskpp/methods/implementation_hho/curl.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

#include "mumps.hpp"
#include "paramloader.hpp"

#define EPS_0   (8.8541878128e-12)
#define MU_0    (4e-7*M_PI)

#define NODE_NAME_CONST     "const"
#define NODE_NAME_SIM       "sim"
#define NODE_NAME_MATERIALS "materials"
#define NODE_NAME_SILO      "silo"
#define NODE_NAME_DOMAIN    "domain"
#define NODE_NAME_BOUNDARY  "boundary"



enum class boundary_type {
    UNDEFINED,
    DIRICHLET,
    NEUMANN,
    IMPEDANCE,
    TFSF_INTERFACE,
};

enum class source_type {
    UNDEFINED,
    HOMOGENEOUS,
    NON_HOMOGENEOUS,
};

struct bc_descriptor
{
    boundary_type       condition;
    source_type         source;
    bool                has_parameter;
};

boundary_type strtobnd(const std::string& str)
{
    if (str == "dirichlet")
        return boundary_type::DIRICHLET;

    if (str == "neumann")
        return boundary_type::NEUMANN;

    if (str == "impedance")
        return boundary_type::IMPEDANCE;

    if (str == "tfsf")
        return boundary_type::TFSF_INTERFACE; 

    return boundary_type::UNDEFINED;
}

template<typename T>
struct complex3
{
    using complex_type = std::complex<T>;

    complex_type x, y, z;

    complex3()
    {}

    complex3(const complex_type& px, const complex_type& py, const complex_type& pz)
        : x(px), y(py), z(pz)
    {}
};

template<typename T>
std::ostream&
operator<<(std::ostream& os, const complex3<T>& c3)
{
    os << "[" << c3.x << ", " << c3.y << ", " << c3.z << "]";
    return os;
}

template<typename T>
class config_loader
{
    sol::state lua;

public:
    using point_type = disk::point<T,3>;
    using complex_type = std::complex<T>;
    using complex3_type = complex3<T>;
    using real_type = T;

private:
    void register_lua_usertypes()
    {
        sol::usertype<complex_type> complex_type_ut =
            lua.new_usertype<complex_type>("complex",
                sol::constructors<
                    complex_type(),
                    complex_type(const real_type&, const real_type&)
                    >()
                );
        
        sol::usertype<complex3_type> complex3_type_ut =
            lua.new_usertype<complex3_type>("complex3",
                sol::constructors<
                    complex3_type(),
                    complex3_type(const complex_type&, const complex_type&, const complex_type&)
                    >()
                );

        auto yr =[](const complex3_type& c3) -> auto { return real(c3.y); };
        auto yi =[](const complex3_type& c3) -> auto { return imag(c3.y); };
        lua["yr"] = yr;
        lua["yi"] = yi;
    }

public:
    config_loader()
    {
        lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::io);
        register_lua_usertypes();

        lua[NODE_NAME_CONST] = lua.create_table();
        lua[NODE_NAME_SIM] = lua.create_table();
        lua[NODE_NAME_MATERIALS] = lua.create_table();
        lua[NODE_NAME_SILO] = lua.create_table();
        
        lua[NODE_NAME_CONST]["eps0"] = EPS_0;
        lua[NODE_NAME_CONST]["mu0"] = MU_0;
        lua[NODE_NAME_DOMAIN] = lua.create_table();
        lua[NODE_NAME_BOUNDARY] = lua.create_table();
    }
    sol::state& lua_state() { return lua; }

    std::string mesh_filename() const {
        auto mfn = lua[NODE_NAME_SIM]["mesh_filename"];
        if (not mfn.valid())
        {
            std::cout << "[CONFIG]: Mesh file name not specified ";
            std::cout << "(sim.mesh_filename)" << std::endl;
            throw std::invalid_argument("sim.mesh_filename");
        }

        return mfn;
    }

    bool load(const std::string& fn)
    {
        bool success = true;
        lua.script_file(fn);

        if ( not lua[NODE_NAME_MATERIALS]["epsilon"].valid() )
        {
            std::cout << "[CONFIG]: Function epsilon() not defined" << std::endl;
            success = false;
        }

        if ( not lua[NODE_NAME_MATERIALS]["mu"].valid() )
        {
            std::cout << "[CONFIG]: Function mu() not defined" << std::endl;
            success = false;
        }

        if ( not lua[NODE_NAME_MATERIALS]["sigma"].valid() )
        {
            std::cout << "[CONFIG]: Function sigma() not defined" << std::endl;
            success = false;
        }

        if ( not lua[NODE_NAME_SIM]["frequency"].valid() )
        {
            std::cout << "[CONFIG]: Parameter frequency not defined" << std::endl;
            success = false;
        }

        return success;
    }

    /*
    bool load(Mesh& msh, const std::string& fn)
    {
        bool success = true;
        lua.script_file(fn);
        
        if ( not lua["epsilon"].valid() )
        {
            std::cout << "CONFIG PROBLEM: Function epsilon() not defined" << std::endl;
            success = false;
        }

        if ( not lua["mu"].valid() )
        {
            std::cout << "CONFIG PROBLEM: Function mu() not defined" << std::endl;
            success = false;
        }

        if ( not lua["sigma"].valid() )
        {
            std::cout << "CONFIG PROBLEM: Function sigma() not defined" << std::endl;
            success = false;
        }

        if ( not lua["frequency"].valid() )
        {
            std::cout << "CONFIG PROBLEM: Parameter frequency not defined" << std::endl;
            success = false;
        }

        if ( lua["mesh_transform"].valid() )
        {
            auto f = [&](const point_type& pt) -> point_type {
                typename point_type::value_type new_x, new_y, new_z;

                sol::tie(new_x, new_y, new_z) =
                    lua["mesh_transform"](pt.x(), pt.y(), pt.z());

                return point_type(new_x, new_y, new_z);
            };

            msh.transform(f);
        }

        return success;
    }
    */

    complex_type epsilon(size_t tag)
    {
        double re, im;
        sol::tie(re, im) = lua[NODE_NAME_MATERIALS]["epsilon"](tag);
        return std::complex<double>(re, im);
    }

    complex_type mu(size_t tag)
    {
        double re, im;
        sol::tie(re, im) = lua[NODE_NAME_MATERIALS]["mu"](tag);
        return std::complex<double>(re, im);
    }

    real_type sigma(size_t tag)
    {
        return lua[NODE_NAME_MATERIALS]["sigma"](tag);
    }

    complex_type volume_source(size_t tag, const point_type& pt)
    {
        return complex_type(0.0);
    }

    real_type frequency()
    {
        return lua[NODE_NAME_SIM]["frequency"];
    }

    size_t order()
    {
        return lua[NODE_NAME_SIM]["order"];
    }

    #if 0
    Eigen::Matrix<complex_type,3,1>
    dirichlet_data(size_t tag, const point_type& pt)
    {
        Eigen::Matrix<complex_type,3,1> ret = Eigen::Matrix<complex_type,3,1>::Zero();
        auto bnd = lua[NODE_NAME_BOUNDARY][tag];
        if (not bnd.valid())
            return ret;

        auto src = bnd["source"];
        if (not src.valid())
            return ret;

        priv::E_complex_field<real_type> field =
            bnd["source"](tag, pt.x(), pt.y(), pt.z());

        ret(0) = complex_type(field.Ex_re, field.Ex_im);
        ret(1) = complex_type(field.Ey_re, field.Ey_im);
        ret(2) = complex_type(field.Ez_re, field.Ez_im);

        return ret;
    }
    #endif


    bool is_dirichlet(size_t bndnum)
    {
        auto bnd_data = lua[NODE_NAME_BOUNDARY][bndnum];
        if (not bnd_data.valid())
            return true;

        std::string bndtype = bnd_data["kind"];
        if (bndtype == "e_field")
            return true;

        return false;
    }

    bc_descriptor
    boundary_type(size_t bndnum)
    {
        bc_descriptor ret;
        ret.condition = boundary_type::UNDEFINED;
        ret.source = source_type::UNDEFINED;
        ret.has_parameter = false;

        auto bnd_data = lua[NODE_NAME_BOUNDARY][bndnum];
        if ( bnd_data.valid() )
        {
            auto kind = bnd_data["kind"];
            if ( kind.valid() )
            {
                ret.condition = strtobnd(kind);
            
                auto source = bnd_data["source"];
                if ( source.valid() )
                    ret.source = source_type::NON_HOMOGENEOUS;
                else
                    ret.source = source_type::HOMOGENEOUS;
                
                auto value = bnd_data["value"];
                if ( value.valid() )
                    ret.has_parameter = true;
            }
        }

        return ret;
    }

    bool is_scattered_field_region(size_t di)
    {
        auto dom_data = lua[NODE_NAME_DOMAIN][di];
        if ( dom_data.valid() )
        {
            auto scattered_field = dom_data["scattered_field"];
            if (scattered_field.valid() and scattered_field == true)
                return true;
        }

        return false;
    }

    real_type
    bnd_param_value(size_t bndnum)
    {
        real_type z_data = lua[NODE_NAME_BOUNDARY][bndnum]["value"];
        return z_data;    
    }

    Eigen::Matrix<complex_type,3,1>
    eval_volume_source(size_t di, const point_type& pt)
    {
        Eigen::Matrix<complex_type,3,1> ret = Eigen::Matrix<complex_type,3,1>::Zero();

        auto dom_node = lua[NODE_NAME_DOMAIN][di];
        if (not dom_node.valid())
            return ret;
        
        auto source = dom_node["source"];
        if (not source.valid())
            return ret;

        complex3_type src = source(di, pt.x(), pt.y(), pt.z());

        ret(0) = src.x;
        ret(1) = src.y;
        ret(2) = src.z;
        return ret;
    }

    Eigen::Matrix<complex_type,3,1>
    eval_boundary_source(size_t bndnum, const point_type& pt)
    {
        sol::function source = lua[NODE_NAME_BOUNDARY][bndnum]["source"];
        complex3_type src = source(bndnum, pt.x(), pt.y(), pt.z());

        Eigen::Matrix<complex_type,3,1> ret;
        ret(0) = src.x;
        ret(1) = src.y;
        ret(2) = src.z;
        return ret;
    }
};
















template<typename Mesh, typename ScalT = typename Mesh::coordinate_type>
class maxwell_hho_assembler
{
    using coordinate_type = typename Mesh::coordinate_type;
    using scalar_type = ScalT;
    using mesh_type = Mesh;
    using face_type = typename Mesh::face_type;
    using cell_type = typename Mesh::cell_type;

    //Mesh                                msh;
    disk::hho_degree_info               hdi;

    std::vector<Triplet<scalar_type>>   triplets;
    std::vector<bool>                   is_dirichlet;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    const size_t INVALID_OFFSET = (size_t) ~0;

    size_t face_basis_size(void) const {
        return disk::vector_basis_size(hdi.face_degree(), 2, 2);
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
    Matrix<scalar_type, Dynamic, 1>     RHS;

    size_t                              syssz;
    size_t                              sysfcs;


    maxwell_hho_assembler()
    {}

    maxwell_hho_assembler(Mesh& msh, const disk::hho_degree_info& hdi,
                          const std::vector<bool>& is_dirichlet)
    {
        initialize(msh, hdi, is_dirichlet);
    }

    void clear()
    {
        triplets.clear();
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

        make_tables(msh);

        LHS = SparseMatrix<ScalT>(syssz, syssz);
        RHS = Matrix<ScalT, Dynamic, 1>::Zero(syssz);

        std::cout << "Assembler initialized: " << sysfcs << " faces in system, ";
        std::cout << syssz << " DoFs" << std::endl;
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<ScalT, Dynamic, Dynamic>& lhsc,
             const Matrix<ScalT, Dynamic, 1>& rhs,
             const Matrix<ScalT, Dynamic, 1>& dirichlet_data)
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
                {
                    RHS.segment(cofsi, fbs) += -lhsc.block(lofsi, lofsj, fbs, fbs)*dirichlet_data.segment(lofsj, fbs);
                    continue;
                }

                auto cofsj = get_system_offset(msh, fcs[fj]);
                for (size_t i = 0; i < fbs; i++)
                    for(size_t j = 0; j < fbs; j++)
                        triplets.push_back( Triplet<ScalT>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
            }

            RHS.segment(cofsi, fbs) += rhs.segment(fi*fbs, fbs);
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        std::cout << "Matrix has " << LHS.nonZeros() << " nonzeros." << std::endl; 
    }

    disk::dynamic_vector<ScalT>
    get_expanded_solution(const Mesh& msh, disk::dynamic_vector<ScalT>& sol)
    {
        auto fbs = face_basis_size();

        disk::dynamic_vector<ScalT> ret = disk::dynamic_vector<ScalT>::Zero( fbs*msh.faces_size() );

        for (size_t i = 0; i < sysfcs; i++)
        {
            auto in_offset = i*fbs;
            auto out_offset = expand_table.at(i)*fbs;
            ret.segment(out_offset, fbs) = sol.segment(in_offset, fbs);
        }

        return ret;
    }

    disk::dynamic_vector<ScalT>
    get_element_dofs(const Mesh& msh, const typename Mesh::cell& cl, disk::dynamic_vector<ScalT>& sol)
    {
        auto fbs = face_basis_size();
        auto fcs = faces(msh, cl);
        disk::dynamic_vector<ScalT> ret = disk::dynamic_vector<ScalT>::Zero( fbs*fcs.size() );

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            if ( not is_in_system(msh, fc) )
                continue;

            auto ofs = get_system_offset(msh, fc);
            ret.segment(i*fbs, fbs) = sol.segment(ofs, fbs);
        }

        return ret;
    }
};


template<typename Mesh>
struct hho_maxwell_solver_state
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = std::complex<double>;

    mesh_type                                       msh;
    disk::hho_degree_info                           hdi;
    maxwell_hho_assembler<mesh_type, scalar_type>   assm;
    disk::dynamic_vector<scalar_type>               sol;
    disk::dynamic_vector<scalar_type>               sol_full;
    disk::dynamic_vector<scalar_type>               reco;
};



template<typename Mesh, typename clT>
auto
compute_element_contribution(hho_maxwell_solver_state<Mesh>& state,
                             config_loader<clT>& cfg,
                             typename Mesh::cell_type& cl)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename hho_maxwell_solver_state<Mesh>::scalar_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& hdi = state.hdi;

    double eps0 = EPS_0, mu0 = MU_0; // XXX

    auto omega = 2*M_PI*cfg.frequency();
    auto jw = scalar_type(0,omega);
    auto jwmu0 = scalar_type(0,omega*mu0);
    auto k0sq = (eps0*omega)*(mu0*omega);

    auto rd = hdi.reconstruction_degree();
    auto cd = hdi.cell_degree();
    auto fd = hdi.face_degree();

    auto di = msh.domain_info(cl);


    auto epsr = cfg.epsilon( di.tag() );
    auto mur = cfg.mu( di.tag() );
    auto sigma = cfg.sigma( di.tag() );
    auto Z = std::sqrt( (jw*mur*mu0) / (sigma + jw*epsr*eps0) );
    auto CR = disk::curl_reconstruction_pk(msh, cl, hdi);
    auto ST = disk::curl_hdg_stabilization(msh, cl, hdi);
    auto MM = disk::make_vector_mass_oper(msh, cl, hdi);

    auto stabparam = omega*std::sqrt(real(epsr*eps0)/real(mur));
    //auto stabparam = omega*std::sqrt(epsr*eps0/mur);

    Matrix<scalar_type, Dynamic, Dynamic> lhst =
        (1./mur) * CR.second + stabparam*ST;

    Matrix<scalar_type, Dynamic, Dynamic> lhs =
        lhst - (k0sq*epsr - jwmu0*sigma)*MM;

    Matrix<scalar_type, Dynamic, 1> rhs =
        Matrix<scalar_type, Dynamic, 1>::Zero(lhs.rows());

    const auto cb = make_vector_monomial_basis(msh, cl, rd);
    auto cbs = cb.size();

    auto qps = disk::integrate(msh, cl, 2*cb.degree()+2);
    for (auto& qp : qps)
    {
        auto phi = cb.eval_functions(qp.point());
        rhs.segment(0, cbs) +=
            qp.weight() * phi * cfg.eval_volume_source(di.tag(), qp.point());
    }

    auto fcs = faces(msh, cl);
    auto fbs = disk::vector_basis_size(hdi.face_degree(), 2, 2);

    Matrix<scalar_type, Dynamic, 1> dirichlet_data =
        Matrix<scalar_type, Dynamic, 1>::Zero(fcs.size() * fbs);
    
    for (size_t iF = 0; iF < fcs.size(); iF++)
    {
        auto& fc = fcs[iF];
        auto bi = msh.boundary_info(fc);
        if (not bi.is_boundary())
            continue;

        auto face_tag = bi.tag();

        auto bcd = cfg.boundary_type(face_tag);

        if ( bcd.condition == boundary_type::IMPEDANCE )
        {
            auto bnd_Z = Z;
            if ( bcd.has_parameter )
                bnd_Z = cfg.bnd_param_value(face_tag);

            if ( bcd.source == source_type::NON_HOMOGENEOUS )
            {
                auto f = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
                    return cfg.eval_boundary_source(face_tag, pt);
                };
                auto jw = scalar_type(0,omega);
                auto [Y, y] = disk::make_plane_wave_term<std::complex<double>>(msh, fc, fd, f);
                lhs.block(cbs+iF*fbs, cbs+iF*fbs, fbs, fbs) += (jwmu0/bnd_Z)*Y;
                rhs.segment(cbs+iF*fbs, fbs) += 2.0*(jwmu0/bnd_Z)*y;
            }
            else
            {
                auto Y = disk::make_impedance_term(msh, fc, fd);
                lhs.block(cbs+iF*fbs, cbs+iF*fbs, fbs, fbs) += (jwmu0/bnd_Z)*Y;
            }
        }

        if ( bi.is_internal() and bcd.condition == boundary_type::TFSF_INTERFACE )
        {
            if ( cfg.is_scattered_field_region(di.tag()) )
            {
                auto bnd_Z = Z;
                if ( bcd.has_parameter )
                    bnd_Z = cfg.bnd_param_value(face_tag);

                auto n = normal(msh, cl, fc);

                auto f = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
                    return cfg.eval_boundary_source(face_tag, pt);
                };

                const auto fb = disk::make_vector_monomial_tangential_basis<
                    mesh_type, face_type, scalar_type>(
                        msh, fc, hdi.face_degree()); 

                disk::dynamic_matrix<scalar_type> tfsf_mass =
                    disk::dynamic_matrix<scalar_type>::Zero(fbs, fbs);

                disk::dynamic_vector<scalar_type> tfsf_rhs =
                    disk::dynamic_vector<scalar_type>::Zero(fbs);

                auto qps = disk::integrate(msh, fc, 2*hdi.face_degree());
                for (const auto& qp : qps)
                {
                    auto f_phi = fb.eval_functions(qp.point());
                    tfsf_mass += qp.weight() * f_phi * f_phi.transpose();
                    tfsf_rhs += qp.weight() * f_phi * f(qp.point());
                }

                disk::dynamic_vector<scalar_type> s =
                    disk::dynamic_vector<scalar_type>::Zero(lhs.rows()); 
                s.segment(cbs+iF*fbs, fbs) = tfsf_mass.ldlt().solve(tfsf_rhs);
                rhs += lhst*s;
                rhs.segment(cbs+iF*fbs, fbs) += (jwmu0/bnd_Z)*tfsf_rhs;
            }

        }

    }

    return std::make_tuple(lhs, rhs, dirichlet_data);
}



template<typename Mesh, typename clT>
void
assemble(hho_maxwell_solver_state<Mesh>& state, config_loader<clT>& cfg)
{
    auto& msh = state.msh;
    auto& hdi = state.hdi;
    auto& assm = state.assm;

    std::vector<bool> is_dirichlet;
    is_dirichlet.resize( msh.faces_size() );
    for (auto& fc : faces(msh))
    {
        auto bi = msh.boundary_info(fc);
        if (not bi.is_boundary())
            continue;

        if (bi.is_internal())
            continue;
        
        auto face_tag = bi.tag();

        auto bcd = cfg.boundary_type(face_tag);

        if ( bcd.condition == boundary_type::IMPEDANCE )
            continue;
        
        if ( bcd.condition == boundary_type::NEUMANN )
            continue;
        
        is_dirichlet[ offset(msh, fc) ] = true;
    }

    state.assm.initialize(state.msh, state.hdi, is_dirichlet);

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        if (cell_i%100 == 0)
        {
            std::cout << "\rAssemblying system matrix : " << cell_i << "/";
            std::cout << msh.cells_size() << std::flush;
        }
        auto cbs = disk::vector_basis_size(hdi.cell_degree(), 3, 3);
        auto [lhs, rhs, dirichlet_data] = compute_element_contribution(state, cfg, cl);
        auto [LC, bC] = disk::static_condensation(lhs, rhs, cbs);
        assm.assemble(msh, cl, LC, bC, dirichlet_data);
        cell_i++;
    }
    std::cout << "\rAssemblying system matrix : " << cell_i << "/";
    std::cout << msh.cells_size() << std::endl;
    assm.finalize();

    std::cout << "Assembly finished." << std::endl;
}

template<typename Mesh, typename clT>
void
solve(hho_maxwell_solver_state<Mesh>& state, config_loader<clT>& cfg)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename hho_maxwell_solver_state<Mesh>::scalar_type;

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol = state.sol;
    auto& hdi = state.hdi;

    sol = disk::dynamic_vector<scalar_type>::Zero(assm.syssz);
    std::cout << "Running MUMPS" << std::endl;
    sol = mumps_lu(assm.LHS, assm.RHS);

    std::cout << "Expanding solution" << std::endl;
    auto cd = hdi.cell_degree();
    auto fd = hdi.face_degree();
    auto cbs = disk::vector_basis_size(cd, 3, 3);
    auto fbs = disk::vector_basis_size(fd, 2, 2);
    auto fullsz = cbs*msh.cells_size() + fbs*msh.faces_size();

    state.sol_full = disk::dynamic_vector<scalar_type>::Zero(fullsz);
    state.reco = disk::dynamic_vector<scalar_type>::Zero(cbs*msh.cells_size());

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto [lhs, rhs, dd] = compute_element_contribution(state, cfg, cl);
        auto cb = disk::make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, state.hdi.cell_degree());

        auto edofs = assm.get_element_dofs(msh, cl, sol);
        edofs += dd;

        Matrix<scalar_type, Dynamic, 1> esol = disk::static_decondensation(lhs, rhs, edofs);
    
        state.sol_full.segment(cbs*cell_i, cbs) = esol.head(cbs);

        auto CR = disk::curl_reconstruction_pk(msh, cl, hdi);
        state.reco.segment(cbs*cell_i, cbs) = CR.first*esol;

        auto fcs = faces(msh, cl);
        for (size_t iF = 0; iF < fcs.size(); iF++)
        {
            auto& fc = fcs[iF];
            auto lofs = cbs + fbs*iF;
            auto gofs = cbs*msh.cells_size() + fbs*offset(msh, fc);
            state.sol_full.segment(gofs, fbs) = esol.segment(lofs, fbs);
        }
        cell_i++;
    }
    std::cout << "Done" << std::endl;
}

template<typename Mesh, typename clT>
void
save_to_silo(hho_maxwell_solver_state<Mesh>& state, config_loader<clT>& cfg)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = std::complex<double>;

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol = state.sol;
    auto& hdi = state.hdi;

    auto& lua = cfg.lua_state();

    std::vector<scalar_type> data_ex( msh.cells_size() );
    std::vector<scalar_type> data_ey( msh.cells_size() );
    std::vector<scalar_type> data_ez( msh.cells_size() );

    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto cb = disk::make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, state.hdi.cell_degree());
        auto cbs = cb.size();
        Matrix<scalar_type, Dynamic, 1> esolT = state.sol_full.segment(cbs*cell_i, cbs);

        auto bar = barycenter(msh, cl);

        auto phi = cb.eval_functions(bar);
        auto ls = phi.transpose()*esolT;

        data_ex[cell_i] = ls(0);
        data_ey[cell_i] = ls(1);
        data_ez[cell_i] = ls(2);

        cell_i++;
    }

    std::vector<double> data_mag_e( msh.cells_size() );

    for (size_t i = 0; i < msh.cells_size(); i++)
    {
        auto ex = data_ex[i];
        auto ey = data_ey[i];
        auto ez = data_ez[i];

        data_mag_e[i] = real(std::sqrt( ex*conj(ex) + ey*conj(ey) + ez*conj(ez) ));
    }

    std::string silo_fn = lua[NODE_NAME_SILO]["filename"];

    disk::silo_database silo_db;
    silo_db.create(silo_fn);
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<scalar_type> ex("ex", data_ex);
    silo_db.add_variable("mesh", ex);

    disk::silo_zonal_variable<scalar_type> ey("ey", data_ey);
    silo_db.add_variable("mesh", ey);

    disk::silo_zonal_variable<scalar_type> ez("ez", data_ez);
    silo_db.add_variable("mesh", ez);

    disk::silo_zonal_variable<double> mag_e("mag_e", data_mag_e);
    silo_db.add_variable("mesh", mag_e);
}

template<typename Mesh, typename clT>
clT
compute_return_loss(hho_maxwell_solver_state<Mesh>& state,
    config_loader<clT>& cfg, size_t face_tag)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = std::complex<double>;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol_full = state.sol_full;
    auto& hdi = state.hdi;

    auto cbs = disk::vector_basis_size(hdi.cell_degree(), 3, 3);
    auto num_cells = msh.cells_size();
    auto c_ofs = cbs*num_cells;

    size_t elem_count = 0;
    scalar_type P_incident = 0.0;
    scalar_type P_reflected_1 = 0.0;
    scalar_type P_reflected_2 = 0.0;
    scalar_type P_total = 0.0;

    for ( auto& fc : faces(msh) )
    {
        auto bi = msh.boundary_info(fc);
        auto tag = bi.tag();

        if (bi.tag() == face_tag)
        {
            const auto fb = disk::make_vector_monomial_tangential_basis<
                mesh_type, face_type, scalar_type>(
                    msh, fc, hdi.face_degree());
            auto fbs = fb.size();
            
            auto Y = disk::make_impedance_term(msh, fc, hdi.face_degree());

            auto ofs = c_ofs + fbs*offset(msh,fc);

            Matrix<scalar_type, Dynamic, 1> fdofs = sol_full.segment(ofs, fbs);
            const auto qps_f = integrate(msh, fc, 2*hdi.face_degree() );
            for (auto& qp : qps_f)
            {
                Matrix<scalar_type, Dynamic, 3> fphi = fb.eval_functions(qp.point());
                Matrix<scalar_type, 3, 1> Ecalc = Matrix<scalar_type, 3, 1>::Zero();

                for (size_t kk = 0; kk < fdofs.size(); kk++)
                    Ecalc += fdofs(kk)*fphi.row(kk);

                auto f = [&](const point_type& pt) -> Matrix<std::complex<double>,3,1> {
                    return cfg.eval_boundary_source(face_tag, pt);
                };

                Matrix<scalar_type, 3, 1> Einc = f(qp.point());
                P_reflected_1 += qp.weight() * (Ecalc - Einc).dot(Einc);
                P_reflected_2 += qp.weight() * (Ecalc - Einc).dot(Ecalc-Einc);
                P_incident += qp.weight() * (Einc).dot(Einc);
                P_total += qp.weight() * (Ecalc).dot(Ecalc);
            }
            elem_count++;
        }
    }

    std::cout << elem_count << std::endl;
    std::cout << "Ref1:  " << P_reflected_1 << std::endl;
    std::cout << "Ref2:  " << P_reflected_2 << std::endl;
    std::cout << "Inc:   " << P_incident << std::endl;
    std::cout << "Total: " << P_total << std::endl;

    std::cout << "S11 (ref1) = " << 20.0*log10(  abs(P_reflected_1/P_incident) ) << std::endl;
    std::cout << "S11 (ref2) = " << 10.0*log10(  abs(P_reflected_2/P_incident) ) << std::endl;

    return 20.0*log10(  abs(P_reflected_1/P_incident) );
}

template<typename Mesh, typename clT>
clT
compute_error(hho_maxwell_solver_state<Mesh>& state, config_loader<clT>& cfg)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename hho_maxwell_solver_state<Mesh>::scalar_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol_full = state.sol_full;
    auto& hdi = state.hdi;

    auto& lua = cfg.lua_state();
    sol::function lua_ref_fun = lua["reference_solution"];

    size_t cell_i = 0;
    scalar_type err = 0.0;
    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);
        auto enabled = lua["err_enabled"];
        if ( enabled.valid() )
            if ( not enabled( di.tag() ) )
                continue;

        auto ref_fun = [&](const point_type& pt) -> Matrix<scalar_type, 3, 1> {
            complex3<clT> field = lua_ref_fun(di.tag(), pt.x(), pt.y(), pt.z());
            Matrix<scalar_type, 3, 1> ret;
            ret(0) = field.x;
            ret(1) = field.y;
            ret(2) = field.z;
            return ret;
        };

        auto cb = disk::make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, state.hdi.cell_degree());
        auto cbs = cb.size();
        Matrix<scalar_type, Dynamic, 1> num_sol = state.sol_full.segment(cbs*cell_i, cbs);

        Matrix<scalar_type, Dynamic, Dynamic> MM = 
            Matrix<scalar_type, Dynamic, Dynamic>::Zero(cbs, cbs);
    
        Matrix<scalar_type, Dynamic, 1> rhs =
            Matrix<scalar_type, Dynamic, 1>::Zero(cbs);

        auto qps = disk::integrate(msh, cl, 2*hdi.cell_degree()+1);
        for (auto& qp : qps)
        {
            auto phi = cb.eval_functions(qp.point());
            MM += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * phi * ref_fun( qp.point() );

            Matrix<scalar_type, 3, 1> ls = phi.transpose()*num_sol;
            Matrix<scalar_type, 3, 1> vdiff = ls - ref_fun(qp.point());
            err += qp.weight() * vdiff.dot(vdiff);
        }
        
        //Matrix<scalar_type, Dynamic, 1> ref_sol = MM.ldlt().solve(rhs);
        //Matrix<scalar_type, Dynamic, 1> diff = ref_sol - num_sol;
        //err += diff.dot(MM*diff);
        cell_i++;
    }

    return std::sqrt( real(err) );
}

template<typename Mesh, typename clT>
clT
compute_reconstruction_error(hho_maxwell_solver_state<Mesh>& state, config_loader<clT>& cfg)
{
    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using scalar_type = typename hho_maxwell_solver_state<Mesh>::scalar_type;
    using point_type = typename Mesh::point_type;

    auto& msh = state.msh;
    auto& assm = state.assm;
    auto& sol_full = state.sol_full;
    auto& hdi = state.hdi;

    auto& lua = cfg.lua_state();
    sol::function lua_ref_fun = lua["reference_reconstruction"];

    size_t cell_i = 0;
    scalar_type err_int = 0.0;
    scalar_type err_MM = 0.0;

    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);

        auto ref_fun = [&](const point_type& pt) -> Matrix<scalar_type, 3, 1> {
            complex3<clT> field = lua_ref_fun(di.tag(), pt.x(), pt.y(), pt.z());
            Matrix<scalar_type, 3, 1> ret;
            ret(0) = field.x;
            ret(1) = field.y;
            ret(2) = field.z;
            return ret;
        };

        auto cb = disk::make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, state.hdi.cell_degree());
        auto cbs = cb.size();
        Matrix<scalar_type, Dynamic, 1> num_sol = state.reco.segment(cbs*cell_i, cbs);

        Matrix<scalar_type, Dynamic, Dynamic> MM = 
            Matrix<scalar_type, Dynamic, Dynamic>::Zero(cbs, cbs);
    
        Matrix<scalar_type, Dynamic, 1> rhs =
            Matrix<scalar_type, Dynamic, 1>::Zero(cbs);

        auto qps = disk::integrate(msh, cl, 2*hdi.cell_degree()+1);
        for (auto& qp : qps)
        {
            auto phi = cb.eval_functions(qp.point());
            MM += qp.weight() * phi * phi.transpose();
            rhs += qp.weight() * phi * ref_fun( qp.point() );

            Matrix<scalar_type, 3, 1> ls = phi.transpose()*num_sol;
            Matrix<scalar_type, 3, 1> vdiff = ls - ref_fun(qp.point());
            err_int += qp.weight() * vdiff.dot(vdiff);
        }
        
        Matrix<scalar_type, Dynamic, 1> ref_sol = MM.ldlt().solve(rhs);
        Matrix<scalar_type, Dynamic, 1> diff = ref_sol - num_sol;
        err_MM += diff.dot(MM*diff);
        cell_i++;
    }

    std::cout << err_int << " " << err_MM << std::endl;

    return std::sqrt( real(err_int) );
}

template<typename Mesh, typename clT>
void
register_lua_functions(hho_maxwell_solver_state<Mesh>& state, config_loader<clT>& cfg)
{
    auto lua_assemble = [&]() {
        assemble(state, cfg);
    };

    auto lua_solve = [&]() {
        solve(state, cfg);
    };

    auto lua_save_to_silo = [&]() {
        save_to_silo(state, cfg);
    };

    auto lua_compute_error = [&]() {
        return compute_error(state, cfg);
    };

    auto lua_compute_reconstruction_error = [&]() {
        return compute_reconstruction_error(state, cfg);
    };

    auto lua_compute_return_loss = [&](size_t face_tag) {
        return compute_return_loss(state, cfg, face_tag);
    };

    auto& lua = cfg.lua_state();

    lua["assemble"] = lua_assemble;
    lua["solve"] = lua_solve;
    lua["save_to_silo"] = lua_save_to_silo;
    lua["compute_error"] = lua_compute_error;
    lua["compute_reconstruction_error"] = lua_compute_reconstruction_error;
    lua["compute_return_loss"] = lua_compute_return_loss;
}

template<typename Mesh, typename clT>
void
run_maxwell_solver(hho_maxwell_solver_state<Mesh>& state, config_loader<clT>& cfg)
{
    rusage_monitor rm;

    size_t order = cfg.order();
    disk::hho_degree_info chdi( disk::priv::hdi_named_args{
                                  .rd = (size_t) order,
                                  .cd = (size_t) order,
                                  .fd = (size_t) order
                                } );

    state.hdi = chdi;

    auto num_elems = state.msh.cells_size();
    std::cout << "Mesh has " << num_elems << " elements" << std::endl;

    std::cout << "Registering lua stuff" << std::endl;
    register_lua_functions(state, cfg);

    auto& lua = cfg.lua_state();

    std::cout << "Calling user code" << std::endl;
    lua["xxx"]();
}

int main(int argc, const char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <param filename>" << std::endl;
        return 1;
    }

    using real_type = double;
    using complex_type = std::complex<real_type>;

    config_loader<real_type> cfg;
    cfg.load(argv[1]);

    auto mesh_filename = cfg.mesh_filename();

    /* Netgen 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.mesh$") ))
    {
        //std::cout << "Guessed mesh format: Netgen 3D" << std::endl;
        //auto msh = disk::load_netgen_3d_mesh<real_type>(mesh_filename);

        //run_maxwell_solver(msh, cfg);

        //return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3g$") ))
    {
        std::cout << "Guessed mesh format: GMSH polyhedral" << std::endl;

        using mesh_type = disk::generic_mesh<real_type,3>;
        hho_maxwell_solver_state<mesh_type> state;

        disk::gmsh_geometry_loader< mesh_type > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(state.msh);

        run_maxwell_solver(state, cfg);

        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<real_type,3>;
        hho_maxwell_solver_state<mesh_type> state;
        
        disk::gmsh_geometry_loader< mesh_type > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(state.msh);

        run_maxwell_solver(state, cfg);

        return 0;
    }

    std::cout << "Unknown mesh format, can't proceed." << std::endl;
    return 1;
}
