/*
 * contact_report.hpp -- the one who informs.
 *
 * Two layers, deliberately split:
 *
 *   evaluate()   model-agnostic RAW quantities, one FaceQ per face of the
 *                contact view: sigma_n, ||sigma_t||, augmented tractions
 *                P_n, ||P_t||, u_n, u_t.  This is the shared machinery that
 *                was duplicated across the three drivers.
 *
 *   write_*()    model-specific derived fields (ratios, zones, thresholds)
 *                + write_contact_vtu.  Tresca compares against a constant
 *                s; Coulomb computes s(u) = -F [P_n]_- per face and
 *                classifies separation / slip / stick.
 *
 * The augmented tractions use the theta = 1 form (P = sigma - gamma u), as
 * in the constraint itself, matching the original drivers: they are what
 * the formulation actually bounds, and the gap against the discrete stress
 * measures consistency of the discretisation.
 */

#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include "diskpp/methods/hho"

#include "spatial_operator.hpp"
#include "newton_solver.hpp"
#include "write_vtu.hpp"

namespace hho_contact {

template<typename Mesh>
class ContactReport
{
public:
    using T           = typename Mesh::coordinate_type;
    using matrix_type = disk::dynamic_matrix<T>;
    using vector_type = disk::dynamic_vector<T>;
    using tensor_vec  = disk::static_vector<T, Mesh::dimension>;
    using bc_type     = disk::vector_boundary_conditions<Mesh>;

    struct FaceQ
    {
        typename Mesh::face_type fc;
        T sigma_n;       // discrete normal stress
        T sigma_t_norm;  // ||sigma_t(u)||
        T P_n;           // augmented normal traction  sigma_n - gamma_n u_n
        T P_t_norm;      // ||sigma_t - gamma_t u_t||
        T u_n;           // normal displacement
        T u_t0;          // first tangential component (signed; the 2D scalar)
        T u_t_norm;      // ||u_t||  (what to look at in 3D)
        T bx, by, bz;    // face barycentre (bz = 0 in 2D)
    };

private:
    const NewtonSolver<Mesh>&     solver;
    const spatial_operator<Mesh>& op;
    const bc_type&                bnd_C;    // which faces count as contact
    T gn0, gt0;

public:
    ContactReport(const NewtonSolver<Mesh>& sv, const spatial_operator<Mesh>& op_,
                  const bc_type& contact_view, T gamma_n_0, T gamma_t_0)
        : solver(sv), op(op_), bnd_C(contact_view), gn0(gamma_n_0), gt0(gamma_t_0)
    {}

    // ---- raw quantities, one per contact face ---------------------------
    std::vector<FaceQ> evaluate(const hho_state<T>& st) const
    {
        std::vector<FaceQ> out;

        size_t ci = 0;
        for (auto& cl : op.msh)
        {
            auto fcs = disk::faces(op.msh, cl);

            bool touches = false;
            for (auto& fc : fcs)
                if (bnd_C.is_contact_face(fc)) touches = true;
            if (!touches) { ci++; continue; }

            const auto cell_infos = op.di.cellDegreeInfo(op.msh, cl);
            auto grad_basis = disk::make_sym_matrix_monomial_basis(
                                  op.msh, cl, cell_infos.grad_degree());
            vector_type u_local = solver.gather(cl, ci, st);

            size_t face_offset = op.cbs;
            for (auto& fc : fcs)
            {
                const auto   face_info = op.di.degreeInfo(op.msh, fc);
                const size_t face_dofs = disk::vector_face_dofs(op.msh, face_info);

                if (bnd_C.is_contact_face(fc))
                {
                    /* same penalty the solve uses, from the single
                     * definition in contact_terms.hpp -- otherwise the
                     * reported P_n is not the P_n being solved for */
                    const T gamma_n_h = disk::nitsche_gamma(op.msh, fc, gn0, op.mu);
                    const T gamma_t_h = disk::nitsche_gamma(op.msh, fc, gt0, op.mu);

                    const auto n          = disk::normal(op.msh, cl, fc);
                    const auto barycentre = disk::barycenter(op.msh, fc);
                    auto face_basis = disk::make_vector_monomial_basis(
                                          op.msh, fc, face_info.degree());

                    auto traction_op = disk::eval_traction_operator(
                            op.msh, grad_basis, op.G[ci],
                            barycentre, n, op.mu, op.lam);
                    /* The augmented tractions below must read the SAME
                     * displacement trace the solve used, or the report
                     * measures a different quantity than the one solved
                     * for: cell trace u_T|_F in the cell variant, face
                     * trace u_F otherwise.                              */
                    auto trace_op    = disk::eval_displacement_trace_operator(
                            op.msh, cl, face_basis, op.cell_deg, barycentre,
                            traction_op.rows(), face_offset, op.variant);

                    vector_type sigma_n_op = traction_op * n;
                    vector_type trace_n_op = trace_op    * n;

                    // augmented normal traction (theta = 1 form)
                    const T P_n = (sigma_n_op - gamma_n_h*trace_n_op)
                                  .dot(u_local);

                    // augmented tangential traction
                    matrix_type P_t_op =
                          (traction_op - sigma_n_op * n.transpose())
                        - gamma_t_h * (trace_op - trace_n_op * n.transpose());
                    tensor_vec P_t = P_t_op.transpose() * u_local;

                    // discrete stress and displacement
                    tensor_vec traction = traction_op.transpose() * u_local;
                    const T    sigma_n  = traction.dot(n);
                    tensor_vec sigma_t  = traction - sigma_n * n;

                    tensor_vec u_face = trace_op.transpose() * u_local;
                    const T    u_n    = u_face.dot(n);
                    tensor_vec u_t    = u_face - u_n * n;

                    T bz = T(0);
                    if constexpr (Mesh::dimension == 3) bz = barycentre.z();
                    out.push_back({fc, sigma_n, sigma_t.norm(),
                                   P_n, P_t.norm(),
                                   u_n, u_t(0), u_t.norm(),
                                   barycentre.x(), barycentre.y(), bz});
                }
                face_offset += face_dofs;
            }
            ci++;
        }
        return out;
    }

    // ---- bulk field -----------------------------------------------------
    void write_bulk(const std::string& base, const hho_state<T>& st,
                    double warp = 0.0) const
    {
        disk::write_vtu(base + ".vtu", op.msh, st.cells, op.cell_deg, warp);
    }

    // ---- CSV: all barycentre quantities, one row per contact face -------
    // friction < 0  -> Tresca (coeff is the constant threshold s)
    // friction >= 0 -> Coulomb (coeff is F, per-face s = -F [P_n]_-)
    void write_csv(const std::string& path, const hho_state<T>& st,
                   T coeff, bool coulomb) const
    {
        auto q = evaluate(st);

        // sort by barycentre so the file is a clean profile along Gamma_C
        std::sort(q.begin(), q.end(), [](const FaceQ& a, const FaceQ& b) {
            if (a.bx != b.bx) return a.bx < b.bx;
            if (a.by != b.by) return a.by < b.by;
            return a.bz < b.bz;
        });

        std::ofstream ofs(path);
        ofs << std::scientific << std::setprecision(10);
        ofs << "bx,by,bz,sigma_n,sigma_t_norm,P_n,P_t_norm,"
               "u_n,u_t0,u_t_norm,s,slip_ratio,zone\n";

        for (const auto& f : q)
        {
            // per-face friction bound and zone classification
            T s, slip_ratio; int zone;
            if (coulomb)
            {
                const T Pn_neg = std::min(f.P_n, T(0));
                s = -coeff * Pn_neg;
                if (f.P_n > T(0))          { zone = 0; slip_ratio = T(0); }   // separation
                else {
                    slip_ratio = (s > T(1e-30)) ? f.P_t_norm / s : T(0);
                    zone = (slip_ratio > T(1)) ? 1 : 2;                        // slip : stick
                }
            }
            else   // Tresca: constant threshold, symmetry_condition (no separation)
            {
                s = coeff;
                slip_ratio = (s > T(1e-30)) ? f.P_t_norm / s : T(0);
                zone = (slip_ratio >= T(1)) ? 1 : 2;
            }

            ofs << f.bx << ',' << f.by << ',' << f.bz << ','
                << f.sigma_n << ',' << f.sigma_t_norm << ','
                << f.P_n << ',' << f.P_t_norm << ','
                << f.u_n << ',' << f.u_t0 << ',' << f.u_t_norm << ','
                << s << ',' << slip_ratio << ',' << zone << '\n';
        }
        std::cout << "  wrote " << path << " (" << q.size()
                  << " contact faces)\n";
    }

    // ---- Tresca: ratios against the constant threshold s ---------------
    void write_tresca(const std::string& base, const hho_state<T>& st,
                      T s) const
    {
        const auto q = evaluate(st);
        const size_t n = q.size();

        std::vector<typename Mesh::face_type> faces(n);
        std::vector<double> stress_ratio(n), traction_ratio(n),
                            augmented_ratio(n), normal_stress(n),
                            normal_displacement(n), tangential_displacement(n);

        for (size_t i = 0; i < n; i++)
        {
            faces[i]                   = q[i].fc;
            stress_ratio[i]            = q[i].sigma_t_norm / s;
            traction_ratio[i]          = std::min(q[i].P_t_norm, s) / s;
            augmented_ratio[i]         = q[i].P_t_norm / s;
            normal_stress[i]           = q[i].sigma_n;
            normal_displacement[i]     = q[i].u_n;
            tangential_displacement[i] = q[i].u_t0;
        }

        disk::write_contact_vtu(base + "_contact.vtu", op.msh, faces,
            {{"stress_ratio",            stress_ratio},
             {"traction_ratio",          traction_ratio},
             {"augmented_ratio",         augmented_ratio},
             {"normal_stress",           normal_stress},
             {"normal_displacement",     normal_displacement},
             {"tangential_displacement", tangential_displacement}});
    }

    // ---- Coulomb: per-face threshold s(u) = -F [P_n]_-  + zoning --------
    void write_coulomb(const std::string& base, const hho_state<T>& st,
                       T F) const
    {
        const auto q = evaluate(st);
        const size_t n = q.size();

        std::vector<typename Mesh::face_type> faces(n);
        std::vector<double> zone(n),               // 0 separation, 1 slip, 2 stick
                            friction_bound(n),     // s = -F [P_n]_-
                            slip_ratio(n),         // ||P_t|| / s
                            normal_stress(n),
                            contact_pressure(n),   // [P_n]_-
                            friction_traction(n),  // ||[P_t]_s||
                            normal_displacement(n),
                            tangential_displacement(n);

        for (size_t i = 0; i < n; i++)
        {
            const T P_n_neg   = std::min(q[i].P_n, T(0));
            const T s_coulomb = -F * P_n_neg;

            T face_zone, face_ratio;
            if (q[i].P_n > T(0))                       // no contact
            {
                face_zone  = 0.0;
                face_ratio = 0.0;
            }
            else
            {
                face_ratio = (s_coulomb > T(1e-30))
                             ? q[i].P_t_norm / s_coulomb : 0.0;
                face_zone  = (face_ratio > 1.0) ? 1.0 : 2.0;
            }

            faces[i]                   = q[i].fc;
            zone[i]                    = face_zone;
            friction_bound[i]          = s_coulomb;
            slip_ratio[i]              = face_ratio;
            normal_stress[i]           = q[i].sigma_n;
            contact_pressure[i]        = P_n_neg;
            friction_traction[i]       = std::min(q[i].P_t_norm, s_coulomb);
            normal_displacement[i]     = q[i].u_n;
            tangential_displacement[i] = q[i].u_t0;
        }

        disk::write_contact_vtu(base + "_contact.vtu", op.msh, faces,
            {{"zone",                    zone},
             {"friction_bound",          friction_bound},
             {"slip_ratio",              slip_ratio},
             {"normal_stress",           normal_stress},
             {"contact_pressure",        contact_pressure},
             {"friction_traction",       friction_traction},
             {"normal_displacement",     normal_displacement},
             {"tangential_displacement", tangential_displacement}});
    }
};

} // namespace hho_contact