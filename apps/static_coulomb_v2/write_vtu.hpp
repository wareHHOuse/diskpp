/*
 * write_vtu.hpp -- VTK UnstructuredGrid output, any mesh type, 2D or 3D.
 *
 * Two writers:
 *   write_vtu          the bulk field, with NODAL displacement (vertex
 *                      averaged from the HHO cell polynomials) so that a
 *                      viewer can warp the mesh
 *   write_contact_vtu  a boundary as its own mesh, one value per face
 *
 * Cell types are chosen for generality rather than specificity:
 *   2D  ->  VTK_POLYGON    (7)   valid for triangles, quads, hexagons, ...
 *   3D  ->  VTK_POLYHEDRON (42)  valid for tets, hexes, Voronoi cells, ...
 * Both accept the connectivity DiSk++ provides directly, so no per-element
 * node-ordering table is needed.  VTK_POLYHEDRON additionally wants the
 * face-based description, which is exactly what a polytopal mesh stores:
 * faces(msh, cl), and fc.point_ids() for each face.
 */

#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <array>
#include <cmath>

#include "diskpp/bases/bases.hpp"

namespace disk
{

/* ==========================================================================
 * Bulk field with nodal displacement
 * ========================================================================*/
template<typename Mesh>
void
write_vtu(const std::string&  path,
          const Mesh&         msh,
          const std::vector<Matrix<typename Mesh::coordinate_type,
                                   Dynamic, 1>>&  u_cells,
          const size_t        cell_deg,
          const double        warp_factor = 0.0)
{
    using T = typename Mesh::coordinate_type;
    constexpr size_t DIM = Mesh::dimension;

    const size_t n_pts   = msh.points_size();
    const size_t n_cells = msh.cells_size();

    std::vector<point<T,DIM>> pt_of_id;
    pt_of_id.reserve(n_pts);
    for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        pt_of_id.push_back(*itor);

    // ---- vertex-averaged displacement -----------------------------------
    // The HHO cell unknown is a polynomial per cell, hence discontinuous.
    // Evaluate it at the vertices of its own cell and average over the
    // cells meeting at each vertex.
    std::vector<std::array<T,3>> u_node(n_pts, {0.0, 0.0, 0.0});
    std::vector<int> count(n_pts, 0);

    // ---- connectivity ---------------------------------------------------
    std::vector<std::vector<size_t>> conn(n_cells);
    std::vector<std::vector<std::vector<size_t>>> cell_faces;   // 3D only

    {
        size_t ci = 0;
        for (auto& cl : msh)
        {
            auto cb = make_vector_monomial_basis(msh, cl, cell_deg);

            // vertices of this cell.  A generic 3D cell stores its faces,
            // so the vertex set is gathered from them.
            std::vector<size_t> vids;
            if constexpr (DIM == 2)
            {
                for (auto pid : cl.point_ids()) vids.push_back(pid);
            }
            else
            {
                std::vector<std::vector<size_t>> face_list;
                std::set<size_t> vset;
                for (auto& fc : faces(msh, cl))
                {
                    std::vector<size_t> fpts;
                    for (auto pid : fc.point_ids())
                    {
                        fpts.push_back(pid);
                        vset.insert(pid);
                    }
                    face_list.push_back(fpts);
                }
                cell_faces.push_back(face_list);
                vids.assign(vset.begin(), vset.end());
            }

            conn[ci] = vids;

            for (auto pid : vids)
            {
                auto ub = eval(u_cells[ci], cb.eval_functions(pt_of_id[pid]));
                for (size_t d = 0; d < DIM; d++)
                    u_node[pid][d] += ub(d);
                count[pid]++;
            }
            ci++;
        }
    }

    for (size_t i = 0; i < n_pts; i++)
        if (count[i] > 0)
            for (size_t d = 0; d < 3; d++) u_node[i][d] /= count[i];

    // ---- write ----------------------------------------------------------
    std::ofstream ofs(path);
    ofs << std::scientific << std::setprecision(8);

    ofs << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << n_pts
        << "\" NumberOfCells=\"" << n_cells << "\">\n";

    ofs << "<Points>\n"
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < n_pts; i++)
    {
        const auto& p = pt_of_id[i];
        ofs << p.x() + warp_factor*u_node[i][0] << " "
            << p.y() + warp_factor*u_node[i][1] << " ";
        if constexpr (DIM == 3) ofs << p.z() + warp_factor*u_node[i][2] << "\n";
        else                    ofs << "0.0\n";
    }
    ofs << "</DataArray>\n</Points>\n";

    ofs << "<Cells>\n"
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto& c : conn) { for (auto id : c) ofs << id << " "; ofs << "\n"; }
    ofs << "</DataArray>\n"
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    {
        size_t off = 0;
        for (auto& c : conn) { off += c.size(); ofs << off << " "; }
    }
    ofs << "\n</DataArray>\n"
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < n_cells; i++)
        ofs << (DIM == 2 ? 7 : 42) << " ";       // POLYGON / POLYHEDRON
    ofs << "\n</DataArray>\n";

    if constexpr (DIM == 3)
    {
        // VTK_POLYHEDRON face stream, per cell:
        //   nFaces, nPts_f0, p..., nPts_f1, p..., ...
        ofs << "<DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\">\n";
        std::vector<size_t> faceoffsets;
        size_t running = 0;
        for (auto& fl : cell_faces)
        {
            ofs << fl.size() << " ";
            running += 1;
            for (auto& f : fl)
            {
                ofs << f.size() << " ";
                for (auto id : f) ofs << id << " ";
                running += 1 + f.size();
            }
            ofs << "\n";
            faceoffsets.push_back(running);
        }
        ofs << "</DataArray>\n"
            << "<DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\">\n";
        for (auto o : faceoffsets) ofs << o << " ";
        ofs << "\n</DataArray>\n";
    }

    ofs << "</Cells>\n";

    ofs << "<PointData Vectors=\"displacement\" Scalars=\"u_mag\">\n"
        << "<DataArray type=\"Float64\" Name=\"displacement\" "
        << "NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < n_pts; i++)
        ofs << u_node[i][0] << " " << u_node[i][1] << " " << u_node[i][2] << "\n";
    ofs << "</DataArray>\n"
        << "<DataArray type=\"Float64\" Name=\"u_mag\" format=\"ascii\">\n";
    for (size_t i = 0; i < n_pts; i++)
        ofs << std::sqrt(u_node[i][0]*u_node[i][0]
                       + u_node[i][1]*u_node[i][1]
                       + u_node[i][2]*u_node[i][2]) << "\n";
    ofs << "</DataArray>\n</PointData>\n";

    ofs << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";

    std::cout << "  wrote " << path << " (" << n_pts << " points, "
              << n_cells << " cells)\n";
}

/* ==========================================================================
 * A boundary as its own mesh, one value per face (CELL data)
 *
 * In 2D each face is a segment (VTK_LINE); in 3D a polygon (VTK_POLYGON).
 * Data stays piecewise constant per face, exactly as computed -- apply
 * "Cell Data to Point Data" in the viewer for a smooth curve or surface.
 *
 * Taking the faces themselves, rather than precomputed endpoints, keeps the
 * caller dimension-agnostic: it only collects the faces it visited.
 * ========================================================================*/
template<typename Mesh>
void
write_contact_vtu(const std::string& path,
                  const Mesh&        msh,
                  const std::vector<typename Mesh::face_type>& face_list,
                  const std::vector<std::pair<std::string,
                                              std::vector<double>>>& fields)
{
    using T = typename Mesh::coordinate_type;
    constexpr size_t DIM = Mesh::dimension;

    std::vector<point<T,DIM>> pt_of_id;
    pt_of_id.reserve(msh.points_size());
    for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        pt_of_id.push_back(*itor);

    // renumber: only the vertices actually used
    std::map<size_t, size_t> local_of_global;
    std::vector<size_t>      global_of_local;
    std::vector<std::vector<size_t>> conn;

    for (auto& fc : face_list)
    {
        std::vector<size_t> c;
        for (auto pid : fc.point_ids())
        {
            auto it = local_of_global.find(pid);
            if (it == local_of_global.end())
            {
                const size_t lid = global_of_local.size();
                local_of_global[pid] = lid;
                global_of_local.push_back(pid);
                c.push_back(lid);
            }
            else
                c.push_back(it->second);
        }
        conn.push_back(c);
    }

    const size_t n_pts  = global_of_local.size();
    const size_t n_cell = conn.size();

    std::ofstream ofs(path);
    ofs << std::scientific << std::setprecision(8);

    ofs << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << n_pts
        << "\" NumberOfCells=\"" << n_cell << "\">\n";

    ofs << "<Points>\n"
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < n_pts; i++)
    {
        const auto& p = pt_of_id[global_of_local[i]];
        ofs << p.x() << " " << p.y() << " ";
        if constexpr (DIM == 3) ofs << p.z() << "\n";
        else                    ofs << "0.0\n";
    }
    ofs << "</DataArray>\n</Points>\n";

    ofs << "<Cells>\n"
        << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto& c : conn) { for (auto id : c) ofs << id << " "; ofs << "\n"; }
    ofs << "</DataArray>\n"
        << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    {
        size_t off = 0;
        for (auto& c : conn) { off += c.size(); ofs << off << " "; }
    }
    ofs << "\n</DataArray>\n"
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < n_cell; i++)
        ofs << (DIM == 2 ? 3 : 7) << " ";         // LINE / POLYGON
    ofs << "\n</DataArray>\n</Cells>\n";

    ofs << "<CellData>\n";
    for (auto& [name, vals] : fields)
    {
        ofs << "<DataArray type=\"Float64\" Name=\"" << name
            << "\" format=\"ascii\">\n";
        for (auto v : vals) ofs << v << "\n";
        ofs << "</DataArray>\n";
    }
    ofs << "</CellData>\n";

    ofs << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";

    std::cout << "  wrote " << path << " (" << n_cell << " boundary faces)\n";
}

} // namespace disk