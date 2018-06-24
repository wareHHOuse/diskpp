/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <algorithm>
#include <array>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>
#include <string>
#include <unistd.h>
#include <vector>

class Nonconforming2D
{
 private:
   std::vector<std::array<double, 3>> m_vertices;
   std::vector<std::array<size_t, 3>> m_edges;
   std::vector<std::array<size_t, 4>> m_triangles;
   std::vector<std::array<size_t, 5>> m_quads;
   std::vector<std::array<size_t, 7>> m_hexagons;
   std::vector<bool>                  m_notfuse;

   std::pair<bool, size_t>
   find(const std::array<size_t, 5>& quad, const size_t vert)
   {
      for (size_t i = 0; i < 4; i++) {
         if (quad[i] == vert) return std::make_pair(true, i);
      }

      return std::make_pair(false, 5);
   }

   std::array<size_t, 4>
   circularpermut(const std::array<size_t, 5>& quad, const size_t vert)
   {
      std::array<size_t, 4> permut = {0, 0, 0, 0};
      if (quad[0] == vert) {
         permut[0] = quad[0];
         permut[1] = quad[1];
         permut[2] = quad[2];
         permut[3] = quad[3];

         return permut;
      } else if (quad[1] == vert) {
         permut[0] = quad[1];
         permut[1] = quad[2];
         permut[2] = quad[3];
         permut[3] = quad[0];

         return permut;
      } else if (quad[2] == vert) {
         permut[0] = quad[2];
         permut[1] = quad[3];
         permut[2] = quad[0];
         permut[3] = quad[1];

         return permut;
      } else if (quad[3] == vert) {
         permut[0] = quad[3];
         permut[1] = quad[0];
         permut[2] = quad[1];
         permut[3] = quad[2];

         return permut;
      } else
         return permut;
   }

 public:
   Nonconforming2D() = default;

   void
   readMeditmesh(const std::string& filename)
   {
      std::cout << "Guessed mesh format: MEDIT 2D" << std::endl;

      std::ifstream ifs(filename);
      std::string   keyword;

      if (!ifs.is_open()) {
         std::cout << "Error opening " << filename << std::endl;
      }

      ifs >> keyword;
      if (keyword != "MeshVersionFormatted") {
         std::invalid_argument("Expected keyword \"MeshVersionFormatted\"");
      }

      size_t format;
      ifs >> format;

      if (format != 2) {
         std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
         std::invalid_argument("Expected format");
      }

      ifs >> keyword;
      if (keyword != "Dimension") {
         std::invalid_argument("Expected keyword \"Dimension\"");
      }

      size_t dim;
      ifs >> dim;

      if (dim != 3) {
         std::cout << "Expected dimension = 3 (here: " << dim << ")" << std::endl;
         std::invalid_argument("Wrong dimension");
      }

      size_t elements_to_read = 0;

      ifs >> keyword;
      while (keyword != "End") {
         std::cout << keyword << std::endl;
         if (keyword == "Vertices") {
            ifs >> elements_to_read;
            m_vertices.reserve(elements_to_read);

            std::array<double, 3> vertice = {0.0, 0.0, 0.0};

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> vertice[0] >> vertice[1] >> vertice[2] >> keyword;
               m_vertices.push_back(vertice);
            }

         } else if (keyword == "Triangles") {
            ifs >> elements_to_read;
            m_triangles.reserve(elements_to_read);

            std::array<size_t, 4> tri = {0, 0, 0, 0};

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> tri[0] >> tri[1] >> tri[2] >> tri[3];
               m_triangles.push_back(tri);
            }
         } else if (keyword == "Quadrilaterals") {
            ifs >> elements_to_read;
            m_quads.reserve(elements_to_read);

            std::array<size_t, 5> quad = {0, 0, 0, 0, 0};

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> quad[0] >> quad[1] >> quad[2] >> quad[3] >> quad[4];
               m_quads.push_back(quad);
            }
         } else if (keyword == "Edges") {
            ifs >> elements_to_read;
            m_edges.reserve(elements_to_read);

            std::array<size_t, 3> edge = {0, 0, 0};

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> edge[0] >> edge[1] >> edge[2];
               m_edges.push_back(edge);
            }
         } else if (keyword == "Tetrahedra") {
            throw std::invalid_argument("mesh is not tetrahedral");
         } else if (keyword == "Hexahedra") {
            throw std::invalid_argument("mesh is not hexahedral");
         } else {
            std::invalid_argument("Error parsing Medit file");
         }

         ifs >> keyword;
      }

      ifs.close();
   }

   void
   readDiskmesh(const std::string& filename)
   {
      std::cout << "Guessed mesh format: DISK 2D" << std::endl;

      std::ifstream ifs(filename);
      std::string   keyword;

      if (!ifs.is_open()) {
         std::cout << "Error opening " << filename << std::endl;
      }

      size_t elements_to_read = 0;

      ifs >> elements_to_read;
      m_vertices.reserve(elements_to_read);

      std::array<double, 3> vertice = {0.0, 0.0, 0.0};

      for (size_t i = 0; i < elements_to_read; i++) {
         ifs >> vertice[0] >> vertice[1];
         m_vertices.push_back(vertice);
      }

      ifs >> elements_to_read;
      m_quads.reserve(elements_to_read);

      std::array<size_t, 5> quad = {0, 0, 0, 0, 0};

      for (size_t i = 0; i < elements_to_read; i++) {
         // Warning: We need to invert vertice 2 and 3
         ifs >> quad[0] >> quad[1] >> quad[3] >> quad[2];
         quad[0] += 1;
         quad[1] += 1;
         quad[2] += 1;
         quad[3] += 1;
         m_quads.push_back(quad);
      }

      ifs >> elements_to_read;
      m_edges.reserve(elements_to_read);

      std::array<size_t, 3> edge = {0, 0, 0};

      for (size_t i = 0; i < elements_to_read; i++) {
         ifs >> edge[0] >> edge[1];
         edge[0] += 1;
         edge[1] += 1;
         m_edges.push_back(edge);
      }

      ifs.close();
   }

   void
   writeMeditmesh(const std::string& filename)
   {
      std::ofstream mesh_file(filename);
      if (!mesh_file.is_open()) {
         std::cout << "Unable to open file " << filename << std::endl;
         abort();
      }

      mesh_file << "MeshVersionFormatted 2" << std::endl;
      mesh_file << "Dimension 3" << std::endl;
      mesh_file << " " << std::endl;

      mesh_file << "Vertices" << std::endl;
      mesh_file << m_vertices.size() << std::endl;

      for (auto& vert : m_vertices) {
         mesh_file << vert[0] << " " << vert[1] << " " << vert[2] << " " << 1 << std::endl;
      }

      mesh_file << " " << std::endl;
      mesh_file << "Edges" << std::endl;
      mesh_file << m_edges.size() << std::endl;

      for (auto& edge : m_edges) {
         mesh_file << edge[0] << " " << edge[1] << " " << edge[2] << std::endl;
      }

      mesh_file << " " << std::endl;
      mesh_file << "Triangles" << std::endl;
      mesh_file << m_triangles.size() << std::endl;

      for (auto& tri : m_triangles) {
         mesh_file << tri[0] << " " << tri[1] << " " << tri[2] << " " << tri[3] << std::endl;
      }

      mesh_file << " " << std::endl;
      mesh_file << "Quadrilaterals" << std::endl;
      mesh_file << m_quads.size() << std::endl;

      for (auto& quad : m_quads) {
         mesh_file << quad[0] << " " << quad[1] << " " << quad[2] << " " << quad[3] << " "
                   << quad[4] << std::endl;
      }

      mesh_file << " " << std::endl;
      mesh_file << "Hexagons" << std::endl;
      mesh_file << m_hexagons.size() << std::endl;

      for (auto& hexa : m_hexagons) {
         mesh_file << hexa[0] << " " << hexa[1] << " " << hexa[2] << " " << hexa[3] << " "
                   << hexa[4] << " " << hexa[5] << " " << hexa[6] << std::endl;
      }

      mesh_file << " " << std::endl;
      mesh_file << "End" << std::endl;

      mesh_file.close();
   }

   void
   modifyMesh(const size_t& nb_hexa)
   {
      if (nb_hexa > 0) {
         std::random_device rd;        // Will be used to obtain a seed for the random number engine
         std::mt19937       gen(rd()); // Standard mersenne_twister_engine seeded with rd()
         std::uniform_int_distribution<> dis(0, 3);

         m_notfuse.assign(m_quads.size(), true);

         std::vector<size_t> quad_to_fuse;

         const size_t incr = m_quads.size() / nb_hexa;
         for (size_t i = 0; i < m_quads.size(); i += incr)
            quad_to_fuse.push_back(i);

         for (size_t num_cell : quad_to_fuse) {
            if (m_notfuse[num_cell]) {
               const std::array<size_t, 5> quadf = m_quads[num_cell];
               const size_t                ind0  = dis(gen);
               const size_t                ind1  = (ind0 + 1) % 4;
               const size_t                ind2  = (ind0 + 2) % 4;
               const size_t                ind3  = (ind0 + 3) % 4;

               const size_t p0 = quadf[ind0];
               const size_t p1 = quadf[ind1];
               const size_t p2 = quadf[ind2];
               const size_t p3 = quadf[ind3];

               // check if is not a boundary edge
               bool bnd = false;
               for (auto& edge : m_edges) {
                  if (p0 == edge[0] && p1 == edge[1]) {
                     bnd = true;
                     break;
                  } else if (p0 == edge[1] && p1 == edge[0]) {
                     bnd = true;
                     break;
                  }
               }

               if (!bnd) {
                  for (size_t i = 0; i < m_quads.size(); i++) {
                     if (m_notfuse[i] && i != num_cell) {
                        std::array<size_t, 5>         quad    = m_quads[i];
                        const std::pair<bool, size_t> find_p0 = find(quad, p0);
                        if (find_p0.first) {
                           const std::pair<bool, size_t> find_p1 = find(quad, p1);
                           if (find_p1.first) {
                              // we have found the neighbourd
                              m_notfuse[i]                       = false;
                              m_notfuse[num_cell]                = false;
                              const std::array<size_t, 4> permut = circularpermut(quad, p0);

                              const std::array<size_t, 7> hexagon = {
                                p0, permut[1], permut[2], p1, p2, p3, 1};
                              m_hexagons.push_back(hexagon);
                           }
                        }
                     }
                  }
               }
            }
         }

         std::cout << "We have created " << m_hexagons.size() << " hexagons" << std::endl;

         // remove old quad
         std::vector<std::array<size_t, 5>> new_quads;
         new_quads.reserve(m_quads.size());
         size_t nb_quad = 0;

         for (size_t i = 0; i < m_quads.size(); i++) {
            if (m_notfuse[i]) {
               nb_quad++;
               new_quads.push_back(m_quads[i]);
            }
         }

         m_quads.clear();
         m_quads.reserve(nb_quad);
         for (size_t i = 0; i < nb_quad; i++) {
            m_quads.push_back(new_quads[i]);
         }
      }
   }

   void
   random(const double h)
   {
      std::vector<bool> point_bnd;
      point_bnd.assign(m_vertices.size(), false);

      for (auto& edge : m_edges) {
         point_bnd[edge[0] - 1] = true;
         point_bnd[edge[1] - 1] = true;
      }

      std::random_device rd;        // Will be used to obtain a seed for the random number engine
      std::mt19937       gen(rd()); // Standard mersenne_twister_engine seeded with rd()
      std::uniform_real_distribution<> dis(-h, h);

      for (size_t i = 0; i < m_vertices.size(); i++) {
         if (!point_bnd[i]) {
            m_vertices[i][0] += dis(gen);
            m_vertices[i][1] += dis(gen);
         }
      }
   }
};

namespace priv {
class Edge
{
 public:
   Edge() {}
   Edge(size_t ap0, size_t ap1)
   {
      if (ap0 < ap1) {
         p0 = ap0;
         p1 = ap1;
      } else {
         p0 = ap1;
         p1 = ap0;
      }
   }

   size_t              p0, p1;
   std::vector<size_t> faces;

   inline friend bool
   operator==(const Edge& lhs, const Edge& rhs)
   {
      if (lhs.p0 == rhs.p0 && lhs.p1 == rhs.p1) {
         return true;
      }

      return false;
   }
};

class Face
{
 public:
   std::vector<size_t> vertices;
   std::vector<size_t> edges;
   size_t              id;
   int                 v0, v1;

   Face() : id(0), v0(-1), v1(-1) {}
   Face(const std::vector<size_t>& vert, const size_t aid) : id(aid), v0(-1), v1(-1)
   {
      vertices.clear();
      vertices = vert;
      std::vector<size_t>::iterator min =
        std::min_element(std::begin(vertices), std::end(vertices));
      std::rotate(vertices.begin(), min, vertices.end());
   }

   void
   sort()
   {
      std::vector<size_t>::iterator min = std::min_element(std::begin(edges), std::end(edges));
      std::rotate(edges.begin(), min, edges.end());
   }

   inline friend bool
   operator==(const Face& lhs, const Face& rhs)
   {
      if (lhs.vertices.size() == rhs.vertices.size()) {
         for (size_t i = 0; i < lhs.vertices.size(); i++) {
            if (std::find(lhs.vertices.begin(), lhs.vertices.end(), rhs.vertices[i]) ==
                lhs.vertices.end()) {
               return false;
            }
         }

         return true;
      }

      return false;
   }
};

class Volume
{
 public:
   std::vector<size_t> vertices;
   std::vector<size_t> faces;

   Volume() {}
   Volume(const std::vector<size_t>& vert)
   {
      vertices.clear();
      vertices = vert;
      std::vector<size_t>::iterator min =
        std::min_element(std::begin(vertices), std::end(vertices));
      std::rotate(vertices.begin(), min, vertices.end());
   }

   void
   sort()
   {
      std::vector<size_t>::iterator min = std::min_element(std::begin(faces), std::end(faces));
      std::rotate(faces.begin(), min, faces.end());
   }

   inline friend bool
   operator==(const Volume& lhs, const Volume& rhs)
   {
      if (lhs.vertices.size() == rhs.vertices.size()) {
         for (size_t i = 0; i < lhs.vertices.size(); i++) {
            if (lhs.vertices[i] != rhs.vertices[i]) {
               return false;
            }
         }

         return true;
      }

      return false;
   }
};

} // end priv

class Nonconforming3D
{
   std::vector<std::array<double, 3>> m_vertices;
   std::vector<priv::Edge>            m_edges;
   std::vector<priv::Face>            m_faces;
   std::vector<priv::Volume>          m_tetra;
   std::vector<priv::Volume>          m_hexa;
   std::vector<priv::Volume>          m_poly;
   std::vector<bool>                  m_notfuse;
   std::vector<size_t>                m_bnd_faces;

   size_t
   find(priv::Edge& ed)
   {
      for (size_t i = 0; i < m_edges.size(); i++) {
         if (ed == m_edges[i]) {
            m_edges[i].faces.push_back(ed.faces[0]);
            return i;
         }
      }

      m_edges.push_back(ed);

      return m_edges.size() - 1;
   }

   size_t
   find(priv::Face& face)
   {
      for (size_t i = 0; i < m_faces.size(); i++) {
         if (face == m_faces[i]) {
            if (m_faces[i].v0 == -1) {
               m_faces[i].v0 = face.v0;
            } else {
               m_faces[i].v1 = face.v0;
            }
            return i;
         }
      }

      // find edge
      size_t num_vert = face.vertices.size();
      size_t face_id  = m_faces.size();
      for (size_t i = 0; i < num_vert; i++) {
         priv::Edge e0(face.vertices[i], face.vertices[(i + 1) % num_vert]);
         e0.faces.push_back(face_id);
         face.edges.push_back(find(e0));
      }

      m_faces.push_back(face);

      return face_id;
   }

   bool
   is_inside(const priv::Face fc, const size_t pt)
   {
      for (auto& vert : fc.vertices) {
         if (vert == pt) {
            return true;
         }
      }

      return false;
   }

   std::pair<bool, std::pair<std::array<size_t, 2>, size_t>>
   has_bnd_face(const priv::Face fc)
   {
      std::pair<bool, std::pair<std::array<size_t, 2>, size_t>> ret;
      size_t                                                    face_bnd = 0;
      std::array<size_t, 2>                                     faces({0, 0});

      for (auto& edge_id : fc.edges) {
         for (auto& face_id : m_edges[edge_id].faces) {
            const priv::Face face = m_faces[face_id];
            if (face.v0 >= 0 && face.v1 == -1) {
               faces[face_bnd++] = face_id;
            }
         }
         if (face_bnd == 2) {
            ret.first  = true;
            ret.second = std::make_pair(faces, edge_id);
            return ret;
         } else {
            face_bnd = 0;
         }
      }
      ret.first  = false;
      ret.second = std::make_pair(faces, 0);

      return ret;
   }

   size_t
   fuse_face(const std::pair<std::array<size_t, 2>, size_t> to_fuse)
   {
      const auto edge = m_edges[to_fuse.second];
      const auto f0   = m_faces[to_fuse.first[0]];
      const auto f1   = m_faces[to_fuse.first[1]];

      auto f0_pts = f0.vertices;
      auto f1_pts = f1.vertices;
      // std::cout << edge.p0 << "-" << edge.p1 << std::endl;
      // std::cout << f0_pts[0] << "-" << f0_pts[1] << "-" << f0_pts[2] << "-" << f0_pts[3]
      //           << std::endl;
      // std::cout << f1_pts[0] << "-" << f1_pts[1] << "-" << f1_pts[2] << "-" << f1_pts[3]
      //           << std::endl;

      const auto p0_it = std::find(f0_pts.begin(), f0_pts.end(), edge.p0);
      std::rotate(f0_pts.begin(), p0_it, f0_pts.end());

      if (f0_pts[1] == edge.p1) {
         const auto p1_it = std::find(f0_pts.begin(), f0_pts.end(), edge.p1);
         std::rotate(f0_pts.begin(), p1_it, f0_pts.end());
      }

      const auto it = std::find(f1_pts.begin(), f1_pts.end(), f0_pts[0]);
      std::rotate(f1_pts.begin(), it, f1_pts.end());

      // std::cout << f0_pts[0] << "-" << f0_pts[1] << "-" << f0_pts[2] << "-" << f0_pts[3]
      //           << std::endl;
      // std::cout << f1_pts[0] << "-" << f1_pts[1] << "-" << f1_pts[2] << "-" << f1_pts[3]
      //           << std::endl;

      std::vector<size_t> pts = f0_pts;

      for (size_t i = 2; i < f1_pts.size(); i++) {
         pts.push_back(f1_pts[i]);
      }

      //    std::cout << pts[0] << "-" << pts[1] << "-" << pts[2] << "-" << pts[3] << "-" << pts[4]
      //    << "-"
      // << pts[5] << std::endl;

      priv::Face new_fc(pts, f0.id);

      for (auto& e : f0.edges) {
         if (e != to_fuse.second) new_fc.edges.push_back(e);
      }

      for (auto& e : f1.edges) {
         if (e != to_fuse.second) new_fc.edges.push_back(e);
      }

      new_fc.v1      = -1;
      size_t face_id = m_faces.size();
      m_faces.push_back(new_fc);

      return face_id;
   }

 public:
   Nonconforming3D() = default;

   bool
   readMeditmesh(const std::string& filename)
   {
      std::cout << "Guessed mesh format: MEDIT 3D" << std::endl;

      std::ifstream ifs(filename);
      std::string   keyword;

      if (!ifs.is_open()) {
         std::cout << "Error opening " << filename << std::endl;
      }

      ifs >> keyword;
      if (keyword != "MeshVersionFormatted") {
         std::invalid_argument("Expected keyword \"MeshVersionFormatted\"");
      }

      size_t format;
      ifs >> format;

      if (format != 2) {
         std::cout << "Expected format 2 (here: " << format << ")" << std::endl;
         std::invalid_argument("Expected format");
      }

      ifs >> keyword;
      if (keyword != "Dimension") {
         std::invalid_argument("Expected keyword \"Dimension\"");
      }

      size_t dim;
      ifs >> dim;

      if (dim != 3) {
         std::cout << "Expected dimension = 3 (here: " << dim << ")" << std::endl;
         std::invalid_argument("Wrong dimension");
      }

      size_t elements_to_read = 0;

      size_t num_vol = 0;

      ifs >> keyword;
      while (keyword != "End") {
         std::cout << keyword << std::endl;
         if (keyword == "Vertices") {
            ifs >> elements_to_read;
            m_vertices.reserve(elements_to_read);

            std::array<double, 3> vertice = {0.0, 0.0, 0.0};

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> vertice[0] >> vertice[1] >> vertice[2] >> keyword;
               m_vertices.push_back(vertice);
            }

         } else if (keyword == "Edges") {
            ifs >> elements_to_read;

            size_t p0, p1, id;

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> p0 >> p1 >> id;
            }
         } else if (keyword == "Triangles") {
            ifs >> elements_to_read;

            size_t p0, p1, p2, id;

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> p0 >> p1 >> p2 >> id;
               priv::Face tri(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1}), id);
               tri.v0 = -1;
               m_bnd_faces.push_back(find(tri));
            }
         } else if (keyword == "Quadrilaterals") {
            ifs >> elements_to_read;

            size_t p0, p1, p2, p3, id;

            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> p0 >> p1 >> p2 >> p3 >> id;
               priv::Face q(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1}), id);
               q.v0 = -1;
               m_bnd_faces.push_back(find(q));
            }
         } else if (keyword == "Tetrahedra") {
            ifs >> elements_to_read;

            size_t p0, p1, p2, p3, id;
            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> p0 >> p1 >> p2 >> p3 >> id;
               priv::Volume tet(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1}));
               size_t       tet_id = m_tetra.size();

               priv::Face tri1(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1}), 0);
               tri1.v0 = tet_id;
               tet.faces.push_back(find(tri1));

               priv::Face tri2(std::vector<size_t>({p0 - 1, p3 - 1, p2 - 1}), 0);
               tri2.v0 = tet_id;
               tet.faces.push_back(find(tri2));

               priv::Face tri3(std::vector<size_t>({p1 - 1, p2 - 1, p3 - 1}), 0);
               tri3.v0 = tet_id;
               tet.faces.push_back(find(tri3));

               priv::Face tri4(std::vector<size_t>({p0 - 1, p3 - 1, p1 - 1}), 0);
               tri4.v0 = tet_id;
               tet.faces.push_back(find(tri4));

               m_tetra.push_back(tet);
            }
         } else if (keyword == "Hexahedra") {
            ifs >> elements_to_read;

            size_t p0, p1, p2, p3, p4, p5, p6, p7, id;
            for (size_t i = 0; i < elements_to_read; i++) {
               ifs >> p0 >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> id;
               priv::Volume hex(std::vector<size_t>(
                 {p0 - 1, p1 - 1, p2 - 1, p3 - 1, p4 - 1, p5 - 1, p6 - 1, p7 - 1}));
               size_t       hexid = m_hexa.size();

               priv::Face q1(std::vector<size_t>({p0 - 1, p1 - 1, p2 - 1, p3 - 1}), 0);
               q1.v0 = hexid;
               hex.faces.push_back(find(q1));

               priv::Face q2(std::vector<size_t>({p0 - 1, p4 - 1, p7 - 1, p3 - 1}), 0);
               q2.v0 = hexid;
               hex.faces.push_back(find(q2));

               priv::Face q3(std::vector<size_t>({p4 - 1, p5 - 1, p6 - 1, p7 - 1}), 0);
               q3.v0 = hexid;
               hex.faces.push_back(find(q3));

               priv::Face q4(std::vector<size_t>({p5 - 1, p1 - 1, p2 - 1, p6 - 1}), 0);
               q4.v0 = hexid;
               hex.faces.push_back(find(q4));

               priv::Face q5(std::vector<size_t>({p3 - 1, p7 - 1, p6 - 1, p2 - 1}), 0);
               q5.v0 = hexid;
               hex.faces.push_back(find(q5));

               priv::Face q6(std::vector<size_t>({p0 - 1, p4 - 1, p5 - 1, p1 - 1}), 0);
               q6.v0 = hexid;
               hex.faces.push_back(find(q6));

               m_hexa.push_back(hex);
            }
         } else {
            std::invalid_argument("Error parsing Medit file");
         }

         ifs >> keyword;
      }
      ifs.close();

      // order element
      for (auto& fc : m_faces) {
         fc.sort();
      }

      for (auto& tet : m_tetra) {
         tet.sort();
      }

      for (auto& hex : m_hexa) {
         hex.sort();
      }

      return true;
   }

   void
   writeMeditmesh(const std::string& filename)
   {
      std::ofstream mesh_file(filename);
      if (!mesh_file.is_open()) {
         std::cout << "Unable to open file " << filename << std::endl;
         abort();
      }

      mesh_file << "MeshVersionFormatted 2" << std::endl;
      mesh_file << "Dimension 3" << std::endl;
      mesh_file << " " << std::endl;

      mesh_file << "Vertices"
                << "\t" << m_vertices.size() << std::endl;

      for (auto& vert : m_vertices) {
         mesh_file << vert[0] << " " << vert[1] << " " << vert[2] << std::endl;
      }

      mesh_file << "Volumes->Faces"
                << "\t" << m_tetra.size() + m_hexa.size() + m_poly.size() << std::endl;

      for (auto& tet : m_tetra) {
         mesh_file << 4 << " " << tet.faces[0] << " " << tet.faces[1] << " " << tet.faces[2] << " "
                   << tet.faces[3] << std::endl;
      }

      for (auto& hex : m_hexa) {
         mesh_file << 6 << " " << hex.faces[0] << " " << hex.faces[1] << " " << hex.faces[2] << " "
                   << hex.faces[3] << " " << hex.faces[4] << " " << hex.faces[5] << std::endl;
      }

      for (auto& vol : m_poly) {
         mesh_file << vol.faces.size();
         for (auto& fc : vol.faces) {
            mesh_file << " " << fc;
         }
         mesh_file << std::endl;
      }

      mesh_file << "Volumes->Vertices"
                << "\t" << m_tetra.size() + m_hexa.size() + m_poly.size() << std::endl;

      for (auto& tet : m_tetra) {
         mesh_file << 4 << " " << tet.vertices[0] << " " << tet.vertices[1] << " "
                   << tet.vertices[2] << " " << tet.vertices[3] << std::endl;
      }

      for (auto& hex : m_hexa) {
         mesh_file << 8 << " " << hex.vertices[0] << " " << hex.vertices[1] << " "
                   << hex.vertices[2] << " " << hex.vertices[3] << " " << hex.vertices[4] << " "
                   << hex.vertices[5] << " " << hex.vertices[6] << " " << hex.vertices[7]
                   << std::endl;
      }

      for (auto& vol : m_poly) {
         mesh_file << vol.vertices.size();
         for (auto& fc : vol.vertices) {
            mesh_file << " " << fc;
         }
         mesh_file << std::endl;
      }

      mesh_file << "Faces->Edges"
                << "\t" << m_faces.size() << std::endl;

      for (auto& fc : m_faces) {
         mesh_file << fc.edges.size();
         for (auto& e : fc.edges) {
            mesh_file << " " << e;
         }
         mesh_file << std::endl;
      }

      mesh_file << "Faces->Vertices"
                << "\t" << m_faces.size() << std::endl;

      for (auto& fc : m_faces) {
         mesh_file << fc.vertices.size();
         for (auto& e : fc.vertices) {
            mesh_file << " " << e;
         }
         mesh_file << std::endl;
      }

      mesh_file << "Faces->Control volumes"
                << "\t" << m_bnd_faces.size() << std::endl;

      for (size_t i = 0; i < m_bnd_faces.size(); i++) {
         mesh_file << m_bnd_faces[i] << " " << m_faces[m_bnd_faces[i]].id << std::endl;
      }

      mesh_file << "Edges"
                << "\t" << m_edges.size() << std::endl;

      for (auto& edge : m_edges) {
         mesh_file << edge.p0 << " " << edge.p1 << std::endl;
      }

      mesh_file << "End" << std::endl;

      mesh_file.close();
   }

   void
   modifyMesh(const size_t& nb_hexa)
   {
      if (nb_hexa > 0 && nb_hexa <= m_hexa.size() / 2) {
         std::random_device rd;        // Will be used to obtain a seed for the random number
         std::mt19937       gen(rd()); // Standard mersenne_twister_engine seeded with rd()
         std::uniform_int_distribution<> dis(0, 4);

         m_notfuse.assign(m_hexa.size(), true);
         std::vector<bool> save_edges(m_edges.size(), true);
         std::vector<bool> save_faces(m_faces.size(), true);

         std::vector<size_t> hexa_to_fuse;

         const size_t incr = m_hexa.size() / nb_hexa;
         for (size_t i = 0; i < m_hexa.size(); i += incr)
            hexa_to_fuse.push_back(i);

         for (size_t num_cell : hexa_to_fuse) {
            if (m_notfuse[num_cell]) {
               const auto   hexa   = m_hexa[num_cell];
               const size_t ind0   = dis(gen);
               const size_t f0_ind = hexa.faces[ind0];
               const auto   f0     = m_faces[f0_ind];

               // check if is not a boundary edge
               bool bnd = false;
               if (f0_ind < m_bnd_faces.size()) {
                  bnd = true;
               }

               if (!bnd) {
                  // find the other volume
                  int vf = f0.v0;
                  if (vf == num_cell) {
                     vf = f0.v1;
                  }

                  if (m_notfuse[vf]) {
                     m_notfuse[vf]       = false;
                     m_notfuse[num_cell] = false;
                     save_faces[f0_ind]  = false;

                     const auto hexa2 = m_hexa[vf];

                     // add points
                     priv::Volume vol(hexa.vertices);
                     for (auto vert : hexa2.vertices) {
                        if (!is_inside(f0, vert)) vol.vertices.push_back(vert);
                     }

                     // if one of the edges if on a boundary we supress it and fuse two surfaces
                     const auto bnd_face = has_bnd_face(f0);
                     if (bnd_face.first) {
                        const auto to_fuse           = bnd_face.second;
                        save_faces[to_fuse.first[0]] = false;
                        save_faces[to_fuse.first[1]] = false;
                        save_edges[to_fuse.second]   = false;

                        size_t face_id = fuse_face(to_fuse);

                        save_faces.push_back(true);
                        vol.faces.push_back(face_id);
                        m_bnd_faces.push_back(face_id);

                        // add faces
                        for (auto& fc : hexa.faces) {
                           if (fc != f0_ind && fc != to_fuse.first[0] && fc != to_fuse.first[1])
                              vol.faces.push_back(fc);
                        }

                        for (auto& fc : hexa2.faces) {
                           if (fc != f0_ind && fc != to_fuse.first[0] && fc != to_fuse.first[1])
                              vol.faces.push_back(fc);
                        }
                     } else {
                        // add faces

                        for (auto& fc : hexa.faces) {
                           if (fc != f0_ind) vol.faces.push_back(fc);
                        }

                        for (auto& fc : hexa2.faces) {
                           if (fc != f0_ind) vol.faces.push_back(fc);
                        }
                     }

                     m_poly.push_back(vol);
                  }
               }
            }
         }

         std::cout << "We have created " << m_poly.size() << " polyedres" << std::endl;

         // remove old edges
         std::vector<priv::Edge> new_edge;
         std::vector<int>        num_edge(m_edges.size(), -1);
         size_t                  nb_edge = 0;

         for (size_t i = 0; i < m_edges.size(); i++) {
            if (save_edges[i]) {
               num_edge[i] = nb_edge++;
               new_edge.push_back(m_edges[i]);
            }
         }

         m_edges.clear();
         m_edges.reserve(nb_edge);
         for (size_t i = 0; i < nb_edge; i++) {
            m_edges.push_back(new_edge[i]);
         }

         // remove old faces
         std::vector<priv::Face> new_face;
         std::vector<int>        num_face(m_faces.size(), -1);
         size_t                  nb_face = 0;

         for (size_t i = 0; i < m_faces.size(); i++) {
            if (save_faces[i]) {
               num_face[i] = nb_face++;
               new_face.push_back(m_faces[i]);
            }
         }

         m_faces.clear();
         m_faces.reserve(nb_face);
         for (size_t i = 0; i < nb_face; i++) {
            // new num
            for (size_t j = 0; j < new_face[i].edges.size(); j++) {
               new_face[i].edges[j] = num_edge[new_face[i].edges[j]];
            }
            new_face[i].sort();
            m_faces.push_back(new_face[i]);
         }

         // remove old bnd faces
         std::vector<size_t> new_bnd;
         size_t              nb_bnd = 0;

         for (size_t i = 0; i < m_bnd_faces.size(); i++) {
            if (save_faces[m_bnd_faces[i]]) {
               new_bnd.push_back(num_face[m_bnd_faces[i]]);
               nb_bnd++;
            }
         }

         m_bnd_faces.clear();
         m_bnd_faces.reserve(nb_bnd);
         for (size_t i = 0; i < nb_bnd; i++) {
            m_bnd_faces.push_back(new_bnd[i]);
         }

         // remove old hexa
         std::vector<priv::Volume> new_hexa;
         new_hexa.reserve(m_hexa.size());
         size_t nb_hexa = 0;

         for (size_t i = 0; i < m_hexa.size(); i++) {
            if (m_notfuse[i]) {
               nb_hexa++;
               new_hexa.push_back(m_hexa[i]);
            }
         }

         m_hexa.clear();
         m_hexa.reserve(nb_hexa);
         for (size_t i = 0; i < nb_hexa; i++) {
            for (size_t j = 0; j < new_hexa[i].faces.size(); j++) {
               new_hexa[i].faces[j] = num_face[new_hexa[i].faces[j]];
            }
            m_hexa.push_back(new_hexa[i]);
         }

         // update tetra et poly
         for (size_t i = 0; i < m_tetra.size(); i++) {
            for (size_t j = 0; j < m_tetra[i].faces.size(); j++) {
               m_tetra[i].faces[j] = num_face[m_tetra[i].faces[j]];
            }
         }

         for (size_t i = 0; i < m_poly.size(); i++) {
            for (size_t j = 0; j < m_poly[i].faces.size(); j++) {
               m_poly[i].faces[j] = num_face[m_poly[i].faces[j]];
            }
            m_poly[i].sort();
         }
      }
   }

   //  void
   //  random(const double h)
   //  {
   //     std::vector<bool> point_bnd;
   //     point_bnd.assign(m_vertices.size(), false);

   //     for (auto& edge : m_edges) {
   //        point_bnd[edge[0] - 1] = true;
   //        point_bnd[edge[1] - 1] = true;
   //     }

   //     std::random_device rd;        // Will be used to obtain a seed for the random number
   //     engine std::mt19937       gen(rd()); // Standard mersenne_twister_engine seeded with rd()
   //     std::uniform_real_distribution<> dis(-h, h);

   //     for (size_t i = 0; i < m_vertices.size(); i++) {
   //        if (!point_bnd[i]) {
   //           m_vertices[i][0] += dis(gen);
   //           m_vertices[i][1] += dis(gen);
   //        }
   //     }
   //  }
};

/* To use:
 * ./utils/createnonconforming.cpp -n number_hexagons input_mesh.*  output_mesh.medit2d
 */

int
main(int argc, char** argv)
{
   size_t nb_hexa = 0;
   bool   random  = false;
   double h       = 0;
   int    ch;
   while ((ch = getopt(argc, argv, "n:r:")) != -1) {
      switch (ch) {
         case 'n': nb_hexa = atoi(optarg); break;
         case 'r':
            h      = atof(optarg);
            random = true;
            break;
         default: std::cout << "wrong arguments" << std::endl; exit(1);
      }
   }

   argc -= optind;
   argv += optind;

   if (argc == 0) {
      std::cout << "Error" << std::endl;
      return 0;
   }

   char* filename = argv[0];
   char* output   = argv[1];

   if (std::regex_match(filename, std::regex(".*\\.medit2d$"))) {
      Nonconforming2D mesh;

      mesh.readMeditmesh(filename);
      mesh.modifyMesh(nb_hexa);
      if (random) mesh.random(h);
      mesh.writeMeditmesh(output);
   } else if (std::regex_match(filename, std::regex(".*\\.quad$"))) {
      Nonconforming2D mesh;

      mesh.readDiskmesh(filename);
      mesh.modifyMesh(nb_hexa);
      if (random) mesh.random(h);
      mesh.writeMeditmesh(output);
   } else if (std::regex_match(filename, std::regex(".*\\.mesh$"))) {
      Nonconforming3D mesh;

      mesh.readMeditmesh(filename);
      mesh.modifyMesh(nb_hexa);
      mesh.writeMeditmesh(output);
   } else {
      std::cout << "Unknown format" << std::endl;
   }
}

// lire maillage
// modifier maillage
// ecrire maillage