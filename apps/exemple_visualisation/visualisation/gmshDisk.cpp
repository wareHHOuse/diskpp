/*
 *       /\
 *      /__\       Nicolas Pignet (C) 2017
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    Gmsh tools
 *  /_\/_\/_\/_\
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

 #include "gmshMesh.h"
 #include "gmshData.h"
 #include "gmshElement.h"
 #include "gmshDisk.h"
 #include<vector>
 #include<string>
 #include<utility>
 #include<fstream>
 #include<iostream>

 namespace visu{
 
// Edge convertEdge( const auto& edge, const size_t index)
// {
//    auto nodes = edge.point_ids();
//    
//    assert(nodes.size() == 2);
//    
//    std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1 };
//       
//    Edge tmpedge(tmpnodes, index, 0, 0);
//       
//    return tmpedge;
// }
// 
// 
// Triangle convertTriangle( const auto& triangle, const size_t index)
// {
//    auto nodes = triangle.point_ids();
//    
//    assert(nodes.size() == 3);
//    
//    std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1 };
//       
//    Triangle tri(tmpnodes, index, 0, 0);
//       
//    return tri;
// }
// 
// Quadrangle convertQuadrangle( const auto& quadrangle, const size_t index)
// {
//    auto nodes = quadrangle.point_ids();
//    
//    assert(nodes.size() == 4);
//    
//    std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1, nodes[3] + 1 };
//       
//    Quadrangle quad(tmpnodes, index, 0, 0);
//       
//    return quad;
// }
// 
// Tetrahedron convertTetrahedron( const auto& tetrahedron, const size_t index)
// {
//    auto nodes = tetrahedron.point_ids();
//    
//    assert(nodes.size() == 4);
//    
//    std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1, nodes[3] + 1 };
//       
//    Tetrahedron tetra(tmpnodes, index, 0, 0);
//       
//    return tetra;
// }
// 
// Hexahedron convertHexahedron( const auto& hexahedron, const size_t index)
// {
//    auto nodes = hexahedron.point_ids();
//    
//    assert(nodes.size() == 8);
//    
//    std::vector<size_t> tmpnodes = { nodes[0] + 1, nodes[1] + 1, nodes[2] + 1, nodes[3] + 1,
//                                     nodes[4] + 1, nodes[5] + 1, nodes[6] + 1, nodes[7] + 1, };
//       
//    Hexahedron hexa(tmpnodes, index, 0, 0);
//       
//    return hexa;
// }
    
   

} //visu