#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
#       /\        Matteo Cicuttin (C) 2016-2021
#      /__\       matteo.cicuttin@enpc.fr
#     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
#    /\    /\
#   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
#  /_\/_\/_\/_\   methods.
#
# This file is copyright of the following authors:
# Nicolas Pignet  (C) 2021                     nicolas.pignet@enpc.fr
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# If you use this code or parts of it for scientific publications, you
# are required to cite it as following:
#
# Hybrid High-Order methods for finite elastoplastic deformations
# within a logarithmic strain framework.
# M. Abbas, A. Ern, N. Pignet.
# International Journal of Numerical Methods in Engineering (2019)
# 120(3), 303-327
# DOI: 10.1002/nme.6137
#

import logging
import medcoupling
import argparse

class Edge:
    def __init__(self):
        self.nodes = []

class Face:
    def __init__(self):
        self.nodes = []
        self.edges = []

class Volume:
    def __init__(self):
        self.nodes = []
        self.edges = []
        self.faces = []


def rotate(liste):
    minpos = liste.index(min(liste))
    list_tmp = liste[minpos:] + liste[:minpos]
    size = len(list_tmp)
    if size > 2 and list_tmp[-1] < list_tmp[1]:
        list_tmp[1:] = list_tmp[::-1][0:size-1]

    return tuple(list_tmp)

class CellConverter:
    def __init__(self):
        _med_types = 'POINT1 SEG2 TRI3 QUAD4 TETRA4 HEXA8 PYRA5 PENTA6 SEG3 TRI6 QUAD8 TETRA10 HEXA20 PYRA13 PENTA15 SEG4 TRI7 QUAD9 PENTA18 HEXA27'.split()
        self.mdata = {}
        for i in range(len(_med_types)):
            self.mdata[getattr(medcoupling, 'NORM_%s'%_med_types[i])] = _med_types[i]
        self.mdim = {'POINT1' : 0, 'SEG2': 1, 'TRI3': 2, "QUAD4":2, "TETRA4":3, "HEXA8":3, "PYRA5":3, \
             "PENTA6":3, "SEG3":1, "TRI6":2, "QUAD8":2, "TETRA10":3, "HEXA20":3, "PYRA13":3, \
                  "PENTA15":3, "SEG4":1, "TRI7":2, "QUAD9":2, "PENTA18":3, "HEXA27":3 }
        self.mnodes = {'POINT1' : 1, 'SEG2': 2, 'TRI3': 3, "QUAD4":4, "TETRA4":4, "HEXA8":8, "PYRA5":5, \
             "PENTA6":6, "SEG3":2, "TRI6":3, "QUAD8":4, "TETRA10":4, "HEXA20":8, "PYRA13":5, \
                  "PENTA15":6, "SEG4":2, "TRI7":3, "QUAD9":4, "PENTA18":6, "HEXA27":8 }

        self.list_edges = {}
        self.list_faces = {}

        self.mstype = {'SEG2': 'SEG', 'TRI3': "TRI", "QUAD4":"QUAD", "TETRA4":"TETRA", "HEXA8":"HEXA", "PYRA5":"PYRA", \
             "PENTA6":"PENTA", "SEG3":'SEG', "TRI6":"TRI", "QUAD8":"QUAD", "TETRA10":"TETRA", "HEXA20":"HEXA", "PYRA13":"PYRA", \
                  "PENTA15":"PENTA", "SEG4":'SEG', "TRI7":"TRI", "QUAD9":"QUAD", "PENTA18":"PENTA", "HEXA27":"HEXA" }
        self.mfaces = {"TRI" : ((0,1),(1,2),(2,0)), "QUAD": ((0,1),(1,2),(2,3),(3,0)),
                       "TETRA" : [["TRI3", [0,1,2]], ["TRI3", [1,2,3]], ["TRI3", [0,1,3]], ["TRI3", [0,2,3]]],
                       "HEXA" : [["QUAD4", [0,1,2,3]], ["QUAD4", [0,3,7,4]], ["QUAD4", [2,6,7,3]], \
                           ["QUAD4", [1,2,6,5]],["QUAD4", [0,1,5,4]], ["QUAD4", [4,7,6,5]]],
                        "PENTA" : [["TRI3", [0,1,2]], ["QUAD4", [1,4,5,2]], ["QUAD4", [0,1,4,3]], \
                           ["QUAD4", [0,2,5,3]], ["TRI3", [3,4,5]]],
                        "PYRA" : [["QUAD4", [0,1,2,3]], ["TRI3", [4,3,2]], ["TRI3", [4,2,1]], \
                           ["TRI3", [4,1,0]], ["TRI3", [4,0,3]]]}

    def getType(self, typeC):
        return self.mdata[typeC]

    def getDim(self, typeC):
        return self.mdim[typeC]

    def addCell(self, mtype, nodes, edges, faces, volumes):
        if self.mdim[mtype] == 0:
            logging.info("0D-cell are not supported. We skip them...")
        elif self.mdim[mtype] == 1:
            mnodes = rotate(nodes[0:self.mnodes[mtype]])
            if mnodes in self.list_edges:
                return self.list_edges[mnodes]
            else:
                edge = Edge()
                edge.nodes = mnodes
                edges.append(edge)
                self.list_edges[edge.nodes] = len(edges)-1
                return len(edges)-1
        elif self.mdim[mtype] == 2:
            mnodes = rotate(nodes[0:self.mnodes[mtype]])
            if mnodes in self.list_faces:
                return self.list_faces[mnodes]
            else:
                face = Face()
                face.nodes = mnodes
                list_egdes = self.mfaces[self.mstype[mtype]]
                for edge in list_egdes:
                    enodes = [nodes[edge[0]], nodes[edge[1]]]
                    eid = self.addCell("SEG2", enodes, edges, faces, volumes)
                    face.edges.append(eid)
                faces.append(face)
                self.list_faces[face.nodes] = len(faces)-1
                return len(faces)-1
        elif self.mdim[mtype] == 3:
            mnodes = rotate(nodes[0:self.mnodes[mtype]])
            volu = Volume()
            volu.nodes = mnodes
            list_faces = self.mfaces[self.mstype[mtype]]
            for face in list_faces:
                fnodes = []
                for nid in face[1]:
                    fnodes.append(nodes[nid])
                fid = self.addCell(face[0], fnodes, edges, faces, volumes)
                volu.faces.append(fid)
                volu.edges += faces[fid].edges

            volu.edges = sorted(set(tuple(volu.edges)))
            volumes.append(volu)
            return len(volumes)-1
        else:
            raise RuntimeError("Should not be here")


class PolyMesh:

    def __init__(self, level=logging.INFO):
        logging.basicConfig(format='%(levelname)s:%(message)s', level=level)
        self.edges = []
        self.faces = []
        self.volumes = []
        self.edges_grp = []
        self.faces_grp = []
        self.volumes_grp = []
        self.nodes_grp = []


    def createMesh(self, args):
        medmesh = medcoupling.MEDFileUMesh(args.fileNameInput)
        self.dimension = medmesh.getSpaceDimension()

        nb_nodes = medmesh.getNumberOfNodes()
        logging.info("Reading %d nodes..."%(nb_nodes))
        self.nodes = medmesh.getCoords().getValuesAsTuple()

        cc = CellConverter()

        non_empty_levs = medmesh.getNonEmptyLevels()
        cells_shift = 0
        logging.info("Reading %d levels..."%(len(non_empty_levs)))
        for lev in non_empty_levs:
            mesh_lev = medmesh[lev]
            types_at_level = mesh_lev.getAllGeoTypesSorted()
            loc_glo_ind = {}
            cell_off = cells_shift
            mdim = 0
            for medcoupling_cell_type in types_at_level :
                cells_by_type = mesh_lev.giveCellsWithType(medcoupling_cell_type).getValues()
                mtype = cc.getType(medcoupling_cell_type)
                logging.info("**Reading %d %s..."%(len(cells_by_type), mtype))
                mdim = cc.getDim(mtype)
                for cell in cells_by_type :
                    element_nodes_med = mesh_lev.getNodeIdsOfCell(cell)
                    cell_id = cc.addCell(mtype, element_nodes_med, self.edges, self.faces, self.volumes)
                    loc_glo_ind[cell_off] = cell_id
                    cell_off += 1

            for group in medmesh.getGroupsOnSpecifiedLev(lev):
                ids = cells_shift + medmesh.getGroupArr(lev, group)
                val = [loc_glo_ind[cid] for cid in ids.getValues()]
                if mdim == 0:
                    logging.info("0D-groups are not supported. We skip them...")
                if mdim == 1:
                    self.edges_grp.append([group, val])
                elif mdim == 2:
                    self.faces_grp.append([group, val])
                elif mdim == 3:
                    self.volumes_grp.append([group, val])
                else:
                    raise RuntimeError("Should not be here")

            cells_shift+=mesh_lev.getNumberOfCells()


        #clean nodes
        nodes_to_keep = {}
        for edge in self.edges:
            for node in edge.nodes:
                if node not in nodes_to_keep:
                    nodes_to_keep[node] = True
        for face in self.faces:
            for node in face.nodes:
                if node not in nodes_to_keep:
                    nodes_to_keep[node] = True
        for vol in self.volumes:
            for node in vol.nodes:
                if node not in nodes_to_keep:
                    nodes_to_keep[node] = True

        node_id = 0
        new_nodes = []
        for node in nodes_to_keep:
            new_nodes.append(self.nodes[node])
            nodes_to_keep[node] = node_id
            node_id += 1

        self.nodes = new_nodes

        for edge in self.edges:
            edge.nodes = [nodes_to_keep[edge.nodes[i]] for i in range(len(edge.nodes))]
        for face in self.faces:
            face.nodes = [nodes_to_keep[face.nodes[i]] for i in range(len(face.nodes))]
        for vol in self.volumes:
            vol.nodes = [nodes_to_keep[vol.nodes[i]] for i in range(len(vol.nodes))]

        # add grp of nodes
        for group in medmesh.getGroupsOnSpecifiedLev(1):
            ids = medmesh.getGroupArr(1, group).getValues()
            val = []
            for nid in ids:
                if nid in nodes_to_keep:
                    val.append(nodes_to_keep[nid])
            self.nodes_grp.append([group, val])

    def writeMesh(self, args):
        f = open(args.fileNameOutput, 'w')

        f.write("**BeginMesh\n")
        f.write("*Version %d\n"%1)
        f.write("*Dimension %d\n"%self.dimension)
        f.write("*Nodes %d\n"%len(self.nodes))
        i_node = 0
        for node in self.nodes:
            if len(node) == 2:
                f.write("%d %f %f %f \n"%(i_node, node[0], node[1], 0.0))
            else:
                f.write("%d %f %f %f \n"%(i_node, node[0], node[1], node[2]))
            i_node += 1

        if len(self.edges) > 0:
            f.write("*Edges->Nodes %d\n"%len(self.edges))
            i_edge = 0
            for edge in self.edges:
                f.write("%d %d %d %d \n"%(i_edge, 2, edge.nodes[0], edge.nodes[1]))
                i_edge += 1

        if len(self.faces) > 0:
            f.write("*Faces->Nodes %d\n"%len(self.faces))
            i_face = 0
            for face in self.faces:
                f.write("%d %d "%(i_face, len(face.nodes)))
                for node in face.nodes:
                    f.write("%d "%(node))
                f.write("\n")
                i_face += 1
            f.write("*Faces->Edges %d\n"%len(self.faces))
            i_face = 0
            for face in self.faces:
                f.write("%d %d "%(i_face, len(face.edges)))
                for edge in face.edges:
                    f.write("%d "%(edge))
                f.write("\n")
                i_face += 1

        if len(self.volumes) > 0:
            f.write("*Volumes->Nodes %d\n"%len(self.volumes))
            i_vol = 0
            for vol in self.volumes:
                f.write("%d %d "%(i_vol, len(vol.nodes)))
                for node in vol.nodes:
                    f.write("%d "%(node))
                f.write("\n")
                i_vol += 1
            f.write("*Volumes->Edges %d\n"%len(self.volumes))
            i_vol = 0
            for vol in self.volumes:
                f.write("%d %d "%(i_vol, len(vol.edges)))
                for edge in vol.edges:
                    f.write("%d "%(edge))
                f.write("\n")
                i_vol += 1
            i_vol = 0
            f.write("*Volumes->Faces %d\n"%len(self.volumes))
            for vol in self.volumes:
                f.write("%d %d "%(i_vol, len(vol.faces)))
                for face in vol.faces:
                    f.write("%d "%(face))
                f.write("\n")
                i_vol += 1

        i_grp = 0
        if len(self.nodes_grp) > 0:
            f.write("*Nodes->Groups %d\n"%len(self.nodes_grp))
            for grp in self.nodes_grp:
                logging.info("Group naming: %s (med name) -> %d (internal id)"%(grp[0], i_grp))
                f.write("%s %d %d "%(grp[0], i_grp, len(grp[1])))
                for cell in grp[1]:
                    f.write("%d "%cell)
                i_grp+=1
                f.write("\n")

        if len(self.edges_grp) > 0:
            f.write("*Edges->Groups %d\n"%len(self.edges_grp))
            for grp in self.edges_grp:
                logging.info("Group naming: %s (med name) -> %d (internal id)"%(grp[0], i_grp))
                f.write("%s %d %d "%(grp[0], i_grp, len(grp[1])))
                for cell in grp[1]:
                    f.write("%d "%cell)
                i_grp+=1
                f.write("\n")

        if len(self.faces_grp) > 0:
            f.write("*Faces->Groups %d\n"%len(self.faces_grp))
            for grp in self.faces_grp:
                logging.info("Group naming: %s (med name) -> %d (internal id)"%(grp[0], i_grp))
                f.write("%s %d %d "%(grp[0], i_grp, len(grp[1])))
                for cell in grp[1]:
                    f.write("%d "%cell)
                i_grp+=1
                f.write("\n")

        if len(self.volumes_grp) > 0:
            f.write("*Volumes->Groups %d\n"%len(self.volumes_grp))
            for grp in self.volumes_grp:
                logging.info("Group naming: %s (med name) -> %d (internal id)"%(grp[0], i_grp))
                f.write("%s %d %d "%(grp[0], i_grp, len(grp[1])))
                for cell in grp[1]:
                    f.write("%d "%cell)
                i_grp+=1
                f.write("\n")


        f.write("**EndMesh")
        f.close()

# This a python script to convert a mesh from MED format to poly2d or poly3d format
# This script need the medcoupling python module (present in salome)
# For each group, we add an unique internal id. This internal id has to be used ins DiSk++

# Inside a terminal
# If you have not medcoupling but salome
# ./salome shell
# Then use the script (or directly)
# python med2poly.py -fi myInputMesh.med -fo myOutputMesh.poly2d

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-fi','--filenameIn', dest="fileNameInput", required=True, help="med file path")
    parser.add_argument('-fo','--filenameOut', dest="fileNameOutput", default=True, help="mesh contained in filepath to be considered")
    parser.add_argument('-v','--verbose', dest = "verbosity", type=int , default=1, help="verbosity. Default 1 (print info). 0 only errors reported. 2 and higher all debug messages")

    args = parser.parse_args()

    pmesh = PolyMesh()
    pmesh.createMesh(args)
    pmesh.writeMesh(args)
