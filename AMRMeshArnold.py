import numpy as np
import operator
import matplotlib.pyplot as plt
import copy
import TopFileTool
from Utility import *

'''
Implement the algorithm in
Arnold, Douglas N., Arup Mukherjee, and Luc Pouly.
"Locally adapted tetrahedral meshes using bisection."
SIAM Journal on Scientific Computing 22.2 (2000): 431-448.
'''


class MarkedTet():
    def __init__(self, nodes,  faceEdges, flag):
        '''
        :param nodes: n0,n1,n2,n3
        :param refEdge: MAKE SURE the refine edge is (n0,n1)
        :param faceEdges: (a0,b0),(a1,b1),(a2,b2),(a3,b3) and  (ai,bi) corresponding to face not contain ni
        :param flag: element flag 0 or 1
        MAKE SURE the edge pair (a,b) must satisfy a < b
        :return:
        '''
        self.nodes = nodes
        self.faceEdges = faceEdges
        self.flag = flag
    def get_nodes(self):
        return self.nodes
    def get_ref_edge(self):
        return self.nodes[0],self.nodes[1]
    def get_opposite_edge(self):
        return self.nodes[2],self.nodes[3]


    def bisect_tet(self,n):
        '''
        bisect the tet
        :param n: new node number
        :return: return two marked tet
        '''
        type,nf = self.get_type()
        nodes,faceEdges,flag =[],[],[]
        for i in range(2):
            #tet i is the tet,  ni, n, n2, n3

            a,b = self.faceEdges[1-i]
            #find the other node on the inherited face
            for val in [self.nodes[i],self.nodes[2],self.nodes[3]]:
                if val not in self.faceEdges[1 - i]:
                    c = val
                    break

            nodes.append([a, b, c, n])


            # the edge marked on the inherited face is (a,b), it is also refinement edge -> faceEdges[3] = (a,b)
            # cut faces, : opposite of the new node
            # new faces, : opposite of the new node(default), if type= 'Pf' its marked edge is the edge
            # connecting the new vertex to
            # by default : the other 3 marked edges, we first put the opposite of the new node
            # and if the type = 'Pf', the new face is self.nodes[i] opposite face, the marked edge is (n, nf)
            faceEdges.append([pair_sort(b,c),pair_sort(a,c),pair_sort(a,b),(a,b)])

            if type == 'Pf':
                # for coplanar, we know the new node is either a nor b
                faceEdges[i][(0 if self.nodes[i] == a else 1)] = pair_sort(n,nf)

            flag.append( 1 if type == 'Pu' else 0)



        return MarkedTet(nodes[0],faceEdges[0],flag[0]),MarkedTet(nodes[1],faceEdges[1],flag[1])

    def get_type(self):
        '''
        :return: the type of the element, if the type is Pf, also return the node one the the
        coplanar face, but not on the refinement edge
        '''


        faceEdgeNodes = set()
        for e in self.faceEdges:
            for i in range(2):
                faceEdgeNodes.add(e[i])
        coplanar = False if len(faceEdgeNodes) == 4 else True

        nf = -1
        if(coplanar):
            for n in faceEdgeNodes:
                if n != self.nodes[0] and n != self.nodes[1]:
                    nf = n
                    break
            #it is Pf if flag is set, otherwise Pu
            return ('Pf', nf)  if self.flag == 1 else ('Pu',nf)
        else:
            return 'AOM',nf








class AmrMeshArnold():

    def __init__(self, mshfile):
        '''
        nodes: a list of float[3], node coordinates
        midNode: map from edge pair (a,b) a<b,  to its midpoint id
                 If it does not have midpoint node, the midpoint id is -1.
                 This array is used to detect hanging nodes
        markedElems: a list of markedTet class
        initElems: a list of int[4], original mesh element nodes


        :param mshfile:
        '''
        nodes,initElems,_ = TopFileTool.read_tet(mshfile)

        self.nodes = nodes
        self.initElems = initElems

        #build dictionary edge->midNode

        self._build_mid_node()
        self._build_marked_elems()





    def _build_mid_node(self):
        '''
        This function should only be called in the class constructor, it uses initElems
        to initialize midNode mesh
        '''

        initElems = self.initElems
        midNode = {}
        for e in range(len(initElems)):
            nodes = initElems[e][0:4]
            for i in range(4):
                for j in range(i + 1, 4):
                    midNode[pair_sort(nodes[i], nodes[j])] = -1

        self.midNode = midNode


    def _max_pair(self,nodes,edgeSort):
        '''
        find the edge that has maximal edgeSort number
        :param nodes: 3 nodes forming triangle or 4 nodes forming tetrahedron
        :param edgeSort: a dictionary maps edge to its value
        :return: the edge that has maximal value
        '''
        maxLenId = -1
        nNodes = len(nodes)
        n0,n1 = -1,-1
        for i in range(nNodes):
            for j in range(i + 1, nNodes):
                if(edgeSort[pair_sort(nodes[i], nodes[j])] > maxLenId):
                    n0,n1 = pair_sort(nodes[i], nodes[j])
                    maxLenId = edgeSort[pair_sort(nodes[i], nodes[j])]
        return n0,n1,maxLenId

    def _build_marked_elems(self):
        '''
        This function should only be called in the class constructor after the _buildMidNode,
        it uses initElems to build marked elements

        Strictly order the edges of the mesh in an arbitrary but fixed manner,
        for example, by length with a well-defined tie-breaking rule. Then choose the maximal
        edge of each tetrahedron as its refinement edge and the maximal edge of each face as
        its marked edge. Unset the flag on all tetrahedra.
        '''
        nodes = self.nodes
        edgeSort = copy.deepcopy(self.midNode)
        nEdge = len(edgeSort)
        for edgePair in edgeSort:
            n0,n1 = np.array(nodes[edgePair[0]]),np.array(nodes[edgePair[1]])
            edgeSort[edgePair] = np.linalg.norm(n0 - n1)

        edgeSortList = sorted(edgeSort.items(), key=operator.itemgetter(1)) #now edge Sort is a list

        for i in range(nEdge):
            edgeSort[(edgeSortList[i][0])] = i

        markedElems = []
        initElems = self.initElems
        for ele in initElems:
            #find the max edge in the elem

            n0,n1,_ = self._max_pair(ele,edgeSort)
            n2,n3 = [n for n in ele if n not in [n0,n1]]
            a0,b0,_ = self._max_pair([n1,n2,n3],edgeSort)
            a1,b1,_ = self._max_pair([n0,n2,n3],edgeSort)

            #find the max edge on the faces
            markedElems.append(MarkedTet([n0,n1,n2,n3],  [(a0,b0),(a1,b1),(n0,n1),(n0,n1)], 0))


        self.markedElems = markedElems

    def bisection(self, rElems):
        '''
        rElems: element list, need to be refined
        update: nodes, markedElems, midNode
        '''
        nodes = self.nodes
        markedElems = self.markedElems
        midNode = self.midNode

        for e in rElems:
            elem = markedElems[e]
            newNode = midNode[elem.get_ref_edge()]
            if(newNode == -1):
                #The edge is cut first time, add new nodes,
                newNode = len(nodes)
                n0,n1 = elem.get_ref_edge()
                y = [(nodes[n0][0] + nodes[n1][0]) / 2.0, (nodes[n0][1] + nodes[n1][1]) / 2.0,
                     (nodes[n0][2] + nodes[n1][2]) / 2.0]
                nodes.append(y)
                #mark the refined edge's midpoint
                midNode[elem.get_ref_edge()] = newNode

                #add new edge, for there are 4 new edges
                for n in elem.nodes:
                    midNode[pair_sort(n, newNode)] = -1

            else:
                #add new edge, for there are potentially 2 new edges on the cut faces
                for n in  elem.get_opposite_edge():
                    if pair_sort(n, newNode) not in midNode:
                        midNode[pair_sort(n, newNode)] = -1

            # add new elements
            elem0, elem1 = elem.bisect_tet(newNode)
            markedElems[e] = elem0
            markedElems.append(elem1)


    def _nonconformity(self,e):
        '''
        loop the element's all edges, check it has hanging nodes or not
        :param e: element id
        :return: True: has hanging nodes, otherwise False
        '''
        elem = self.markedElems[e]
        midNode = self.midNode
        nodes = elem.get_nodes()
        for i in range(4):
            for j in range(i+1,4):
                if midNode[pair_sort(nodes[i],nodes[j])] >= 0:
                    return True


    def check_conformity(self):
        '''
        check whether current mesh conforming or not
        :return: rElems: a list of element id. These elements have hanging nodes, need to be cut
        '''
        markedElems = self.markedElems
        nElems = len(markedElems)
        rElems = []
        for e in range(nElems):
            if self._nonconformity(e):
                rElems.append(e)
        return rElems


    def refine(self, rElems):
        '''
        Refine mesh
        :param rElems: a list of int, for marked element ids to be refined
        :return:
        '''
        print('Refining mesh...')
        while rElems:
            print('In Refining mesh loop, for %d tetrahedrons...' %(len(rElems)))
            self.bisection(rElems)
            rElems = self.check_conformity()


    def coarsen(self):
        #todo
        print('todo, have not implemented yet')

    def output_to_top(self,mshfile):
        '''
        output the results to top file
        :param mshfile: output top file name
        :return:
        '''
        nodes = self.nodes
        markedElems = self.markedElems
        elems = []
        for i in range(len(markedElems)):
            elems.append(markedElems[i].nodes)
        #todo no boundaries
        TopFileTool.write_tet(nodes,elems, [], mshfile)



    def find(self,point):
        '''
        find elements that contains this point
        :param point: float[3]
        :return: a list of element ids
        '''
        markedElems = self.markedElems
        nodes = self.nodes
        nElems = len(markedElems)
        rElems = []
        for e in range(nElems):
            ele = markedElems[e]
            xyz = np.array([nodes[ele.nodes[0]],nodes[ele.nodes[1]],nodes[ele.nodes[2]],nodes[ele.nodes[3]]])
            if inside_tet(xyz, point):
                rElems.append(e)
        return rElems




if __name__ == '__main__':
    test = 2
    if(test == 0):
        amrMesh = AmrMeshArnold('Trirectangular_tet_mesh.top')
        amrMesh.refine([0])
        amrMesh.output_to_top('Refined_mesh.top')

    elif(test == 1):
        #Test: read a cube[-1,1][-1,1][-1,1], refine the element containing (0.2., 0.3, 0.4)
        amrMesh = AmrMeshArnold('cube.top')
        refineLevel = 50
        singularity = np.array([0.2,0.3,0.4])
        for i in range(refineLevel):
            rElems = amrMesh.find(singularity)
            amrMesh.refine(rElems)
        amrMesh.output_to_top('Refined_mesh.top')

        from MeshQual import Mesh
        tet_mesh = Mesh('Refined_mesh.top')
        AR = tet_mesh.AR
        tet_mesh.plot()

        tet_mesh = Mesh('cube.top')
        AR = tet_mesh.AR
        tet_mesh.plot()
    elif(test == 2):
        print('AMR Mesh Maubach')
        amrMesh = AmrMeshArnold('domain.top')
        refineLevel = 50
        singularity = np.array([2.3,0.3,0.4])
        for i in range(refineLevel):
            rElems = amrMesh.find(singularity)
            amrMesh.refine(rElems)
        amrMesh.output_to_top('Refined_mesh.top')

        from MeshQual import Mesh
        tet_mesh = Mesh('Refined_mesh.top')
        AR = tet_mesh.AR
        tet_mesh.plot()







