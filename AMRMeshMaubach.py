import TopFileTool
from Utility import *


'''
Implement the algorithm in

Bartels, Sören, and Patrick Schreier. "Local coarsening of simplicial finite element meshes generated by bisections."
BIT Numerical Mathematics 52.3 (2012): 559-569.

Stevenson, Rob. "The completion of locally refined simplicial partitions created by bisection."
Mathematics of computation 77.261 (2008): 227-241.
APA


Slightly modified from
Maubach, Joseph M. "Local bisection refinement for n-simplicial grids generated by reflection."
SIAM Journal on Scientific Computing 16.1 (1995): 210-227.
'''





class MarkedTet():
    def __init__(self, nn, type, bFlag = False):
        '''
        :param nn: n0,n1,n2,n3
        :param type: element flag 0 , 1,...,d-1

        constructor for markedTet
        nn = [n0, n1, n2, n3]
        type = gam
        bFlag: boundary flag, False: if it is not on the boundary; True: it is potentially on the boundary,
        it can also not be on the boundary
        '''
        self.nn = nn
        self.type = type
        self.bFlag = bFlag

    def get_ref_edge(self):
        nn = self.nn
        return [nn[0], nn[-1]]

    def get_cut_face(self):
        nn = self.nn
        return [[nn[0],nn[1],nn[-1]],[nn[0],nn[2],nn[-1]]]

    def get_opposite_edge(self):
        nn = self.nn
        return [nn[1],nn[2]]


    def bisect_tet(self,n, Flag = False):
        '''
        bisect the tet
        :param n: new node number
        :return: return two marked tet
        bisect the tet into 2
        T = (z0,z1 ... zd)_g
        s1(T) = (z0, y, z1,..., z_g, z_{g+1},..., z_{d-1})_{g+1} mod d
        s2(T) = (z0, y, z1,..., z_g, z_{d-1},..., z_{g+1})_{g+1} mod d
        '''
        type = self.type
        nn = self.nn
        d = 3 # the dimension of the problem

        nodes =[[],[]]
        nodes[0] = [nn[0] , n, *nn[1:d]]
        nodes[1] = [nn[-1], n, *nn[1:type + 1], *nn[d - 1 : type: -1]]

        return MarkedTet(nodes[0],(type+1)%d, Flag),MarkedTet(nodes[1],(type+1)%d, Flag)










class AmrMeshMaubach():

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
        nodes,initElems, initBoundaries = TopFileTool.read_tet(mshfile)

        self.nodes = nodes
        self.initElems = initElems

        #build dictionary edge->midNode
        self._build_mid_node()

        self._build_marked_elems()

        self._build_boundary_tris(initBoundaries)


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


    def _build_marked_elems(self):
        '''
        This function should only be called in the class constructor, it uses initElems
        to initialize marked elements
        '''
        initElems = self.initElems
        markedElems = []

        for ele in initElems:
            markedElems.append(MarkedTet([ele[0], ele[1], ele[2], ele[3]], 0))

        self.markedElems = markedElems

    def _build_boundary_tris(self, initBoundaries):
        '''
        This function should only be called in the class constructor, it uses initElems
        to initialize marked elements
        '''
        boundaries = set()
        for i,j,k in initBoundaries:
            boundaries.add(triplet_sort(i,j,k))
        self.boundaries = boundaries

    def _update_boundary_eles(self):
        '''
        This function should only be called in the class constructor, it uses initElems
        and boundaries to update the boundary flag in marked elements
        '''
        markedElems = self.markedElems
        boundaries = self.boundaries
        for e in markedElems:

            nn = e.nn
            if((triplet_sort(nn[0],nn[1],nn[2]) in boundaries) or (triplet_sort(nn[0],nn[1],nn[3]) in boundaries) or
               (triplet_sort(nn[0],nn[2],nn[3]) in boundaries) or (triplet_sort(nn[1],nn[2],nn[3]) in boundaries)):
                e.bFlag = True


    def bisection(self, rElems):
        '''
        rElems: element list, need to be refined
        update: nodes, markedElems, midNode
        '''
        nodes = self.nodes
        markedElems = self.markedElems
        midNode = self.midNode
        boundaries = self.boundaries

        for e in rElems:
            elem = markedElems[e]
            newNode = midNode[pair_sort(*elem.get_ref_edge())]
            if(newNode == -1):
                #The edge is cut first time, add a new node(the midpoint),
                newNode = len(nodes)
                n0,nd = elem.get_ref_edge()
                y = [(nodes[n0][0] + nodes[nd][0]) / 2.0, (nodes[n0][1] + nodes[nd][1]) / 2.0,
                     (nodes[n0][2] + nodes[nd][2]) / 2.0]
                nodes.append(y)
                #mark the refined edge's midpoint
                midNode[pair_sort(*elem.get_ref_edge())] = newNode

                #add new edge, for there are 4 new edges
                for n in elem.nn:
                    midNode[pair_sort(n, newNode)] = -1

            else:
                #add new edges, for there are potentially 2 new edges on the cut faces
                for n in  elem.get_opposite_edge():
                    if pair_sort(n, newNode) not in midNode:
                        midNode[pair_sort(n, newNode)] = -1

            # update boundary faces, for boundary elements
            if(elem.bFlag):
                #the element is potentially on the boundary
                cutFace = elem.cut_face()
                for i,j,k in cutFace:
                    if triplet_sort(i,j,k) in boundaries:
                        boundaries.remove(triplet_sort(i,j,k))
                        boundaries.add(triplet_sort(i, j, newNode))
                        boundaries.add(triplet_sort(j, k, newNode))


            # add new elements
            elem0, elem1 = elem.bisect_tet(newNode, elem.bFlag)
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
        nn = elem.nn
        for i in range(4):
            for j in range(i+1,4):
                if midNode[pair_sort(nn[i],nn[j])] >= 0:
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
            elems.append(markedElems[i].nn)
        boundaries = list(self.boundaries)
        TopFileTool.write_tet(nodes,elems, boundaries, mshfile)



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
            xyz = np.array([nodes[ele.nn[0]],nodes[ele.nn[1]],nodes[ele.nn[2]],nodes[ele.nn[3]]])
            if inside_tet(xyz, point):
                rElems.append(e)
        return rElems




if __name__ == '__main__':
    print('AMR Mesh Maubach')
    amrMesh = AmrMeshMaubach('domain.top')
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










