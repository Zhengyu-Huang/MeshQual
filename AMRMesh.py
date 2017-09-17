import sys
import numpy as np
import matplotlib.pyplot as plt

def pair_sort(a,b):
    return (a,b) if a < b else (b,a)

def terahedron():
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

    def bisect_tet(self,n):
        type = self.get_type()
        nodes1 = [a0, b0, c0, n]
        nodes1 = [a1, b1, c1, n]
        faceEdges1 =
        faceEdges2 =
        flag1 = 1 if type == 'Pu' else 0
        flag2 = 1 if type == 'Pu' else 0


        return tet1, tet2





class AMRmesh( ):

    def __init__(self, mshfile):
        '''
        nodes: float[:,3] node coordinates
        node_parents : int[:,2], its two parents, -1 means no parents
        elems: int[:,5] its 4 nodes and its type(for refinement)
        edge_t: int[:] if 1 it is cut else 0


        :param mshfile:
        '''
        self.mshfile = mshfile
        try:
            fid = open(mshfile, "r")
        except IOError:
            print("File '%s' not found." % mshfile)
            sys.exit()

        nodes = []
        elems = []

        print('Reading mesh ...')

        line = fid.readline()  # should be "Nodes FluidNodes"

        while line:
            line = fid.readline()
            data = line.split()
            if data[0] == 'Elements':
                break;
            nodes.append(list(map(float, data[1:4])))

        while line:
            line = fid.readline()
            data = line.split()
            if data[0] == 'Elements':
                break;
            elem = list(map(int, data[2:6]))
            elem = [x - 1 for x in elem]  # fixed top file 1-based index
            elem.append(0) # the generation
            elems.append(elem)

        self.nodes = nodes
        self.elems = elems

        #build dictionary edge->midNode

        self._buildMidNode()

        self.nNodes = len(nodes)
        self.nElems = len(elems)

    def bisect_tet(self,i):
        '''
        bisect tetrahedron e_i
        :return: marked tet1 and marked tet2
        '''

    def _buildMidNode(self):
        elems = self.elems
        midNode = {}
        for e in range(len(elems)):
            nodes = elems[e][0:4]
            for i in range(4):
                for j in range(i + 1, 4):
                    midNode[pair_sort(nodes[i], nodes[j])] = -1
        self.midNode = midNode


    def bisection(self, rElems):
        '''
        rElems: element list, need to be refined
        :return:
        '''
        nodes = self.nodes
        elems = self.elems
        midNode = self.midNode

        for e in rElems:

            newNode = midNode[e.get_refEdge()]
            if(newNode == -1):
                #The edge is cut first time, add new nodes,
                newNode = len(nodes)
                n0,n1 = e.get_refEdge()
                y = [(nodes[n0][0] + nodes[n1][0]) / 2.0, (nodes[n0][1] + nodes[n1][1]) / 2.0,
                     (nodes[n0][2] + nodes[n1][2]) / 2.0]
                nodes.append(y)
                #mark the edge
                midNode[e.get_refEdge()] = newNode

                #add new edge
                n2,n3 = e.get_#todo should check first?
                midNode[pair_sort(n2, newNode)] = -1
                midNode[pair_sort(n3, newNode)] = -1



            # add new elements
            tet1, tet2 = elems.bisect_tet(newNode)
            elems[e] = tet1
            elems.append(tet2)


    def _nonconformity(self,e):
        elems = self.elems
        midNode = self.midNode
        nodes = elems[e][0:4]
        for i in range(4):
            for j in range(i+1,4):
                if midNode[pair_sort(nodes[i],nodes[j])] >= 0:
                    return True


    def check_conformity(self):
        elems = self.elems
        rElems = []
        for e in elems:
            if self._nonconformity(e):
                rElems.append(e)
        return rElems


    def refine(self, rElems):
        while not rElems:
            self.bisection(rElems)
            rElems = self.check_conformity()

    def coarsen(self):
        #todo

    def output_to_top(self):
        #todo







if __name__ == '__main__':





