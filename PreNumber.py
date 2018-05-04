'''
Alkämper, Martin, Fernando Gaspoz, and Robert Klöfkorn.
"A Weak Compatibility Condition for Newest Vertex Bisection in any dimension."
arXiv preprint arXiv:1711.03141 (2017).
'''
from TopFileTool import write_tet, read_tet
from Utility import pair_sort,triplet_sort
import numpy as np
import operator



def edge_in_face(edge, face):
   '''
   :param edge: [n1, n2]
   :param face: [m1,m2,m3]
   :return: the edge is on the face or not
   '''

   if edge[0] in set(face) and edge[1] in set(face):
       return True
   else:
       return False





class PreNumber:
    def __init__(self,inputfile,outputfile):
        self.nodes, self.eles, self.boundaryNames, self.boundaries = read_tet(inputfile)
        self._build_faces()
        self._mark_edge()
        self._global_ordering()
        self._renumber()
        self.write_topfile(outputfile)



    def _build_faces(self):
        '''
        :return: build a face map faces. which map the face to its neighbor tets
        '''
        eles = self.eles
        faces = {}
        for ei in range(len(eles)):
            e = eles[ei]
            for i in range(4):
                enen = [e[j] for j in range(4) if (j != i)]

                if triplet_sort(enen[0], enen[1], enen[2]) in faces:
                    faces[triplet_sort(enen[0], enen[1], enen[2])].append(ei)
                else:
                    faces[triplet_sort(enen[0], enen[1], enen[2])] = [ei]

        self.faces = faces

    def _check_on_boundary(self, fp):
        '''
        :param fp: [m1,m2,m3]
        :return: check whether the face is on the computational domain boundary
        '''
        faces = self.faces
        es = faces[fp]
        return True if len(es) == 1 else False

    def _mark_edge(self):
        '''
        for each element,
        mark its longest edge, if it has multiple longest edges,
        mark the one with largest value, i + j,
        here i, j are the two node ids
        :return: refineEdge
        '''

        eles = self.eles
        nodes = self.nodes
        refineEdges = []

        for e in eles: #loop all elements
            maxLen, maxSum = -1.0, 0
            for i in range(4):
                for j in range(i):
                    ni, nj = e[i], e[j]
                    len = np.sqrt((nodes[ni][0] - nodes[nj][0])**2 + (nodes[ni][1] - nodes[nj][1])**2 + (nodes[ni][2] - nodes[nj][2])**2)
                    sum = ni + nj
                    if(len > maxLen or (len == maxLen and sum > maxSum)):
                        refineEdge = pair_sort(ni, nj)
                        maxLen = len
                        maxSum = sum

            refineEdges.append(refineEdge)

        self.refineEdges = refineEdges

    def _refine_edge(self, ele, nodeMap):
        '''
        :param ele: int
        :param nodeMap: a map
        :return: the refinement edge [n0, nd] of the element
        '''
        min, max = np.inf, -1
        for n in ele:
            if nodeMap[n] < min:
                min = nodeMap[n]
                n0 = n
            if nodeMap[n] > max:
                max = nodeMap[n]
                nd = n
        return [n0,nd]

    def _putIntoNonBoundaryFaceSet(self, elem_id, nonBoundaryFaceSet):
        '''
        :param elem_id: element id
        :param nonBoundaryFaceSet: a list F
        :return: put the faces of the element which are not on the fluid domain boundary to the set
        '''
        eles = self.eles
        refineEdges = self.refineEdges

        for fp in [triplet_sort(eles[elem_id][1], eles[elem_id][2], eles[elem_id][3]),
                   triplet_sort(eles[elem_id][0], eles[elem_id][2], eles[elem_id][3]),
                   triplet_sort(eles[elem_id][0], eles[elem_id][1], eles[elem_id][3]),
                   triplet_sort(eles[elem_id][0], eles[elem_id][1], eles[elem_id][2])]:
            if self._check_on_boundary(fp):
                continue
            else:
                if edge_in_face(refineEdges[elem_id], fp):
                    nonBoundaryFaceSet.insert(0, fp)
                else:
                    nonBoundaryFaceSet.append(fp)



    def _global_ordering(self):
        '''
        Order all nodes, we build the map nodeMap, the order is
        (nodeMap[n], n) < = > (nodeMap[m], m)
        todo: different from the one in the paper, to make it O(N)
        the nodeMap number can be equal for different nodes
        :return:
        '''
        faces = self.faces
        eles = self.eles
        refineEdges = self.refineEdges

        elementMarker = np.zeros(len(eles))
        nodeMap = {}
        nonBoundaryFaceSet = []

        ### Put the first element
        elementMarker[0] = 1
        self._putIntoNonBoundaryFaceSet(0, nonBoundaryFaceSet)
        # order these nodes of the first element
        n0, nd = refineEdges[0]
        nodeMap[n0] = 0
        nodeMap[nd] = 3
        temp = 1
        for ni in eles[0]:
            if ni != n0 and ni != nd:
                nodeMap[ni] = temp
                temp = temp + 1


        while nonBoundaryFaceSet: # the element
            f = nonBoundaryFaceSet.pop(0)
            es = faces[f]
            vf1, vf2, vf3 = f
            for i in range(2):
                #prepare to add the nodes of eles[es[1-i]]
                if(elementMarker[es[1 - i]]):
                     continue

                v = eles[es[i]][0] + eles[es[i]][1] + eles[es[i]][2] + \
                    eles[es[i]][3] - vf1 - vf2 - vf3

                vp = eles[es[1 - i]][0] + eles[es[1 - i]][1] + eles[es[1 - i]][2] + \
                    eles[es[1 - i]][3] - vf1 - vf2 - vf3
                if vp not in nodeMap: #vp is not in the nodeMap
                    if((nodeMap[v],v) > (nodeMap[vf1],vf1) and (nodeMap[v],v) > (nodeMap[vf2],vf2) and (nodeMap[v],v) > (nodeMap[vf3],vf3)) or \
                            ((nodeMap[v],v) < (nodeMap[vf1],vf1) and (nodeMap[v],v) < (nodeMap[vf2],vf2) and (nodeMap[v],v) < (nodeMap[vf3],vf3)):
                        #no RefEdge insert before v or after v
                        nodeMap[vp] = nodeMap[v]
                    else:
                        #RefEdge insert after v
                        nodeMap[vp] = nodeMap[v]
                elementMarker[es[1 - i]] = 1

                #mark T'
                #loop face of T'
                self._putIntoNonBoundaryFaceSet(es[1 - i], nonBoundaryFaceSet)


        self.nodeMap = nodeMap



    def _renumber(self):
        '''
        :return: renumber the element nodes, and put them in newElems
        '''
        nodeMap = self.nodeMap
        eles = self.eles
        newElems = []

        for e in eles: #loop all elements
            #the tetrahedron center is P
            n0,n1,n2,n3 = e
            tetMap = {n0:nodeMap[n0], n1:nodeMap[n1], n2:nodeMap[n2], n3:nodeMap[n3],}
            sortedTetMap = sorted(tetMap.items(), key=operator.itemgetter(1, 0))
            newElems.append([sortedTetMap[0][0], sortedTetMap[1][0], sortedTetMap[2][0], sortedTetMap[3][0]])
        self.newElems = newElems



    def write_topfile(self, outputfile):
        nodes = self.nodes
        eles = self.newElems
        boundaries = self.boundaries
        boundaryNames = self.boundaryNames
        write_tet(nodes, eles, boundaryNames, boundaries,outputfile)
        print('Be careful, all these elements have tag = 0')







if __name__ == '__main__':
    #naca2D()
    naca2D = PreNumber('domain2.top','domain2_iso.top')

