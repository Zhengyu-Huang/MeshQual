__author__ = 'zhengyuh'
'''
This file implement the pre refine algorithm in Rob Stevenson , The completion of locally refined simplicial partitions created by bisection
'''
from TopFileTool import write_tet, read_tet
from Utility import pair_sort,triplet_sort

class PreRefine:
    def __init__(self,inputfile = 'domain.top',outputfile = 'domain1.top'):
        self.nodes, self.eles, self.boundaryNames, self.boundaries = read_tet(inputfile)
        self._build_faces()
        self.refine()
        self.write_topfile(outputfile)

    def _build_faces(self):
        eles = self.eles
        faces = {}
        for e in eles:
            for i in range(4):
                enen = [e[j] for j in range(4) if (j != i)]

                faces[triplet_sort(enen[0],enen[1],enen[2])] = -1
        self.faces = faces

    def refine(self):
        '''
        insert a node P on the tetrahedron center, to divide it into three new tetrahedra
        Each child tetrahedron t are divided into 3 tetrahedron by adding a node Q at the original face center

        The P ----> Z1  Q ---> Z2
        Original edge A,B if A < B: A ---> Z0 B ---> Z3
                          if A > B: A ---> Z3 B ---> Z0

        :return:
        '''

        nodes = self.nodes
        boundaries = self.boundaries
        eles = self.eles

        newNodes = self.nodes

        newElems = []
        newBoundaries = [[] for i in range(len(boundaries))]

        for e in eles: #loop all elements
            #the tetrahedron center is P
            n0,n1,n2,n3 = e
            P = [(nodes[n0][0] + nodes[n1][0] + nodes[n2][0] + nodes[n3][0]) / 4.0,
                 (nodes[n0][1] + nodes[n1][1] + nodes[n2][1] + nodes[n3][1]) / 4.0,
                 (nodes[n0][2] + nodes[n1][2] + nodes[n2][2] + nodes[n3][2]) / 4.0]
            newNodes.append(P)
            np = len(newNodes) - 1
            for i in range(4): # loop each face, consider the sub-element containing the face
                #The face, nodes enen, not include ni
                enen = [e[j] for j in range(4) if (j != i)]
                nq = self.faces[triplet_sort(enen[0],enen[1],enen[2])]
                if nq < 0 :

                    Q = [(nodes[enen[0]][0] + nodes[enen[1]][0] + nodes[enen[2]][0]) / 3.0,
                         (nodes[enen[0]][1] + nodes[enen[1]][1] + nodes[enen[2]][1]) / 3.0,
                         (nodes[enen[0]][2] + nodes[enen[1]][2] + nodes[enen[2]][2]) / 3.0]
                    newNodes.append(Q)
                    nq = self.faces[triplet_sort(enen[0], enen[1], enen[2])] = len(newNodes) - 1

                # The face center is Q nq
                for j in range(3): # loop these three edges of the face
                    fnfn = [enen[k] for k in range(3) if (k != j)]
                    ni,nj = pair_sort(fnfn[0], fnfn[1])
                    newElems.append([ni, np, nq, nj])


        for i in range(len(boundaries)):

            for n0,n1,n2 in boundaries[i]:
                nq = self.faces[triplet_sort(n0,n1,n2)]
                newBoundaries[i].append([n0, n1, nq])
                newBoundaries[i].append([n0, n2, nq])
                newBoundaries[i].append([n1, n2, nq])





        self.newElems = newElems
        self.newBoundaries = newBoundaries



    def write_topfile(self, outputfile='domain.top'):
        nodes = self.nodes
        eles = self.newElems
        boundaries = self.newBoundaries
        boundaryNames = self.boundaryNames
        write_tet(nodes, eles, boundaryNames, boundaries)







if __name__ == '__main__':
    #naca2D()
    naca2D = PreRefine('domain.top','Naca.top')

