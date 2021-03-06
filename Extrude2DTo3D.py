import copy as copy
import numpy as np
import TopFileTool



def Extrude2DTo3D(nlayers, zmin, zmax, extrudeBoundaryName = "SymmetrySurface", inputfile = "domain.top", outputfile="domain.3D.top"):


    nodes2D, elems2D, boundaryNames2D, boundaries2D = TopFileTool.read_tri(inputfile)
    nnodes2D, nelems2D, nboundaries2D = len(nodes2D), len(elems2D), len(boundaryNames2D)

    dz = (zmax - zmin) / nlayers


    # build boundary map
    boundaryNodes2D = [set() for i in range(nboundaries2D)]

    for i in range(nboundaries2D):
        boundary = boundaries2D[i]
        for n0, n1 in boundary:
            print(n0, n1)
            boundaryNodes2D[i].add(n0)
            boundaryNodes2D[i].add(n1)

    # bulid new 3D mesh, extude to prism and cut each to 3
    nodes3D = np.zeros((nnodes2D * (nlayers+1), 3), dtype=float)
    for i in range(nlayers+1):
        nodes3D[i * nnodes2D:(i + 1) * nnodes2D, 0:2] = nodes2D
        nodes3D[i * nnodes2D:(i + 1) * nnodes2D, 2] = i*dz + zmin

    elems3D = np.zeros((nelems2D * (nlayers) * 3, 4), dtype=int)

    boundaries3D = [[] for i in range(nboundaries2D+1)] # last one is for top and bottom boundaries
    boundaryNames3D = copy.copy(boundaryNames2D)
    boundaryNames3D.append(extrudeBoundaryName )

    nelems3D = 0
    v1, v2, v3 = np.zeros(4, dtype=int),np.zeros(4, dtype=int),np.zeros(4, dtype=int)
    for elem in elems2D:

        if (elem[0] < elem[1] and elem[0] < elem[2]):
            VI = [elem[0], elem[1], elem[2], elem[0]+nnodes2D, elem[1]+nnodes2D, elem[2]+nnodes2D]
        if (elem[1] < elem[0] and elem[1] < elem[2]):
            VI = [elem[1], elem[2], elem[0], elem[1]+nnodes2D, elem[2]+nnodes2D, elem[0]+nnodes2D]
        if (elem[2] < elem[0] and elem[2] < elem[1]):
            VI = [elem[2], elem[0], elem[1], elem[2]+nnodes2D, elem[0]+nnodes2D, elem[1]+nnodes2D]


        if (min(VI[1], VI[5]) < min(VI[2], VI[4])):
            v1[:] = [VI[0], VI[1], VI[2], VI[5]]
            v2[:] = [VI[0], VI[1], VI[5], VI[4]]
            v3[:] = [VI[0], VI[4], VI[5], VI[3]]
        else:
            v1[:] = [VI[0], VI[1], VI[2], VI[4]]
            v2[:] = [VI[0], VI[4], VI[2], VI[5]]
            v3[:] = [VI[0], VI[4], VI[5], VI[3]]

        e1, e2, e3 = nelems3D + 0, nelems3D + 1, nelems3D + 2

        for ilayer in range(nlayers):
            elems3D[e1 + 3*nelems2D*ilayer,:] = v1 + ilayer*nnodes2D
            elems3D[e2 + 3*nelems2D*ilayer,:] = v2 + ilayer*nnodes2D
            elems3D[e3 + 3*nelems2D*ilayer,:] = v3 + ilayer*nnodes2D

        nelems3D = nelems3D + 3



    # construct boundary face
    for elem in elems3D:
        for ibc in range(len(boundaryNames2D)):
            for iface in range(4):
                n0, n1, n2 = [elem[i] for i in range(4) if i != iface]
                if(n0%nnodes2D in boundaryNodes2D[ibc]  and  n1%nnodes2D in boundaryNodes2D[ibc]  and n2%nnodes2D in boundaryNodes2D[ibc]):
                    boundaries3D[ibc].append([n0,n1,n2])

    for elem in elems2D:
        boundaries3D[-1].append(elem)
        boundaries3D[-1].append([elem[0]+nnodes2D*nlayers, elem[1]+nnodes2D*nlayers, elem[2]+nnodes2D*nlayers])




    TopFileTool.write_tet(nodes3D,elems3D, boundaryNames3D, boundaries3D, outputfile)



if __name__ == "__main__":
    nlayers = 2

    zmin, zmax = 0.0, 1.0

    extrudeBoundaryName = "SymmetrySurface"

    inputfile, outputfile = "domain.refined.top", "domain.3D.top"

    Extrude2DTo3D(nlayers, zmin, zmax, extrudeBoundaryName, inputfile, outputfile)


