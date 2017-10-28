__author__ = 'zhengyuh'
import sys

def read_tet(mshfile):
    '''
    This function reads top file, mshfile, extract the node coordinates in the list nodes,
    and element node numbers in elems(the node number is 0-based index instead of top file'
    1-based index)
    :param mshfile: top file name
    :return: nodes, elems
    '''

    try:
        file = open(mshfile, "r")
    except IOError:
        print("File '%s' not found." % mshfile)
        sys.exit()

    nodes = []
    elems = []

    print('Reading mesh ...')

    # Read nodes
    line = file.readline()  # should be "Nodes FluidNodes"
    while line:
        line = file.readline()
        data = line.split()
        if data[0] == 'Elements':
            break;
        nodes.append(list(map(float, data[1:4])))

    # Read Elements
    line = file.readline()
    while line:
        data = line.split()
        if data[0] == 'Elements':
            break;
        elem = list(map(int, data[2:6]))
        elem = [x - 1 for x in elem]  # fixed top file 1-based index
        elems.append(elem)
        line = file.readline()

    file.close()

    return nodes, elems

def write_tet(nodes,elems, boundaries, mshfile = 'domain.top'):
    '''
    This function writes top file, mshfile, node coordinates are  in the list nodes,
    and element node numbers are in elems(the node number is 0-based index instead of top file'
    1-based index)
    :param mshfile: top file name
    :param nodes: a list of node coordinates
    :param elems: a list of elems node number
    '''
    file = open(mshfile, 'w')
    nNodes,nElems = len(nodes), len(elems)

    file.write('Nodes FluidNodes\n')
    for nN in range(nNodes):
        file.write('%d  %.12f  %.12f  %.12f\n'%(nN + 1, nodes[nN][0],nodes[nN][1],nodes[nN][2]))

    file.write('Elements FluidMesh using FluidNodes\n')
    for nE in range(nElems):
        file.write('%d  %d  %d  %d  %d  %d\n'%(nE + 1, 5, elems[nE][0] + 1, elems[nE][1] + 1, elems[nE][2] + 1,elems[nE][3] + 1))

    nBound = len(boundaries)
    for nB in range(nBound):
        file.write('Elements InletFiexedSurface using FluidNodes\n')
        boundary = boundaries[nB]
        nTris = len(boundary)

        for nT in range(nTris):
            nE += 1
            file.write('%d  %d  %d  %d  %d\n' % (nE, 4, boundary[nT][0] + 1, boundary[nT][1] + 1, boundary[nT][2] + 1))

    file.close()

