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
    boundaryNames = []
    boundaries = []

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

    # Read Boundaries
    while line:
        boundary = []
        boundaryName = data[1]
        line = file.readline()
        while line:
            data = line.split()
            if data[0] == 'Elements':
                break;
            tri = list(map(int, data[2:5]))
            tri = [x - 1 for x in tri]  # fixed top file 1-based index
            boundary.append(tri)
            line = file.readline()
        boundaries.append(boundary)
        boundaryNames.append(boundaryName)

    file.close()

    return nodes, elems, boundaryNames, boundaries

def write_tet(nodes,elems, boundaryNames, boundaries, mshfile = 'domain.top'):
    '''
    This function writes top file, mshfile, node coordinates are  in the list nodes,
    and element node numbers are in elems(the node number is 0-based index instead of top file'
    1-based index)
    :param mshfile: top file name
    :param nodes: a list of node coordinates
    :param elems: a list of elems node number
    :param boundaries: a list of several lists, each sublist is a list of boundary triangle node numbers
    '''

    file = open(mshfile, 'w')
    nNodes,nElems = len(nodes), len(elems)

    file.write('Nodes FluidNodes\n')
    for nN in range(nNodes):
        file.write('%d  %.12f  %.12f  %.12f\n'%(nN + 1, nodes[nN][0],nodes[nN][1],nodes[nN][2]))

    file.write('Elements Volume_0 using FluidNodes\n')
    for nE in range(nElems):
        file.write('%d  %d  %d  %d  %d  %d\n'%(nE + 1, 5, elems[nE][0] + 1, elems[nE][1] + 1, elems[nE][2] + 1,elems[nE][3] + 1))

    nBounds = len(boundaries)
    for i in range(nBounds):

        file.write('Elements %s using FluidNodes\n' %(boundaryNames[i] + '_' + str(i+1)))
        boundary = boundaries[i]
        nTris = len(boundary)
        for nT in range(nTris):
            nE += 1
            file.write('%d  %d  %d  %d  %d\n' % (nE, 4, boundary[nT][0] + 1, boundary[nT][1] + 1, boundary[nT][2] + 1))

    file.close()

