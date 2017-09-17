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

def write_tet(mshfile,nodes,elems):
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
    for i in range(nNodes):
        file.write('%d  %.12f  %.12f  %.12f\n'%(i+1, nodes[i][0],nodes[i][1],nodes[i][2]))

    file.write('Elements FluidMesh using FluidNodes\n')
    for i in range(nElems):
        file.write('%d  %d  %d  %d  %d  %d\n'%(i+1, 5, elems[i][0] + 1, elems[i][1] + 1, elems[i][2] + 1,elems[i][3] + 1))

    file.close()

