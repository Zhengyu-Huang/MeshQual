__author__ = 'zhengyuh'
import numpy as np

def KuhnSimplex(xx,yy,zz):
    '''
    :param xx:
    :param yy:
    :param zz:
    cut each small cube into 6 tetrahedrons
    T(1 2 3） (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
    T(1 3 2） (0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1)
    T(2 1 3) (0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)
    T(3 1 2) (0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)
    T(3 2 1) (0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)
    T(2 3 1) (0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)
    :return:
    '''

    nx,ny,nz = len(xx),len(yy),len(zz)

    #number the nodes
    #number the elems
    #number the 

