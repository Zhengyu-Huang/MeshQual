__author__ = 'zhengyuh'
import numpy as np

def pair_sort(a,b):
    '''
    :param a: int
    :param b: int
    :return: sorted pair
    '''
    return (a,b) if a < b else (b,a)

def triplet_sort(a,b,c):
    '''
    :param a: int
    :param b: int
    :param c: int
    :return: sorted triplet
    '''
    a,b = pair_sort(a,b)
    b,c = pair_sort(b,c)
    a,b = pair_sort(a,b)
    return (a,b,c)

def inside_tet(xyz,point):
    '''
    :param xyz:   numpy array, double[4][3]
    :param point: numpy array, double[3]
    :return: whether the point is in the tetrahedron or not
    '''

    for i in range(4):
        n0,n1,n2 = [n for n in range(4) if n != i]
        alpha = np.dot(xyz[n0, :] - point, np.cross(xyz[n1, :] - point, xyz[n2, :] - point))/\
                np.dot(xyz[n0, :] - xyz[i, :], np.cross(xyz[n1, :] - xyz[i, :], xyz[n2, :] - xyz[i, :]))

        if(alpha < 0):
            return False

    return True