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


def inside_tri(xy, point):
    '''
    :param xy:   numpy array, double[3][2]
    :param point: numpy array, double[2]
    :return: whether the point is in the triangle or not
    '''

    for i in range(3):
        n0,n1 = [n for n in range(3) if n != i]
        alpha = np.cross(xy[n0, :] - point,  xy[n1, :] - point) / np.cross(xy[n0, :] - xy[i, :], xy[n1, :] - xy[i, :])

        if(alpha < 0):
            return False

    return True


def intersect_tri(xy, line):
    '''
    :param xy:   numpy array, double[3][2]
    :param point: numpy array, double[4] x0,y0 --- x1,y1
    :return: whether the linesegment interesects the triangle or not
    '''

    p0, p1 = line[0:2], line[2:4]

    for i in range(3):
        n0,n1 = [n for n in range(3) if n != i]
        alpha = np.cross(p0 - xy[n0,:], p1 - p0)/ np.cross(xy[n1,:] - xy[n0,:], p1 - p0)
        beta = np.cross(xy[n0, :] - p0, xy[n1,:] - xy[n0,:]) / np.cross(p1 - p0, xy[n1, :] - xy[n0, :])

        if(alpha >= 0 and alpha <= 1.0 and beta >= 0 and beta <= 1.0):
            return True

    return False
