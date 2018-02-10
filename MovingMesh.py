'''
This file move the mesh toward the interface

First test 2D case
'''

import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Mesh2D:
    def __init__(self, top_file = None):
        '''
        nodes_:
        edges_:
        elem_:
        boundary_:
        '''

        edges_ = []
        elems_ = []

        if(top_file):
            print('have not implemented yet')
        else:
            '''
            structure mesh (-a,-a),(-a,a),(a,-a),(a,a)
            n nodes on each direction
            '''
            ##Parameters
            a = 2.0
            n = 21
            nodes_ = np.empty([n*n,2])
            xx, yy= np.meshgrid(np.linspace(-a,a,n), np.linspace(-a,a,n))
            nodes_[:, 0] = np.reshape(xx, (1, -1), order='F')
            nodes_[:, 1] = np.reshape(yy, (1, -1), order='F')

            # y direction
            for i in range(n-1):
                for j in range(n):
                    edges_.append([n*i+j, n*(i+1)+j])

            # x direction
            for j in range(n - 1):
                for i in range(n):
                    edges_.append([n * i + j, n * i + j + 1])

            for i in range(n-1):
                for j in range(n-1):
                    edges_.append([n*i+j, n*(i+1) + j + 1])
                    elems_.extend([[n*(i+1) + j+1, n*(i+1)+j, n*i+j],[n*(i+1) + (j+1), n*i+j, n*i+j+1]])

            boundary_ = np.zeros(n*n)
            boundary_[0:n] = 1
            boundary_[0::n] = 1
            boundary_[n-1::n] = 1
            boundary_[n*(n-1):] = 1

        self.edges_ = np.array(edges_, dtype=int)
        self.nodes_ = np.array(nodes_)
        self.elems_ = np.array(elems_, dtype=int)
        self.boundary_ = boundary_
        self.n_edges_ = len(edges_)
        self.n_nodes_ = len(nodes_)


    def visual(self, nodes = None):
        if nodes is None:
            nodes = self.nodes_

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        ax1.triplot(nodes[:,0], nodes[:,1], self.elems_, 'go-', lw=0.1)
        ax1.set_title('Mesh')


        ax1.add_patch(
            patches.Circle(
                (0.0, 0.0),
                0.81,
                fill=False  # remove background
            ))

        plt.show()



def DistanceToWall(nodes, case):
    n_nodes,dim = nodes.shape
    dis_to_wall, grand_dis_to_wall = np.zeros(n_nodes), np.zeros([n_nodes,2])
    if(case == 1):
        '''
        The structure is a circle centered as (0,0) with radius r
        '''
        r = 0.81
        for i in range(n_nodes):
            x = nodes[i,:]
            d = np.linalg.norm(x) - r

            if(np.linalg.norm(x) > 0.0):
                grand_dis_to_wall[i,:] = -x*np.sign(d)/np.linalg.norm(x)
            dis_to_wall[i] = np.fabs(d)



    return dis_to_wall, grand_dis_to_wall



def MovingMesh(mesh,case):
    '''
    The governing equation is
    :param mesh:
    :return:
    '''

    n_edges = mesh.n_edges_

    nodes0 = mesh.nodes_

    edges = mesh.edges_


    edges_len = np.empty(n_edges)

    boundary = mesh.boundary_

    inner = boundary == 0

    n_nodes, dim = mesh.nodes_.shape

    k = 1.0
    eps = 0.1
    dt = 0.1
    TOL = 1e-3
    MAXITE = 10

    for i in range(n_edges):
        n1,n2 = edges[i]
        edges_len[i] = np.linalg.norm(nodes0[n1,:] - nodes0[n2,:])

    f_int,f_ext = np.empty([n_nodes,dim]), np.empty([n_nodes,dim])

    nodes = copy.copy(nodes0)

    for ite in range(MAXITE):

        dis_to_wall, grad_dis_to_wall = DistanceToWall(nodes, case)
        f_int[:],f_ext[:] = 0,0


        for i in range(n_edges):
            n1,n2 = edges[i]
            x1,x2 = nodes[n1,:], nodes[n2,:]
            df = k*(np.linalg.norm(x1-x2) - edges_len[i])*(x2 - x1)

            f_int[n1,:] = f_int[n1,:] + df
            f_int[n2,:] = f_int[n2,:] - df


        f_ext = eps*grad_dis_to_wall

        nodes[inner] = nodes[inner] + dt*(f_ext[inner] + f_int[inner])


        if(np.linalg.norm(f_ext[inner] + f_int[inner]) < TOL):
            break

    return nodes



if __name__ == '__main__':
    mesh = Mesh2D()
    #mesh.visual()
    nodes = MovingMesh(mesh,1)



    mesh.visual(nodes)
