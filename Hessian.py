import sys
import numpy as np
import matplotlib.pyplot as plt



'''
The mesh has two dual cells, like
      0*----1*------2*
    e0  e1 e2  e3    e4
 3*      4*      5*       6*
    e5  e6  e7  e8  e9
    7*     8*       9*
dual cells are at node 5 and node 6



'''

class simpleMesh:
    def __init__(self,u,xy):
        self.dim = 2
        self.nElems = 10
        self.nNodes = 10
        self.elems = np.array([[0,3,4],[0,4,1],[1,4,5],[1,5,2],[2,5,6],[3,7,4],[4,7,8],[4,8,5],[5,8,9],[5,9,6]],dtype=int)
        self.nToElems = np.array([[0,1,2,7,6,5],[2,3,4,9,8,7]],dtype=int)
        self.nToNodes = np.array([[0,1,5,8,7,3],[1,2,6,9,8,4]],dtype=int)
        self.nc = np.array([4,5],dtype=int)


        self.u = u
        self.xy = xy
        self.xy_c = (xy[self.nc[0],:]  + xy[self.nc[1],:])/2.0
        self.e_c = xy[self.nc[0], :] - xy[self.nc[1], :]

        self._compute_volume()




    def _volume(self,e):
        n1, n2, n3 = self.elems[e, :]
        xy = self.xy
        x1, y1 = xy[n1, :]
        x2, y2 = xy[n2, :]
        x3, y3 = xy[n3, :]
        return np.fabs(0.5*((x1 - x3)*(y2 - y3) - (x2 - x3)*(y1 - y3)))

    def _compute_volume(self):
        nElems = self.nElems
        vol = np.empty(nElems)
        for e in range(nElems):
            vol[e] = self._volume(e)

        self.volumes = vol

        dual_vol = np.zeros(2)
        for i in range(2):
            for e in self.nToElems[i]:
                dual_vol[i] += vol[e]/3.0
        self.dual_volumes = dual_vol

    def plot_mesh(self):
        xy = self.xy
        xy_c = self.xy_c

        plt.figure('mesh')
        plt.plot(xy[:,0], xy[:,1],'ro')
        plt.plot(xy_c[0],xy_c[1],'bo')
        for n1,n2,n3 in self.elems:
            plt.plot([xy[n1,0],xy[n2,0]],[xy[n1,1],xy[n2,1]],'-')
            plt.plot([xy[n2,0],xy[n3,0]],[xy[n2,1],xy[n3,1]],'-')
            plt.plot([xy[n3,0],xy[n1,0]],[xy[n3,1],xy[n1,1]],'-')
        plt.show()




    def least_square_gradient_hessian(self):
        dim = self.dim
        nToNodes = self.nToNodes
        gradient = np.zeros((2, dim))
        xy = self.xy
        nc = self.nc
        u = self.u
        for i in range(2):
            A = np.zeros((dim,dim))
            b = np.zeros(dim)

            for n in nToNodes[i]:
                d_xy = xy[n,:] - xy[nc[i],:]
                d_u = u[n] - u[nc[i]]

                A += np.outer(d_xy, d_xy)
                b += d_u*d_xy


            gradient[i,:] = np.linalg.solve(A,b)

        print('gradient from least square is \n',gradient)

        return np.dot(gradient[0,:] - gradient[1,:],xy[nc[0],:] - xy[nc[1],:])


    def _ele_jacobian(self, e):
        '''
        The baseis for iso parametric map are N1 = s, N2 = t, N3 = 1 - s - t
        x = N1 x1 + N2 x2 + N3 x3
        y = N1 y1 + N2 y2 + N3 y3

        jac = d(x,y)/d(s,t) = x1 - x3, x2 - x3
                              y1 - y3, y2 - y3
        :param e: int
        :return: jac of element e
        '''

        n1, n2, n3 = self.elems[e,:]
        xy = self.xy
        x1, y1 = xy[n1,:]
        x2, y2 = xy[n2, :]
        x3, y3 = xy[n3, :]
        return np.array([[x1-x3,x2-x3],[y1-y3,y2-y3]])

    def _dN_ds(self):
        '''
        The baseis for iso parametric map are N1 = s, N2 = t, N3 = 1 - s - t
        dN/ds =  1  0
                 0  1
                -1 -1
        :param e:
        :return:
        '''
        return np.array([[1.0, 0.0],[0.0, 1.0],[-1.0, -1.0]])


    def _gradient(self, e):
        '''

        :param e: int
        :return: the gradient of element e
        '''

        u  = self.u[self.elems[e,:]]
        dN_ds = self._dN_ds()
        ds_dx = np.linalg.inv(self._ele_jacobian(e))
        return np.dot(np.dot(u, dN_ds),ds_dx)




    def nodal_galerkin_gradient_hessian(self):
        nToElems = self.nToElems
        dim = self.dim
        xy = self.xy
        nc = self.nc
        vol = self.volumes
        dual_vol = self.dual_volumes

        gradient = np.zeros((2,dim))
        for i in range(2):
            for e in nToElems[i]:
                gradient[i,:] += (self._gradient(e))*vol[e]
            #weighted-average by volume
            gradient[i,:] /= (3.0*dual_vol[i])

        print('gradient from nodal galerkin is \n', gradient)

        return np.dot(gradient[0, :] - gradient[1, :], xy[nc[0], :] - xy[nc[1], :])

    def _hessian(self,i):
        '''
        :param i: int 0 or 1
        :return: the hessian matrix at node self.nc[i]
        '''
        dim = self.dim

        dual_vol = self.dual_volumes[i]

        vol = self.volumes

        H = np.zeros((dim,dim))

        for e in self.nToElems[i]:
            du = self._gradient(e)


            dN_ds = self._dN_ds()
            dphi_ds = np.zeros(2)
            for k in range(3):
                if self.nc[i] == self.elems[e,k]:
                    dphi_ds = dN_ds[k,:]

            ds_dx = np.linalg.inv(self._ele_jacobian(e))
            dphi = np.dot(dphi_ds,ds_dx)

            H += vol[e]*np.outer(du,dphi)
        H /= -dual_vol

        return H


    def fem_hessian(self):
        hessian = []
        for i in range(2):
            hessian.append(self._hessian(i))

        e_c = self.e_c

        hessian_c = 0.5*(hessian[0] + hessian[1])

        print('fem hessian at edge center is \n', hessian_c)

        return np.dot(e_c, np.dot(hessian_c, e_c))



def fun(xy):
    return xy[:,0]**3 + xy[:,1]**2

def dfun(xy):
    return 2*xy[:,0]**1

def ddfun(xy):
    x,y = xy
    return np.array([[2.0, 0.0],[0.0, 2.0]])



if __name__ == "__main__":

    xy = np.array([[-0.5,0.5],[0.0, 0.5],[0.5,0.5],
                   [-1.0,0.0],[-0.3,0.0],[0.2,0.0],[1.0,0.0],
                   [-0.5,-0.5],[0.0,-0.5],[0.6,-0.5]])

    u = fun(xy)

    mesh = simpleMesh(u,xy)

    mesh.plot_mesh()
    hessian_lsq = mesh.least_square_gradient_hessian()
    hessian_ng = mesh.nodal_galerkin_gradient_hessian()
    hessian_fem = mesh.fem_hessian()
    print('hessian_lsq: ', hessian_lsq, 'hessian_ng: ', hessian_ng, 'hessian_fem: ', hessian_fem)

    xy_c = mesh.xy_c
    e_c = mesh.e_c
    hessian = np.dot(e_c, np.dot(ddfun(xy_c),e_c))

    print('exact hessian is ', hessian)

    gradient = dfun(xy[mesh.nc])

    print('exact gradient is ', gradient)




















