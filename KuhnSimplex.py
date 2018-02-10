__author__ = 'zhengyuh'
import numpy as np
from TopFileTool import write_tet


class KuhnSimplex:
    def __init__(self,x,y,z, boundaryNames):
        self.x,  self.y,  self.z = x,y,z
        self.nx, self.ny, self.nz = len(x), len(y), len(z)
        self.nodes = self.create_nodes()
        self.eles = self.create_tet()
        self.boundaries = self.create_boundaries()
        self.boundaryNames = boundaryNames


    def _node_id(self,ix,iy,iz):
        nx, ny, nz = self.nx, self.ny, self.nz
        return ix + iy*nx + iz*nx*ny

    def create_nodes(self):
        x,y,z = self.x, self.y, self.z
        nx, ny, nz = self.nx, self.ny, self.nz


        nodes = np.empty((nx*ny*nz, 3))

        xx,yy,zz = np.meshgrid(x,y,z)

        nodes[:, 0] = np.reshape(np.reshape(xx, (-1, nz)), (1, -1), order='F')
        nodes[:, 1] = np.reshape(np.reshape(yy, (-1, nz)), (1, -1), order='F')
        nodes[:, 2] = np.reshape(np.reshape(zz, (-1, nz)), (1, -1), order='F')

        return nodes





    def create_tet(self):
        '''
        :param xx:
        :param yy:
        :param zz:
        cut each small cube into 6 tetrahedrons

        nodes = [n0,n1,n2,n3,n4,n5,n6,n7]

        n0: (0,0,0)
        n1: (1,0,0)
        n2: (0,1,0)
        n3: (1,1,0)
        n4: (0,0,1)
        n5: (1,0,1)
        n6: (0,1,1)
        n7: (1,1,1)

        T(1 2 3) n0(0, 0, 0), n1(1, 0, 0), n3(1, 1, 0), n7(1, 1, 1)
        T(1 3 2) n0(0, 0, 0), n1(1, 0, 0), n5(1, 0, 1), n7(1, 1, 1)
        T(2 1 3) n0(0, 0, 0), n2(0, 1, 0), n3(1, 1, 0), n7(1, 1, 1)
        T(3 1 2) n0(0, 0, 0), n2(0, 1, 0), n6(0, 1, 1), n7(1, 1, 1)
        T(3 2 1) n0(0, 0, 0), n4(0, 0, 1), n6(0, 1, 1), n7(1, 1, 1)
        T(2 3 1) n0(0, 0, 0), n4(0, 0, 1), n5(1, 0, 1), n7(1, 1, 1)

        :return:
        '''
        nx, ny, nz = self.nx, self.ny, self.nz

        eles = np.empty((6*(nx - 1) * (ny - 1) * (nz - 1), 4), dtype=int)
        for k in range(nz - 1):
            for j in range(ny - 1):
                for i in range(nx - 1):
                    cubeId = i + j*(nx-1) + k*(nx-1)*(ny-1)
                    nn = [self._node_id(i, j, k), self._node_id(i + 1, j, k), self._node_id(i, j + 1, k), self._node_id(i + 1, j + 1, k),
                          self._node_id(i, j, k + 1), self._node_id(i + 1, j, k + 1), self._node_id(i, j + 1, k + 1), self._node_id(i + 1, j + 1, k + 1)]
                    eles[6 * cubeId:6 * (cubeId + 1), :] = [[nn[0], nn[1], nn[3], nn[7]],
                                                  [nn[0], nn[1], nn[5], nn[7]],
                                                  [nn[0], nn[2], nn[3], nn[7]],
                                                  [nn[0], nn[2], nn[6], nn[7]],
                                                  [nn[0], nn[4], nn[6], nn[7]],
                                                  [nn[0], nn[4], nn[5], nn[7]]]

                    # eles[6 * cubeId:6 * (cubeId + 1), :] = [[nn[0], nn[1], nn[3], nn[7]],
                    #                               [nn[0], nn[5], nn[1], nn[7]],
                    #                               [nn[0], nn[3], nn[2], nn[7]],
                    #                               [nn[0], nn[2], nn[6], nn[7]],
                    #                               [nn[0], nn[6], nn[4], nn[7]],
                    #                               [nn[0], nn[4], nn[5], nn[7]]]

        return eles

    def create_boundaries(self):
        nx, ny, nz = self.nx, self.ny, self.nz
        tri = [[],[],[],[],[],[]]
        # bottom z = zmin
        for j in range(ny - 1):
            for i in range(nx - 1):
                k = 0
                tri[0].append([self._node_id(i, j, k), self._node_id(i + 1, j, k), self._node_id(i + 1, j + 1, k)])
                tri[0].append([self._node_id(i, j, k), self._node_id(i + 1, j + 1, k), self._node_id(i, j + 1, k)])


        # top  z = zmax
        for j in range(ny - 1):
            for i in range(nx - 1):
                k = nz - 1
                tri[1].append([self._node_id(i, j, k), self._node_id(i + 1, j, k), self._node_id(i + 1, j + 1, k)])
                tri[1].append([self._node_id(i, j, k), self._node_id(i + 1, j + 1, k), self._node_id(i, j + 1, k)])

        # left x = xmin
        for k in range(nz - 1):
            for j in range(ny - 1):
                i = 0
                tri[2].append([self._node_id(i, j, k), self._node_id(i , j + 1, k),     self._node_id(i, j + 1, k + 1)])
                tri[2].append([self._node_id(i, j, k), self._node_id(i , j + 1, k + 1), self._node_id(i, j , k + 1)])
        # right x = xmax
        for k in range(nz - 1):
            for j in range(ny - 1):
                i = nx - 1
                tri[3].append([self._node_id(i, j, k), self._node_id(i, j + 1, k), self._node_id(i, j + 1, k + 1)])
                tri[3].append([self._node_id(i, j, k), self._node_id(i, j + 1, k + 1), self._node_id(i, j, k + 1)])
        # front y = ymin
        for k in range(nz - 1):
            for i in range(nx - 1):
                j = 0
                tri[4].append([self._node_id(i, j, k), self._node_id(i, j , k + 1), self._node_id(i + 1, j, k + 1)])
                tri[4].append([self._node_id(i, j, k), self._node_id(i + 1, j, k + 1), self._node_id(i + 1, j, k)])
        # back y = ymax
        for k in range(nz - 1):
            for i in range(nx - 1):
                j = ny - 1
                tri[5].append([self._node_id(i, j, k), self._node_id(i, j, k + 1), self._node_id(i + 1, j, k + 1)])
                tri[5].append([self._node_id(i, j, k), self._node_id(i + 1, j, k + 1), self._node_id(i + 1, j, k)])


        return tri

    def write_topfile(self, outputfile = 'domain.top', volFunc = lambda x: 0):
        nodes = self.nodes
        eles = self.eles
        boundaries = self.boundaries
        boundaryNames = self.boundaryNames
        write_tet(nodes, eles, boundaryNames, boundaries, outputfile, volFunc)


    def plot_mesh(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        eles  = self.eles
        nodes = self.nodes
        for e in eles:
            for i in range(4):
                for j in range(i+1,4):
                    x = [nodes[e[i],0],nodes[e[j],0]]
                    y = [nodes[e[i],1],nodes[e[j],1]]
                    z = [nodes[e[i],2],nodes[e[j],2]]
                    ax.plot(x,y,z)
        plt.show()






def geomspace(x0, xn, dx0, ratio, includeX0 = True):
    xx = [x0] if includeX0 else []
    x, dx = x0, dx0

    while (x0 < xn and x < xn) or (x0 > xn and x > xn):
        dx *= ratio
        x = x + dx
        xx.append(x)



    return xx


def parachute3D():
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    xcl,xcr = -15.0, 15.0
    xc = np.linspace(xcl,xcr, 61)
    dxc = (xcr - xcl)/60.0
    xRatio = 1.4
    xl = geomspace(xcl, -150, -dxc, xRatio, False)
    xr = geomspace(xcr, 150, dxc, xRatio, False)


    x = [*reversed(xl),*xc,*xr]


    y = x
    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-100 200]
    zcl, zcr = 25.0, 50.0
    zc = np.linspace(zcl, zcr, 51)
    dzc = (zcr - zcl)/50.0
    zRatio = 1.4
    zl = geomspace(zcl, -100, -dzc, zRatio, False)
    zr = geomspace(zcr, 150, dzc, zRatio, False)
    z = [*reversed(zl), *zc, *zr]

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    simpleKuhnSimplex = KuhnSimplex(x,y,z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()

def naca2D():

    xcl,xcr = -1.0/3.0, 1.5-1.0/3.0
    xc = np.linspace(xcl, xcr, 151)
    dxc = xc[1] - xc[0]
    xRatio = 1.5
    xl = geomspace(xcl, -20, -dxc, xRatio, False)
    xr = geomspace(xcr, 20, dxc, xRatio, False)
    x = [*reversed(xl), *xc, *xr]

    ycl, ycr = -1./3., 1./3.
    yc = np.linspace(ycl, ycr, 68)
    dyc = yc[1] - yc[0]
    yRatio = 1.5
    yl = geomspace(ycl, -10, -dyc, yRatio, False)
    yr = geomspace(ycr,  10, dyc, yRatio, False)
    y = [*reversed(yl), *yc, *yr]

    z = [0.0, 0.005]
    boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()














'''
def parachute3D():
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    xcl,xcr = -15.2, 15.2
    xc = np.linspace(xcl,xcr, 10)
    dxc = (xcr - xcl)/60.0
    xRatio = 1.5
    xl = geomspace(xcl, -80, -dxc, xRatio, False)
    xr = geomspace(xcr, 80, dxc, xRatio, False)


    x = [*reversed(xl),*xc,*xr]


    y = x
    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-100 200]
    zcl, zcr = 25.2, 50.2
    zc = np.linspace(zcl, zcr, 6)
    dzc = (zcr - zcl)/50.0
    zRatio = 1.5
    zl = geomspace(zcl, -80, -dzc, zRatio, False)
    zr = geomspace(zcr, 80, dzc, zRatio, False)
    z = [*reversed(zl), *zc, *zr]

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    simpleKuhnSimplex = KuhnSimplex(x,y,z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()
'''



def Tube():
    #Important the name of these boundaries are important
    x = np.linspace(-5.0, 5.0, 100)
    y = np.linspace(-0.2, 0.2, 4)
    z = np.linspace(-0.2, 0.2, 4)
    boundaryNames = ['SymmetrySurface', 'SymmetrySurface', 'InletFixedSurface',
                     'InletFixedSurface', 'SymmetrySurface', 'SymmetrySurface']
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()

def shockTube():
    #Important the name of these boundaries are important
    x = np.linspace(-5.0, 5.0, 400)
    y = np.linspace(-0.2, 0.2, 4)
    z = np.linspace(-0.2, 0.2, 4)
    boundaryNames = ['SymmetrySurface', 'SymmetrySurface', 'InletFixedSurface',
                     'OutletFixedSurface', 'SymmetrySurface', 'SymmetrySurface']
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile(volFunc=lambda xc: 0 if xc[0] < 1.0e-6 else 1)



def uniform2D():
    # for circle
    # x = np.linspace(-2.0, 8.0, 100*2)
    # y = np.linspace(-4.0, 4.0, 80*2)

    # for NACA
    x = np.linspace(-25.0, 25.0, 100)
    y = np.linspace(-25.0, 25.0, 100)



    z = [0.0, 0.001]
    boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


def uniform3D():
    # for sphere
    x = np.linspace(-40.0, 80.0, 3*23)
    y = np.linspace(-40.0, 40.0, 2*23)
    z = np.linspace(-40.0, 40.0, 2*23)



    boundaryNames=['InletFixedSurface','InletFixedSurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


def sphere3D_trial1():
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    xcl, xcr, xcn = -2.0, 4.0, 52
    xc = np.linspace(xcl, xcr, xcn)
    dxc = (xcr - xcl) /(xcn - 1)
    xRatio = 1.5
    xl = geomspace(xcl, -40, -dxc, xRatio, False)
    xr = geomspace(xcr, 80, dxc, xRatio, False)

    x = [*reversed(xl), *xc, *xr]


    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-100 200]
    ycl, ycr, ycn = -2.0, 2.0, 35
    yc = np.linspace(ycl, ycr, ycn)
    dyc = (ycr - ycl) / (ycn - 1)
    yRatio = 1.5
    yl = geomspace(ycl, -40, -dyc, yRatio, False)
    yr = geomspace(ycr, 40, dyc, yRatio, False)
    y = [*reversed(yl), *yc, *yr]

    z = y

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()

def sphere3D():
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    xcl, xcr, xcn = -10.0, 20.0, 52
    xc = np.linspace(xcl, xcr, xcn)
    dxc = (xcr - xcl) /(xcn - 1)
    xRatio = 1.2
    xl = geomspace(xcl, -40, -dxc, xRatio, False)
    xr = geomspace(xcr, 80, dxc, xRatio, False)

    x = [*reversed(xl), *xc, *xr]


    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-100 200]
    ycl, ycr, ycn = -10.0, 10.0,35
    yc = np.linspace(ycl, ycr, ycn)
    dyc = (ycr - ycl) / (ycn - 1)
    yRatio = 1.2
    yl = geomspace(ycl, -40, -dyc, yRatio, False)
    yr = geomspace(ycr, 40, dyc, yRatio, False)
    y = [*reversed(yl), *yc, *yr]

    z = y

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()

if __name__ == '__main__':
    #uniform2D()
    #shockTube()
    #parachute3D()
    #uniform2D()
    shockTube()
    # x = [1,2]
    # y = [4,5]
    # z = [7,8]
    #
    # simpleKuhnSimplex = KuhnSimplex(x,y,z)
    #
    # print('Writing to top file')
    # simpleKuhnSimplex.plot_mesh()
