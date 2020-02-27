__author__ = 'zhengyuh'
import numpy as np
from TopFileTool import write_tet, write_tri


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




class KuhnSimplex2D:
    def __init__(self,x,y, boundaryNames):
        self.x,  self.y = x,y
        self.nx, self.ny = len(x), len(y)
        self.nodes = self.create_nodes()
        self.eles = self.create_tri()
        self.boundaries = self.create_boundaries()
        self.boundaryNames = boundaryNames


    def _node_id(self,ix,iy):
        nx, ny = self.nx, self.ny
        return ix + iy*nx

    def create_nodes(self):
        x,y = self.x, self.y
        nx, ny = self.nx, self.ny


        nodes = np.empty((nx*ny, 2))

        xx,yy = np.meshgrid(x,y)

        nodes[:, 0] = np.reshape(xx, (1, -1))
        nodes[:, 1] = np.reshape(yy, (1, -1))


        return nodes





    def create_tri(self):
        '''
        cut each small quad into 2 triangles

        nodes = [n0,n1,n2,n3]

        n0: (0,0)
        n1: (1,0)
        n2: (0,1)
        n3: (1,1)

        T(0, 1, 2)
        T(1, 3, 2)


        :return:
        '''
        nx, ny = self.nx, self.ny

        eles = np.empty((2*(nx - 1) * (ny - 1), 3), dtype=int)

        for j in range(ny - 1):
            for i in range(nx - 1):
                quadId = i + j*(nx-1)
                nn = [self._node_id(i, j), self._node_id(i + 1, j), self._node_id(i, j + 1), self._node_id(i + 1, j + 1)]
                eles[2 * quadId:2 * (quadId + 1), :] = [[nn[0], nn[1], nn[2]],
                                                        [nn[3], nn[2], nn[1]]]


        return eles

    def create_boundaries(self):
        nx, ny = self.nx, self.ny
        bc = [[],[],[],[]]

        # left x = xmin

        for j in range(ny - 1):
            i = 0
            bc[0].append([self._node_id(i, j), self._node_id(i , j + 1)])
        # right x = xmax
        for j in range(ny - 1):
            i = nx - 1
            bc[1].append([self._node_id(i, j), self._node_id(i, j + 1)])
        # front y = ymin
        for i in range(nx - 1):
            j = 0
            bc[2].append([self._node_id(i, j), self._node_id(i + 1, j)])

        # back y = ymax
        for i in range(nx - 1):
            j = ny - 1
            bc[3].append([self._node_id(i, j), self._node_id(i + 1, j)])
        return bc


    def write_topfile(self, outputfile = 'domain.top', volFunc = lambda x: 0):
        nodes = self.nodes
        eles = self.eles
        boundaries = self.boundaries
        boundaryNames = self.boundaryNames
        write_tri(nodes, eles, boundaryNames, boundaries, outputfile, volFunc)



# def geomspace(x0, xn, dx0, ratio, includeX0 = True):
#     xx = [x0] if includeX0 else []
#     x, dx = x0, dx0
#
#     while (x0 < xn and x < xn) or (x0 > xn and x > xn):
#         x = x + dx
#         xx.append(x)
#         dx *= ratio
#
#
#
#     return xx


def geomspace(x0, xn, dx0, includeX0, type, ratio):
    '''
    :param x0: float, start point
    :param xn: float, end point (end point is beyond xn)
    :param dx0: float, increment at x0
    :param includeX0, bool, should the array include x0 or not
    :param type: string 'num' or 'dist'
    :param ratio: list of float
           [float, int; float, int, float ...]
           [float float; float, float; , float ...]
    :return: xx
           ratio has n pairs(ci, yi), y of the last pair is optional
           xx array has n segments, in each segment the increment dx
           is a geometric sequence of ratio ci
           if type = 'num', each segment has yi points,
           the ratio is interpreted as c1, n1, c2, n2 ...
           the incremental array is dx, dx*c1, dx*c1^2 ... dx*c1^(n1-1),
           dx*c1^(n1-1)*c2, dx*c1^(n1-1)*c2^2 ...dx*c1^(n1-1)*c2^(n2-1)...
           if type = 'dist', each segment has length > yi - yi-1
           the ratio is interpreted as c1, n1, c2, n2 ...
           the incremental array is dx, dx*c1, dx*c1^2 ... dx*c1^(n1-1),
           dx*c1^(n1-1)*c2, dx*c1^(n1-1)*c2^2 ...dx*c1^(n1-1)*c2^(n2-1)...
           ni in this case is computed based on yi
    '''

    n = (len(ratio) + 1)//2
    xx = [x0] if includeX0 else []

    if(type == 'num'):
        x, dx = x0, dx0
        if len(ratio) % 2:
            ratio.append(int(1e12))#append any large number
        for i in range(n):
            for j in range(ratio[2*i + 1]):
                x = x + dx
                xx.append(x)
                dx *= ratio[2 * i]
                if(x0 < xn and x > xn) or (x0 > xn and x < xn):
                    break

    if(type == 'dist'):
        x, dx = x0, dx0
        if len(ratio) % 2:
            ratio.append(xn)
        for i in range(n):
            xc = x0 + ratio[2*i + 1]
            while (x0 < xc and x < xc) or (x0 > xc and x > xc):
                x = x + dx
                xx.append(x)
                dx *= ratio[2 * i]



    return xx


def symmetry(xx, type):
    '''
    :param xx: array {x0, x1, x2}
    :param type: string, 'left' or 'right'
    :return: mirror the array to left or right
             if type = 'left' xx = {2x0 - x2, 2x0 - x1, x0, x1, x2}
             if type = 'right' xx = {x0, x1, x2, 2x2 - x1, 2x2 - x0}
    '''

    n = len(xx)

    if(type == 'left'):
        xx_temp = []
        xc = xx[0]
        for i in range(n - 1, 0, -1):
            xx_temp.append(2*xc - xx[n - 1 - i])
        xx = xx_temp + xx
    if(type == 'right'):
        xc = xx[-1]
        for i in range(n - 1):
            xx.append(2*xc - xx[n - 1 - i - 1])
    return xx



def parachute3D(type=0):
    '''
    To generate initial mesh for 3D parachute
    :param type: 0: just for canopy, 1: for both the whole parachute
    :return:
    '''
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    xcl, xcr, nx = -15.0, 15.0, 61
    xc = np.linspace(xcl, xcr, nx)
    dxc = (xcr - xcl) / (nx - 1)
    xRatio = 1.1
    xl = geomspace(xcl, -80, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 80, dxc, False, 'num', [xRatio])


    x = [*reversed(xl),*xc,*xr]


    y = x
    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-20 200]

    if(type == 0):
        zcl, zcr, nz = 25.0, 50.0, 51

    elif(type == 1):
        #zcl, zcr, nz = -10.0, 60.0, 141
        zcl, zcr, nz = -15.0, 65.0, 161

    zc = np.linspace(zcl, zcr, nz)
    dzc = (zcr - zcl)/(nz - 1)
    zRatio = 1.1
    zl = geomspace(zcl, -20, -dzc, False, 'num', [zRatio])
    zr = geomspace(zcr, 180, dzc, False, 'num', [zRatio])
    z = [*reversed(zl), *zc, *zr]

    #shift the domain
    #z = np.array(z) - 4.7487995381941190e+01

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    print('len(x) ', len(x), ' len(y) ', len(y), ' len(z) ', len(z))
    simpleKuhnSimplex = KuhnSimplex(x,y,z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()



def cylinderStretched3D():
    '''
    To generate initial mesh for 3D parachute
    :param type: 0: just for canopy, 1: for both the whole parachute
    :return:
    '''
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    xcl, xcr, nx = -1.0, 1.0, 100
    xc = np.linspace(xcl, xcr, nx)
    dxc = (xcr - xcl) / (nx - 1)
    xRatio = 1.1
    xl = geomspace(xcl, -10, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 50, dxc, False, 'num', [xRatio])
    x = [*reversed(xl),*xc,*xr]


    ycl, ycr, ny = -1.5, 1.5, 32
    yc = np.linspace(ycl, ycr, ny)
    dyc = (ycr - ycl) / (ny - 1)
    yRatio = 1.1
    yl = geomspace(ycl, -10, -dyc, False, 'num', [yRatio])
    yr = geomspace(ycr, 10, dyc, False, 'num', [yRatio])
    y = [*reversed(yl),*yc,*yr]

    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-20 200]

   
   
    zcl, zcr, nz = -1.0, 1.0, 100
    zc = np.linspace(zcl, zcr, nz)
    dzc = (zcr - zcl)/(nz - 1)
    zRatio = 1.1
    zl = geomspace(zcl, -10, -dzc, False, 'num', [zRatio])
    zr = geomspace(zcr, 10, dzc, False, 'num', [zRatio])
    z = [*reversed(zl), *zc, *zr]

    #shift the domain
    #z = np.array(z) - 4.7487995381941190e+01

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    print('len(x) ', len(x), ' len(y) ', len(y), ' len(z) ', len(z))
    simpleKuhnSimplex = KuhnSimplex(x,y,z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()

def subscaleParachute3D():
    '''
    To generate initial mesh for 3D parachute
    :param type: 0: just for canopy, 1: for both the whole parachute
    :return:
    '''
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    scale = 0.04
    xcl, xcr, nx = -15.0*scale, 15.0*scale, 61
    xc = np.linspace(xcl, xcr, nx)
    dxc = (xcr - xcl) / (nx - 1)
    xRatio = 1.1
    xl = geomspace(xcl, -80*scale, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 80*scale, dxc, False, 'num', [xRatio])


    x = [*reversed(xl),*xc,*xr]


    y = x
    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-20 200]

    
    zcl, zcr, nz = -15.0*scale, 65.0*scale, 161

    zc = np.linspace(zcl, zcr, nz)
    dzc = (zcr - zcl)/(nz - 1)
    zRatio = 1.1
    zl = geomspace(zcl, -20*scale, -dzc, False, 'num', [zRatio])
    zr = geomspace(zcr, 180*scale, dzc, False, 'num', [zRatio])
    z = [*reversed(zl), *zc, *zr]

    #shift the domain
    #z = np.array(z) - 4.7487995381941190e+01

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    print('len(x) ', len(x), ' len(y) ', len(y), ' len(z) ', len(z))
    simpleKuhnSimplex = KuhnSimplex(x,y,z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()





# def RAE2822():
#
#     xcl,xcr,xcn = -1.0/9.0, 1.5 - 1.0/9.0, 1501
#     xc = np.linspace(xcl, xcr, xcn)
#     dxc = xc[1] - xc[0]
#     xRatio = 1.1
#     xl = geomspace(xcl, -30, -dxc, xRatio, False)
#     xr = geomspace(xcr, 30, dxc, xRatio, False)
#     x = [*reversed(xl), *xc, *xr]
#
#
#     ycl, ycr, ycn = -0.75, 0.75, 151
#     yc = np.linspace(ycl, ycr, ycn)
#     dyc = yc[1] - yc[0]
#     yRatio = 1.1
#     yl = geomspace(ycl, -30, -dyc, yRatio, False)
#     yr = geomspace(ycr,  30, dyc, yRatio, False)
#     y = [*reversed(yl), *yc, *yr]
#
#     z = [0.0, 0.005]
#     boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
#                    'InletFixedSurface','InletFixedSurface','InletFixedSurface']
#     simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)
#
#     print('Writing to top file')
#     simpleKuhnSimplex.write_topfile()




def naca2D():
    '''
    This is the domain file used for AMR heaving NACA
    :return:
    '''

    #xcl,xcr,xcn = -1.0/9.0, 1.25 - 1.0/9.0, 1251
    eps = 5.e-4
    xcl, xcr, xcn = -0.5 + eps, 5. - 0.5 + eps, 1001
    xc = np.linspace(xcl, xcr, xcn)
    dxc = xc[1] - xc[0]
    xRatio = 1.1
    xl = geomspace(xcl, -50, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 50, dxc, False, 'num', [xRatio])
    x = [*reversed(xl), *xc, *xr]

    # ycl, ycr, ycn = -0.15, 0.15, 301
    ycl, ycr, ycn = -1.0, 1.0, 401
    yc = np.linspace(ycl, ycr, ycn)
    dyc = yc[1] - yc[0]
    yRatio = 1.1
    yl = geomspace(ycl, -50, -dyc, False, 'num', [yRatio])
    yr = geomspace(ycr,  50, dyc, False, 'num', [yRatio])
    y = [*reversed(yl), *yc, *yr]

    z = [0.0, 0.005]
    boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


def naca2D():
    '''
    This is the domain file used for AMR heaving NACA
    :return:
    '''

    #xcl,xcr,xcn = -1.0/9.0, 1.25 - 1.0/9.0, 1251
    eps = 5.e-4
    xcl, xcr, xcn = -0.5 + eps, 5. - 0.5 + eps, 1001
    xc = np.linspace(xcl, xcr, xcn)
    dxc = xc[1] - xc[0]
    xRatio = 1.1
    xl = geomspace(xcl, -50, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 50, dxc, False, 'num', [xRatio])
    x = [*reversed(xl), *xc, *xr]

    # ycl, ycr, ycn = -0.15, 0.15, 301
    ycl, ycr, ycn = -1.0, 1.0, 401
    yc = np.linspace(ycl, ycr, ycn)
    dyc = yc[1] - yc[0]
    yRatio = 1.1
    yl = geomspace(ycl, -50, -dyc, False, 'num', [yRatio])
    yr = geomspace(ycr,  50, dyc, False, 'num', [yRatio])
    y = [*reversed(yl), *yc, *yr]

    z = [0.0, 0.005]
    boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


def RAE2822():
    '''
    This is the domain file used for AMR Wall Model RAE2822
    :return:
    '''

    xcl,xcr,xcn = -1.0/9.0, 1.5 - 1.0/9.0, 1501
    xc = np.linspace(xcl, xcr, xcn)
    dxc = xc[1] - xc[0]
    xRatio = 1.1
    xl = geomspace(xcl, -50, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 50, dxc, False, 'num', [xRatio])
    x = [*reversed(xl), *xc, *xr]

    # ycl, ycr, ycn = -0.15, 0.15, 301
    ycl, ycr, ycn = -0.075, 0.075, 151
    yc = np.linspace(ycl, ycr, ycn)
    dyc = yc[1] - yc[0]
    yRatio = 1.1
    yl = geomspace(ycl, -50, -dyc, False, 'num', [yRatio])
    yr = geomspace(ycr,  50, dyc, False, 'num', [yRatio])
    y = [*reversed(yl), *yc, *yr]

    z = [0.0, 0.002]
    boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()





def heavingNaca2D():

    xcl, dxc1 = -0.1, 0.2
    xcr, dx = 1 - xcl, 1.e-3

    xRatio = [1, dxc1, 1.1, 0.5 - xcl]
    xc = geomspace(xcl, xcr, dx, True, 'dist', xRatio)

    dxc = xc[1] - xc[0]
    xRatio = 1.1
    xl = geomspace(xc[0], -50, -dxc,  False, 'num', [xRatio])
    x = [*reversed(xl), *xc]
    x = symmetry(x, 'right')

    print('length x is ', len(x))

    ycl, ycr, ycn = -0.3, 0.0, 301
    yc = np.linspace(ycl, ycr, ycn)
    dyc = yc[1] - yc[0]

    yRatio = [1.1]
    yl = geomspace(ycl, -40, -dyc,  False, 'num', yRatio)
    y = [*reversed(yl), *yc]
    y = symmetry(y, 'right')
    print('length y is ', len(y))

    z = [0.0, 0.005]
    boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()








def Tube():
    #Important the name of these boundaries are important
    #x = np.linspace(-5.0,  5.0, 401)
    x = np.linspace(-5.0, 5.0, 8001)
    y = np.linspace(-0.05, 0.05, 5)
    z = np.linspace(-0.05, 0.05, 5)
    boundaryNames = ['SymmetrySurface', 'SymmetrySurface', 'InletFixedSurface',
                     'InletFixedSurface', 'SymmetrySurface', 'SymmetrySurface']
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()

def shockTube():
    #Important the name of these boundaries are important
    #x = np.linspace(-5.0, 5.0, 101)
    #y = np.linspace(-0.25, 0.25, 6)
    #z = np.linspace(-0.25, 0.25, 6)

    x = np.linspace(-2.5, 2.5, 101)
    y = np.linspace(-0.05, 0.05, 3)
    z = np.linspace(-0.05, 0.05, 3)
    boundaryNames = ['SymmetrySurface', 'SymmetrySurface', 'InletFixedSurface',
                     'OutletFixedSurface', 'SymmetrySurface', 'SymmetrySurface']
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile(volFunc=lambda xc: 0 if xc[0] < 1.0e-6 else 1)



def uniform2D():
    # for circle
    x = np.linspace(-4.0, 8.0, 120)
    y = np.linspace(-4.0, 4.0, 80)

    # for NACA
    # x = np.linspace(-1.0, 2.0, 20)
    # y = np.linspace(-0.3, 0.3, 20)



    z = [0, 0.1]
    boundaryNames=['SymmetrySurface','SymmetrySurface','InletFixedSurface',
                   'InletFixedSurface','InletFixedSurface','InletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x,y,z,boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


def uniform3D():
    # for sphere

    # x = np.linspace(-1.0,      4.0,       51)
    # y = np.linspace(-1.0,      1.0,       21)
    # z = np.linspace(-1.0,     16.0,      171)

    x = np.linspace(-2.0,      4.0,       31)
    y = np.linspace(-2.0,      2.0,       21)
    z = np.linspace(-2.0,      8.0,       51)

    # x = np.linspace(-0.1 ,0.1,   26)
    # y = np.linspace(-0.1, 0.1,   26)
    # z = np.linspace(-0.2, 1.0,   51)



    boundaryNames=['InletFixedSurface','SymmetrySurface','InletFixedSurface',
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
    xl = geomspace(xcl, -40, -dxc,  False, 'num', [xRatio])
    xr = geomspace(xcr, 80, dxc,  False, 'num', [xRatio])

    x = [*reversed(xl), *xc, *xr]


    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-100 200]
    ycl, ycr, ycn = -2.0, 2.0, 35
    yc = np.linspace(ycl, ycr, ycn)
    dyc = (ycr - ycl) / (ycn - 1)
    yRatio = 1.5
    yl = geomspace(ycl, -40, -dyc,  False, 'num', [yRatio])
    yr = geomspace(ycr, 40, dyc,  False, 'num', [yRatio])
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
    xl = geomspace(xcl, -40, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 80, dxc, False, 'num', [xRatio])

    x = [*reversed(xl), *xc, *xr]


    # z direction, the parachute is in [35.7358, 39.2198]
    # the fluid domain is [-100 200]
    ycl, ycr, ycn = -10.0, 10.0,35
    yc = np.linspace(ycl, ycr, ycn)
    dyc = (ycr - ycl) / (ycn - 1)
    yRatio = 1.2
    yl = geomspace(ycl, -40, -dyc,False, 'num', [yRatio])
    yr = geomspace(ycr, 40, dyc, False, 'num', [yRatio])
    y = [*reversed(yl), *yc, *yr]

    z = y

    boundaryNames = ['InletFixedSurface' for i in range(6)]
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()

def deformedShape():
    x = np.linspace( -5.0,   15.0, 41)
    y = np.linspace( -5.0,    5.0, 21)
    z = np.linspace(-10.0,   10.0, 41)


    boundaryNames = ['OutletFixedSurface', 'OutletFixedSurface', 'InletFixedSurface',
                     'OutletFixedSurface', 'OutletFixedSurface', 'OutletFixedSurface']
    simpleKuhnSimplex = KuhnSimplex(x, y, z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


def F16():
    '''
    To generate initial mesh for F16
    :return:
    '''
    # x direction, the parachute is in [-7.7235, 7.7235]
    # the fluid domain is [-150 150]
    #body = np.array([[-7, 566.797],[-186.325, 186.325],[38.8403, 217.66]])
    #body = np.array([[-7, 566.797], [-186.325, 186.325], [38.8403, 135.373]])

    xcl, xcr, nx = -20.0, 780.0, 41
    xc = np.linspace(xcl, xcr, nx)
    dxc = (xcr - xcl) / (nx - 1)
    xRatio = 1.1
    xl = geomspace(xcl, -3000, -dxc, False, 'num', [xRatio])
    xr = geomspace(xcr, 3000, dxc, False, 'num', [xRatio])
    x = [*reversed(xl),*xc,*xr]

    ycl, ycr, ny = -300.0, 300.0, 31
    yc = np.linspace(ycl, ycr, ny)
    dyc = (ycr - ycl) / (ny - 1)
    yRatio = 1.1
    yl = geomspace(ycl, -1500, -dyc, False, 'num', [yRatio])
    yr = geomspace(ycr, 1500, dyc, False, 'num', [yRatio])
    y = [*reversed(yl), *yc, *yr]

    zcl, zcr, nz = 0.0, 240.0, 49
    zc = np.linspace(zcl, zcr, nz)
    dzc = (zcr - zcl)/(nz - 1)
    zRatio = 1.1
    zl = geomspace(zcl, -3000, -dzc, False, 'num', [zRatio])
    zr = geomspace(zcr, 3000, dzc, False, 'num', [zRatio])
    z = [*reversed(zl), *zc, *zr]



    boundaryNames = ['InletFixedSurface' for i in range(6)]
    print('len(x) ', len(x), ' len(y) ', len(y), ' len(z) ', len(z))
    simpleKuhnSimplex = KuhnSimplex(x,y,z, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


def quad2D():


    '''
    To generate 2D quad to test
    :return:
    '''

    xl, xr, nx = -1.0, 1.0, 10
    x = np.linspace(xl, xr, nx)


    yl, yr, ny = 0.0, 1.0, 10
    y = np.linspace(yl, yr, ny)


    boundaryNames = ['InletFixedSurface',  'InletFixedSurface',  'SymmetrySurface',  'SymmetrySurface']
    print('len(x) ', len(x), ' len(y) ', len(y))
    simpleKuhnSimplex = KuhnSimplex2D(x, y, boundaryNames)

    print('Writing to top file')
    simpleKuhnSimplex.write_topfile()


if __name__ == '__main__':
    #uniform2D()
    #RAE2822()
    #naca2D();
    #RAE2822();
    #uniform2D();
    #heavingNaca2D()
    #Tube()
    #shockTube()
    #parachute3D(1)
    #subscaleParachute3D()
    #deformedShape()
    #F16()
    quad2D()
    #cylinderStretched3D()
    # x = [1,2]
    # y = [4,5]
    # z = [7,8]
    #
    # simpleKuhnSimplex = KuhnSimplex(x,y,z)
    #
    # print('Writing to top file')
    # simpleKuhnSimplex.plot_mesh()
