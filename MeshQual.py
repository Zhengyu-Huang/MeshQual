import sys
import numpy as np
import matplotlib.pyplot as plt

def _tetrahedron_volume(nodes):
    '''
    :param nodes[4,3]:
    :return: tetrahedron volume
    '''
    return np.fabs(np.dot(nodes[0, :] - nodes[3, :],
                          np.cross(nodes[1, :] - nodes[3, :], nodes[2, :] - nodes[3, :]))) / 6.0

def _tetrahedron_surf(nodes):
    '''
    :param nodes[4,3]:
    :return: tetrahedron SA_1^2 + SA_2^2 + SA_3^2 + SA_4^2
    '''
    SA2 = 0.25 * (np.linalg.norm(np.cross(nodes[1, :] - nodes[3, :], nodes[2, :] - nodes[3, :]))**2
               + np.linalg.norm(np.cross(nodes[2, :] - nodes[0, :], nodes[3, :] - nodes[0, :]))**2
               + np.linalg.norm(np.cross(nodes[3, :] - nodes[1, :], nodes[0, :] - nodes[1, :]))**2
               + np.linalg.norm(np.cross(nodes[0, :] - nodes[2, :], nodes[1, :] - nodes[2, :]))**2)


    return SA2

def _tetrahedron_radius(nodes):
    '''
    :param nodes[4,3]:
    :return: tetrahedron circumscribed sphere radius and inscribed sphere radius
    '''

    vol = np.fabs(np.dot(nodes[0, :] - nodes[3, :],
                         np.cross(nodes[1, :] - nodes[3, :], nodes[2, :] - nodes[3, :]))) / 6.0

    S = 0.5 * (np.linalg.norm(np.cross(nodes[1, :] - nodes[3, :], nodes[2, :] - nodes[3, :]))
               + np.linalg.norm(np.cross(nodes[2, :] - nodes[0, :], nodes[3, :] - nodes[0, :]))
               + np.linalg.norm(np.cross(nodes[3, :] - nodes[1, :], nodes[0, :] - nodes[1, :]))
               + np.linalg.norm(np.cross(nodes[0, :] - nodes[2, :], nodes[1, :] - nodes[2, :])))

    IR = 3 * vol / S

    # https: // math.stackexchange.com / questions / 1087011 / calculating -
    # the - radius - of - the - circumscribed - sphere - of - an - arbitrary - tetrahedron

    a = np.linalg.norm(nodes[0, :] - nodes[1, :])
    a1 = np.linalg.norm(nodes[2, :] - nodes[3, :])
    b = np.linalg.norm(nodes[0, :] - nodes[2, :])
    b1 = np.linalg.norm(nodes[1, :] - nodes[3, :])
    c = np.linalg.norm(nodes[0, :] - nodes[3, :])
    c1 = np.linalg.norm(nodes[1, :] - nodes[2, :])
    p = (a * a1 + b * b1 + c * c1) / 2.0
    CR = np.sqrt(p * (p - a * a1) * (p - b * b1) * (p - c * c1)) / (6 * vol)

    return CR, IR


def _tetrahedron_edgelen(nodes):
    '''
    :param nodes[4,3]:
    :return: edgemin, edgemax, edgeave edgerms
    '''
    e = 0
    edgelen = np.empty(6)
    for n1 in range(4):
        for n2 in range(n1 + 1, 4):
            edgelen[e] = np.linalg.norm(nodes[n1, :] - nodes[n2, :])
            e = e + 1
    edge_min = np.min(edgelen)
    edge_max = np.max(edgelen)
    edge_ave = np.mean(edgelen)
    edge_rms = np.sqrt(np.mean(edgelen ** 2))

    return edge_min, edge_max, edge_ave, edge_rms

def _tetrahedron_quantities(nodes):
    '''
    :param nodes[4,3]:
    :return: tetrahedron circumscribed sphere radius and inscribed sphere radius
    '''

    vol = np.fabs(np.dot(nodes[0, :] - nodes[3, :],
                          np.cross(nodes[1, :] - nodes[3, :], nodes[2, :] - nodes[3, :]))) / 6.0

    SA =  np.array([0.5 * np.linalg.norm(np.cross(nodes[1, :] - nodes[3, :], nodes[2, :] - nodes[3, :])),
                    0.5 * np.linalg.norm(np.cross(nodes[2, :] - nodes[0, :], nodes[3, :] - nodes[0, :])),
                    0.5 * np.linalg.norm(np.cross(nodes[3, :] - nodes[1, :], nodes[0, :] - nodes[1, :])),
                    0.5 * np.linalg.norm(np.cross(nodes[0, :] - nodes[2, :], nodes[1, :] - nodes[2, :]))])

    S = SA[0] + SA[1] + SA[2] + SA[3]


    S2 = SA[0]**2 + SA[1]**2 + SA[2]**2 + SA[3]**2

    IR = 3 * vol / S

    # https: // math.stackexchange.com / questions / 1087011 / calculating -
    # the - radius - of - the - circumscribed - sphere - of - an - arbitrary - tetrahedron

    a = np.linalg.norm(nodes[0, :] - nodes[1, :])
    a1 = np.linalg.norm(nodes[2, :] - nodes[3, :])
    b = np.linalg.norm(nodes[0, :] - nodes[2, :])
    b1 = np.linalg.norm(nodes[1, :] - nodes[3, :])
    c = np.linalg.norm(nodes[0, :] - nodes[3, :])
    c1 = np.linalg.norm(nodes[1, :] - nodes[2, :])
    p = (a * a1 + b * b1 + c * c1) / 2.0
    CR = np.sqrt(p * (p - a * a1) * (p - b * b1) * (p - c * c1)) / (6 * vol)

    edgelen = np.array([a,a1,b,b1,c,c1])

    edge_min = np.min(edgelen)

    edge_max = np.max(edgelen)

    edge_ave = np.mean(edgelen)

    edge_rms = np.sqrt(np.mean(edgelen ** 2))


    return vol, edge_min, edge_max, edge_ave, edge_rms, CR, IR, S2

class mesh:
    def __init__(self,mshfile):
        try:
            fid = open(mshfile, "r")
        except IOError:
            print("File '%s' not found." % mshfile)
            sys.exit()

        nodes = []
        elems = []

        print('Reading mesh ...')

        line = fid.readline() # should be "Nodes FluidNodes"

        while line:
            line = fid.readline()
            data = line.split()
            if data[0] == 'Elements':
                break;
            nodes.append(list(map(float, data[1:4])))

        while line:
            line = fid.readline()
            data = line.split()
            if data[0] == 'Elements':
                break;
            elems.append(list(map(int, data[2:6])))


        self.nodes = np.array(nodes, dtype = float)
        self.elems = np.array(elems,dtype = int) - 1 #fixed top file 1-based index



        self.nNodes = len(nodes)
        self.nElems = len(elems)

        ################# Element quantity
        print('Node number is ', self.nNodes, ' Element number is ', self.nElems)
        print('Computing Elements Quantities ...')
        self.compute_elem_quantities()
        print('Computing Elements Aspect Ratios ... ...')
        self.compute_elem_aspect_ratios()


    def num_nodes(self):
        return self.nNodes

    def num_elem(self):
        return self.nElems

    def get_nodes(self):
        return self.nodes

    def get_elements(self):
        return self.elems



    def compute_elem_quantities(self, print_freq=5000):
        '''
        compute element volumes: self.volumes
        min edge length: self.edge_min
        max edge length: self.edge_max
        average edge length: self.edge_ave
        root mean square edge length: self.edge_rms
        radius of the circumscribed sphere: self.CR
        radius of the inscribed sphere: self.IR
        :return:
        '''
        nElems = self.nElems

        self.vol = np.empty(nElems)
        self.edge_min = np.empty(nElems)
        self.edge_max = np.empty(nElems)
        self.edge_ave = np.empty(nElems)
        self.edge_rms = np.empty(nElems)
        self.CR = np.empty(nElems)
        self.IR = np.empty(nElems)
        self.SAsquare = np.empty(nElems)


        for i in range(nElems):
            if(i%print_freq == 0):
                print('Finish ', i, ' Elements')

            nodes = self.nodes[self.elems[i,:]]

            # self.vol[i] = _tetrahedron_volume(nodes)
            #
            # self.edge_min[i],self.edge_max[i],self.edge_ave[i],self.edge_rms[i] = _tetrahedron_edgelen(nodes)
            #
            # self.CR[i], self.IR[i] = _tetrahedron_radius(nodes)
            #
            # self.SAsquare[i] = _tetrahedron_surf(nodes)

            self.vol[i], self.edge_min[i],self.edge_max[i],self.edge_ave[i],self.edge_rms[i], self.CR[i], self.IR[i] ,self.SAsquare[i]  = _tetrahedron_quantities(nodes)


    def compute_elem_aspect_ratios(self):
        '''
        A comparison of tetrahedron quality measures
        V.N.Parthasarathy C.M.Graichen and A.F.Hathaway
        :return: 7 aspect ratios of in this paper
        '''
        nElems = self.nElems

        AR = np.empty((nElems, 7))

        AR[:, 0] = self.CR / self.IR
        AR[:, 1] = self.edge_max / self.IR
        AR[:, 2] = self.CR / self.edge_max
        AR[:, 3] = self.edge_max / self.edge_min
        AR[:, 4] = self.vol**4 / self.SAsquare**3
        AR[:, 5] = self.edge_ave**3 / self.vol
        AR[:, 6] = self.edge_rms**3 / self.vol

        self.AR = AR

    def plot(self, plot_sample = 100):
        #Sort Aspect Ratio
        AR = self.AR
        AR.sort(axis=0)


        vol = self.vol
        vol.sort(axis=0)

        f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8, sharex=True, sharey=False)
        ax1.plot(AR[:,0], '-')
        ax1.set_title(r'CR/IR(3.0)')
        ax2.plot(AR[:,1], '-')
        ax2.set_title(r'$S_{max}/IR(4.9)$')
        ax3.plot(AR[:,2], '-')
        ax3.set_title(r'$CR/S_{max}(0.613)$')
        ax4.plot(AR[:,3], '-')
        ax4.set_title(r'$S_{max}/S_{min}(1.0)$')
        ax5.plot(AR[:,4], '-')
        ax5.set_title(r'$V^4/SA^3(4.58e-4)$')
        ax6.plot(AR[:,5], '-')
        ax6.set_title(r'$S^3_{avg}/V(8.48)$')
        ax7.plot(AR[:,6], '-')
        ax7.set_title(r'$S^3_{rms}/V(8.48)$')
        ax8.plot(vol, '-')
        ax8.set_title(r'Vol')
        # Fine-tune figure; make subplots close to each other and hide x ticks for
        # all but bottom plot.
        # f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()



if __name__ == '__main__':
    ##############Test elem quality functions
    equilateral_tet = np.array([[-0.5,0,0],
                                [0.5,0.0,0.],
                                [0, np.sqrt(3.0)/2.0, 0.0],
                                [0.0, 1.0/(2*np.sqrt(3.0)), np.sqrt(2.0/3.0)]])
    right_tet = np.array([[0.0,0,0],
                          [1.0, 0.0, 0.],
                          [0.0, 1.0, 0.0],
                          [0.0, 0.0, 1.0]])

    trirectangular_tet = np.array([[1.0,0,0],
                          [0.0, 0.0, 0.],
                          [0.0, 2.0, 0.0],
                          [0.0, 0.0, 3.0]])

    for tet in [equilateral_tet, right_tet, trirectangular_tet]:
        vol = _tetrahedron_volume(tet)
        CR,IR = _tetrahedron_radius(tet)
        e_min, e_max, e_ave, e_rms = _tetrahedron_edgelen(tet)
        SAsquare = _tetrahedron_surf(tet)
        print('vol is ',vol)
        print('SAsquare is ', SAsquare)
        print('e_min is ', e_min, ' e_max is ', e_max, ' e_ave is ', e_ave, ' e_rms is ', e_rms)
        print('CR is ', CR, ' IR is ', IR, 'CR/IR', CR/IR)
        print('\n\n')

    ##############Test class
    tet_mesh = mesh('Equilateral_tet_mesh.top')
    AR = tet_mesh.AR
    print(AR)
    tet_mesh.plot()




