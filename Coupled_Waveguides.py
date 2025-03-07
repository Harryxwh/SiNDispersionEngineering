import numpy as np
from Waveguide import *
class Coupled_Waveguides():

    c           =   299792458
    u0          =   4 * np.pi * 1e-7
    component_name_list = ['Ex','Ey','Ez','Hx','Hy','Hz']

    def __init__(self, n1, n0, gap_x, gap_y, wavelength,
                 name1,name2, ModeIdx1, ModeIdx2, param_file_name,Plot_field,Plot_index):
        self.load_param(param_file_name)
        self.n1 = n1
        self.n0 = n0
        self.gap_x = gap_x                  # unit :um
        self.gap_y = gap_y                  # unit :um
        self.bend_radius_outer = self.bend_radius_inner +\
                                 self.gap_x + (self.WG1_width + self.WG2_width)/2
        self.wavelength = wavelength        # unit :um
        Field_shape_padded = np.array([self.Convert_to_num_of_cells(self.FDE_height_padded,axis='y'),
                                       self.Convert_to_num_of_cells(self.FDE_width_padded,axis='x')])
        self.WG1 = Waveguide(self.WG1_width,self.WG1_height,
                             bendradius=self.bend_radius_inner,PML_len=self.PML_len,
                             foldername=name1,ModeIdx=1,
                             Field_shape_padded=Field_shape_padded,Plot=Plot_field)
        self.WG2 = Waveguide(self.WG2_width,self.WG2_height,
                             bendradius=self.bend_radius_outer,PML_len=self.PML_len,
                             foldername=name2,ModeIdx=1,
                             Field_shape_padded=Field_shape_padded,Plot=Plot_field)
        # self.neff_1 = self.WG1.Modes_info_list[self.WG1.ModeIdx-1].neff
        # self.neff_2 = self.WG2.Modes_info_list[self.WG2.ModeIdx-1].neff

        # unit: rad/rad.
        self.beta_ang_1 = self.WG1.Modes_info_list[self.WG1.ModeIdx-1].beta_ang
        self.beta_ang_2 = self.WG2.Modes_info_list[self.WG2.ModeIdx-1].beta_ang
        self.beta_ave   = (self.beta_ang_2 + self.beta_ang_1)/2
        self.delta_beta = (self.beta_ang_2 - self.beta_ang_1)/2

        self.Calculate_index_profile()
        if Plot_index:
            self.Plot_index_profile()
        self.Field_dict_uncoupled = self.Shifted_field()

    # Load the class variables from Param.csv
    def load_param(self,param_file_name):
        dtype = [('key', str), ('value', str),('dtype', str), ('comments', str)]
        data = np.loadtxt(param_file_name, delimiter=',', dtype=str,skiprows=1)
        for param in data:
            exec("self."+param[0]+" = "+param[2]+"("+param[1]+")")
            # print("self."+param[0]+" = "+str(eval("self."+param[0])))
            # print(param[3])
            # print("---------------------------")

    # Convert the unit from "um" to "num of cells"
    def Convert_to_num_of_cells(self,len_in_um,axis):
        if axis=='x':
            return int(len_in_um * self.Num_of_cells_x/self.FDE_width)
        else:
            return int(len_in_um * self.Num_of_cells_y/self.FDE_height)

    # Average bend radius. (unit: m)
    @property
    def Bend_radius_ave(self):
        return (self.bend_radius_outer + self.bend_radius_inner) /2 *1e-6

    # Wavenumber in vacuum. (unit: rad/m)
    @property
    def k0(self):
        return 2*np.pi / self.wavelength  * 1e6

    #shift_x: (unit: # of data points) Distance between the centers of two waveguides in x axis
    #shift_y: (unit: # of data points) Distance between the centers of two waveguides in y axis
    #shift_x * shift_y = 0
    @property
    def Calculate_shift(self):
        if self.gap_x == 0:
            shift_x = 0
        else:
            shift_x = self.Convert_to_num_of_cells(self.gap_x +
                                                   (self.WG1.width+self.WG2.width)/2,axis='x')

        if self.gap_y == 0:
            shift_y = 0
        else:
            shift_y = self.Convert_to_num_of_cells(self.gap_y +
                                                   (self.WG1.height+self.WG2.height)/2,axis='y')
        return shift_x,shift_y

    #Compressed_shape=(len_y,len_x)
    #The largest possible shape of shifted field profile w/o zero padding
    @property
    def Compressed_shape(self):
        #Field_shape=(len_y,len_x)
        Field_shape = np.concatenate((np.array(self.WG1.Field_shape).reshape((1,2)),
                                    np.array(self.WG2.Field_shape).reshape((1,2))),
                                     axis = 0)
        Field_shape = np.min(Field_shape,axis=0)
        len_y = Field_shape[0]
        len_x = Field_shape[1]
        shift_x,shift_y = self.Calculate_shift
        # The maxmimum possible size of the region where data of both modes exist
        compressed_len_y_max = len_y - shift_y
        compressed_len_x_max = len_x - shift_x
        assert compressed_len_y_max > 0
        assert compressed_len_x_max > 0
        compressed_len_y = int(compressed_len_y_max * self.Shrink_ratio)
        compressed_len_x = int(compressed_len_x_max * self.Shrink_ratio)
        return np.array([compressed_len_y,compressed_len_x])

    #Return N, N1, N2
    def Calculate_index_profile(self):
        compressed_len_y = self.Compressed_shape[0]
        compressed_len_x = self.Compressed_shape[1]
        shift_x,shift_y = self.Calculate_shift

        WG1_width = self.Convert_to_num_of_cells(self.WG1.width,axis='x')
        WG1_height= self.Convert_to_num_of_cells(self.WG1.height,axis='y')
        WG2_width = self.Convert_to_num_of_cells(self.WG2.width,axis='x')
        WG2_height= self.Convert_to_num_of_cells(self.WG2.height,axis='y')

        assert compressed_len_y > max(WG1_height+shift_y, WG2_height+shift_y)
        assert compressed_len_x > max(WG1_width+shift_x, WG2_width+shift_x)

        # N(x,y):The overall index matrix
        self.N = np.ones(self.Compressed_shape) * self.n0
        self.N[int(compressed_len_y/2-WG1_height/2-shift_y/2) :
                int(compressed_len_y/2+WG1_height/2-shift_y/2),
                int(compressed_len_x/2-WG1_width/2-shift_x/2) :
                int(compressed_len_x/2+WG1_width/2-shift_x/2)] = self.n1
        self.N[int(compressed_len_y/2-WG2_height/2+shift_y/2) :
                int(compressed_len_y/2+WG2_height/2+shift_y/2),
                int(compressed_len_x/2-WG2_width/2+shift_x/2) :
                int(compressed_len_x/2+WG2_width/2+shift_x/2)] = self.n1
        # N1(x,y):The index matrix of WG1 in the absence of WG2
        self.N1 = np.ones(self.Compressed_shape) * self.n0
        self.N1[int(compressed_len_y/2-WG1_height/2-shift_y/2) :
                int(compressed_len_y/2+WG1_height/2-shift_y/2),
                int(compressed_len_x/2-WG1_width/2-shift_x/2) :
                int(compressed_len_x/2+WG1_width/2-shift_x/2)] = self.n1
        # N2(x,y):The index matrix of WG2 in the absence of WG1
        self.N2 = np.ones(self.Compressed_shape) * self.n0
        self.N2[int(compressed_len_y/2-WG2_height/2+shift_y/2) :
                int(compressed_len_y/2+WG2_height/2+shift_y/2),
                int(compressed_len_x/2-WG2_width/2+shift_x/2) :
                int(compressed_len_x/2+WG2_width/2+shift_x/2)] = self.n1

    def Plot_index_profile(self,label='N',
                           save_name='./results/index_profile.png',dpi=150):
        plt.figure(figsize=(5,4))
        if label == 'N':
            plt.imshow(self.N)
        elif label == 'N1':
            plt.imshow(self.N1)
        else:
            plt.imshow(self.N2)
        plt.gca().invert_yaxis()
        plt.savefig(save_name,dpi=dpi)
        plt.show()

    # Calculate the shifted fields; Return Field_dict_uncoupled
    def Shifted_field(self):
        Field_dict_uncoupled = dict.fromkeys(self.component_name_list)
        shift_x,shift_y = self.Calculate_shift
        shift_y = int(shift_y + (1-self.Shrink_ratio) * (self.Compressed_shape[0]/2))
        shift_x = int(shift_x + (1-self.Shrink_ratio) * (self.Compressed_shape[1]/2))
        for component in self.component_name_list:
            if shift_y == 0:
                Field_profile_1 = self.WG1.Field_dict[component][:,shift_x:]   #Left part of E_1y
                Field_profile_2 = self.WG2.Field_dict[component][:,:-shift_x]  #Right part of E_2y
            elif shift_x == 0:
                Field_profile_1 = self.WG1.Field_dict[component][shift_y:,:]   #Left part of E_1y
                Field_profile_2 = self.WG2.Field_dict[component][:-shift_y,:]  #Right part of E_2y
            else:
                Field_profile_1 = self.WG1.Field_dict[component][shift_y:,shift_x:]   #Left part of E_1y
                Field_profile_2 = self.WG2.Field_dict[component][:-shift_y,:-shift_x]  #Right part of E_2y
            assert np.all(np.shape(Field_profile_1) == self.Compressed_shape)
            assert np.all(np.shape(Field_profile_2) == self.Compressed_shape)
            Field_dict_uncoupled[component] = [Field_profile_1,Field_profile_2]

        return Field_dict_uncoupled

    def Kappa(self,p_,q_,Nq):
        p = p_ - 1
        q = q_ - 1
        K = np.sum((self.N*self.N-Nq*Nq) * (np.conj(self.Field_dict_uncoupled['Ex'][p]) * \
                                                    self.Field_dict_uncoupled['Ex'][q] + \
                                            np.conj(self.Field_dict_uncoupled['Ey'][p]) * \
                                                    self.Field_dict_uncoupled['Ey'][q] + \
                                            np.conj(self.Field_dict_uncoupled['Ez'][p]) * \
                                                    self.Field_dict_uncoupled['Ez'][q] ))

        K  = K * self.k0 / (self.c * self.u0)

        return K

    def C(self,p_,q_):
        p = p_ - 1
        q = q_ - 1
        C = np.sum (np.conj(self.Field_dict_uncoupled['Ex'][p]) *\
                            self.Field_dict_uncoupled['Hy'][q] -\
                    np.conj(self.Field_dict_uncoupled['Ey'][p]) *\
                            self.Field_dict_uncoupled['Hx'][q] +\
                            self.Field_dict_uncoupled['Ex'][q] *\
                    np.conj(self.Field_dict_uncoupled['Hy'][p] -
                            self.Field_dict_uncoupled['Ey'][q])*\
                    np.conj(self.Field_dict_uncoupled['Hx'][p]))
        return C

    def Chi(self,p_,Np):
        p = p_ - 1
        Chi = np.sum((self.N*self.N-Np*Np)*(np.conj(self.Field_dict_uncoupled['Ex'][p] * \
                                                    self.Field_dict_uncoupled['Ex'][p] + \
                                            np.conj(self.Field_dict_uncoupled['Ey'][p]) * \
                                                    self.Field_dict_uncoupled['Ey'][p] + \
                                            np.conj(self.Field_dict_uncoupled['Ez'][p]) * \
                                                    self.Field_dict_uncoupled['Ez'][p]) ))
        print(self.k0 / (self.c * self.u0))
        Chi  = Chi * self.k0 / (self.c * self.u0)
        return Chi

    def P_factor(self,p_, q_):
        p = p_ - 1
        q = q_ - 1
        P = 2*np.real(np.sum(np.conj(self.Field_dict_uncoupled['Ex'][p]) *\
                            self.Field_dict_uncoupled['Hy'][q] -
                            np.conj(self.Field_dict_uncoupled['Ey'][p]) *\
                            self.Field_dict_uncoupled['Hx'][q]))
        return P

    # C(p,q) = int_0^Inf dx dy [ Ep * np.conj(Hq) ]
    def C_factor(self,p_, q_, N, Np):
        p = p_ - 1
        q = q_ - 1
        C = np.sum((self.N*self.N-Np*Np)*(
                                np.conj(self.Field_dict_uncoupled['Ex'][p] * \
                                        self.Field_dict_uncoupled['Ex'][q] + \
                                np.conj(self.Field_dict_uncoupled['Ey'][p]) * \
                                        self.Field_dict_uncoupled['Ey'][q] + \
                                np.conj(self.Field_dict_uncoupled['Ez'][p]) * \
                                        self.Field_dict_uncoupled['Ez'][q]) ))
        return np.real(C)

    # D/Dz [A B] = F [A B]
    # F: Differential Euqation Matrix
    @property
    # def F_Matrix(self):
    #     M = np.array([[self.P_factor(1,1), self.P_factor(1,2)],
    #                   [self.P_factor(2,1), self.P_factor(2,2)]])
    #     assert not np.linalg.det(M) == 0
    #     M_inv = np.linalg.inv(M)
    #     N = self.N
    #     N1 = self.N1
    #     N2 = self.N2
    #     C = np.array([[self.C_factor(1,1,N,N1),self.C_factor(2,1,N,N2)],
    #                   [self.C_factor(1,2,N,N1),self.C_factor(2,2,N,N2)]])

    #     F = M_inv * C
    #     return F
    def F_Matrix(self):
        C_12 = self.C(1,2)
        C_21 = self.C(2,1)
        K_12 = self.Kappa(1,2,self.N2)
        K_21 = self.Kappa(2,1,self.N1)
        Chi_1 = self.Chi(1,self.N1)
        Chi_2 = self.Chi(2,self.N2)

        # delta_beta is in the unit of rad/rad
        # The unit of the rest of the elements should also be converted to rad/rad
        Alpha_a = (K_21 * C_12 - Chi_1)/(1-np.abs(C_12)**2) * self.Bend_radius_ave
        Alpha_b = (K_12 * C_21 - Chi_2)/(1-np.abs(C_12)**2) * self.Bend_radius_ave
        Kappa_a = (K_12 - C_12 * Chi_2)/(1-np.abs(C_12)**2) * self.Bend_radius_ave
        Kappa_b = (K_21 - C_21 * Chi_1)/(1-np.abs(C_12)**2) * self.Bend_radius_ave
        F = np.array([[Alpha_a + self.delta_beta, -Kappa_a],
                      [-Kappa_b, Alpha_b - self.delta_beta]])

        F = F
        return F

    def Find_supermodes(self):
        F = self.F_Matrix
        F_sym =  np.diag(np.ones(2))*((F[0][0] + F[1][1])) /2
        F_antisym = F - F_sym

        Eigenvalues = np.array(np.linalg.eig(F_antisym)[0])
        Eigenvectors = np.array(np.linalg.eig(F_antisym)[1])

        Rabi_freq = Eigenvalues.reshape((1,2))
        beta_minus_ave = Rabi_freq #+ np.array([1,-1]) * self.delta_beta

        print("wavelength = ",self.wavelength)
        print("F = ",F)
        print("F_sym = ",F_sym)
        print("F_antisym = ",F_antisym)
        print("beta_ave = ",self.beta_ave)
        print("beta_minus_ave = ",beta_minus_ave)
        print("delta_beta = ",self.delta_beta)
        print("S = ",Rabi_freq)
        print("EigenVec = ",Eigenvectors)
        print("\n")
        return Rabi_freq , beta_minus_ave