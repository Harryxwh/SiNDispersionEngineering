'''
class Coupled_Waveguides: Represent two coupled waveguides.

Parameters:
n_core              : refractive index of the waveguide core
n_cladding              : refractive index of the waveguide cladding
gap_x           : gap between to waveguides in x direction
gap_y           : gap between to waveguides in y direction
wavelength      : wavelength in vacuum of light
name1           : foldername of data of the first WG
name2           : foldername of data of the second WG
ModeIdx1        : eigenmode index of the first WG
ModeIdx2        : eigenmode index of the second WG
param_file_name : file name of the parameters
Plot_field      : True or False. whether to plot all the mode profiles.
Plot_index      : True or False. Whether to plot all the index profiles.
'''
import numpy as np
from Waveguide import *
class Coupled_Waveguides():

    c           =   299792458
    u0          =   4 * np.pi * 1e-7
    n_air       =   1
    # n_si        =   3.48
    component_name_list = ['Ex','Ey','Ez','Hx','Hy','Hz']

    def __init__(self, n_core, n_cladding,
                 gap_x, gap_y, wavelength,
                 name1,name2, ModeIdx1, ModeIdx2,
                 param_file_name,Plot_field,Plot_index):
        self.load_param(param_file_name)

        self.n_core = n_core
        self.n_cladding = n_cladding
        self.gap_x = gap_x                  # unit :um
        self.gap_y = gap_y                  # unit :um
        if gap_x > 0:
            self.bend_radius_outer = self.bend_radius_inner +\
                self.gap_x + (self.WG1_width + self.WG2_width)/2       # unit :um
        else:
            self.bend_radius_outer = self.bend_radius_inner
        self.wavelength = wavelength        # unit :nm
        self.FDE_width = self.FDE_x_max - self.FDE_x_min
        self.FDE_height = self.FDE_y_max - self.FDE_y_min

        FDE_shape_padded = (self.FDE_width_padded,self.FDE_height_padded)
        Num_of_cells     = (self.Num_of_cells_x,self.Num_of_cells_y)
        self.WG1 = Waveguide(self.WG1_width,self.WG1_height,self.bend_radius_inner,
                             (self.FDE_x_min,self.FDE_x_max),
                             (self.FDE_y_min,self.FDE_y_max),
                             PML_len=self.PML_len,
                             Num_of_cells= Num_of_cells,
                             foldername=name1,ModeIdx=1,
                             FDE_shape_padded=FDE_shape_padded,Plot=Plot_field)
        self.WG2 = Waveguide(self.WG2_width,self.WG2_height,self.bend_radius_outer,
                             (self.FDE_x_min,self.FDE_x_max),
                             (self.FDE_y_min,self.FDE_y_max),
                             PML_len=self.PML_len,
                             Num_of_cells= Num_of_cells,
                             foldername=name2,ModeIdx=1,
                             FDE_shape_padded=FDE_shape_padded,Plot=Plot_field)

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
            return int(len_in_um * self.Num_of_cells_x/ self.FDE_width)
        else:
            return int(len_in_um * self.Num_of_cells_y/ self.FDE_height)
    # Convert the unit from "num of cells" to "um"
    def Convert_to_um(self,num_cells,axis):
        if axis=='x':
            return self.FDE_width * num_cells / self.Num_of_cells_x
        else:
            return self.FDE_height * num_cells / self.Num_of_cells_y
    # Convert the unit of ticks of plots from "num of cells" to "um"
    def Convert_ticks(self,xticks_prev,yticks_prev):
        shift_x,shift_y = self.Calculate_shift
        yticks      = self.Convert_to_um(yticks_prev,axis='y') + self.FDE_y_min \
                        - (self.FDE_height_padded - self.FDE_height)/2 \
                        + self.Convert_to_um(shift_y,axis='y')
        yticks      = np.round(yticks, 2)
        xticks      = self.Convert_to_um(xticks_prev - np.max(xticks_prev)/2
                                         ,axis='x') -\
                                        (self.WG1.width - self.WG2.width)/4
        xticks      = np.round(xticks, 2)
        return xticks,yticks

    # Average bend radius. (unit: m)
    @property
    def Bend_radius_ave(self):
        return (self.bend_radius_outer + self.bend_radius_inner) /2 *1e-6

    # Wavenumber in vacuum. (unit: rad/m)
    @property
    def k0(self):
        return 2*np.pi / self.wavelength  * 1e9

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
        return (shift_x,shift_y)
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

    # The coordinate of WGs (unit: num of cells)
    @property
    def Structure_coordinates(self):
        x_min_padded = self.FDE_x_min - (self.FDE_width_padded - self.FDE_width)/2
        y_min_padded = self.FDE_y_min - (self.FDE_height_padded - self.FDE_height)/2
        WG1_y = self.Convert_to_num_of_cells(self.WG1_y - y_min_padded,axis='y')
        WG2_y = self.Convert_to_num_of_cells(self.WG2_y - y_min_padded,axis='y')
        WG1_x = self.Convert_to_num_of_cells(self.WG1_x - x_min_padded,axis='x')
        WG2_x = self.Convert_to_num_of_cells(self.WG2_x - x_min_padded,axis='x')
        Ly_cladding = self.Convert_to_num_of_cells(self.Ly_cladding - y_min_padded, axis='y')
        return (WG1_x,WG1_y,WG2_x,WG2_y,Ly_cladding)

    @property
    def WG_geometry_in_num_of_cells(self):
        WG1_width = self.Convert_to_num_of_cells(self.WG1.width,axis='x')
        WG1_height= self.Convert_to_num_of_cells(self.WG1.height,axis='y')
        WG2_width = self.Convert_to_num_of_cells(self.WG2.width,axis='x')
        WG2_height= self.Convert_to_num_of_cells(self.WG2.height,axis='y')
        return WG1_width,WG1_height,WG2_width,WG2_height

    def Calculate_index_profile(self):
        compressed_len_y = self.Compressed_shape[0]
        compressed_len_x = self.Compressed_shape[1]
        shift_x,shift_y = self.Calculate_shift

        WG1_width,WG1_height,WG2_width,WG2_height = self.WG_geometry_in_num_of_cells

        assert compressed_len_y > max(WG1_height+shift_y, WG2_height+shift_y)
        assert compressed_len_x > max(WG1_width+shift_x, WG2_width+shift_x)

        # The coordinate of WGs (unit: num of cells)
        (WG1_x,WG1_y,WG2_x,WG2_y,Ly_cladding)  = self.Structure_coordinates

        # N(x,y):The overall index matrix
        self.N = np.ones(self.Compressed_shape) * self.n_cladding
        self.N[Ly_cladding:,:] = self.n_air
        # self.N[:Ly_BOX,:] = self.n_si
        self.N[ int(WG1_y - WG1_height/2 -shift_y) :
                int(WG1_y + WG1_height/2 -shift_y),
                int(WG1_x - WG1_width/2  -shift_x) :
                int(WG1_x + WG1_width/2  -shift_x)] = self.n_core
        self.N[ int(WG2_y - WG2_height/2):
                int(WG2_y + WG2_height/2),
                int(WG2_x - WG2_width/2) :
                int(WG2_x + WG2_width/2)] = self.n_core

        # N1(x,y):The index matrix of WG1 in the absence of WG2
        self.N1 = np.ones(self.Compressed_shape) * self.n_cladding
        self.N1[Ly_cladding:,:] = self.n_air
        self.N1[int(WG1_y - WG1_height/2 -shift_y) :
                int(WG1_y + WG1_height/2 -shift_y),
                int(WG1_x - WG1_width/2  -shift_x) :
                int(WG1_x + WG1_width/2  -shift_x)] = self.n_core

        # N2(x,y):The index matrix of WG2 in the absence of WG1
        self.N2 = np.ones(self.Compressed_shape) * self.n_cladding
        self.N2[Ly_cladding:,:] = self.n_air
        self.N2[int(WG2_y - WG2_height/2):
                int(WG2_y + WG2_height/2),
                int(WG2_x - WG2_width/2) :
                int(WG2_x + WG2_width/2)] = self.n_core

    def Plot_index_profile(self,label='N',
                           save_name='./results/index_profile.pdf',dpi=300):

        colormap = "YlOrRd"
        WG1_x,WG1_y,WG2_x,WG2_y,Ly_cladding  = self.Structure_coordinates
        WG1_width,WG1_height,WG2_width,WG2_height = self.WG_geometry_in_num_of_cells
        shift_x,shift_y = self.Calculate_shift

        figsize = (10,6)
        fig = plt.figure(figsize=figsize)
        if label == 'N':
            image = self.N
        elif label == 'N1':
            image = self.N1
        else:
            image = self.N2

        im = plt.imshow(image,cmap= colormap)
        plt.gca().invert_yaxis()
        # yticks_prev = np.linspace(0,np.shape(image)[0],figsize[1])
        # xticks_prev = np.linspace(0,np.shape(image)[1],figsize[0])
        xticks_prev = np.array([0,
                                int(WG1_x - WG1_width/2-shift_x),
                                int(WG1_x + WG1_width/2-shift_x),
                                int(WG2_x - WG2_width/2),
                                int(WG2_x + WG2_width/2),
                                np.shape(image)[1]])
        yticks_prev = np.linspace(0,np.shape(image)[0],figsize[1]*2)
        xticks,yticks = self.Convert_ticks(xticks_prev,yticks_prev)
        plt.xticks(xticks_prev, xticks, rotation=0)
        plt.yticks(yticks_prev, yticks, rotation=0)
        plt.ylabel(r'Y($\mu m$)')
        plt.xlabel(r'X($\mu m$)')
        plt.title("Refractive index distribution "+label+"(x,y)")
        cbar = plt.colorbar(im, orientation='vertical', shrink=0.4, pad=0.05)
        # cbar.set_label('Refrative Index', fontsize=10)
        cbar.set_ticks(np.linspace(np.min(image),np.max(image),5))  # 自定义刻度
        cbar.ax.tick_params(labelsize=10)  # 设置刻度字体大小
        plt.savefig(save_name,dpi=dpi)
        plt.close()
        # plt.show()

    # field_coefficients = [A,B]     E_supermode = A * E_1 + B * E_2
    # E_1,E_2 are the eigenmodes of separate WGs. A,B are the coefficients
    def Plot_field_profile(self,field_coefficients,
                           field_name, title,
                           Plot_log = False,
                           save_name='./results/field_profile_',
                           dpi=400):
        fonttype = "Helvetica"
        fontsize = 18
        linewidth = 0.3
        # colormap = "jet"
        # colormap = "turbo"
        colormap = "bwr"
        cbar_num_of_pts = 5
        num_of_plots = len(field_coefficients)
        figsize =  (40, 6*num_of_plots)

        name_list = ['Abs','Re','Im']
        # name_list = ['Re']

        fig, ax = plt.subplots(num_of_plots,len(name_list),figsize=figsize,dpi=dpi)
        plt.subplots_adjust(left=0.05, right=0.95, wspace =0.1, hspace =0.2)   #调整子图间距

        for plot_idx in range(num_of_plots):
            coeffis =   field_coefficients[plot_idx]
            field   =   coeffis[0]*self.Field_dict_uncoupled[field_name][0] +\
                        coeffis[1]*self.Field_dict_uncoupled[field_name][1]
            if np.real(coeffis[0]) * np.real(coeffis[1]) < 0:
                mode_name = "Asymmetric Mode "
            else:
                mode_name = "Symmetric Mode "
            field_list = [np.abs(field),np.real(field),np.imag(field)]
            # field_list = [np.real(field)]
            if Plot_log:
                # Only calculate log for nonzero terms

                field_abs = field_list[0]
                non_zero_mask = field_abs != 0
                log_arr = np.full_like(field_abs,
                                       fill_value=np.min(field_abs[non_zero_mask]),
                                       dtype=float)
                log_arr[non_zero_mask] = np.log(field_abs[non_zero_mask])
                print(np.min(field_abs[non_zero_mask]))
                field_list[0] = log_arr
            # Converting the unit of ticks to um
            yticks_prev = np.linspace(0,np.shape(field)[0],6)
            xticks_prev = np.linspace(0,np.shape(field)[1],10)
            xticks,yticks = self.Convert_ticks(xticks_prev,yticks_prev)

            shift_x,shift_y = self.Calculate_shift
            (WG1_x,WG1_y,WG2_x,WG2_y,Ly_cladding) = self.Structure_coordinates
            WG1_width,WG1_height,WG2_width,WG2_height = self.WG_geometry_in_num_of_cells

            # coordinates of WGs (unit: num of cells)
            WG1_x_arr = [WG1_x - 0.5*WG1_width - shift_x,
                        WG1_x + 0.5*WG1_width - shift_x,
                        WG1_x + 0.5*WG1_width - shift_x,
                        WG1_x - 0.5*WG1_width - shift_x,
                        WG1_x - 0.5*WG1_width - shift_x,]
            WG1_y_arr = [WG1_y - 0.5*WG1_height - shift_y,
                        WG1_y - 0.5*WG1_height - shift_y,
                        WG1_y + 0.5*WG1_height - shift_y,
                        WG1_y + 0.5*WG1_height - shift_y,
                        WG1_y - 0.5*WG1_height - shift_y,]
            WG2_x_arr = [WG2_x - 0.5*WG2_width ,
                        WG2_x + 0.5*WG2_width ,
                        WG2_x + 0.5*WG2_width ,
                        WG2_x - 0.5*WG2_width ,
                        WG2_x - 0.5*WG2_width ]
            WG2_y_arr = [WG2_y - 0.5*WG2_height ,
                        WG2_y - 0.5*WG2_height ,
                        WG2_y + 0.5*WG2_height ,
                        WG2_y + 0.5*WG2_height ,
                        WG2_y - 0.5*WG2_height ]

            for idx in range(0,3):
                # Plot the field profile
                im = ax[plot_idx,idx].imshow(field_list[idx], cmap=colormap)
                # Plot the boundaries of the WGs
                ax[plot_idx,idx].plot(WG1_x_arr,WG1_y_arr,color='black', linewidth=linewidth)
                ax[plot_idx,idx].plot(WG2_x_arr,WG2_y_arr,color='black', linewidth=linewidth)

                ax[plot_idx,idx].set_title(mode_name+name_list[idx]+'('+field_name+')',fontsize=fontsize*1.5)
                cbar = fig.colorbar(im, ax=ax[plot_idx,idx], orientation='vertical',
                                    label='', shrink=0.8, pad=0.02)
                cbar.set_ticks(np.linspace(np.max(field_list[idx]),
                                        np.min(field_list[idx]),
                                        cbar_num_of_pts))
                cbar.ax.tick_params(labelsize=fontsize)
                ax[plot_idx,idx].set_xticks(xticks_prev)
                ax[plot_idx,idx].set_xticklabels(xticks,fontsize=fontsize)
                ax[plot_idx,idx].set_yticks(yticks_prev)
                ax[plot_idx,idx].set_yticklabels(yticks,fontsize=fontsize)
                ax[plot_idx,idx].set_xlabel(r'X($\mu m$)',fontsize=fontsize)
                ax[plot_idx,idx].set_ylabel(r'Y($\mu m$)',fontsize=fontsize)
                ax[plot_idx,idx].invert_yaxis()
                ax[plot_idx,idx].tick_params(axis='both',labelsize=fontsize)

        plt.title(title)
        plt.savefig(save_name+title+".pdf",dpi=dpi)
        plt.close()
        # plt.show()
        return

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
                Field_profile_1 = self.WG1.Field_dict[component][shift_y:,:]   #Higher part of E_1y
                Field_profile_2 = self.WG2.Field_dict[component][:-shift_y,:]  #Lower part of E_2y
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
        # plt.imshow(self.N*self.N-Nq*Nq)
        plt.gca().invert_yaxis()
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
        Chi = np.sum((self.N*self.N-Np*Np)*
                    (np.conj(self.Field_dict_uncoupled['Ex'][p]) * \
                            self.Field_dict_uncoupled['Ex'][p] + \
                    np.conj(self.Field_dict_uncoupled['Ey'][p]) * \
                            self.Field_dict_uncoupled['Ey'][p] + \
                    np.conj(self.Field_dict_uncoupled['Ez'][p]) * \
                            self.Field_dict_uncoupled['Ez'][p] ))
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

        return F

    def Find_supermodes(self):
        F = self.F_Matrix
        F_sym =  np.diag(np.ones(2))*((F[0][0] + F[1][1])) /2
        F_antisym = F - F_sym

        Eigenvalues = np.array(np.linalg.eig(F_antisym)[0]).reshape((1,2))
        Eigenvectors = np.array(np.linalg.eig(F_antisym)[1])

        print("wavelength = {:.1f} nm".format(self.wavelength))
        print("F = ",F)
        print("F_sym = ",F_sym)
        print("F_antisym = ",F_antisym)
        print("beta_ave = ",self.beta_ave)
        print("delta_beta = ",self.delta_beta)
        print("S = ",Eigenvalues)
        print("EigenVec = ",Eigenvectors)
        print("\n")
        return Eigenvalues, Eigenvectors

    def Export_kappa(self,filename_kappa="./results/straight_WG_Kappa.txt"):
        K_12 = self.Kappa(1,2,self.N2)
        K_21 = self.Kappa(2,1,self.N1)
        print("K_12 = {:.6f}".format(K_12))
        print("K_21 = {:.6f}".format(K_21))
        with open(filename_kappa,"a") as f:
            f.write("{:.1f}".format(self.wavelength) + ",{:.6f}".format(K_12) + ",{:.6f}\n".format(K_21))
