'''
class Waveguide: Represent a single waveguide.

Parameters:
width               : width of the waveguide. (unit:um)
height              : height the waveguide. (unit:um)
bendradius          : bend radius of the waveguide. (unit:um)
FDE_x               : geometry of the FDE region. format:(x_min,x_max), unit: um
FDE_y               : geometry of the FDE region. format:(y_min,y_max), unit: um
PML_len             : Width of the boundary PML layer (unit: num of data points)
Num_of_cells        : num of data points in the FDE region. format:(Num_of_cells_x,Num_of_cells_y)
name                : foldername of the mode profile
ModeIdx             : Eigen mode index
FDE_shape_padded    : Shape of the field matrix afer zero padding (unit: um)
Plot                : True or False. Whether to plot all the mode profiles.

Methods:

'''

import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
MIN_VALUE = np.complex128(10**(-20),0)      # value used as zero when padding the field profile


class Waveguide():

    #constants
    c       =   299792458
    u0      =   4 * np.pi * 1e-7
    component_name_list = ['Ex','Ey','Ez','Hx','Hy','Hz']
    PML_len = 23        # num of points for the PML layer

    # Struct for mode info
    Mode = namedtuple('Mode', ['modeidx','neff','ng', 'loss', 'polarization','beta_ang'])

    #width:  num of data points for |x|<a
    #height: num of data points for |y|<b
    def __init__(self, width, height, bendradius,           # geometry of the WG
                 FDE_x, FDE_y, PML_len, Num_of_cells,       # geometry of the FDE region
                 foldername, ModeIdx, FDE_shape_padded,     # parameters of the mode profile
                 Plot):
        # geometry of the WG
        self.width  = width                 # unit :um
        self.height = height                # unit :um
        self.bend_radius = bendradius       # unit :um

        # geometry of the FDE region
        self.FDE_x_min,self.FDE_x_max = FDE_x   # format:(x_min,x_max), unit: um
        self.FDE_y_min,self.FDE_y_max = FDE_y   # format:(y_min,y_max), unit: um
        self.FDE_width = self.FDE_x_max - self.FDE_x_min    # unit: um
        self.FDE_height = self.FDE_y_max - self.FDE_y_min   # unit: um
        self.PML_len = PML_len                              # unit: num of data points
        self.Num_of_cells_x,self.Num_of_cells_y = Num_of_cells  # unit: num of data points

        # parameters of the mode profile
        self.name = foldername
        self.ModeIdx = ModeIdx
        self.Modes_info_list = self.Load_modes_info(foldername)
        self.FDE_width_padded,self.FDE_height_padded = FDE_shape_padded # unit: um

        # unit: num of cells
        Field_shape_padded = np.array([self.Convert_to_num_of_cells(
                                        self.FDE_height_padded,axis='y'),
                                        self.Convert_to_num_of_cells(
                                        self.FDE_width_padded,axis='x')])
        self.Load_all_field_profile(
            foldername = foldername, ModeIdx = ModeIdx,
            Field_shape_padded=Field_shape_padded,plot= Plot,PML_len=self.PML_len)
        # print("shape of loaded field matrix:",self.Field_shape)

    def str2complex(self,s):
        str = s.replace(" ","").replace("i","j")
        return complex(str)

    '''
    Note: All field profile matrixes use the first index for Y, the second index for X
    E.g. Ex[100:][:] select the upper part of Ex, Ey[:][:-100] select the left part of Ey
    '''
    def Load_modes_info(self,foldername):
        info_file_name = foldername + "/Mode_info.txt"
        Modes_info = []
        # dtype = [('modeidx', 'i2'), ('neff', str),
        #          ('ng', str), ('loss', 'f8'), ('polarization','f8')]
        # data_read = np.loadtxt(info_file_name, delimiter=',', dtype=dtype,skiprows=1)
        with open(info_file_name,'r') as f:
            data_read = f.readlines()
            for line in data_read[1:]:
                modeidx = int(line.split(',')[0])
                neff = self.str2complex(line.split(',')[1])
                ng = self.str2complex(line.split(',')[2])
                loss = float(line.split(',')[3])
                polarization = float(line.split(',')[4])
                beta_ang = float(line.split(',')[5])
                mode = self.Mode(modeidx, neff, ng, loss, polarization, beta_ang)
                Modes_info.append(mode)
        return Modes_info

    def Load_field_profile(self,foldername,ModeIdx,component,
                           Field_shape_padded,PML_len = 23):
        filename = './' + foldername + '/Mode' + str(ModeIdx) + '_' + component + '.txt'
        Field = np.array([[]])
        with open(filename) as f:
            lines = f.readlines()
            Field = np.array([self.str2complex(s) for s in lines[0].split('\t')])
            for line in lines[1:]:
                line_arr = line.split('\t')
                Field = np.c_[Field,np.array([self.str2complex(s) for s in line_arr])]

        Field = Field[int((PML_len)/2):np.shape(Field)[0]-int((PML_len)/2),
                      int((PML_len)/2):np.shape(Field)[1]-int((PML_len)/2)]
        Field = self.Zero_Padding_of_Field_Profile(Field, Field_shape_padded)
        self.Field_shape = np.shape(Field)
        return Field

    # Pad zero. Lx_padded and Ly_padded in unit of num of points
    # shape_padded : format:(Ly_padded, Lx_padded)
    def Zero_Padding_of_Field_Profile(self,Field,shape_padded):
        Ly_padded, Lx_padded = shape_padded
        Ly_unpad, Lx_unpad = np.shape(Field)

        if Ly_padded <= Ly_unpad:
            # No need to pad zero
            idx_y_start = 0
            idx_y_end   = Ly_unpad
            Ly_padded   = Ly_unpad
        else:
            idx_y_start = int(Ly_padded/2-Ly_unpad/2)
            idx_y_end   = int(Ly_padded/2+Ly_unpad/2)
        if Lx_padded <= Lx_unpad:
            # No need to pad zero
            idx_x_start = 0
            idx_x_end   = Lx_unpad
            Lx_padded   = Lx_unpad
        else:
            idx_x_start = int(Lx_padded/2-Lx_unpad/2)
            idx_x_end   = int(Lx_padded/2+Lx_unpad/2)

        Field_padded = np.ones((Ly_padded, Lx_padded),
                               dtype=np.complex128) * np.max([MIN_VALUE,np.min(Field)])

        if not idx_y_end-idx_y_start == Ly_unpad:
            idx_y_end   = idx_y_end + 1
            print("Y Warning when padding zero!")
        if not idx_x_end-idx_x_start == Lx_unpad:
            idx_x_end   = idx_x_end + 1
            print("X Warning when padding zero!")
        assert idx_y_end-idx_y_start == Ly_unpad
        assert idx_x_end-idx_x_start == Lx_unpad

        Field_padded[idx_y_start:idx_y_end,
                     idx_x_start:idx_x_end] = Field
        return Field_padded

    def Load_all_field_profile(self,foldername,ModeIdx,
                               Field_shape_padded,plot=True,PML_len = 23):

        Ex = self.Load_field_profile(foldername,ModeIdx,'Ex',Field_shape_padded,PML_len)
        Ey = self.Load_field_profile(foldername,ModeIdx,'Ey',Field_shape_padded,PML_len)
        Ez = self.Load_field_profile(foldername,ModeIdx,'Ez',Field_shape_padded,PML_len)
        Hx = self.Load_field_profile(foldername,ModeIdx,'Hx',Field_shape_padded,PML_len)
        Hy = self.Load_field_profile(foldername,ModeIdx,'Hy',Field_shape_padded,PML_len)
        Hz = self.Load_field_profile(foldername,ModeIdx,'Hz',Field_shape_padded,PML_len)

        component_list = [Ex,Ey,Ez,Hx,Hy,Hz]
        component_name_list = ['Ex','Ey','Ez','Hx','Hy','Hz']
        self.Field_dict = dict(zip(component_name_list,component_list))
        self.Field_dict_normalize()     # normalize the fields to let P_total = 1

        if plot ==True:
            component_list = list(self.Field_dict.values())
            self.Plot_all_field_profile(component_list,(self.FDE_width*0.8,self.FDE_height*0.8),
                                        save_name='./results/'+foldername[19:25]+'_all_field_profile.jpg')
        return

    # plot_Log        : Whether to plot Log(E)
    # plot_field_size : Size of the region to be plotted (unit: um)
    def Plot_all_field_profile(self,component_list, plot_field_size ,plot_log = False,
                               save_name='./results/all_field_profile.jpg',dpi=300):
        #Plot parameters
        # save_name='./results/All_field_profile_'+self.name+'.jpg'
        font = {'family': 'serif',
                'serif': 'Helvetica',
                'weight': 'normal',
                'size': 10}
        plt.rc('font', **font)
        grid_linewidth = 1
        colormap = "jet"
        figsize = (30,15)

        fig, ax = plt.subplots(4,3,figsize=figsize,dpi=dpi)
        ax = ax.flatten()
        plt.subplots_adjust(left=0.05, right=0.95,
                            bottom=0.05,top=0.95,
                            wspace =0.1, hspace =0.2)   #调整子图间距

        plt.grid(linewidth=grid_linewidth, alpha=0.3)

        for i in range(4):
            for j in range(3):
                x_size = self.Convert_to_num_of_cells(plot_field_size[0],axis='x')
                y_size = self.Convert_to_num_of_cells(plot_field_size[1],axis='y')

                field = component_list[j if i<2 else j+3]
                field = field[int(np.shape(field)[0]/2-y_size/2) : int(np.shape(field)[0]/2+y_size/2),
                              int(np.shape(field)[1]/2-x_size/2) : int(np.shape(field)[1]/2+x_size/2)]
                if plot_log:
                    field = np.log(np.abs(field))

                yticks_prev = np.arange(0,np.shape(field)[0],
                                        int(np.shape(field)[0]/figsize[1]))
                xticks_prev = np.arange(0,np.shape(field)[1],
                                        int(np.shape(field)[1]/figsize[0]))
                xticks,yticks = self.Convert_ticks(xticks_prev,yticks_prev,plot_field_size)

                if i%2 == 0:
                    im = ax[i*3+j].imshow(np.real(field),cmap=colormap)
                    ax[i*3+j].set_title('Re('+self.component_name_list[j if i<2 else j+3]+')')
                else:
                    im = ax[i*3+j].imshow(np.imag(field),cmap=colormap)
                    ax[i*3+j].set_title('Im('+self.component_name_list[j if i<2 else j+3]+')')
                cbar = fig.colorbar(im, ax=ax[i*3+j], orientation='vertical',
                            label='', shrink=0.6, pad=0.02)

                ax[i*3+j].set_xticks(xticks_prev)
                ax[i*3+j].set_xticklabels(xticks)
                ax[i*3+j].set_yticks(yticks_prev)
                ax[i*3+j].set_yticklabels(yticks)
                ax[i*3+j].set_xlabel(r'X($\mu m$)',fontsize=4)
                ax[i*3+j].set_ylabel(r'Y($\mu m$)',fontsize=4)
                ax[i*3+j].invert_yaxis()
                ax[i*3+j].tick_params(axis='both',labelsize=5)

        plt.savefig(save_name, dpi=dpi)
        plt.close()
        # plt.show()
        return

    # Converting the unit of ticks to um
    def Convert_ticks(self,xticks_prev,yticks_prev,plot_field_size,shift_tuple=(0,0)):

        shift_x,shift_y = shift_tuple
        yticks      = self.Convert_to_um(yticks_prev+ shift_y,axis='y') - plot_field_size[1]/2
                        #self.FDE_y_min - (self.FDE_height_padded-self.FDE_height)/2
        yticks      = np.round(yticks, 2)
        xticks      = self.Convert_to_um(xticks_prev+ shift_x,axis='x') - plot_field_size[0]/2
                        #self.FDE_x_min - (self.FDE_width_padded-self.FDE_width)/2
        xticks      = np.round(xticks, 2)
        return xticks,yticks

    def Convert_to_num_of_cells(self,len_in_um,axis):
        if axis=='x':
            return int(len_in_um * self.Num_of_cells_x / self.FDE_width)
        else:
            return int(len_in_um * self.Num_of_cells_y / self.FDE_height)

    def Convert_to_um(self,num_cells,axis):
        if axis=='x':
            return self.FDE_width * num_cells / self.Num_of_cells_x
        else:
            return self.FDE_height * num_cells / self.Num_of_cells_y

    # The optical power carried by the eigen mode in the waveguide
    def P_total(self):
        Field_dict = self.Field_dict
        P_total = 2*np.real(np.sum(np.conj(Field_dict['Ex']) * Field_dict['Hy'] -
                           np.conj(Field_dict['Ey']) * Field_dict['Hx']))
        pass
        return P_total

    # Normalize fields to make P_total = 1
    def Field_dict_normalize(self):
        P_total = self.P_total()
        component_list = list(self.Field_dict.values())
        component_list_normalized = np.copy(component_list)
        for idx in range(np.shape(component_list)[0]):
            component_list_normalized[idx] = component_list[idx] / np.sqrt(P_total)
        self.Field_dict = dict(zip(self.component_name_list,component_list_normalized))
        return

