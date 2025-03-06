import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple


'''
class Waveguide()

Parameters:
width       : width of the waveguide. (unit:um)
height      : height the waveguide. (unit:um)
bendradius  : bend radius of the waveguide. (unit:um)
PML_len     : Width of the boundary PML layer (unit: num of data points)
name        : foldername of the mode profile
ModeIdx     : Eigen mode index
Plot        : True or False. Whether to plot all the mode profiles.

'''

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
    def __init__(self, width, height, bendradius, PML_len,
                 foldername, ModeIdx, Field_shape_padded, Plot):
        self.width  = width                 # unit :um
        self.height = height                # unit :um
        self.bend_radius = bendradius      # unit :um
        self.PML_len = PML_len              # unit: num of data points
        self.ModeIdx = ModeIdx
        self.Modes_info_list = self.Load_modes_info(foldername)

        self.Field_dict= self.Load_all_field_profile(
            foldername = foldername,ModeIdx = ModeIdx,
            Field_shape_padded=Field_shape_padded,plot= Plot,PML_len=self.PML_len)
        # print("shape of loaded field matrix:",self.Field_shape)
        self.Field_dict_normalize()     # normalize the fields to let P_total = 1


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
            #print("# of lines:",len(lines))
            #print("# of columns:",len(lines[0].split('\t')))

        Field = Field[int((PML_len)/2):np.shape(Field)[0]-int((PML_len)/2),
                      int((PML_len)/2):np.shape(Field)[1]-int((PML_len)/2)]
        Field = self.Zero_Padding_of_Field_Profile(Field, Field_shape_padded)
        self.Field_shape = np.shape(Field)
        return Field

    # Pad zero. Lx_padded and Ly_padded in unit of num of points
    def Zero_Padding_of_Field_Profile(self,Field,shape_padded):
        Ly_padded, Lx_padded = shape_padded
        Field_padded = np.zeros(shape_padded,dtype=np.complex64)
        Ly_unpad, Lx_unpad = np.shape(Field)
        idx_y_start = int(Ly_padded/2-Ly_unpad/2)
        idx_y_end   = int(Ly_padded/2+Ly_unpad/2)
        idx_x_start = int(Lx_padded/2-Lx_unpad/2)
        idx_x_end   = int(Lx_padded/2+Lx_unpad/2)
        if not idx_y_end-idx_y_start == Ly_unpad:
            idx_y_end   = idx_y_end + 1
        if not idx_x_end-idx_x_start == Lx_unpad:
            idx_x_end   = idx_x_end + 1
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
        Field_dict = dict(zip(component_name_list,component_list))
        if plot ==True:
            self.Plot_all_field_profile(component_list)

        return Field_dict


    def Plot_all_field_profile(self,component_list,
                               save_name='./results/all_field_profile.png',dpi=150):
        #Plot parameters
        font = {'family': 'serif',
                'serif': 'Helvetica',
                'weight': 'normal',
                'size': 4}
        plt.rc('font', **font)
        grid_linewidth = 1

        fig, ax = plt.subplots(4,3,figsize=(4,5),dpi=dpi)

        ax = ax.flatten()
        plt.subplots_adjust(wspace =1, hspace =0.8)   #调整子图间距
        plt.ylabel('Y')
        plt.xlabel('X')
        plt.grid(linewidth=grid_linewidth, alpha=0.3)
        # plt.yticks(fontproperties = font["type"], size = font["size"])
        # plt.xticks(fontproperties = font["type"], size = font["size"])
        # plt.rcParams["font.family"] = font["type"]
        # plt.rcParams.update({'font.size': fontsize})
        # plt.ylabel('Y', fontdict={'family' : fonttype, 'size' : fontsize})
        # plt.xlabel('X', fontdict={'family' : fonttype, 'size' : fontsize})

        for i in range(4):
            for j in range(3):
                if i%2 == 0:
                    im = ax[i*3+j].imshow(np.real(component_list[j if i<2 else j+3]))
                    fig.colorbar(im,orientation='vertical')
                    ax[i*3+j].set_title('Re('+self.component_name_list[j if i<2 else j+3]+')')

                else:
                    im = ax[i*3+j].imshow(np.imag(component_list[j if i<2 else j+3]))
                    fig.colorbar(im,orientation='vertical')
                    ax[i*3+j].set_title('Im('+self.component_name_list[j if i<2 else j+3]+')')

                ax[i*3+j].set_xlabel("X")
                ax[i*3+j].set_ylabel("Y")

                ax[i*3+j].invert_yaxis()
                ax[i*3+j].tick_params(axis='both',labelsize=5)


        plt.savefig(save_name, dpi=dpi)
        plt.show()
        return

    # The optical power carried by the eigen mode in the waveguide
    def P_total(self):
        Field_dict = self.Field_dict
        P_total = 2*np.real(np.sum(np.conj(Field_dict['Ex']) * Field_dict['Hy'] -
                           np.conj(Field_dict['Ey']) * Field_dict['Hx']))
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

