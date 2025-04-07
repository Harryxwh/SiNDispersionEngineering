'''
Main program for the simulation of two coupled microrings or bent waveguides
'''
def Plot_field_profile(component,component_name,
                       save_name='./results/field_profile.png',dpi=150):
    fonttype = "Helvetica"
    fontsize = 4
    grid_linewidth = 1

    fig, ax = plt.subplots(1,3,figsize=(8, 4),dpi=150)
    plt.subplots_adjust(wspace =0.6, hspace =0.5)   #调整子图间距

    plt.yticks(fontproperties = fonttype, size = fontsize)
    plt.xticks(fontproperties = fonttype, size = fontsize)
    plt.rcParams["font.family"] = fonttype
    plt.rcParams.update({'font.size': fontsize})
    plt.ylabel('Y', fontdict={'family' : fonttype, 'size' : fontsize})
    plt.xlabel('X', fontdict={'family' : fonttype, 'size' : fontsize})
    plt.legend()

    name_list = ['Abs','Re','Im']
    component_list = [np.abs(component),np.real(component),np.imag(component)]
    for idx in range(3):
        im = ax[idx].imshow(component_list[idx])
        ax[idx].set_title(name_list[idx]+'('+component_name+')')
        fig.colorbar(im,fraction=0.08, pad=0.05,orientation='vertical')
        ax[idx].set_xlabel("X")
        ax[idx].set_ylabel("Y")
        ax[idx].invert_yaxis()
        ax[idx].tick_params(axis='both',labelsize=5)
    plt.savefig(save_name,dpi=dpi)
    plt.show()
    return

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

'''
class Coupled_Waveguides()

Parameters:
n1              : refractive index of the waveguide core
n0              : refractive index of the waveguide cladding
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

from Coupled_Waveguides import *

# sweepable parameters include: wavelength, gap
class Parameter_sweeper():
    #The material dispersion of waveguide and cladding

    def __init__(self,wavelengths,gap,foldername1,foldername2):
        index_file1 = "Si3N4_index.txt"
        index_file2 = "SiO2_index.txt"
        self.n1_arr, self.n0_arr = self.Load_index_data(index_file1,index_file2)
        self.wavl_vec = wavelengths     # Wavelength in Vacuum (unit:um)
        self.gap = gap                  # gap = [gap_x,gap_y] shape=(2,x) (unit:um)
        self.foldername_1 = foldername1
        self.foldername_2 = foldername2
        self.beta_array = np.array([[],[]])

    def Load_index_data(self,filename_n1,filename_n0):
        n1_arr = []
        n0_arr = []
        with open(filename_n1,'r') as f_n1:
            data = f_n1.readlines()
            for line in data:
                line = line.replace('\n','')
                wavl = int(line.split(' ')[0].strip())
                n_re = float(line.split(' ')[1].strip())
                n_im = float(line.split(' ')[2].strip())
                #n1_arr.append(complex(n_re,n_im))
                n1_arr.append(n_re)
        with open(filename_n0,'r') as f_n0:
            data = f_n0.readlines()
            for line in data:
                wavl = int(line.split()[0])
                n_re = float(line.split()[1])
                n_im = float(line.split()[2])
                #n0_arr.append(complex(n_re,n_im))
                n0_arr.append(n_re)
        n1_arr = np.array(n1_arr)
        n0_arr = np.array(n0_arr)
        return n1_arr, n0_arr

    def Create_instance(self,wavl_idx,gap_idx,Plot_field = False,Plot_index = False):
        assert wavl_idx < min(len(self.wavl_vec),len(self.n0_arr),len(self.n1_arr))
        assert gap_idx < len(self.gap)
        wavl = self.wavl_vec[wavl_idx]
        n0 = self.n0_arr[wavl_idx]
        n1 = self.n1_arr[wavl_idx]
        gap_x, gap_y = self.gap[gap_idx]
        path_1 = self.foldername_1 +'/'+ str(int(wavl*1000))
        path_2 = self.foldername_2 +'/'+ str(int(wavl*1000))
        Coupled_WG = Coupled_Waveguides(n1, n0, gap_x = gap_x, gap_y = gap_y,
                                        wavelength = wavl,
                                        name1 = path_1,
                                        name2 = path_2,
                                        ModeIdx1 = 1, ModeIdx2 = 1,
                                        param_file_name="./Param.csv",
                                        Plot_field=Plot_field, Plot_index=Plot_index)
        return Coupled_WG


    def Plot_curve(self,X,Y_arr,Y_legends,
                   X_label,Y_label,
                   dpi=150):

        #Plot parameters
        fonttype = "Helvetica"
        fontsize = 10
        grid_linewidth = 0.8
        plot_linewidth = 1.5
        colors = ['tab:blue', 'tab:orange', 'tab:red', 'tab:green']
        plt.figure(figsize=(5,4))
        for idx in range(np.shape(Y_arr)[1]):
            plt.plot(X,Y_arr[:,idx],label=Y_legends[idx],
                     color=colors[idx], marker='o',
                     linestyle='-', linewidth=plot_linewidth)
        plt.yticks(fontproperties = fonttype, size = fontsize)
        plt.xticks(fontproperties = fonttype, size = fontsize)
        plt.rcParams["font.family"] = fonttype
        plt.rcParams.update({'font.size': fontsize})
        plt.ylabel(Y_label, fontdict={'family' : fonttype, 'size' : fontsize})
        plt.xlabel(X_label, fontdict={'family' : fonttype, 'size' : fontsize})
        plt.legend()
        plt.grid(linewidth=grid_linewidth, alpha=0.3)
        savename = "results/"+str(Y_label)+".png"
        plt.savefig(savename,dpi=dpi)
        plt.tight_layout()
        plt.show()

    def Scan_wavl(self,gap_idx):
        beta_uncoupled_arr = []
        beta_coupled_arr   = []
        beta_list = np.array([[]])
        for idx in range(len(self.wavl_vec)):
            Coupled_WG = self.Create_instance(wavl_idx=idx,gap_idx= gap_idx,
                                              Plot_field= True if idx==0 else False,
                                              Plot_index= True if idx==0 else False)
            beta_uncoupled = np.array([[Coupled_WG.beta_1,
                                        Coupled_WG.beta_2]]) - Coupled_WG.beta_ave
            Rabi_freq, beta_coupled = Coupled_WG.Find_supermodes()
            if idx == 0:
                beta_uncoupled_arr = np.array(beta_uncoupled)
                beta_coupled_arr = np.array(beta_coupled)
            else:
                beta_uncoupled_arr = np.append(beta_uncoupled_arr,
                                               beta_uncoupled,axis=0)
                beta_coupled_arr = np.append(beta_coupled_arr,
                                               beta_coupled ,axis=0)

        return beta_uncoupled_arr, beta_coupled_arr

    def Dispersion(self):
        fre_vec = self.c / self.wavl_vec
        Diff_fre = np.diff(fre_vec,axis=0)
        Diff_beta = np.diff(self.beta_array,axis=1)
        FOD = np.multiply(Diff_beta, 1/(2 * np.pi) * np.reciprocal(Diff_fre))
        self.Plot_curve(self.wavl_vec, FOD, "FOD")

    #Load results from Di Yu
    def Save_results(self,plot=True):
        # filename = '../FDE/results/dispersion_ring_resonator.txt'
        # n_eff_ring = np.loadtxt(filename,delimiter=',',dtype=float,skiprows=1)
        # filename = '../FDE/results/dispersion_bus_waveguide.txt'
        # n_eff_bus  = np.loadtxt(filename,delimiter=',',dtype=float,skiprows=1)
        beta_uncoupled_arr, beta_coupled_arr  = self.Scan_wavl(gap_idx=0)
        filename_uncoupled = "results/beta_uncoupled.txt"
        filename_coupled = "results/beta_coupled.txt"
        with open(filename_uncoupled,'w') as f:
            for wavl_idx in range(len(self.wavl_vec)):
                beta_str = str(np.real(beta_uncoupled_arr[wavl_idx,0])) + "," + \
                            str(np.real(beta_uncoupled_arr[wavl_idx,1])) + "\n"
                f.write(beta_str)

        with open(filename_coupled,'w') as f:
            for wavl_idx in range(len(self.wavl_vec)):
                beta_str = str(np.real(beta_coupled_arr[wavl_idx,0])) + "," + \
                            str(np.real(beta_coupled_arr[wavl_idx,1])) + "\n"
                f.write(beta_str)

        if plot == True:
            self.Plot_curve(self.wavl_vec,
                            np.c_[beta_uncoupled_arr,beta_coupled_arr],
                            ['beta_uncoupled_1','beta_uncoupled_2',
                            'beta_coupled_1','beta_coupled_2'],
                            'wavelength(um)','beta - beta_ave(rad\\rad)')

