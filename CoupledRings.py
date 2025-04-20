'''
Main program for the simulation of two coupled microrings or bent waveguides
'''

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



