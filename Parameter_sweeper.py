'''
class Parameter_sweeper
Description:
This class is used to sweep the parameters of the coupled waveguides.
Parameters that can be swept include the wavelength and the gap between two WGs.

Parameters:
wavelengths             : format: [wavl1,wavl2,...]
gap_arr                 : format: [[gap_x1,gap_y1],[gap_x2,gap_y2],...]
foldername_WG1          : foldername of mode profile data of WG1
foldername_WG2          : foldername of mode profile data of WG2
save_foldername         : foldername to save the results
param_filename          : filename of parameters
save_D_in_csv           : True if you want to save the dispersion curve in csv format. Default is False.
plot_profile            : True if you want to plot the super modes profile. Default is False.
plot_log_scale          : True if you want to plot the super modes profile in log scale. Default is False.

Author: Weihao Xu
Date: May. 12th, 2025
'''
import numpy as np
from Coupled_Waveguides import *
from Data_analyzer import *
from Functions import *

# sweepable parameters include: wavelength, gap
class Parameter_sweeper():
    def __init__(self, wavelengths, gap_arr,
                 foldername_WG1, foldername_WG2,
                 save_foldername, param_filename,
                 save_D_in_csv = False,
                 plot_profile = False, plot_log_scale = False):

        self.wavl_arr = wavelengths             # Wavelength in Vacuum (unit:um)
        self.gap_arr = gap_arr                  # gap = [gap_x,gap_y] shape=(2,x) (unit:um)
        self.foldername_WG1 = foldername_WG1
        self.foldername_WG2 = foldername_WG2
        self.save_foldername = save_foldername
        self.param_filename = param_filename
        self.save_D_in_csv = save_D_in_csv
        self.plot_profile = plot_profile
        self.plot_log_scale = plot_log_scale

    # Sweep wavlength whne the gap between two WGs is fixed
    def Scan_wavl(self, gap, filename_uncoupled, filename_coupled):

        beta_uncoupled_arr     = [] # beta of eigenmodes of two separate WGs
        beta_ave_uncoupled_arr = [] # beta_ave of eigenmodes of two separate WGs
        beta_coupled_arr       = [] # beta of supermodes calculated with CMT
        beta_list = np.array([[]])

        with open(filename_uncoupled,'w') as f:
            f.write("wavelength,beta_1_uncoupled,beta_2_uncoupled,beta_ave_uncoupled\n")

        with open(filename_coupled,'w') as f:
            f.write("wavelength,beta_1_coupled,beta_2_coupled,coeff_supermode_1,coeff_supermode_2\n")

        gap_x, gap_y = gap
        print("gap_x = {:.3f}".format(gap_x)+" um, gap_y = {:.3f}".format(gap_y)+" um \n")

        for idx in range(len(self.wavl_arr)):
            wavl_in_nm = self.wavl_arr[idx]
            wavl_in_um = wavl_in_nm / 1000

            path_1 = self.foldername_WG1 +'/'+ "{:.0f}".format(wavl_in_nm)
            path_2 = self.foldername_WG2 +'/'+ "{:.0f}".format(wavl_in_nm)

            Coupled_WG = Coupled_Waveguides(wavelength = wavl_in_um,
                                            gap_x = gap_x, gap_y = gap_y,
                                            name1 = path_1,name2 = path_2,
                                            ModeIdx1 = 1, ModeIdx2 = 1,
                                            param_file_name=self.param_filename,
                                            Plot_field= True if idx==0 else False,
                                            Plot_index= True if idx==0 else False)

            #shape:(1,2)
            beta_uncoupled = np.array([[Coupled_WG.beta_ang_1,
                                        Coupled_WG.beta_ang_2]]) - Coupled_WG.beta_ave
            beta_coupled, coeff_of_supermodes = Coupled_WG.Find_supermodes()
            # Coupled_WG.Export_kappa()

            # For anti-sym modes, delta beta is always positive
            if beta_coupled[0,0] < 0:
                beta_coupled[0] = np.flip(beta_coupled[0])
                coeff_of_supermodes = np.flip(coeff_of_supermodes,axis=1)
            if idx == 0:
                beta_uncoupled_arr = np.array(beta_uncoupled)
                beta_coupled_arr = np.array(beta_coupled)
                beta_ave_uncoupled_arr = np.array([Coupled_WG.beta_ave])
            else:
                beta_uncoupled_arr = np.append(beta_uncoupled_arr,
                                               beta_uncoupled,axis=0)
                beta_ave_uncoupled_arr = np.append(beta_ave_uncoupled_arr,
                                                np.array([Coupled_WG.beta_ave]),axis=0)
                beta_coupled_arr = np.append(beta_coupled_arr,
                                             beta_coupled ,axis=0)

            with open(filename_uncoupled,'a') as f:
                beta_str =  str(wavl_in_um) + "," + \
                            str(np.real(beta_uncoupled[0,0])) + "," + \
                            str(np.real(beta_uncoupled[0,1])) + "," + \
                            str(np.real(Coupled_WG.beta_ave)) + "\n"
                f.write(beta_str)
            with open(filename_coupled,'a') as f:
                beta_str =  str(wavl_in_um) + "," + \
                            str(np.real(beta_coupled[0,0])) + "," + \
                            str(np.real(beta_coupled[0,1])) + "," + \
                            str(coeff_of_supermodes[0,0]) + "," +\
                            str(coeff_of_supermodes[1,0]) + "," +\
                            str(coeff_of_supermodes[0,1]) + "," +\
                            str(coeff_of_supermodes[1,1]) + "\n"
                f.write(beta_str)

        return beta_uncoupled_arr,beta_coupled_arr,beta_ave_uncoupled_arr

    # Sweep gap betwen two WGs
    # Parameters:
    # calc_needed       : must be True for initial run, can be set to False if beta files already exist.
    # foldername        : foldername to store the results.
    # num_of_wavl_pts   : number of wavelength points to calculate dispersion curve
    def Scan_gap(self, calc_needed = True, num_of_wavl_pts = 100):
        foldername = self.save_foldername
        for gap_idx in range(len(self.gap_arr)):
            gap_x, gap_y = self.gap_arr[gap_idx]

            filename_uncoupled = "beta_uncoupled_gap_"\
                                + "{:.3f}".format(gap_x) +"_" + "{:.3f}".format(gap_y)
            filename_uncoupled  = foldername + "Supermodes attributes using CMT/" + filename_uncoupled.replace(".","_") + ".txt"
            filename_coupled   = "beta_coupled_gap_" \
                                + "{:.3f}".format(gap_x) +"_" + "{:.3f}".format(gap_y)
            filename_coupled    = foldername + "Supermodes attributes using CMT/" + filename_coupled.replace(".","_")+ ".txt"

            if calc_needed:
                beta_uncoupled_arr,beta_coupled_arr,beta_ave_uncoupled_arr = self.Scan_wavl(gap = (gap_x, gap_y),
                            filename_uncoupled=filename_uncoupled,
                            filename_coupled=filename_coupled)

            # Calc dispersion curve
            Analyzer = Data_analyzer(self.wavl_arr, (gap_x, gap_y), self.param_filename,
                                     filename_uncoupled, filename_coupled,
                                     self.foldername_WG1, self.foldername_WG2,
                                     num_of_pts = num_of_wavl_pts, save_D_in_csv = self.save_D_in_csv,
                                     plot_profile = self.plot_profile, plot_log_scale = self.plot_log_scale,
                                     save_foldername = foldername)

