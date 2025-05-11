import numpy as np
from Coupled_Waveguides import *
from Data_analyzer import *
from Functions import *
'''
class Waveguide()
Parameters:

wavelengths             : format: [wavl1,wavl2,...]
gap_arr                 : format: [[gap_x1,gap_y1],[gap_x2,gap_y2],...]
foldername1             : foldername of mode profile data of WG1
foldername2             : foldername of mode profile data of WG2
param_filename          : filename of parameters

Methods:

'''
# sweepable parameters include: wavelength, gap
class Parameter_sweeper():
    def __init__(self,wavelengths,gap_arr,
                 foldername1,foldername2,
                 param_filename="./Param.csv"):

        self.wavl_arr = wavelengths     # Wavelength in Vacuum (unit:um)
        self.gap_arr = gap_arr                  # gap = [gap_x,gap_y] shape=(2,x) (unit:um)
        self.foldername_1 = foldername1
        self.foldername_2 = foldername2
        self.param_filename = param_filename
        self.beta_array = np.array([[],[]])

    # Sweep wavlength whne the gap between two WGs is fixed
    def Scan_wavl(self, gap,
                  filename_uncoupled = "../data/beta_uncoupled.txt",
                  filename_coupled   = "../data/beta_coupled.txt"):
                #   plot_field_profile = False, plot_log = False):

        beta_uncoupled_arr     = [] # beta of eigenmodes of two separate WGs
        beta_ave_uncoupled_arr = [] # beta_ave of eigenmodes of two separate WGs
        beta_coupled_arr       = [] # beta of supermodes calculated with CMT
        beta_list = np.array([[]])

        with open(filename_uncoupled,'w') as f:
            f.write("wavelength,beta_1_uncoupled,beta_2_uncoupled,beta_ave_uncoupled\n")

        with open(filename_coupled,'w') as f:
            f.write("wavelength,beta_1_coupled,beta_2_coupled,coeff_supermode_1,coeff_supermode_2\n")

        for idx in range(len(self.wavl_arr)):
            wavl_in_nm = self.wavl_arr[idx]
            wavl_in_um = wavl_in_nm / 1000
            gap_x, gap_y = gap
            path_1 = self.foldername_1 +'/'+ "{:.0f}".format(wavl_in_nm)
            path_2 = self.foldername_2 +'/'+ "{:.0f}".format(wavl_in_nm)
            print("####################")
            print("gap_x = {:.3f}".format(gap_x)+" um, gap_y = {:.3f}".format(gap_y)+" um")
            print("####################\n")
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
            # Plot field profile of supermodes
            # if plot_field_profile:
            #     CoupledWG.Plot_field_profile(coeff_of_supermodes,field_name = 'Ex',
            #                         title=r"Electric field profile when wavl = "+"{:.0f}".format(wavl_in_nm)+" nm", Plot_log=plot_log)
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
    # foldername        : foldername to store the beta files.
    # num_of_wavl_pts   : number of wavelength points to calculate dispersion curve
    def Scan_gap(self, calc_needed = True, foldername="./Mode number conserved coupling/results/", num_of_wavl_pts = 100):
        for gap_idx in range(len(self.gap_arr)):
            gap_x, gap_y = self.gap_arr[gap_idx]

            filename_uncoupled = "beta_uncoupled_gap_"\
                                + "{:.3f}".format(gap_x) +"_" + "{:.3f}".format(gap_y)
            filename_uncoupled  = foldername + filename_uncoupled.replace(".","_") + ".txt"
            filename_coupled   = "beta_coupled_gap_" \
                                + "{:.3f}".format(gap_x) +"_" + "{:.3f}".format(gap_y)
            filename_coupled    = foldername + filename_coupled.replace(".","_")+ ".txt"

            if calc_needed:
                beta_uncoupled_arr,beta_coupled_arr,beta_ave_uncoupled_arr = self.Scan_wavl(gap = (gap_x, gap_y),
                            filename_uncoupled=filename_uncoupled,
                            filename_coupled=filename_coupled,
                            plot = False)
            # Calc dispersion curve
            Analyzer = Data_analyzer(self.wavl_arr, (gap_x, gap_y), self.param_filename,
                                     filename_uncoupled, filename_coupled,
                                     self.foldername_1, self.foldername_2,
                                     num_of_pts = num_of_wavl_pts,
                                     plot_profile = True, plot_log_scale = False)

