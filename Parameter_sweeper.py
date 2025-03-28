import numpy as np
from Coupled_Waveguides import *
from Data_analyzer import *
'''
class Waveguide()
Parameters:

wavelengths             : format: [wavl1,wavl2,...]
gap                     : format: [[gap_x1,gap_y1],[gap_x2,gap_y2],...]
index_file1             : filename of refractive index of core material (SiN)
index_file2             : filename of refractive index of cladding material (SiO2)
foldername1             : foldername of mode profile data of WG1
foldername2             : foldername of mode profile data of WG2
param_filename          : filename of parameters

Methods:

'''
# sweepable parameters include: wavelength, gap
class Parameter_sweeper():
    def __init__(self,wavelengths,gap,
                 index_file1, index_file2,
                 foldername1,foldername2,
                 param_filename="./Param.csv"):
        #The material dispersion of waveguide and cladding
        index_file1 = index_file1
        index_file2 = index_file2
        self.n1_dict, self.n0_dict = self.Load_index_data(index_file1,index_file2)

        self.wavl_arr = wavelengths     # Wavelength in Vacuum (unit:um)
        self.gap = gap                  # gap = [gap_x,gap_y] shape=(2,x) (unit:um)
        self.foldername_1 = foldername1
        self.foldername_2 = foldername2
        self.param_filename = param_filename
        self.beta_array = np.array([[],[]])

    def Load_index_data(self,filename_n1,filename_n0):
        n1_dict = {}
        n0_dict = {}
        # with open(filename_n1,'r') as f_n1:
        #     data = f_n1.readlines()
        #     for line in data:
        #         line = line.replace('\n','')
        #         wavl = int(line.split(' ')[0].strip())
        #         n_re = float(line.split(' ')[1].strip())
        #         n_im = float(line.split(' ')[2].strip())
        #         #n1_dict.append(complex(n_re,n_im))
        #         n1_dict[wavl] = n_re
        # with open(filename_n0,'r') as f_n0:
        #     data = f_n0.readlines()
        #     for line in data:
        #         wavl = int(line.split()[0])
        #         n_re = float(line.split()[1])
        #         n_im = float(line.split()[2])
        #         #n0_dict.append(complex(n_re,n_im))
        #         n0_dict[wavl] = n_re
        wavl_nm = np.linspace(1000,2000,1001).astype(int)
        wavl_um = wavl_nm / 1000
        n_Si3N4 = (1+3.0249/(1-(0.1353406/wavl_um)**2)+
                  40314/(1-(1239.842/wavl_um)**2))**.5
        n_SiO2  = (1+0.6961663/(1-(0.0684043/wavl_um)**2)+
                   0.4079426/(1-(0.1162414/wavl_um)**2)+
                   0.8974794/(1-(9.896161/wavl_um)**2))**.5
        n1_dict = dict(zip(wavl_nm,n_Si3N4))
        n0_dict = dict(zip(wavl_nm,n_SiO2))
        return n1_dict, n0_dict

    def Create_instance(self,wavl_idx,gap_idx,param_filename,
                        Plot_field = False,Plot_index = False):
        assert wavl_idx < len(self.wavl_arr)
        assert gap_idx < len(self.gap)
        wavl = self.wavl_arr[wavl_idx]
        assert wavl in self.n0_dict
        assert wavl in self.n1_dict
        n0 = self.n0_dict[wavl]
        n1 = self.n1_dict[wavl]
        gap_x, gap_y = self.gap[gap_idx]
        path_1 = self.foldername_1 +'/'+ "{:.0f}".format(wavl)
        path_2 = self.foldername_2 +'/'+ "{:.0f}".format(wavl)
        print("#################")
        print("gap_x = {:.3f}".format(gap_x)+" um, gap_y = {:.3f}".format(gap_y)+" um")
        print("#################\n")
        Coupled_WG = Coupled_Waveguides(n1, n0, gap_x = gap_x, gap_y = gap_y,
                                        wavelength = wavl,
                                        name1 = path_1,
                                        name2 = path_2,
                                        ModeIdx1 = 1, ModeIdx2 = 1,
                                        param_file_name=param_filename,
                                        Plot_field=Plot_field, Plot_index=Plot_index)
        return Coupled_WG


    def Plot_curve(self,X,Y_arr,Y_legends,
                   X_label,Y_label,
                   dpi=300):

        #Plot parameters
        fonttype = "Helvetica"
        fontsize = 10
        grid_linewidth = 0.8
        plot_linewidth = 1.5
        colors = ['tab:blue', 'tab:orange', 'tab:red', 'tab:green']
        plt.figure(figsize=(8,6))
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
        plt.tight_layout()

        savename = "results/"+str(Y_label)+".png"
        plt.savefig(savename,dpi=dpi)
        # plt.show()

    def Scan_wavl(self, gap_idx,
                  filename_uncoupled = "data/beta_uncoupled.txt",
                  filename_coupled   = "data/beta_coupled.txt",
                  plot = False):

        beta_uncoupled_arr     = [] # beta of eigenmodes of two separate WGs
        beta_ave_uncoupled_arr = [] # beta_ave of eigenmodes of two separate WGs
        beta_coupled_arr       = [] # beta of supermodes calculated with CMT
        beta_list = np.array([[]])

        with open(filename_uncoupled,'w') as f:
            f.write("wavelength,beta_1_uncoupled,beta_2_uncoupled,beta_ave_uncoupled\n")

        with open(filename_coupled,'w') as f:
            f.write("wavelength,beta_1_coupled,beta_2_coupled,coeff_supermode_1,coeff_supermode_2\n")

        for idx in range(len(self.wavl_arr)):
            wavl_in_um = self.wavl_arr[idx] / 1000
            Coupled_WG = self.Create_instance(wavl_idx=idx, gap_idx=gap_idx,
                                              param_filename=self.param_filename,
                                              Plot_field= True if idx==0 else False,
                                              Plot_index= True if idx==0 else False)
            #shape:(1,2)
            beta_uncoupled = np.array([[Coupled_WG.beta_ang_1,
                                        Coupled_WG.beta_ang_2]]) - Coupled_WG.beta_ave
            beta_coupled, coeff_of_supermodes = Coupled_WG.Find_supermodes()
            if beta_coupled[0,0] < 0:
                beta_coupled[0] = np.flip(beta_coupled[0])
                coeff_of_supermodes = np.flip(coeff_of_supermodes,axis=0)
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
                            str(coeff_of_supermodes[0,1]) + "," +\
                            str(coeff_of_supermodes[1,0]) + "," +\
                            str(coeff_of_supermodes[1,1]) + "\n"
                f.write(beta_str)
        if plot == True:
            self.Plot_curve(self.wavl_arr,
                            np.c_[beta_uncoupled_arr,beta_coupled_arr],
                            ['beta_uncoupled_1','beta_uncoupled_2',
                            'beta_supermode_1(CMT)','beta_supermode_1(CMT)'],
                            'wavelength(um)',r'$\beta$ - $\beta_{ave}(rad/rad)$')
        return beta_uncoupled_arr,beta_coupled_arr,beta_ave_uncoupled_arr

    def Scan_gap(self):
        for gap_idx in range(len(self.gap)):
            gap_x, gap_y = self.gap[gap_idx]

            filename_uncoupled = "data/beta_uncoupled_gap_"\
                                + "{:.1f}".format(gap_x) +"_" + "{:.1f}".format(gap_y)
            filename_coupled   = "data/beta_coupled_gap_" \
                                + "{:.1f}".format(gap_x) +"_" + "{:.1f}".format(gap_y)
            filename_coupled = filename_coupled.replace(".","_")+ ".txt"
            filename_uncoupled = filename_uncoupled.replace(".","_") + ".txt"
            beta_uncoupled_arr,beta_coupled_arr,beta_ave_uncoupled_arr = self.Scan_wavl(gap_idx=gap_idx,
                                                    filename_uncoupled=filename_uncoupled,
                                                    filename_coupled=filename_coupled,
                                                    plot = False)
            # Calc dispersion curve
            Analyzer = Data_analyzer(self.wavl_arr, (gap_x, gap_y),
                                     filename_uncoupled, filename_coupled,
                                     self.param_filename)



    # Save supermodes results calculated using coupled mode theory
    # def Calculate_Beta_and_save(self,plot=True,
    #                             filename_uncoupled = "data/beta_uncoupled.txt",
    #                             filename_coupled   = "data/beta_coupled.txt"):
    #     (beta_uncoupled_arr,
    #      beta_ave_uncoupled_arr,
    #      beta_coupled_arr) = self.Scan_wavl(gap_idx=0)
    #     # beta of eigenmodes of two separate WGs
    #     with open(filename_uncoupled,'w') as f:
    #         f.write("wavelength,beta_1_uncoupled,beta_2_uncoupled,beta_ave_uncoupled\n")
    #         for wavl_idx in range(len(self.wavl_arr)):
    #             beta_str =  str(self.wavl_arr[wavl_idx]) + "," + \
    #                         str(np.real(beta_uncoupled_arr[wavl_idx,0])) + "," + \
    #                         str(np.real(beta_uncoupled_arr[wavl_idx,1])) + "," + \
    #                         str(np.real(beta_ave_uncoupled_arr[wavl_idx])) + "\n"
    #             f.write(beta_str)

    #     # beta of supermodes calculated with CMT
    #     with open(filename_coupled,'w') as f:
    #         f.write("wavelength,beta_1_coupled,beta_2_coupled\n")
    #         for wavl_idx in range(len(self.wavl_arr)):
    #             beta_str =  str(self.wavl_arr[wavl_idx]) + "," + \
    #                         str(np.real(beta_coupled_arr[wavl_idx,0])) + "," + \
    #                         str(np.real(beta_coupled_arr[wavl_idx,1])) + "\n"
    #             f.write(beta_str)

    #     if plot == True:
    #         self.Plot_curve(self.wavl_arr,
    #                         np.c_[beta_uncoupled_arr,beta_coupled_arr],
    #                         ['beta_uncoupled_1','beta_uncoupled_2',
    #                         'beta_supermode_1(CMT)','beta_supermode_1(CMT)'],
    #                         'wavelength(um)','beta - beta_ave(rad\\rad)')

    # Load supermodes results from Lumerical
    def Load_results(self,plot=True):
        pass