import numpy as np
from Coupled_Waveguides import *
# sweepable parameters include: wavelength, gap
class Parameter_sweeper():
    def __init__(self,wavelengths,gap,foldername1,
                 foldername2,param_filename="./Param.csv"):
        #The material dispersion of waveguide and cladding
        index_file1 = "Si3N4_index.txt"
        index_file2 = "SiO2_index.txt"
        self.n1_dict, self.n0_dict = self.Load_index_data(index_file1,index_file2)

        self.wavl_vec = wavelengths     # Wavelength in Vacuum (unit:um)
        self.gap = gap                  # gap = [gap_x,gap_y] shape=(2,x) (unit:um)
        self.foldername_1 = foldername1
        self.foldername_2 = foldername2
        self.param_filename = param_filename
        self.beta_array = np.array([[],[]])

    def Load_index_data(self,filename_n1,filename_n0):
        n1_dict = {}
        n0_dict = {}
        with open(filename_n1,'r') as f_n1:
            data = f_n1.readlines()
            for line in data:
                line = line.replace('\n','')
                wavl = int(line.split(' ')[0].strip())
                n_re = float(line.split(' ')[1].strip())
                n_im = float(line.split(' ')[2].strip())
                #n1_dict.append(complex(n_re,n_im))
                n1_dict[wavl] = n_re
        with open(filename_n0,'r') as f_n0:
            data = f_n0.readlines()
            for line in data:
                wavl = int(line.split()[0])
                n_re = float(line.split()[1])
                n_im = float(line.split()[2])
                #n0_dict.append(complex(n_re,n_im))
                n0_dict[wavl] = n_re
        return n1_dict, n0_dict

    def Create_instance(self,wavl_idx,gap_idx,param_filename,
                        Plot_field = False,Plot_index = False):
        assert wavl_idx < len(self.wavl_vec)
        assert gap_idx < len(self.gap)
        wavl = self.wavl_vec[wavl_idx]
        wavl_in_nm = int(wavl*1000)
        assert wavl_in_nm in self.n0_dict
        assert wavl_in_nm in self.n1_dict
        n0 = self.n0_dict[wavl_in_nm]
        n1 = self.n1_dict[wavl_in_nm]
        gap_x, gap_y = self.gap[gap_idx]
        path_1 = self.foldername_1 +'/'+ str(wavl_in_nm)
        path_2 = self.foldername_2 +'/'+ str(wavl_in_nm)
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
        beta_uncoupled_arr     = [] # beta of eigenmodes of two separate WGs
        beta_ave_uncoupled_arr = [] # beta_ave of eigenmodes of two separate WGs
        beta_coupled_arr       = [] # beta of supermodes calculated with CMT
        beta_list = np.array([[]])
        for idx in range(len(self.wavl_vec)):
            Coupled_WG = self.Create_instance(wavl_idx=idx, gap_idx=gap_idx,
                                              param_filename=self.param_filename,
                                              Plot_field= True if idx==0 else False,
                                              Plot_index= True if idx==0 else False)
            #shape:(1,2)
            beta_uncoupled = np.array([[Coupled_WG.beta_ang_1,
                                        Coupled_WG.beta_ang_2]]) - Coupled_WG.beta_ave
            beta_coupled = Coupled_WG.Find_supermodes()
            if idx == 0:
                beta_uncoupled_arr = np.array(beta_uncoupled)
                beta_coupled_arr = np.array(beta_coupled)
            else:
                beta_uncoupled_arr = np.append(beta_uncoupled_arr,
                                               beta_uncoupled,axis=0)
                beta_ave_uncoupled_arr = np.append(beta_ave_uncoupled_arr,
                                                np.array([Coupled_WG.beta_ave]),axis=0)
                if beta_coupled[0,0]<0:
                    beta_coupled = np.flip(beta_coupled)
                beta_coupled_arr = np.append(beta_coupled_arr,
                                             beta_coupled ,axis=0)

        return beta_uncoupled_arr, beta_ave_uncoupled_arr, beta_coupled_arr

    def Scan_gap(self):
        pass

    def Dispersion(self):
        fre_vec = self.c / self.wavl_vec
        Diff_fre = np.diff(fre_vec,axis=0)
        Diff_beta = np.diff(self.beta_array,axis=1)
        FOD = np.multiply(Diff_beta, 1/(2 * np.pi) * np.reciprocal(Diff_fre))
        self.Plot_curve(self.wavl_vec, FOD, "FOD")


    # Save supermodes results calculated using coupled mode theory
    def Save_results(self,plot=True,filename_uncoupled = "results/beta_uncoupled.txt",
                                    filename_coupled   = "results/beta_coupled.txt"):
        (beta_uncoupled_arr,
         beta_ave_uncoupled_arr,
         beta_coupled_arr) = self.Scan_wavl(gap_idx=0)
        # beta of eigenmodes of two separate WGs
        with open(filename_uncoupled,'w') as f:
            f.write("beta_1_uncoupled,beta_2_uncoupled,beta_ave_uncoupled\n")
            for wavl_idx in range(len(self.wavl_vec)):
                beta_str =  str(np.real(beta_uncoupled_arr[wavl_idx,0])) + "," + \
                            str(np.real(beta_uncoupled_arr[wavl_idx,1])) + "," + \
                            str(np.real(beta_ave_uncoupled_arr)) + "\n"
                f.write(beta_str)

        # beta of supermodes calculated with CMT
        with open(filename_coupled,'w') as f:
            f.write("beta_1_coupled,beta_2_coupled\n")
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

    # Load supermodes results from Lumerical
    def Load_results(self,plot=True):

        pass