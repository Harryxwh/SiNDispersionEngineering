from Functions import *
from Coupled_Waveguides import *
class Data_analyzer(Coupled_Waveguides):

    def __init__(self, wavl_arr, gap_arr,
                 filename_uncoupled, filename_coupled,
                 param_filename, plot_curve=True):
        self.load_param(param_filename)
        self.wavl_arr = wavl_arr
        self.gap_arr = gap_arr
        self.plot_curve = plot_curve

        beta_uncoupled_arr = self.Load_uncoupled_data(filename_uncoupled)
        beta_coupled_arr = self.Load_coupled_data_CMT(filename_coupled)
        # beta_coupled_lumerical_arr, beta_ave_arr = self.Load_coupled_data_Lumerical(filename_lumerical)

        gap_x, gap_y = self.gap_arr
        if gap_x > 0:
            self.bend_radius_outer = self.bend_radius_inner +\
                    gap_x + (self.WG1_width + self.WG2_width)/2       # unit :um
        else:
            self.bend_radius_outer = self.bend_radius_inner

        self.Calc_Dispersion_curve(beta_uncoupled_arr,beta_coupled_arr,)
                                #    beta_coupled_lumerical_arr,beta_ave_arr)

    def Load_uncoupled_data(self,filename_uncoupled):
        beta_uncoupled_arr = []
        with open(filename_uncoupled,'r') as f:
            data_uncoupled = f.readlines()
            for line in data_uncoupled[1:]:
                beta_uncoupled_arr.append(np.float64(line.split(",")))
        beta_uncoupled_arr = np.array(beta_uncoupled_arr)
        return beta_uncoupled_arr

    def Load_coupled_data_CMT(self, filename_coupled, wavl_idx = 10, print_coeffi = False):
        beta_coupled_arr = []
        coeff_supermode_1_arr = []
        coeff_supermode_2_arr = []
        with open(filename_coupled,'r') as f:
            data_coupled = f.readlines()
            for line in data_coupled[1:]:
                line = line.split(",")
                line_float = np.float64(line[:3])
                if line_float[1]<0:
                    beta_coupled = [line_float[0],line_float[2],line_float[1]]
                else:
                    beta_coupled = line_float[:3]
                assert len(beta_coupled) == 3
                assert len(line) == 7
                coeff_supermode_1 = [str2complex(line[3]),str2complex(line[4])]
                coeff_supermode_2 = [str2complex(line[5]),str2complex(line[6])]
                beta_coupled_arr.append(beta_coupled)
                coeff_supermode_1_arr.append(coeff_supermode_1)
                coeff_supermode_2_arr.append(coeff_supermode_2)
        beta_coupled_arr = np.array(beta_coupled_arr, dtype = np.float64)
        coeff_supermode_1_arr = np.array(coeff_supermode_1_arr, dtype = np.complex64)
        coeff_supermode_2_arr = np.array(coeff_supermode_2_arr, dtype = np.complex64)

        # coefficient of supermode at 1550nm
        if print_coeffi:
            print("wavelength = {:.2f}".format(beta_coupled_arr[wavl_idx,0]) + " um "\
                + "\nSupermode 1 coefficient : "  \
                + 'A = ({0.real:.6f} + {0.imag:.6f}i)'.format(coeff_supermode_1_arr[wavl_idx,0]) + ', '\
                + 'B = ({0.real:.6f} + {0.imag:.6f}i)'.format(coeff_supermode_1_arr[wavl_idx,1])\
                + "\nSupermode 2 coefficient : "  \
                + 'A = ({0.real:.6f} + {0.imag:.6f}i)'.format(coeff_supermode_2_arr[wavl_idx,0]) + ', '\
                + 'B = ({0.real:.6f} + {0.imag:.6f}i)'.format(coeff_supermode_2_arr[wavl_idx,1]))

        return beta_coupled_arr

    def Load_coupled_data_Lumerical(filename_lumerical):
        beta_coupled_lumerical_arr_ori = []
        with open(filename_lumerical,'r') as f:
            data_lumerical = f.readlines()
            beta_ang_mode1 = 0
            for line in data_lumerical[2:]:
                wavelength  = float(line.split(',')[0])/1000 #unit:um
                modeidx     = int(line.split(',')[1])
                neff        = str2complex(line.split(',')[2])
                ng          = str2complex(line.split(',')[3])
                loss        = float(line.split(',')[4])
                polarization= float(line.split(',')[5])
                beta_ang    = float(line.split(',')[6])

                if modeidx == 1:
                    # Mode 1
                    beta_ang_mode1 = beta_ang
                else:
                    # Mode 2
                    beta_ang_mode2 = beta_ang
                    beta_coupled_lumerical_arr_ori.append([wavelength,beta_ang_mode1,beta_ang_mode2])

        beta_coupled_lumerical_arr_ori = np.array(beta_coupled_lumerical_arr_ori)
        # beta_ave = beta_uncoupled_arr[:,3]
        beta_ave_arr = (beta_coupled_lumerical_arr_ori[:,1] + beta_coupled_lumerical_arr_ori[:,2])/2
        beta_coupled_lumerical_arr = np.copy(beta_coupled_lumerical_arr_ori)
        beta_coupled_lumerical_arr[:,1] = beta_coupled_lumerical_arr_ori[:,1] - beta_ave_arr
        beta_coupled_lumerical_arr[:,2] = beta_coupled_lumerical_arr_ori[:,2] - beta_ave_arr
        return beta_coupled_lumerical_arr, beta_ave_arr

    def Plot_curve(self,X,Y_arr,Y_legends,
                    X_label,Y_label,
                    title,marker_list,linestyle_list,
                    colors_list=['green','mediumblue','tomato','orange']*2,
                    bbox_to_anchor=(),
                    text = "",
                    dpi=400):
        #Plot parameters
        figsize = (8,6)
        fonttype = "Helvetica"
        fontsize = 10
        grid_linewidth = 0.8
        plot_linewidth = 1.5

        plt.figure(figsize=figsize)
        for idx in range(np.shape(Y_arr)[1]):
            plt.plot(X,Y_arr[:,idx],label=Y_legends[idx],
                        color=colors_list[idx], marker=marker_list[idx],
                        linestyle=linestyle_list[idx], linewidth=plot_linewidth)

        plt.rcParams["font.family"] = fonttype
        plt.rcParams.update({'font.size': fontsize})
        plt.yticks(fontproperties = fonttype, size = fontsize)
        plt.xticks(fontproperties = fonttype, size = fontsize)
        plt.ylabel(Y_label, fontdict={'family' : fonttype, 'size' : fontsize})
        plt.xlabel(X_label, fontdict={'family' : fonttype, 'size' : fontsize})
        plt.title(title)
        if len(bbox_to_anchor)>0:
            plt.legend(bbox_to_anchor=bbox_to_anchor, loc='upper left', borderaxespad=0)
        else:
            plt.legend(loc='upper left')
        if not text == "":
            plt.text(np.quantile(X,0.75),np.quantile(Y_arr,0.05),
                    text, bbox=dict(boxstyle="round,pad=0.9", fc="white", alpha=0.9))
        plt.grid(linewidth=grid_linewidth, alpha=0.3)
        savename = "results/"+str(title)+".jpg"
        plt.savefig(savename,dpi=dpi)
        plt.tight_layout()
        plt.show()

    def Interpolation(x,y,x_intp,num_of_pts=100):
        cs = CubicSpline(x, y, bc_type='natural')  # bc_type 可选 'natural', 'clamped', 'periodic' 等
        y_intp = cs(x_intp)
        return y_intp

    def Calculate_dispersion_D(self,Beta):
        fre_arr = 3*10**8 / self.wavl_arr
        Beta_1  = First_derivative_central_diff(Beta, fre_arr)
        Beta_1  = Beta_1 *10**(12) * 10**(3)        # unit: ps/km
        D       = First_derivative_central_diff(Beta_1, self.wavl_arr[1:-1])
        D       = D * 10**(-9)                          # unit: ps/km/nm
        return (D, Beta_1)

    def Anomalous_D_bandwidth(self,wavl_arr,beta_supermode):
        start_idx = int(len(beta_supermode)/2)
        left_zero = wavl_arr[start_idx]
        right_zero = wavl_arr[start_idx]
        left_found = False
        right_found = False
        for delta_idx in range(int(len(beta_supermode)/2)):
            if right_found:
                continue
            idx = start_idx + delta_idx
            if beta_supermode[idx] > 0 and beta_supermode[idx+1] < 0:
                right_zero =  wavl_arr[idx] + (wavl_arr[idx+1]-wavl_arr[idx])*\
                        beta_supermode[idx]/(beta_supermode[idx]-beta_supermode[idx+1])
                right_found = True
            elif  beta_supermode[idx+1] > 0:
                right_zero = wavl_arr[idx+1]
            # print("idx:" + str(idx) + ", beta =" + str(beta_supermode[idx]))
            if left_found:
                continue
            idx = start_idx - delta_idx
            if beta_supermode[idx] > 0 and beta_supermode[idx-1] < 0:
                left_zero =  wavl_arr[idx] + (wavl_arr[idx-1]-wavl_arr[idx])*\
                        beta_supermode[idx]/(beta_supermode[idx]-beta_supermode[idx-1])
                left_found = True
            elif  beta_supermode[idx-1] > 0:
                left_zero = wavl_arr[idx-1]
            # print("idx:" + str(idx) + ", beta =" + str(beta_supermode[idx]))
        return right_zero-left_zero

    def Calc_Dispersion_curve(self,beta_uncoupled_arr,beta_coupled_arr,):
                        #beta_coupled_lumerical_arr,beta_ave_arr):


        # unit: rad/m
        # beta_uncoupled_WG1 = (beta_uncoupled_arr[:,1] + beta_uncoupled_arr[:,3])  / self.bend_radius_inner
        # beta_uncoupled_WG2 = (beta_uncoupled_arr[:,2] + beta_uncoupled_arr[:,3])  / self.bend_radius_outer

        # D_WG1, Beta_1_WG1 = self.Calculate_dispersion_D(beta_uncoupled_WG1)
        # D_WG2, Beta_1_WG2 = self.Calculate_dispersion_D(beta_uncoupled_WG2)

        D_ave, Beta_1_ave = self.Calculate_dispersion_D(beta_uncoupled_arr[:,3]/self.Bend_radius_ave)
        beta_CMT_supermode1 = (beta_coupled_arr[:,1])  / self.Bend_radius_ave     # unit: rad/m
        beta_CMT_supermode2 = (beta_coupled_arr[:,2])  / self.Bend_radius_ave     # unit: rad/m

        D_supermode_1, Beta_1_supermode_1 = self.Calculate_dispersion_D(beta_CMT_supermode1)
        D_supermode_2, Beta_1_supermode_2 = self.Calculate_dispersion_D(beta_CMT_supermode2)
        # beta_lumerical_supermode1 = (beta_coupled_lumerical_arr[:,1] +
        #                             beta_uncoupled_arr[:,3])  / bend_radius_ave    # unit: rad/m

        # beta_coupled_lumerical_arr_ori  = np.c_[beta_coupled_lumerical_arr[:,0],
        #                                         beta_coupled_lumerical_arr[:,1] + beta_ave_arr,
        #                                         beta_coupled_lumerical_arr[:,2] + beta_ave_arr]
        # beta_coupled_lumerical_supermode1_ori  = beta_coupled_lumerical_arr_ori[:,1]  / bend_radius_ave
        # beta_coupled_lumerical_supermode2_ori  = beta_coupled_lumerical_arr_ori[:,2]  / bend_radius_ave

        # (D_lumerical_supermode_1, Beta_1_lumerical_supermode_1) = self.Calculate_dispersion_D(beta_coupled_lumerical_arr[:,1]/bend_radius_ave)
        # (D_lumerical_supermode_2, Beta_1_lumerical_supermode_2) = self.Calculate_dispersion_D(beta_coupled_lumerical_arr[:,2]/bend_radius_ave)

        # Dispersion_data = np.c_[D_WG1, D_WG2, D_supermode_1, D_supermode_2,
        #                         D_lumerical_supermode_1, D_lumerical_supermode_2]
        Dispersion_data = np.c_[D_ave, D_supermode_1, D_supermode_2,]
                                # D_lumerical_supermode_1, D_lumerical_supermode_2]
        Dispersion_data [:,[1,2]] += D_ave.reshape(-1,1)

        AD_range_text = "AD range = "+"{:.4f}".format(self.Anomalous_D_bandwidth(
                        self.wavl_arr, D_supermode_2 + D_ave)) +" nm"
        gap_info = "gap_"+"{:.1f}".format(self.gap_arr[0])+\
                    ","+"{:.1f}".format(self.gap_arr[1])
        gap_info = gap_info.replace(".","_")
        if self.plot_curve:
            self.Plot_curve(self.wavl_arr[2:-2] *10**(6) , Dispersion_data,
                    Y_legends=[r'$D_0$','Supermode 1 (CMT)','Supermode 2 (CMT)',
                                'Supermode 1 (FDE)','Supermode 2 (FDE)'],
                    X_label='wavelength(um)',Y_label=r'$D(ps/nm/km)$',
                    title = "Dispersion of coupled modes_"+gap_info,
                    marker_list=["","v","v","o","o","v","v"],
                    linestyle_list=["--","-","-","-","-","-","-","-"],
                    colors_list = ['green','tomato','orange',
                                  'deepskyblue','lightskyblue','black','red']*2,
                    text=AD_range_text)
        return Dispersion_data

