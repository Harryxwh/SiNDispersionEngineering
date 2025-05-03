import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Functions import *
plt.ion()
c = 3*1e8
ng_1550     = 1.595643967049753
freq_1550   = c/(1550*1e-9)
g0          = 673
Max_M_idx   = 3
num_of_pts  = 10000
D_2o        = -283 *1e3 * 2*np.pi      # unit: Hz

def Control_angle_range_0_2_2pi(angle):
    mask_less = np.where(angle < 0)
    mask_more = np.where(angle > 2*np.pi)
    angle_copy = np.copy(angle)
    angle_copy[mask_less] = angle[mask_less] + (np.astype(-angle[mask_less] / (2*np.pi),int)+1) *2*np.pi
    angle_copy[mask_more] = angle[mask_more] - (np.astype(angle[mask_more] / (2*np.pi),int)) *2*np.pi
    return angle_copy

def Control_angle_range_mpi_2_pi(angle):
    mask_less = np.where(angle < - np.pi)
    mask_more = np.where(angle > np.pi)
    angle_copy = np.copy(angle)
    angle_copy[mask_less] = angle[mask_less] + (np.astype((-angle[mask_less]+np.pi)/ (2*np.pi),int)) *2*np.pi
    angle_copy[mask_more] = angle[mask_more] - (np.astype((angle[mask_more]+np.pi) / (2*np.pi),int)) *2*np.pi
    return angle_copy

def Interpolation(x,y,x_intp,num_of_pts=100):
    cs = CubicSpline(x, y, bc_type='natural')  # bc_type 可选 'natural', 'clamped', 'periodic' 等
    y_intp = cs(x_intp)
    return y_intp

def Eigenvalue_arg(freq,g,Lco,LAS,LAL,LBS,LBL,include_theory=False):
    LA      = 2*(LAS + LAL + Lco)     # unit: m
    LB      = 2*(LBS + LBL + Lco)     # unit: m
    beta    = ng_1550/c * freq * 2*np.pi         # unit: /m
    Sum     = np.exp(2* 1j * beta * LAS) + np.exp(2* 1j * beta * LBL)
    Diff    = np.exp(2* 1j * beta * LAS) - np.exp(2* 1j * beta * LBL)
    # print(beta)
    # print("Sum = ", np.angle(Sum)[:10])
    # print("Sum_theory = ", np.angle(np.exp(1j *2 * beta[:10] * LAS)))
    # print("Diff = ", Diff)
    # print(g)

    DL_uco  = (LBL-LAS)
    DL_tot  = (LB-LA)/2
    DL      = DL_uco - DL_tot
    F       = Sum * np.cos(2 * g * Lco) * np.cos(beta * DL) - 1j * Diff * np.sin(beta * DL)

    Ang_offset  = - np.angle(np.exp(1j * beta * (LAS+LBL)))
    Ang_p       = Control_angle_range_0_2_2pi(np.angle(0.5*(F + 1j * np.sqrt((Sum**2-Diff**2) - F**2))) + Ang_offset) - np.pi
    Ang_m       = Control_angle_range_0_2_2pi(np.angle(0.5*(F - 1j * np.sqrt((Sum**2-Diff**2) - F**2))) + Ang_offset) - np.pi

    flip_mask               = np.where(Ang_p < Ang_m)
    Ang_pre                 = np.c_[Ang_p, Ang_m]
    Ang_pre[flip_mask]      = np.flip(Ang_pre[flip_mask],axis=1)
    if include_theory:
        Ang_pre             = np.c_[Ang_pre,Ang_theory_p,Ang_theory_m]
    return Ang_pre

# Load kappa under different wavls of two coupled WGs width 2.8um
def Load_kappa_data():
    c = 3*1e8
    Kappa_arr = []
    with open ("./results/straight_WG_Kappa.txt",'r') as f:
        read_data = f.readlines()
        for line in read_data[1:]:
            line_strip = line.strip()
            wavl = float(line_strip.split(",")[0])
            freq = c/wavl * 1e9
            K_12 = str2complex(line_strip.split(",")[1])
            K_21 = str2complex(line_strip.split(",")[2])
            Kappa = (K_12+K_21)/2
            Kappa = np.real(Kappa)
            Kappa_arr.append([freq,Kappa])
        Kappa_arr = np.array(Kappa_arr)
    return Kappa_arr

# g : dict of coupling strength
def Reson_freq_2D(m,D_ave,g,L,epsilon):
    return D_ave/(2*np.pi)*np.arccos(np.cos(g*L)*np.cos(2*np.pi * epsilon *m))

def FSR_func_2D(m,D_ave,g,L,epsilon):
    return epsilon * D_ave/(2*np.pi)* np.cos(g*L)*np.sin(2*np.pi * epsilon *m) /        \
            (1-np.cos(g*L)**2 * np.cos(2*np.pi * epsilon *m)**2)**0.5

def Dispersion_2D(m,D_ave,g,L,epsilon):
    return D_ave * 2*np.pi* epsilon**2 *np.cos(g*L) *np.sin(g*L)**2 * np.cos(2*np.pi*epsilon*m) / \
            (1-np.cos(g*L)**2 * np.cos(2*np.pi * epsilon *m)**2)**1.5

def Dispersion_2D_num(m,reson_freq):
    FSR = First_derivative_central_diff(reson_freq,m)
    return First_derivative_central_diff(FSR,m[1:-1])

def Reson_freq_3D_numerical(m,D_ave,g,Lco,LAS,LAL,LBS,LBL):
    LA = 2*(Lco+LAS+LAL)
    LB = 2*(Lco+LBS+LBL)
    L_ave = (LA+LB)/2
    freq = m*c/(ng_1550*L_ave)
    # freq = freq_1550 + m * D_ave/(2*np.pi)
    Eigenvalue_Arg  = Eigenvalue_arg(freq,g,Lco,LAS,LAL,LBS,LBL)
    return D_ave/(2*np.pi) * Eigenvalue_Arg

def Reson_freq_3D(m_arr_intp,D_ave,g,Lco,LAS,LAL,LBS,LBL):
    LA      = 2*(LAS + LAL + Lco)     # unit: m
    LB      = 2*(LBS + LBL + Lco)     # unit: m
    L_ave   = (LA+LB)/2
    DL_uc1  = (LBL-LAS)
    DL_uc2  = (LBS-LAL)
    return D_ave/(2*np.pi)*np.arccos( np.cos(2*np.pi * DL_uc1 * m_arr_intp/ L_ave) *\
            np.cos(2*np.pi * DL_uc2 * m_arr_intp/ L_ave) * np.cos(2*g*Lco) -\
            np.sin(2*np.pi * DL_uc1 * m_arr_intp/ L_ave) * np.sin(2*np.pi * DL_uc2 * m_arr_intp/ L_ave) )

def FSR_func_3D(m,reson_freq):
    FSR = First_derivative_central_diff(reson_freq,m)
    return FSR

def Dispersion_3D(m,reson_freq):
    FSR = FSR_func_3D(m,reson_freq)
    return First_derivative_central_diff(FSR,m[1:-1])

def TOD_3D(m,reson_freq):
    GVD = Dispersion_3D(m,reson_freq)
    return First_derivative_central_diff(GVD,m[2:-2])

# returns the anomalous dispersion range around 1550nm of 2D parallel coupled rings (unit: nm)
def AD_range_func(m_arr,D_arr,M,FSR):
    zero_list,zero_idx_list = find_zero(m_arr, D_arr)
    zero_list               = np.array(zero_list)
    zero_list_m_0           = zero_list[np.where(np.abs(zero_list) < M)]      # zeros around mode numer m = 0
    zero_list_m_0_p         = zero_list_m_0[zero_list_m_0>0]
    zero_list_m_0_m         = zero_list_m_0[zero_list_m_0<0]
    if len(zero_list_m_0_p) > 0 and  len(zero_list_m_0_m) > 0 and D_arr[np.argmin(np.abs(m_arr))] > 0:
        return (np.min(zero_list_m_0_p) - np.max(zero_list_m_0_m))*FSR * 1e-9 /100 * 0.8
    else:
        return 0

def epsilon_func(L1,L2):
    return (L2-L1)/(L1+L2)

def Mode_nonconserved_coupling(Max_M_idx,M,
                               plot_coupled_curves=True,
                               coupled_data_arr=[],coupled_data_label_arr=[],
                               param_dict_={},
                               title = "Mode number non-conservation coupling",
                               num_of_pts=100,ylim=()):
    m_arr = np.linspace(-Max_M_idx*M, Max_M_idx*M, num_of_pts).reshape(-1,1)
    M0 = 0
    Y_legends = []
    Y_data = np.array(D1 * (m_arr - M0))
    Y_legends.append(r"Resonator1 $\omega = \omega_0$+$(D_1-D_{ave})$(m-$M_0$)")
    for m in range(1,Max_M_idx):
        Y = D1 * (m_arr - M0) + m*FSR
        Y_data = np.c_[Y_data,Y]
        Y_legends.append(r"Resonator1 $\omega$ = $\omega_0$+$(D_1-D_{ave})$(m-$M_0$)$\pm$"+str(m)+r"$D_{ave}$")
        Y = D1 * (m_arr - M0) - m*FSR
        Y_data = np.c_[Y_data,Y]
        Y_legends.append("")

    Y = D2 * (m_arr - M0)
    Y_data = np.c_[Y_data,Y]
    Y_legends.append(r"Resonator2 $\omega = \omega_0$+$(D_2-D_{ave})$(m-$M_0$)")
    for m in range(1,Max_M_idx):
        Y = D2 * (m_arr - M0) +  m*FSR
        Y_data = np.c_[Y_data,Y]
        Y_legends.append(r"Resonator2 $\omega$ = $\omega_0$+$(D_2-D_{ave})$(m-$M_0$)$\pm$"+str(m)+r"$D_{ave}$")
        Y = D2 * (m_arr - M0) -  m*FSR
        Y_data = np.c_[Y_data,Y]
        Y_legends.append("")
    Y_data = Y_data - D_ave * (m_arr - M0)

    if plot_coupled_curves:
        Y_data = np.c_[Y_data,coupled_data_arr]
        Y_legends = Y_legends + coupled_data_label_arr

    linestyle_name_list = ["dashed","dotted"]*3
    linestyle_list      = ["-"]
    for i in range(Max_M_idx-1):
        for j in range(2):
            linestyle_list.append(linestyle_name_list[i])
    linestyle_list = linestyle_list + linestyle_list
    linestyle_list = linestyle_list + ["-"]*30
    colors_list = ['lightskyblue']*(Max_M_idx*2-1)+ ['lightcoral']*(Max_M_idx*2-1)\
                    +['orange']*2+['dodgerblue']*2+['black']*30

    xticks       = np.arange(-Max_M_idx,Max_M_idx+1)
    xtickslabels = [("$M_0$+" if xtick>0 else "$M_0$") + str(xtick) + "M"  for xtick in xticks]
    xtickslabels[int(len(xtickslabels)/2)] = "$M_0$"
    xticks       = np.arange(-Max_M_idx,Max_M_idx+1)*M

    param_dict = {  "Y_legends"       : Y_legends,
                    "X_label"         : 'mode number m',
                    "Y_label"         : r"Frequency $\omega$/(2$\pi$) (GHz)",
                    "xticks"          : xticks,
                    "xtickslabel"     : xtickslabels,
                    "title"           : title,
                    "marker_list"     : [""]*30,
                    "linestyle_list"  : linestyle_list,
                    "colors_list"     : colors_list,
                    "xlim"            : (-Max_M_idx*M,Max_M_idx*M),
                    "ylim"            : ylim,
                    "bbox_legend"     : (1.05,0.6)}
    if len(param_dict_)>0:
        for key,value in param_dict_.items():
            param_dict[key] = value

    data_arr = (np.c_[m_arr,Y_data/1e9/(2*np.pi)],)

    Plot_curve(data_arr,**param_dict)

def Optimized_2D_parallel(m_arr_intp,L1,Rtot,gL_arr,g0,D_iso):
    # The optimal AD range if the same rings are put in 2D parallel structure
    epsilon_arr     = epsilon_func(L1, L1*Rtot_arr)
    Lco_arr         = gL_arr  / g0
    L2      = L1 * Rtot
    D1      = c/(ng_1550 * L1) *2* np.pi
    D2      = c/(ng_1550 * L2) *2* np.pi
    D_ave   = c/(ng_1550 * (L1+L2)/2) *2 *np.pi
    # D_ave = D1 * 0.99
    epsilon = (L2-L1)/(L1+L2)
    FSR = (D1-D2)/(2*epsilon)
    M = 1/(2*epsilon)

    AD_range_2D_arr = []
    for Lco in Lco_arr:
        Reson_freq_arr  = Reson_freq_2D(m_arr_intp,D_ave,g_arr_intp,Lco,epsilon)
        # D_coupled_arr  = Dispersion_2D_num(m_arr_intp,Reson_freq_arr)
        D_coupled_arr   = Dispersion_2D(m_arr_intp,D_ave,g_arr_intp,Lco,epsilon)
        AD_range_2D_res = AD_range_func(m_arr_intp[2:-2], D_coupled_arr + D_iso,
                                        M = M, FSR = D_ave/(2*np.pi))
        AD_range_2D_arr.append(AD_range_2D_res)
    AD_range_2D_arr    = np.array(AD_range_2D_arr)

    flat_index = np.argmax(AD_range_2D_arr)
    max_index_2D = np.unravel_index(flat_index, np.shape(AD_range_2D_arr))
    best_AD_range_2D = AD_range_2D_arr[max_index_2D]
    best_gL_2D = gL_arr[max_index_2D[0]]
    return best_AD_range_2D, best_gL_2D

def Optimize_3D_offset(m_arr_intp,L_str_A,L_str_B,L_bend,gL_arr,deltaLs_arr,g0,D_iso):
    # The optimal AD range if the rings are put in 3D offset structure
    AD_range_arr = []
    L1 = (L_str_A + L_bend)*2
    L2 = (L_str_B + L_bend)*2
    D1      = c/(ng_1550 * L1) *2* np.pi
    D2      = c/(ng_1550 * L2) *2* np.pi
    D_ave   = c/(ng_1550 * (L1+L2)/2) *2 *np.pi
    epsilon = (L2-L1)/(L1+L2)
    FSR     = (D1-D2)/(2*epsilon)
    M       = 1/(2*epsilon)

    for deltaLs in deltaLs_arr:
        deltaL_str  = L_str_B-L_str_A
        deltaLL     = deltaL_str - deltaLs
        LBS         = L_bend + max(deltaLs,0)
        LBL         = L_bend + deltaLL
        LAS         = L_bend
        LAL         = L_bend - min(deltaLs,0)
        Lco         = L_str_A + min(deltaLs,0)

        AD_range_arr_given_deltaL = []
        for gL in gL_arr:
            g_arr_eff   = g_arr_intp * gL /(g0*Lco)

            Y_p_2D      = Reson_freq_2D(m_arr_intp,D_ave,g_arr_eff,2*Lco,epsilon)
            Y_m_2D      = -Reson_freq_2D(m_arr_intp,D_ave,g_arr_eff,2*Lco,epsilon)
            Y_p_3D      = Reson_freq_3D(m_arr_intp,D_ave,g_arr_eff,Lco,LAS,LAL,LBS,LBL)
            Y_m_3D      = -Reson_freq_3D(m_arr_intp,D_ave,g_arr_eff,Lco,LAS,LAL,LBS,LBL)

            D_iso       = D_2o * np.ones((np.shape(m_arr_intp)[0]-4,))
            D_2D        = Dispersion_2D(m_arr_intp,D_ave,g_arr_eff,2*Lco,epsilon)[2:-2]
            D_3D        = Dispersion_3D(m_arr_intp,Y_m_3D)
            AD_range_res = AD_range_func(m_arr_intp[2:-2], D_3D+D_iso,
                                        M = M, FSR = D_ave/(2*np.pi))
            AD_range_arr_given_deltaL.append(AD_range_res)

        AD_range_arr.append(AD_range_arr_given_deltaL)

    AD_range_arr    = np.array(AD_range_arr)
    flat_index      = np.argmax(AD_range_arr)
    max_index       = np.unravel_index(flat_index, np.shape(AD_range_arr))
    best_AD_range   = AD_range_arr[max_index]
    assert len(max_index)==2
    best_gL         = gL_arr[max_index[1]]
    best_deltaLs    = deltaLs_arr[max_index[0]]

    # Plot image
    xticks = np.arange(0,len(gL_arr),2)
    yticks = np.arange(0,len(deltaLs_arr),2)
    param_dict = {
            "aspect"        : 0.5,
            "figsize"       : (40,4.5),
            "xlabel"        : r"$g_{co} (m^{-1})$",
            "ylabel"        : r"$\delta L_s (\mu m)$",
            "cbar_label"    : "Anomalous Dispersioin range (nm)",
            "cbar_small_ticks" : True,
            "title"         : "Optimizing Anomalous Dispersion range when "+
                                r"$R_{tot}$="+"{:.4f}".format(Rtot),
            "xticks"        : xticks,
            "yticks"        : yticks,
            "xtickslabel"   : ["{:.0f}".format(gL/Lco) for gL in gL_arr[xticks]],
            "ytickslabel"   : ["{:.2f}".format(dLs*1e6) for dLs in deltaLs_arr[yticks]],
            "fontsize"      : 6,
    }
    Plot_im(AD_range_arr,**param_dict)

    return best_gL, best_deltaLs, best_AD_range


if __name__ == '__main__':

    Kappa_arr       = Load_kappa_data()
    Kappa_arr       = np.flip(Kappa_arr,axis=0)
    freq_arr        = Kappa_arr[:,0]

    Rtot_arr        = np.linspace(1.0010,1.0100,10)
    # Rtot_arr        = np.linspace(1.0041,1.0041,1)

    gL_arr          = np.linspace(0.6,1.6,11)

    # Set the geometry of the two resonators
    # unit: m
    for Rtot in Rtot_arr:

        L_bend      = 3.14  *1e-3
        L_str_A     = 1.61  *1e-3
        L_str_B     = (L_bend+L_str_A)*Rtot - L_bend
        deltaL_str  = L_str_B-L_str_A
        L1          = (L_str_A+L_bend) *2
        L2          = L1 * Rtot

        # -L_strA < deltaLs < 0.5 * DeltaL_str
        deltaLs_arr     = np.linspace(-L_str_A*0.03,0.5*deltaL_str,20)   #unit: m

        D1      = c/(ng_1550 * L1) *2* np.pi
        D2      = c/(ng_1550 * L2) *2* np.pi
        D_ave   = c/(ng_1550 * (L1+L2)/2) *2 *np.pi
        epsilon = (L2-L1)/(L1+L2)
        FSR     = (D1-D2)/(2*epsilon)
        M       = 1/(2*epsilon)
        M_tick  = 0.5       # smallest tick of M when ploting
        m_arr_intp      = np.linspace(-Max_M_idx*M, Max_M_idx*M, num_of_pts)
        freq_arr_intp   = freq_1550 + m_arr_intp * D_ave/(2*np.pi)
        g_arr_intp      = Interpolation(freq_arr,Kappa_arr[:,1],freq_arr_intp)
        D_iso           = D_2o * np.ones((np.shape(m_arr_intp)[0]-4,))

        #################################################################
        # Optimized parameters for 3D offset structures
        #################################################################
        best_gL_3D, best_deltaLs_3D, best_AD_range_3D = Optimize_3D_offset(m_arr_intp,L_str_A,L_str_B,
                                                                        L_bend,gL_arr,deltaLs_arr,g0,D_iso)

        deltaLL     = deltaL_str - best_deltaLs_3D
        LBS         = L_bend + max(best_deltaLs_3D,0)
        LBL         = L_bend + deltaLL
        LAS         = L_bend
        LAL         = L_bend - min(best_deltaLs_3D,0)
        Lco         = L_str_A + min(best_deltaLs_3D,0)
        g_arr_eff   = g_arr_intp * best_gL_3D /(g0*Lco)

        Y_p_3D      = Reson_freq_3D(m_arr_intp,D_ave,g_arr_eff,Lco,LAS,LAL,LBS,LBL)
        Y_m_3D      = -Reson_freq_3D(m_arr_intp,D_ave,g_arr_eff,Lco,LAS,LAL,LBS,LBL)

        #################################################################
        # Optimized parameters for 2D parallel structures
        #################################################################
        best_AD_range_2D, best_gL_2D = Optimized_2D_parallel(m_arr_intp,L1,Rtot,gL_arr,g0,D_iso)

        # Compare Reson freq of 2D and 3D structures
        Y_p_2D      = Reson_freq_2D(m_arr_intp,D_ave,g_arr_intp,best_gL_2D/g0,epsilon)
        Y_m_2D      = -Reson_freq_2D(m_arr_intp,D_ave,g_arr_intp,best_gL_2D/g0,epsilon)

        data_arr = np.c_[Y_p_2D,Y_m_2D,Y_p_3D,Y_m_3D]
        data_label_arr = ["","2D parallel structure",
                        "","3D offset structure"]*3
        text = r"$R_{tot}$ = $L_B/L_A$"+ " = {:.4f}".format(Rtot)+"\n"+\
                r"$D_{ave}/(2\pi$) = "+"{:.2f} GHz".format(D_ave/(2*np.pi)*1e-9)+"\n" +\
                "M = "+"{:.1f} ".format(M)+"\n"\
                "3D offset structure:\n"+ \
                r"$g_{co}$"+"= {:.2f} ".format(best_gL_3D/Lco)+r"$m^{-1}$" +"\n" +\
                r"$\delta L_s$"+ " = {:.4f}".format(best_deltaLs_3D*1e6)+r" $\mu m$"+"\n"\
                "2D parallel structure:\n"+ \
                r"$g_{co}L_{co}$"+"= {:.2f}".format(best_gL_2D)

        param_dict  = {"alpha_list": [1]*14+[0.6,0.6]+[1]*4,
                        "text":text,
                        "loc_text":(1.1,0.2),
                        "linespacing":1.8,
                        "title": "Resonant frequency of 3D coupled rings when "+
                            r"$R_{tot}$="+"{:.4f}".format(Rtot)}
        Mode_nonconserved_coupling(Max_M_idx,M,plot_coupled_curves=True,
                                coupled_data_arr=data_arr,
                                coupled_data_label_arr=data_label_arr,
                                param_dict_ = param_dict,
                                num_of_pts=num_of_pts,ylim=(-20,20))

        # Compare Dispersion of 2D and 3D structures
        D_2D        = Dispersion_2D_num(m_arr_intp,Y_p_2D)
        # D_2D        = Dispersion_2D(m_arr_intp,D_ave,g_arr_intp,best_gL_2D/g0,epsilon)[2:-2]
        D_3D        = Dispersion_3D(m_arr_intp,Y_m_3D)
        data_arr    = np.flip(np.c_[D_2D+D_iso, D_3D+D_iso],axis=0)
        data_arr    = (np.c_[m_arr_intp[2:-2], data_arr/1e3/(2*np.pi)],)

        data_label_arr = [r"$D_2$ 2D parallel rings (Optimized)",r"$D_2$ 3D offset rings (Optimized)"]
        linestyle_list = ["-"]*10
        xticks       = np.arange(-Max_M_idx, Max_M_idx+M_tick, M_tick)*M
        xtickslabels = np.array(["{:.0f}".format(3*1e8 / (freq_1550 + D_ave/(2*np.pi)
                                                          * xtick) * 1e9) for xtick in xticks])
        yticks       = ticks_arr(data_arr)

        text  = "3D offset structure:\n"+ \
                r"$g_{co}$"+"= {:.2f} ".format(best_gL_3D/Lco)+r"$m^{-1}$" +"\n" +\
                r"$\delta L_s$"+ " = {:.4f}".format(best_deltaLs_3D*1e6)+r" $\mu m$"+"\n"\
                "Optimal AD range: {:.2f} nm".format(best_AD_range_3D)+"\n"+\
                "2D parallel structure:\n"+ \
                r"$g_{co}L_{co}$"+"= {:.2f}".format(best_gL_2D) +"\n" +\
                "Optimal AD range: {:.2f} nm".format(best_AD_range_2D)

        param_dict = {
            "Y_legends"     : data_label_arr,
            "X_label"       : 'wavelength (nm)',
            "Y_label"       : r"$D_2$/(2$\pi$) (kHz)",
            "title"         : "Optimized Dispersion of 2D and 3D structures when "+
                                    r"$R_{tot}$="+"{:.4f}".format(Rtot),
            "figsize"       : (10,6),
            "marker_list"   : [""]*15,
            "linestyle_list": linestyle_list,
            "colors_list"   : ["lightcoral"]+['Orange']+['lightskyblue']+['deepskyblue']+['black']*10,
            "xticks"        : xticks,
            "xtickslabel"   : np.flip(xtickslabels),
            "yticks"        : yticks,
            "alpha_list"    :[0.6,1,1],
            "plot_linewidth": [2,2],
            # "xlim"          : (-3*M,3*M+1),
            "xlim"          : (-M,M+1),
            "ylim"          : (-2000,2500),
            "AD_region_color"    : True,
            # "bbox_legend"   : (0.9,0.9),
            "text"          : text,
            # "linespacing"   : 1.5,
            "loc_text"      : (0.1,0.6),
            "linespacing"   : 1.8
        }
        Plot_curve(data_arr,**param_dict)
