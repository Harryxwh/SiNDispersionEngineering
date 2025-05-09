'''
Main program for the simulation of two coupled microrings or bent waveguides
Author : Weihao Xu
Some comments
'''
from Coupled_Waveguides import *
from Parameter_sweeper import *
from Data_analyzer import *

if __name__ == '__main__':

    ########################################
    # Plot the super mode
    ########################################
    foldername1 = "../data/Mode Profiles/"+\
                    "InnerRing_Lx_8um_bendr_1000um_20x10um_800x400cells_beta_ang_21_wavls"+"/1550"
    # foldername1 = "../data/Mode Profiles/"+\
    #                 "OuterRingDesigned_20x10um_800x400cells_21_wavls"+"/1550"
    foldername2 = "../data/Mode Profiles/"+\
                    "OuterRingDesigned_20x10um_800x400cells_21_wavls"+"/1550"

    # foldername1 = "../data/Mode Profiles/"+\
    #                 "InnerRing_L_inner_8um_1400_1700_31wavls"+"/1550"
    # foldername2 = "../data/Mode Profiles/"+\
    #                 "InnerRing_L_inner_8um_1400_1700_31wavls"+"/1550"

    foldername1 = "../data/Mode Profiles/"+\
                    "InnerRing_L_inner_8um_1535_1565_31wavls"+"/1550"
    # foldername2 = "../data/Mode Profiles/"+\
    #                 "OuterRingDesigned_L_inner_8um_gap_5um_L_outer_2436nm_1535_1565nm_31wavls"+"/1550"
    foldername2 = "../data/Mode Profiles/"+\
                    "OuterRingDesigned_L_inner_8um_gap_4um_2589nm_1535_1565_31wavls"+"/1550"

    foldername1 = "../data/Mode Profiles/"+\
                    "Ring_vertical_L_2_8um_1400_1700_31wavls_1500x1000"+"/1550"
    foldername2 = "../data/Mode Profiles/"+\
                    "Ring_vertical_L_2_8um_1400_1700_31wavls_1500x1000"+"/1550"

    foldername1 = "../data/Mode Profiles/"+\
                  "Straight_WG_width_2_8um_1480_1620_15wavls"+"/1550"
    foldername2 = "../data/Mode Profiles/"+\
                  "Straight_WG_width_2_8um_1480_1620_15wavls"+"/1550"

    # CoupledWG = Coupled_Waveguides(1.99,1.45,
    #                             gap_x = 2, gap_y = 0,
    #                             wavelength=1.55e-6,
    #                             name1=foldername1,
    #                             name2=foldername2,
    #                             ModeIdx1=1, ModeIdx2=1,
    #                             param_file_name="./config/Param_straight_2_8.csv",
    #                             Plot_field=True,
    #                             Plot_index=True)
    '''
    L_inner 8um gap 5um :
    wavl_idx =  10
    wavelength = 1.54 um
    Supermode 1 coefficient : A = complex(0.998300,0.000000), B = complex(-0.058282,0.000000)
    Supermode 2 coefficient : A = complex(0.058794,0.000000), B = complex(0.998270,0.000000)
    wavl_idx =  20
    wavelength = 1.55 um
    Supermode 1 coefficient : A = complex(0.083300,-0.000000), B = complex(0.996525,0.000000)
    Supermode 2 coefficient : A = complex(0.996415,0.000000), B = complex(-0.084598,-0.000000)
    wavl_idx =  30
    wavelength = 1.56 um
    Supermode 1 coefficient : A = complex(0.029222,-0.000000), B = complex(0.999573,0.000000)
    Supermode 2 coefficient : A = complex(0.999554,0.000000), B = complex(-0.029876,-0.000000)
    '''
    # CoupledWG.Plot_field_profile([[complex(0.998300,0.000000),
    #                                complex(-0.058282,0.000000)],
    #                               [complex(0.058794,0.000000),
    #                                complex(0.998270,0.000000)]],
    #                              field_name = 'Ex',
    #                              title=r"Electric field profile when wavl = 1540 nm",
    #                              Plot_log=False)
    # CoupledWG.Plot_field_profile([[complex(0.083300,-0.000000),
    #                                complex(0.996525,0.000000)],
    #                               [complex(0.996415,0.000000),
    #                                complex(-0.084598,-0.000000)]],
    #                              field_name = 'Ex',
    #                              title=r"Electric field profile when wavl = 1550 nm",
    #                              Plot_log=False)
    # CoupledWG.Plot_field_profile([[complex(0.029222,-0.000000),
    #                                complex(0.999573,0.000000)],
    #                               [complex(0.999554,0.000000),
    #                                complex(-0.029876,-0.000000)]],
    #                              field_name = 'Ex',
    #                              title=r"Electric field profile when wavl = 1560 nm",
    #                              Plot_log=False)

    '''
    L_inner 8um gap 4um:
    wavl_idx =  5
    wavelength = 1.5400 um
    Supermode 1 coefficient : A = complex(0.991335,0.000000), B = complex(-0.131355,0.000000)
    Supermode 2 coefficient : A = complex(0.132438,0.000000), B = complex(0.991191,0.000000)
    wavl_idx =  15
    wavelength = 1.5500 um
    Supermode 1 coefficient : A = complex(0.494740,-0.000000), B = complex(0.869041,0.000000)
    Supermode 2 coefficient : A = complex(0.865982,0.000000), B = complex(-0.500076,-0.000000)
    wavl_idx =  25
    wavelength = 1.5600 um
    Supermode 1 coefficient : A = complex(0.120213,-0.000000), B = complex(0.992748,0.000000)
    Supermode 2 coefficient : A = complex(0.992455,0.000000), B = complex(-0.122609,-0.000000)
    '''

    # CoupledWG.Plot_field_profile([[complex(0.999579,0.000000),
    #                                complex(-0.029000,0.000000)],
    #                               [complex(0.029154,0.000000),
    #                                complex(0.999575,0.000000)]],
    #                              field_name = 'Ex',
    #                              title=r"Electric field profile when wavl = 1540 nm",
    #                              Plot_log=False)
    # CoupledWG.Plot_field_profile([[complex(0.907773,0.000000),
    #                                complex(-0.419463,0.000001)],
    #                               [complex(0.423688,0.000000),
    #                                complex(0.905808,0.000000)]],
    #                              field_name = 'Ex',
    #                              title=r"Electric field profile when wavl = 1550 nm",
    #                              Plot_log=False)
    # CoupledWG.Plot_field_profile([[complex(-0.042748,-0.000000),
    #                                complex(0.999086,0.000000)],
    #                               [complex(0.999120,0.000000),
    #                                complex(0.041949,-0.000000)]],
    #                              field_name = 'Ex',
    #                              title=r"Electric field profile when wavl = 1560 nm",
    #                              Plot_log=False)

    ############# Parallel double rings L 2.8um ##################
    # foldername1 = "../data/Mode Profiles/"+\
    #               "Straight_WG_width_2_8um_1480_1620_15wavls"
    # foldername2 = "../data/Mode Profiles/"+\
    #               "Straight_WG_width_2_8um_1480_1620_15wavls"
    ##############################################################

    ################## Vertical double rings #####################
    # foldername1 = "../data/Mode Profiles/"+\
    #                 "Ring_vertical_L_2_8um_1400_1700_31wavls_1500x1000"
    # foldername2 = "../data/Mode Profiles/"+\
    #                 "Ring_vertical_L_2_8um_1400_1700_31wavls_1500x1000"
    ##############################################################

    ################## Vertical double rings #####################
    filename_uncoupled_gapy_2um     = "../data/beta_uncoupled_gap_0_000_2_000.txt"
    filename_coupled_gapy_2um       = "../data/beta_coupled_gap_0_000_2_000.txt"
    filename_uncoupled_gapy_0_2um   = "../data/beta_uncoupled_gap_0_000_0_200.txt"
    filename_coupled_gapy_0_2um     = "../data/beta_coupled_gap_0_000_0_200.txt"
    ##############################################################

    # unit: nm
    # wavl_arr = np.linspace(1400,1700,31)
    # wavl_arr = np.linspace(1500,1600,11)
    wavl_arr = np.linspace(1535,1565,31)
    # wavl_arr = np.linspace(1540,1560,21)
    # wavl_arr = np.linspace(1480,1620,15)
    # wavl_arr = np.linspace(1505,1595,19)

    # param_filename = "./config/Param_vertical.csv"
    param_filename = "./config/Param_L_inner_8.csv"
    # param_filename = "./config/Param_straight_2_8.csv"
    # param_filename = "./config/Param_800x400.csv"

    ############################################################################################
    # The folders containing the mode profiles
    foldername1 = "../data/Mode Profiles/"+"InnerRing_L_inner_8um_1535_1565_31wavls"
    foldername2 = "../data/Mode Profiles/"+"OuterRingDesigned_L_inner_8um_gap_5um_L_outer_2436nm_1535_1565nm_31wavls"
    # foldername2 = "../data/Mode Profiles/"+"OuterRingDesigned_L_inner_8um_gap_4um_2589nm_1535_1565_31wavls"

    # foldername1 = "../data/Mode Profiles/"+"InnerRing_L_inner_8um_1500_1600_21wavls_2000x800"
    # foldername2 = "../data/Mode Profiles/"+"OuterRingDesigned_L_inner_8um_gap_3um_L_2762nm_1500_1600_21wavls"

    # foldername1 = "../data/Mode Profiles/"+"InnerRing_Lx_8um_bendr_1000um_20x10um_800x400cells_beta_ang_21_wavls"
    # foldername2 = "../data/Mode Profiles/"+"OuterRingDesigned_20x10um_800x400cells_21_wavls"
    ############################################################################################

    ############################################################################################
    # Sweep gap between 2D concentric rings around a certian value, so the mode profiles can be seen as the same
    # if num_of_gaps is 1, the gap is fixed
    start_gap_x   = 2.5
    end_gap_x     = 3.5
    num_of_gaps_x = 11
    gap_arr_x     = np.linspace(start_gap_x,end_gap_x,num_of_gaps_x)

    start_gap_y   = 0
    end_gap_y     = 0
    num_of_gaps_y = 1
    gap_arr_y     = np.linspace(start_gap_y,end_gap_y,num_of_gaps_y)

    A,B = np.meshgrid(gap_arr_x,gap_arr_y)
    gap_arr = np.column_stack((A.ravel(), B.ravel()))

    sweeper = Parameter_sweeper(wavl_arr,gap_arr,
                                index_file1 = "Si3N4_index_Luke.txt",
                                index_file2 = "SiO2_index.txt",
                                foldername1 = foldername1,
                                foldername2 = foldername2,
                                param_filename = param_filename)
    # sweeper.Scan_gap(calc_needed=False, foldername = "Dispersion_scan_near_3um_CMT/", num_of_wavl_pts=1000)
    # sweeper.Scan_gap(calc_needed=True, num_of_wavl_pts=1000)
    ############################################################

    ############################################################################################
    # Analyze Dispersion of 2D concentric rings
    # Calculate Dispersion using propagation constants calculated by CMT and FDE
    param_filename = "./config/Param_L_inner_8.csv"
    wavl_arr = np.linspace(1540,1560,21)

    filename_uncoupled_gap5um = "./results/2D concentric rings/Propagation Constants using CMT/beta_uncoupled_gap_5_000_0_000.txt"
    filename_coupled_gap5um = "./results/2D concentric rings/Propagation Constants using CMT/beta_coupled_gap_5_000_0_000.txt"
    filename_coeffi_gap5um = "./results/2D concentric rings/Supermode Coefficients using CMT/gapx_5um_coeffi.txt"
    filename_FDE_beta_gap5um = "./results/2D concentric rings/Propagation Constant using FDE/Lumerical_supermodes_results_gapx_5um.txt"
    gap_arr = np.array([[5,0]])
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap5um, filename_coupled_gap5um, num_of_pts=200,
                             filename_FDE_beta=filename_FDE_beta_gap5um,
                             save_D_in_csv=False, filename_coeffi=filename_coeffi_gap5um)

    filename_uncoupled_gap4um = "./results/2D concentric rings/Propagation Constants using CMT/beta_uncoupled_gap_4_000_0_000.txt"
    filename_coupled_gap4um = "./results/2D concentric rings/Propagation Constants using CMT/beta_coupled_gap_4_000_0_000.txt"
    filename_coeffi_gap4um = "./results/2D concentric rings/Supermode Coefficients using CMT/gapx_4um_coeffi.txt"
    filename_FDE_beta_gap4um = "./results/2D concentric rings/Propagation Constant using FDE/Lumerical_supermodes_results_gapx_4um.txt"
    gap_arr = np.array([[4,0]])
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap4um, filename_coupled_gap4um, num_of_pts=100,
                             filename_FDE_beta=filename_FDE_beta_gap4um,
                             save_D_in_csv=False, filename_coeffi=filename_coeffi_gap4um)

    # Compare with Dispersion calculated using FDE
    wavl_arr = np.linspace(1505,1595,19)

    filename_uncoupled_gap3um = "./results/2D concentric rings/Propagation Constants using CMT/beta_uncoupled_gap_3_000_0_000.txt"
    filename_coupled_gap3um = "./results/2D concentric rings/Propagation Constants using CMT/beta_coupled_gap_3_000_0_000.txt"
    filename_coeffi_gap3um = "./results/2D concentric rings/Supermode Coefficients using CMT/gapx_3um_coeffi.txt"
    filename_FDE_D_gap3um = "./results/2D concentric rings/Dispersion using FDE/gapx_3um.txt"
    gap_arr = np.array([[3,0]])
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap3um, filename_coupled_gap3um, num_of_pts=150,
                             filename_FDE_D=filename_FDE_D_gap3um,
                             save_D_in_csv=False, filename_coeffi=filename_coeffi_gap3um)

    filename_uncoupled_gap2_8um = "./results/2D concentric rings/Propagation Constants using CMT/beta_uncoupled_gap_2_800_0_000.txt"
    filename_coupled_gap2_8um = "./results/2D concentric rings/Propagation Constants using CMT/beta_coupled_gap_2_800_0_000.txt"
    filename_coeffi_gap2_8um = "./results/2D concentric rings/Supermode Coefficients using CMT/gapx_2.8um_coeffi.txt"
    filename_FDE_D_gap2_8um = "./results/2D concentric rings/Dispersion using FDE/gapx_2.8um.txt"
    gap_arr = np.array([[2.8,0]])
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap2_8um, filename_coupled_gap2_8um, num_of_pts=100,
                             filename_FDE_D=filename_FDE_D_gap2_8um,
                             save_D_in_csv=False, filename_coeffi=filename_coeffi_gap2_8um)
    ############################################################################################

    wavl_arr = np.linspace(1400,1700,31)
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                            filename_uncoupled_gapy_2um, filename_coupled_gapy_2um, num_of_pts=100,
                            save_D_in_csv=True, save_mode = "AS")
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                            filename_uncoupled_gapy_2um, filename_coupled_gapy_2um, num_of_pts=100,
                            save_D_in_csv=True, save_mode = "S")

