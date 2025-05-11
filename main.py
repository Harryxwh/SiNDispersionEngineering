'''
Main program for the simulation of two coupled microrings or bent waveguides

Author: Weihao Xu
Date: May. 11th, 2025

'''
import numpy as np
from Coupled_Waveguides import Coupled_Waveguides
from Parameter_sweeper import Parameter_sweeper
from Data_analyzer import Data_analyzer

if __name__ == '__main__':

    ############################################################################################
    # SETTINGS for 2D concentric rings

    # The folders containing the mode profiles
    foldername1 = "../data/Mode Profiles/"+"InnerRing_L_inner_8um_1500_1600_21wavls_2000x800"
    # foldername1 = "../data/Mode Profiles/"+"InnerRing_L_inner_8um_1535_1565_31wavls"
    foldername2 = "../data/Mode Profiles/"+"OuterRingDesigned_L_inner_8um_gap_3um_L_2762nm_1500_1600_21wavls"
    # foldername2 = "../data/Mode Profiles/"+"OuterRingDesigned_L_inner_8um_gap_4um_2589nm_1535_1565_31wavls"
    # foldername2 = "../data/Mode Profiles/"+"OuterRingDesigned_L_inner_8um_gap_5um_L_outer_2436nm_1535_1565nm_31wavls"

    # Parameters of the FDE similation
    param_filename = "./config/Param_L_inner_8.csv"

    # Wavlength range for the simulation
    wavl_arr = np.linspace(1505,1595,19)        # for gap = 3um
    #  wavl_arr = np.linspace(1535,1565,31)       # for gap = 4um or 5um

    working_dir = "./Mode number conserved coupling/results/2D concentric rings/"
    ############################################################################################

    ############################################################################################
    # SETTINGS for 2D parallel rings

    # param_filename = "./config/Param_straight_2_8.csv"

    # wavl_arr = np.linspace(1480,1620,15)

    # foldername1 = "../data/Mode Profiles/"+\
    #               "Straight_WG_width_2_8um_1480_1620_15wavls"
    # foldername2 = "../data/Mode Profiles/"+\
    #               "Straight_WG_width_2_8um_1480_1620_15wavls"
    ############################################################################################

    ############################################################################################
    # SETTINGS for 3D concentric rings

    # param_filename = "./config/Param_vertical.csv"

    # wavl_arr = np.linspace(1400,1700,31)

    # foldername1 = "../data/Mode Profiles/"+\
    #                 "Ring_vertical_L_2_8um_1400_1700_31wavls_1500x1000"
    # foldername2 = "../data/Mode Profiles/"+\
    #                 "Ring_vertical_L_2_8um_1400_1700_31wavls_1500x1000"
    ############################################################################################


    ############################################################################################
    # Sweep gap between 2D concentric rings around a certian value, so the mode profiles can be seen as the same
    # if num_of_gaps is 1, the gap is fixed

    start_gap_x   = 3
    end_gap_x     = 3
    num_of_gaps_x = 1
    gap_arr_x     = np.linspace(start_gap_x,end_gap_x,num_of_gaps_x)

    start_gap_y   = 0
    end_gap_y     = 0
    num_of_gaps_y = 1
    gap_arr_y     = np.linspace(start_gap_y,end_gap_y,num_of_gaps_y)

    A,B = np.meshgrid(gap_arr_x,gap_arr_y)
    gap_arr = np.column_stack((A.ravel(), B.ravel()))

    sweeper = Parameter_sweeper(wavl_arr,gap_arr,
                                foldername1 = foldername1,
                                foldername2 = foldername2,
                                param_filename = param_filename)
    # sweeper.Scan_gap(calc_needed=True, num_of_wavl_pts=1000)
    # sweeper.Scan_gap(calc_needed= False,
    #                  foldername = working_dir + "Supermodes attributes using CMT/",
    #                  num_of_wavl_pts=1000)
    ############################################################################################

    ############################################################################################
    # Analyze Dispersion of 2D concentric rings
    # Calculate Dispersion using propagation constants calculated by CMT and FDE

    param_filename = "./config/Param_L_inner_8.csv"

    wavl_arr = np.linspace(1540,1560,21)

    # gap = 5um
    gap_arr = np.array([[5,0]])
    filename_uncoupled_gap5um = working_dir + "Supermodes attributes using CMT/beta_uncoupled_gap_5_000_0_000.txt"
    filename_coupled_gap5um = working_dir + "Supermodes attributes using CMT/beta_coupled_gap_5_000_0_000.txt"
    filename_FDE_beta_gap5um = working_dir + "Supermodes attributes using FDE/Lumerical_supermodes_results_gapx_5um.txt"

    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap5um, filename_coupled_gap5um,
                             foldername_1 = foldername1, foldername_2 = foldername2,
                             filename_FDE_beta=filename_FDE_beta_gap5um,
                             save_D_in_csv=True)

    # gap = 4um
    gap_arr = np.array([[4,0]])
    filename_uncoupled_gap4um = working_dir + "Supermodes attributes using CMT/beta_uncoupled_gap_4_000_0_000.txt"
    filename_coupled_gap4um = working_dir + "Supermodes attributes using CMT/beta_coupled_gap_4_000_0_000.txt"
    filename_FDE_beta_gap4um = working_dir + "Supermodes attributes using FDE/Lumerical_supermodes_results_gapx_4um.txt"
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap4um, filename_coupled_gap4um,
                             foldername_1 = foldername1, foldername_2 = foldername2,
                             filename_FDE_beta=filename_FDE_beta_gap4um,
                             save_D_in_csv=True)
    ############################################################################################

    ############################################################################################
    # Compare with Dispersion calculated using FDE

    wavl_arr = np.linspace(1505,1595,19)

    # gap = 3um
    gap_arr = np.array([[3,0]])
    filename_uncoupled_gap3um = working_dir + "Supermodes attributes using CMT/beta_uncoupled_gap_3_000_0_000.txt"
    filename_coupled_gap3um = working_dir + "Supermodes attributes using CMT/beta_coupled_gap_3_000_0_000.txt"
    filename_FDE_D_gap3um = working_dir + "Dispersion using FDE/gapx_3um.txt"
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap3um, filename_coupled_gap3um,
                             foldername_1 = foldername1, foldername_2 = foldername2,
                             filename_FDE_D=filename_FDE_D_gap3um,
                             save_D_in_csv=False)

    # gap = 2.8um
    gap_arr = np.array([[2.8,0]])
    filename_uncoupled_gap2_8um = working_dir + "Supermodes attributes using CMT/beta_uncoupled_gap_2_800_0_000.txt"
    filename_coupled_gap2_8um = working_dir + "Supermodes attributes using CMT/beta_coupled_gap_2_800_0_000.txt"
    filename_FDE_D_gap2_8um = working_dir + "Dispersion using FDE/gapx_2.8um.txt"
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                             filename_uncoupled_gap2_8um, filename_coupled_gap2_8um,
                             foldername_1 = foldername1, foldername_2 = foldername2,
                             filename_FDE_D=filename_FDE_D_gap2_8um,
                             save_D_in_csv=False)
    ############################################################################################

    ############################################################################################
    # 3D concentric rings

    param_filename = "./config/Param_L_inner_8.csv"

    wavl_arr = np.linspace(1400,1700,31)

    gap_arr = np.array([[0,2]])

    filename_uncoupled_gapy_2um     = "../data/beta_uncoupled_gap_0_000_2_000.txt"
    filename_coupled_gapy_2um       = "../data/beta_coupled_gap_0_000_2_000.txt"
    filename_uncoupled_gapy_0_2um   = "../data/beta_uncoupled_gap_0_000_0_200.txt"
    filename_coupled_gapy_0_2um     = "../data/beta_coupled_gap_0_000_0_200.txt"

    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                            filename_uncoupled_gapy_2um, filename_coupled_gapy_2um,
                            foldername_1 = foldername1, foldername_2 = foldername2,
                            name_structure = "3D concentric rings",
                            num_of_pts=100, save_D_in_csv=True, save_mode = "AS")
    Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], param_filename,
                            filename_uncoupled_gapy_2um, filename_coupled_gapy_2um,
                            foldername_1 = foldername1, foldername_2 = foldername2,
                            name_structure = "3D concentric rings",
                            num_of_pts=100, save_D_in_csv=True, save_mode = "S")
    ############################################################################################
