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
    foldername1 = "data/Mode Profiles/"+\
                    "InnerRing_Lx_8um_bendr_1000um_20x10um_800x400cells_beta_ang_21_wavls"+"/1550"
    # foldername1 = "data/Mode Profiles/"+\
    #                 "OuterRingDesigned_20x10um_800x400cells_21_wavls"+"/1550"
    foldername2 = "data/Mode Profiles/"+\
                    "OuterRingDesigned_20x10um_800x400cells_21_wavls"+"/1550"

    # foldername1 = "data/Mode Profiles/"+\
    #                 "InnerRing_L_inner_8um_1400_1700_31wavls"+"/1550"
    # foldername2 = "data/Mode Profiles/"+\
    #                 "InnerRing_L_inner_8um_1400_1700_31wavls"+"/1550"

    foldername1 = "data/Mode Profiles/"+\
                    "InnerRing_L_inner_8um_1535_1565_31wavls"+"/1550"
    foldername2 = "data/Mode Profiles/"+\
                    "OuterRing_L_inner_8um_gap_5um_L_outer_2436nm_1535_1565nm_31wavls"+"/1550"

    CoupledWG = Coupled_Waveguides(1.99,1.45,
                                gap_x = 5, gap_y = 0,
                                wavelength=1.55e-6,
                                name1=foldername1,
                                name2=foldername2,
                                ModeIdx1=1, ModeIdx2=1,
                                param_file_name="./Param_L_inner_8.csv",
                                Plot_field=False,
                                Plot_index=False)
    CoupledWG.Plot_field_profile([[0.866886,-0.503335],
                                  [0.498506 ,0.864092]],
                                 field_name = 'Ex',
                                 title=r"Electric field profile",
                                 Plot_log=False)

    ########################################
    # Scan the gap between waveguides
    ########################################

    ################## Concentric double rings ###################
    foldername1 = "data/Mode Profiles/"+\
                    "InnerRing_L_inner_8um_1535_1565_31wavls"
    foldername2 = "data/Mode Profiles/"+\
                    "OuterRing_L_inner_8um_gap_5um_L_outer_2436nm_1535_1565nm_31wavls"

    # foldername2 = "data/Mode Profiles/"+\
    #                 "OuterRingDesigned_L_inner_8um_gap_4um_2589nm_1535_1565_31wavls"

    # foldername1 = "data/Mode Profiles/"+\
    #                 "InnerRing_Lx_8um_bendr_1000um_20x10um_800x400cells_beta_ang_21_wavls"
    # foldername2 = "data/Mode Profiles/"+\
    #                 "OuterRingDesigned_20x10um_800x400cells_21_wavls"

    # foldername1 = "data/Mode Profiles/"+\
    #               "InnerRing_L_inner_2_8um_1540_1560_11wavls"
    # foldername2 = "data/Mode Profiles/"+\
    #               "OuterRingDesigned_L_inner_2_8um_gap_8um_1540_1560_21wavls"
    ############################################################

    ################## Vertical double rings ###################
    # foldername1 = "data/Mode Profiles/"+\
    #                 "InnerRing_L_inner_8um_1400_1700_31wavls"
    # foldername2 = "data/Mode Profiles/"+\
    #                 "InnerRing_L_inner_8um_1400_1700_31wavls"
    ############################################################

    # filename_uncoupled_gap4um = "results/L_inner_8um_gapx_5&4um/gapx_4um/beta_uncoupled_L_inner_8um_gapx_4um_using_FDE@gap4um.txt"
    # filename_coupled_gap4um = "results/L_inner_8um_gapx_5&4um/gapx_4um/beta_coupled_L_inner_8um_gapx_4um_using_FDE@gap4um.txt"

    # filename_uncoupled_gap5um = "data/beta_uncoupled_gap_5.txt"
    # filename_coupled_gap5um = "data/beta_coupled_gap_5.txt"

    # filename_lumerical_gap4um = "results/L_inner_8um_gapx_5&4um/gapx_4um/Lumerical_supermodes_results_gapx_4um.txt"
    # filename_lumerical_gap5um = "data/Lumerical_supermodes_results_gap_5.txt"

    # unit: nm
    # wavl_arr = np.linspace(1400,1700,31)
    # wavl_arr = np.linspace(1500,1600,11)
    wavl_arr = np.linspace(1535,1565,31)

    # param_filename = "./Param_vertical.csv"
    param_filename = "./Param_L_inner_8.csv"
    # param_filename = "./Param_L_inner_2_8.csv"
    # param_filename = "./Param_800x400.csv"

    start_gap_x   = 0.5
    end_gap_x     = 2.5
    num_of_gaps_x = 3
    gap_arr_x = np.linspace(start_gap_x,end_gap_x,num_of_gaps_x)

    start_gap_y   = 0
    end_gap_y     = 0
    num_of_gaps_y = 1
    gap_arr_y = np.linspace(start_gap_y,end_gap_y,num_of_gaps_y)

    A,B = np.meshgrid(gap_arr_x,gap_arr_y)
    gap_arr = np.column_stack((A.ravel(), B.ravel()))

    sweeper = Parameter_sweeper(wavl_arr,gap_arr,
                                index_file1 = "Si3N4_index_Luke.txt",
                                index_file2 = "SiO2_index.txt",
                                foldername1=foldername1,foldername2=foldername2,
                                param_filename=param_filename)
    # beta_uncoupled_arr,beta_coupled_arr,beta_ave_uncoupled_arr = sweeper.Scan_wavl(gap_idx=0)
    # sweeper.Scan_gap()

    # beta_data_filename = (filename_uncoupled_gap5um, filename_coupled_gap5um)
    # Analyzer = Data_analyzer(wavl_arr, gap_arr[0,:], beta_data_filename, param_filename)
