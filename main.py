'''
Main program for the simulation of two coupled microrings or bent waveguides
'''
from Coupled_Waveguides import *
from Parameter_sweeper import *
# from Plot_functions import Plot_field_profile

if __name__ == '__main__':

    # foldername1 = "data/Mode Profiles/"+\
    #                 "InnerRingXYZ_Lx_8um_bendr_1000um_20x10um_800x400cells_beta_ang_21_wavls"+"/1550"
    # foldername2 = "data/Mode Profiles/"+\
    #                 "OuterRingDesigned_20x10um_800x400cells_21_wavls"+"/1550"

    # foldername1 = "data/Mode Profiles/"+\
    #                 "InnerRing_Lx_8um_bendr_1000um_30x12um_2000x800cells_1540_1560nm_21wavls"+"/1550"
    # foldername2 = "data/Mode Profiles/"+\
    #                 "OuterRingDesigned_30x12um_2000x800cells_1540_1560nm_21wavls"+"/1550"

    # CoupledWG = Coupled_Waveguides(1.99,1.45,gap_x = 5,gap_y = 0,
    #                             wavelength=1.55e-6,
    #                             name1=foldername1,
    #                             name2=foldername2,
    #                             ModeIdx1=1, ModeIdx2=1,
    #                             param_file_name="./Param.csv",
    #                             Plot_field=False,
    #                             Plot_index=True
    #                             )
    # CoupledWG.Plot_field_profile([0.707,-0.707],'Ex',Plot_log=False)

    foldername1 = "data/Mode Profiles/"+\
                    "InnerRing_Lx_8um_bendr_1000um_30x12um_2000x800cells_1540_1560nm_21wavls"
    foldername2 = "data/Mode Profiles/"+\
                    "OuterRingDesigned_30x12um_2000x800cells_1540_1560nm_21wavls"

    # foldername1 = "data/Mode Profiles/"+"InnerRingXYZ_width_8um_bendr_1000um_400x200cells_beta_ang"
    # foldername2 = "data/Mode Profiles/"+"OuterRingDesigned_400x200cells"
    #unit: um
    #wavl_vec = np.linspace(1.54,1.56,21)

    wavl_vec = np.linspace(1.55,1.55,1).reshape(1,)
    gap_vec = np.array([[5,0],
                        [0,5]])
    sweeper = Parameter_sweeper(wavl_vec,gap_vec,
                                foldername1,foldername2,
                                param_filename="./Param.csv")
    sweeper.Save_results()




