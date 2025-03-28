import numpy as np
import matplotlib.pyplot as plt

def str2complex(s):
        str = s.replace(" ","")\
                .replace("(","")\
                .replace(")","")\
                .replace("i","j")
        return complex(str)

def Plot_field_profile(component,component_name,
                       save_name='./results/field_profile.png',dpi=300):
    fonttype = "Helvetica"
    fontsize = 4
    grid_linewidth = 1
    colormap = "jet"

    fig, ax = plt.subplots(1,3,figsize=(10, 4),dpi=dpi)
    plt.subplots_adjust(left=0.05, right=0.95, wspace =0.3, hspace =0.6)   #调整子图间距

    plt.yticks(np.shape(component[0])[0],fontproperties = fonttype, size = fontsize)
    plt.xticks(fontproperties = fonttype, size = fontsize)
    plt.rcParams["font.family"] = fonttype
    plt.rcParams.update({'font.size': fontsize})
    plt.ylabel('Y', fontdict={'family' : fonttype, 'size' : fontsize})
    plt.xlabel('X', fontdict={'family' : fonttype, 'size' : fontsize})
    plt.legend()

    name_list = ['Abs','Re','Im']
    component_list = [np.abs(component),np.real(component),np.imag(component)]
    for idx in range(3):
        im = ax[idx].imshow(component_list[idx], cmap=colormap)
        ax[idx].set_title(name_list[idx]+'('+component_name+')')
        cbar = fig.colorbar(im, ax=ax[idx], orientation='vertical',
                            label='', shrink=0.3, pad=0.02)
        ax[idx].set_xlabel("X")
        ax[idx].set_ylabel("Y")
        ax[idx].invert_yaxis()
        ax[idx].tick_params(axis='both',labelsize=5)
    plt.savefig(save_name,dpi=dpi)
    plt.show()
    return


    # length of the arr will reduce by 2
def First_derivative_central_diff(y,x):
    # 计算差分
    dy = np.diff(y)  # y[i+1] - y[i]
    dx = np.diff(x)  # x[i+1] - x[i]

    # 中心差分
    dy_central = (y[2:] - y[:-2])  # y[i+1] - y[i-1]
    dx_central = (x[2:] - x[:-2])  # x[i+1] - x[i-1]

    # 计算导数
    derivative = dy_central / dx_central
    return derivative