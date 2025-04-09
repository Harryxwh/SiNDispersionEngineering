import numpy as np
import matplotlib.pyplot as plt

def str2complex(s):
        str = s.replace(" ","")\
                .replace("(","")\
                .replace(")","")\
                .replace("i","j")
        return complex(str)

def find_zero(X,Y):
    zero_idx_arr = []
    zero_arr     = []
    assert len(X) == len(Y)
    for idx in range(len(X)-1):
        if Y[idx] * Y[idx+1] == 0:
            if Y[idx] == 0:
                zero_arr.append(X[idx])
                zero_idx_arr.append(idx)
            else:
                zero_arr.append(X[idx+1])
                zero_idx_arr.append(idx+1)
        elif Y[idx] * Y[idx+1] < 0:
            zero = X[idx] + (X[idx+1]-X[idx])*\
                    np.abs(Y[idx])/(np.abs(Y[idx])+np.abs(Y[idx+1]))
            zero_arr.append(zero)
            zero_idx_arr.append(idx)
    return zero_arr,zero_idx_arr

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

'''
Plot_curve: used to plot multiple curves with different X data

parameters:
data_arr: format: ([x1,y11,y12,y13...],[x2,y21,y22,y23...]) in which y11,y12,y13 share x1 as the X data
Y_legends: format: np.array([label1,label2...]) (if multiple curves use the same label, append a "")
X_label: name of the X axis
Y_label: name of the Y axis
xticks: original xticks to be renamed
xtickslabel: new labels for xticks
yticks: used to control the number of ticks on y axis
title: title of the whole figure
marker_list: format:["",".","o"...] each markers corresp to a y data
linestyle_list: format:["-","--","dotted"...] each linestyle corresp to a y data
colors_list: format:['tab:red','tab:orange'...] each color corresp to a y data
figsize: format: (x,y)
fontsize: size of all the text
ylim: set the y range of the whole figure to (-ylim,ylim)
fill_color: whether to fill color in region of y>0
bbox_to_anchor: coordinates for the legends
text: used to show comments
dpi: 300 by default
'''
def Plot_curve(data_arr,Y_legends,
                X_label,Y_label,
                title,
                marker_list,linestyle_list,colors_list,
                figsize=(10,6),fontsize=10,
                xticks=[],xtickslabel=[],
                yticks = [],
                ylim=-1, fill_color=False,
                bbox_to_anchor=(),text="",
                dpi=300,plot_show=False):
    #Plot parameters
    figsize = figsize
    fontsize = fontsize
    fonttype = "Helvetica"
    grid_linewidth = 0.8
    plot_linewidth = 1.5
    # colors_list = ['tab:blue']*3+['tab:red']*3+['tab:orange']+['tab:green']
    savename = "results/"+str(title)+".jpg"
    idx = 0
    plt.figure(figsize=figsize)
    for data_idx in range(len(data_arr)):
        X       = data_arr[data_idx][:,0]           #X is a 2D array
        Y_arr   = data_arr[data_idx][:,1:]          #Y_arr is a 2D array
        for Y_idx in range(np.shape(Y_arr)[1]):
            if Y_legends[idx] == "":
                plt.plot(X,Y_arr[:,Y_idx],
                        color=colors_list[idx], marker=marker_list[idx],
                        linestyle=linestyle_list[idx], linewidth=plot_linewidth)
            else:
                plt.plot(X,Y_arr[:,Y_idx],label=Y_legends[idx],
                            color=colors_list[idx], marker=marker_list[idx],
                            linestyle=linestyle_list[idx], linewidth=plot_linewidth)
            idx = idx + 1

    plt.rcParams["font.family"] = fonttype
    plt.rcParams.update({'font.size': fontsize})
    plt.title(title)

    ymin = plt.ylim()[0]
    ymax = plt.ylim()[1]

    if ylim>0:
        plt.ylim(-ylim,ylim)
    if fill_color:
        plt.axhspan(0, ymax, color='green', alpha=0.1, label='Anamolous Dispersion Region')
        plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        plt.ylim(ymin,ymax)

    if len(xtickslabel)>1:
        plt.xticks(xticks,xtickslabel,fontproperties = fonttype, size = fontsize)
    if len(yticks)>1:
        yticks = ticks_arr([ymin,ymax])
        plt.yticks(yticks,yticks,fontproperties = fonttype, size = fontsize)

    plt.ylabel(Y_label, fontdict={'family' : fonttype, 'size' : fontsize})
    plt.xlabel(X_label, fontdict={'family' : fonttype, 'size' : fontsize})
    if len(bbox_to_anchor)>0:
                plt.legend(bbox_to_anchor=bbox_to_anchor, borderaxespad=0)
    else:
        plt.legend(loc='best')
    if not text == "":
        plt.text(np.quantile(X,0.5),np.quantile(Y_arr,0.01),
                text, bbox=dict(boxstyle="round,pad=0.9", fc="white", alpha=0.9))
    plt.grid(linewidth=grid_linewidth, alpha=0.3)
    plt.tight_layout()
    plt.savefig(savename,dpi=dpi)
    if plot_show:
        plt.show()


# Gives a better array of ticks than the defaults
def ticks_arr(data_arr):
    # For data>0, the largest base. E.g. if max = 7875, largest base is 1000.
    largest_divider_max = 10**int(np.log10(np.max(data_arr))-0.5)
    ticks_arr = np.arange(largest_divider_max,
                          (int(np.max(data_arr)/largest_divider_max)+1)*largest_divider_max,
                          largest_divider_max)
    # if min<0, than for data<0, the largest base.
    if np.min(data_arr)<0:
        largest_divider_min = 10**int(np.log10(-np.min(data_arr)))
        nega_arr = np.arange(int(np.min(data_arr)/largest_divider_min)*largest_divider_min,
                             0, largest_divider_min)
        ticks_arr = np.r_[nega_arr,np.array([0,]),ticks_arr]

    return ticks_arr