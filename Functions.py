import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import CubicSpline
import json

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
ylim: format : (ymin, ymax). set the y range of the whole figure
AD_region_color: whether to fill color in region of y>0
bbox_to_anchor: coordinates for the legends
text: used to show comments
dpi: 300 by default
'''
def Plot_curve(data_arr,
                default_param_filename = "Param_plot_curve.json",
                *args, **kwargs):
    with open(default_param_filename,'r', encoding='UTF-8') as f:
        param_dict = json.load(f)

    for key,value in kwargs.items():
        param_dict[key] = value

    savename = param_dict["foldername"]+str(param_dict["title"])+".jpg"
    idx = 0
    plt.figure(figsize=param_dict["figsize"])
    for data_idx in range(len(data_arr)):
        X       = data_arr[data_idx][:,0]           #X is a 2D array
        Y_arr   = data_arr[data_idx][:,1:]          #Y_arr is a 2D array
        for Y_idx in range(np.shape(Y_arr)[1]):
            # use the same legend as the last curve
            plt.plot(X,Y_arr[:,Y_idx],label=param_dict["Y_legends"][idx],
                    color = param_dict["colors_list"][idx],
                    marker = param_dict["marker_list"][idx],
                    linestyle = param_dict["linestyle_list"][idx],
                    linewidth = param_dict["plot_linewidth"][idx],
                    alpha = param_dict["alpha_list"][idx])
            idx = idx + 1

    plt.rcParams["font.family"] = param_dict["fonttype"]
    plt.rcParams.update({'font.size': param_dict["fontsize"]})
    plt.title(param_dict["title"])

    if len(param_dict["xlim"])>0:
        xlim = param_dict["xlim"]
        plt.xlim(min(xlim),max(xlim))
    xmin = plt.xlim()[0]
    xmax = plt.xlim()[1]

    if len(param_dict["ylim"])>0:
        ylim = param_dict["ylim"]
        plt.ylim(min(ylim),max(ylim))
    ymin = plt.ylim()[0]
    ymax = plt.ylim()[1]

    if param_dict["AD_region_color"] and ymax>0:
        plt.axhspan(0, ymax, color='green', alpha=0.1,
                    label='Anamolous Dispersion Region')
        plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        plt.ylim(ymin,ymax)

    if len(param_dict["xtickslabel"])>1:
        xticks = np.array(param_dict["xticks"])
        xtickslabel = np.array(param_dict["xtickslabel"])
        xticks_mask = np.where((xticks<=xmax) & (xticks>=xmin))
        xticks = xticks[xticks_mask]
        xtickslabel = xtickslabel[xticks_mask]

        plt.xticks(xticks,xtickslabel,
                   fontproperties = param_dict["fonttype"],
                   size = param_dict["fontsize"])

    if len(param_dict["yticks"])>1:
        yticks = ticks_arr([ymin,ymax])
        plt.yticks(yticks,yticks,
                   fontproperties = param_dict["fonttype"],
                   size = param_dict["fontsize"])

    plt.ylabel(param_dict["Y_label"],
               fontdict={'family' : param_dict["fonttype"],
                        'size' : param_dict["fontsize"]})
    plt.xlabel(param_dict["X_label"],
               fontdict={'family' : param_dict["fonttype"],
                         'size' : param_dict["fontsize"]})
    if len(param_dict["bbox_legend"])>0:
        plt.legend(bbox_to_anchor=param_dict["bbox_legend"],
                   borderaxespad=0)
    else:
        plt.legend(loc='best')

    if not param_dict["text"] == "":
        plt.text(xmin + (xmax-xmin)*param_dict["loc_text"][0],
                 ymin + (ymax-ymin)*param_dict["loc_text"][1],
                param_dict["text"],linespacing = param_dict["linespacing"],
                bbox=dict(boxstyle="round,pad=0.9", fc="white", alpha=0.9))
    plt.grid(linewidth=param_dict["grid_linewidth"], alpha=0.3)
    plt.tight_layout()
    plt.savefig(savename,dpi=param_dict["dpi"])
    if param_dict["plot_show"]:
        plt.show()

def Interpol(x,y,x_intp):
    cs = CubicSpline(x, y, bc_type='natural')
    # bc_type 可选 'natural', 'clamped', 'periodic' 等
    y_intp = cs(x_intp)
    return y_intp

# Gives a better array of ticks than the defaults
def ticks_arr(data_arr, small_ticks = False):
    # For data>0, the largest base. E.g. if max = 7875, largest base is 1000.
    if small_ticks:
        largest_divider_max = 10**int(np.log10(np.max(data_arr))-0.5)
    else:
        largest_divider_max = 10**int(np.log10(np.max(data_arr)))
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



def Plot_im(data_arr, point_arr = [],
            default_param_filename = "Param_plot_image.json",
            *args, **kwargs):

    with open(default_param_filename,'r', encoding='UTF-8') as f:
        param_dict = json.load(f)

    for key,value in kwargs.items():
        param_dict[key] = value

    fig, ax = plt.subplots(figsize=param_dict["figsize"],
                           dpi = param_dict["dpi"])
    if param_dict["norm"] == "zero_divided":
        norm = TwoSlopeNorm(vmin=np.min(data_arr), vcenter=0, vmax=np.max(data_arr))
    else:
        norm = None
    im  = ax.imshow(data_arr,cmap = param_dict["colormap"],
                    norm=norm,aspect=param_dict["aspect"])
    if len(point_arr) > 0:
        ax.scatter(point_arr[:,0],point_arr[:,1],
                   s=param_dict["point_size"], c=param_dict["point_color"], marker='^')

    # shrink: 缩放比例，pad: 间距
    cbar = fig.colorbar(im, ax=ax, orientation=param_dict["cbar_orientation"],
                        shrink=param_dict["shrink"], pad=param_dict["pad"])
    cbar.set_ticks(ticks_arr(data_arr,param_dict["cbar_small_ticks"]))
    cbar.ax.tick_params(labelsize=param_dict["fontsize"]*0.8)
    cbar.set_label(param_dict["cbar_label"])

    if len(param_dict["ytickslabel"]) >0:
        plt.yticks(param_dict["yticks"],
                   param_dict["ytickslabel"],
                   fontproperties = param_dict["fonttype"],
                   size = param_dict["fontsize"])
    if len(param_dict["xtickslabel"]) >0:
        plt.xticks(param_dict["xticks"],
                   param_dict["xtickslabel"],
                   fontproperties = param_dict["fonttype"],
                   size = param_dict["fontsize"])

    plt.xlabel(param_dict["xlabel"])
    plt.ylabel(param_dict["ylabel"])
    plt.title(param_dict["title"])
    plt.savefig("./results/"+param_dict["title"]+".jpg")
    plt.show()