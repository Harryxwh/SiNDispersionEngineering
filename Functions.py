'''
These functions are frequently used in all other files of the project. They are therefore packaged together.

Author: Weihao Xu
Date: May. 5th, 2025
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import CubicSpline
import json

# Constants
c           = 299792458
mm          = 1e-3
um          = 1e-6
nm          = 1e-9
mu0         = 4 * np.pi * 1e-7
epsilon0    = 8.854187817e-12
n_air       =   1
component_name_list = ['Ex','Ey','Ez','Hx','Hy','Hz']


# Convert a string to a complex number
# e.g. str2complex("1+2i") = 1 + 2j
def str2complex(s):
        str = s.replace(" ","")\
                .replace("(","")\
                .replace(")","")\
                .replace("i","j")
        return complex(str)

# Find the zero points of a function
# X: X data, Y: Y data
# The function returns the zero points and their indices
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
    return zero_arr, zero_idx_arr


# Interpolate the data using cubic spline
# x: X data, y: Y data, x_intp: X data for interpolation
def Interpolation(x,y,x_intp):
    cs = CubicSpline(x, y, bc_type='natural')
    # bc_type 可选 'natural', 'clamped', 'periodic' 等
    y_intp = cs(x_intp)
    return y_intp

# Gives a better array of ticks than the defaults
# For example, if the data is between -2000 and 900, the ticks are set to be -2000, -1000, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900.
def Auto_ticks(data_arr, small_ticks = False):
    Auto_ticks = np.array([0,])
    # For data>0, the step is set be log10(np.max). E.g. if max = 7875, the step is set to be 1000.
    if np.max(data_arr)>0:
        # small_ticks: If True, when the max is below 10^3.5, the ticks are set with a step of 100, while when the max is above 10^3.5, the ticks are set with a step of 1000. If False, the ticks are always set with a step of 1000 as long as log10(np.max) > 3.
        if small_ticks:
            largest_divider_max = 10**int(np.log10(np.max(data_arr))-0.5)
        else:
            largest_divider_max = 10**int(np.log10(np.max(data_arr)))
        posi_arr= np.arange(largest_divider_max,
                            (int(np.max(data_arr)/largest_divider_max)+1)*largest_divider_max,
                            largest_divider_max)
        Auto_ticks = np.r_[Auto_ticks,posi_arr]
    # if min<0, than for data<0, the step is set to be log10(-np.min)
    if np.min(data_arr)<0:
        largest_divider_min = 10**int(np.log10(-np.min(data_arr)))
        nega_arr = np.arange(int(np.min(data_arr)/largest_divider_min)*largest_divider_min,
                             0, largest_divider_min)
        Auto_ticks = np.r_[nega_arr,Auto_ticks]

    return Auto_ticks

# Calculate the first order derivative of a function
# Note: length of the arr will be reduced by 2
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


def Load_material_index(wavl_um):
    n_Si3N4 = (1+3.0249/(1-(0.1353406/wavl_um)**2)+
                40314/(1-(1239.842/wavl_um)**2))**.5
    n_SiO2  = (1+0.6961663/(1-(0.0684043/wavl_um)**2)+
                0.4079426/(1-(0.1162414/wavl_um)**2)+
                0.8974794/(1-(9.896161/wavl_um)**2))**.5
    return n_Si3N4, n_SiO2



'''
Plot_curve():   used to plot multiple curves with different X data

Parameters:
data_arr                : format: ([x1,y11,y12,y13...],[x2,y21,y22,y23...]), where y11,y12,y13 share x1 as the X data
default_param_filename  : name of the json file where the default parameters are stored.

The following default parameters are used for plotting and are stored in a json file. (default: Param_plot_curve.json)
They can be changed using a dictionary loaded as **kwargs when necessary:
"title"             : title of the figure. Also used as a part of the filename
"figsize"           : a tuple or list with two elements.
"dpi"               : 300 by default.
"grid_linewidth"    : linewidth of the grids. 0.8 by default
"plot_linewidth"    : a list of the linewidth of the plotted lines (each item for a y).
"marker_list"       : a list of the marker of the plotted lines (each item for a y).
"linestyle_list"    : a list of the linestyle of the plotted lines (each item for a y).
"colors_list"       : a list of the color of the plotted lines (each item for a y).
"alpha_list"        : a list of the transparency of the plotted lines (each item for a y).
"Y_legends"         : a list of the lengends of the plotted lines (each item for a y).
"X_label"           : label for the X axis.
"Y_label"           : label for the Y axis.
"xticks"            : list of the location on X axis where a tick exists.
"xtickslabel"       : list of the labels to replace xticks.
"yticks"            : list of the location on Y axis where a tick exists.
"autoset_yticks"    : if 1, the yticks are set using the function "Auto_ticks". if 0, the default yticks are used.
"xlim"              : a tuple whose max and min values control the plotting range of X axis.
"ylim"              : a tuple whose max and min values control the plotting range of Y axis.
"fonttype"          : font used for texts in the figure.
"fontsize"          : fontsize of the texts in the figure.
"AD_region_color"   : if True, the area above 0 is colored in green.
"bbox_legend"       : a list of the location of the legend. (percentage of the figure size)
"text"              : text to be added to the figure.
"linespacing"       : spacing between lines of the text.
"loc_text"          : a list of the location of the text. (percentage of the figure size)
"foldername"        : folder where the figure is saved.
"comment"           : comment to be added to the filename.
"plot_show"         : if True, the figure is shown after being saved.
'''
def Plot_curve(data_arr,
                default_param_filename = "Param_plot_curve.json",
                *args, **kwargs):
    with open(default_param_filename,'r', encoding='UTF-8') as f:
        param_dict = json.load(f)

    for key,value in kwargs.items():
        param_dict[key] = value

    savename = param_dict["foldername"]+str(param_dict["title"])+param_dict['comment']+".pdf"
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
        yticks = np.array(param_dict["yticks"])
        yticks_mask = np.where((yticks<=ymax) & (yticks>=ymin))
        yticks = yticks[yticks_mask]
        ytickslabels = ["{:.2f}".format(ytick) for ytick in yticks]
        plt.yticks(yticks,ytickslabels,
                    fontproperties = param_dict["fonttype"],
                    size = param_dict["fontsize"])
    elif param_dict["autoset_yticks"] == 1:
        yticks = Auto_ticks([ymin,ymax])
        ytickslabels = ["{:.2f}".format(ytick) for ytick in yticks]
        plt.yticks(yticks,ytickslabels,
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
    plt.savefig(savename,format='pdf',dpi=param_dict["dpi"])

    if param_dict["plot_show"]:
        plt.show()
    plt.close()

'''
Plot_im():  Plot the image of a 2D array

Parameters:
data_arr                : 2D array to be plotted
point_arr               : 2D array of points to be marked on the image. (optional)
default_param_filename  : name of the json file where the default parameters are stored.

The following default parameters are used for plotting and are stored in a json file. (default: Param_plot_image.json)
They can be changed using a dictionary loaded as **kwargs when necessary:
point_arr           : 2D array of points to be marked on the image. (optional)
"title"             : title of the figure. Also used as a part of the filename
"figsize"           : a tuple or list with two elements.
"dpi"               : 300 by default.
"colormap"          : colormap of the image.
"norm"              : if "zero_in_center", the color map is divided into two parts, one for positive values and one for negative values.
"aspect"            : aspect ratio of the image.
"point_size"        : size of the points to be marked on the image.
"point_color"       : color of the points to be marked on the image.
"point_marker"      : marker of the points to be marked on the image.
"cbar_orientation"  : orientation of the color bar.
"shrink"            : shrink in size of the color bar. value between 0 and 1. default: 0.8.
"pad"               : distance between the color bar and the figure. value between 0 and 1. default: 0.05.
"cbar_small_ticks"  : the ticks of the color bar are set using the function "Auto_ticks". If True, small_ticks are enabled.
"cbar_label"        : label of the color bar.
"yticks"            : list of the location on Y axis where a tick exists.
"ytickslabel"       : list of the labels to replace yticks.
"xticks"            : list of the location on X axis where a tick exists.
"xtickslabel"       : list of the labels to replace xticks.
"fonttype"          : font used for texts in the figure.
"fontsize"          : fontsize of the texts in the figure.
"xlabel"            : label for the X axis.
"ylabel"            : label for the Y axis.
"foldername"        : folder where the figure is saved.
"plot_show"         : if True, the figure is shown after being saved.
"comment"           : comment to be added to the filename.
"plot_show"         : if True, the figure is shown after being saved.
'''
def Plot_im(data_arr, point_arr = [],
            default_param_filename = "Param_plot_image.json",
            *args, **kwargs):

    with open(default_param_filename,'r', encoding='UTF-8') as f:
        param_dict = json.load(f)

    for key,value in kwargs.items():
        param_dict[key] = value

    fig, ax = plt.subplots(figsize=(param_dict["figsize"][0],param_dict["figsize"][1]),
                           dpi = param_dict["dpi"])
    # norm = "linear"
    norm = None
    if param_dict["norm"] == "zero_in_center":
        norm = TwoSlopeNorm(vmin=np.min(data_arr), vcenter=0, vmax=np.max(data_arr))
    im  = ax.imshow(data_arr,cmap = param_dict["colormap"],
                    norm=norm,aspect=param_dict["aspect"])
    if len(point_arr) > 0:
        ax.scatter(point_arr[:,0],point_arr[:,1],
                   s=param_dict["point_size"], c=param_dict["point_color"], marker=param_dict["point_marker"])

    # shrink: 缩放比例，pad: 间距
    cbar = fig.colorbar(im, ax=ax, orientation=param_dict["cbar_orientation"],
                        shrink=param_dict["shrink"], pad=param_dict["pad"])
    if param_dict["autoset_ticks"] == 1:
        cbar.set_ticks(Auto_ticks(data_arr,param_dict["cbar_small_ticks"]))
    cbar.ax.tick_params(labelsize=param_dict["fontsize"]*0.8)
    cbar.set_label(param_dict["cbar_label"],size=param_dict["fontsize"])

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

    plt.xlabel(param_dict["xlabel"],size = param_dict["fontsize"])
    plt.ylabel(param_dict["ylabel"],size = param_dict["fontsize"])
    plt.title(param_dict["title"],size = param_dict["fontsize"]*1.5)
    plt.savefig(param_dict["foldername"]+param_dict["title"]+".pdf")
    if param_dict["plot_show"]:
        plt.show()
    plt.close()

'''
Plot_field_profile(): Plot the field profile stored in txt files

Parameters:
field                   : 2D array of the field to be plotted.
field_name              : name of the field to be plotted. E.g. "Ex","Hy".
default_param_filename:  name of the json file where the default parameters are stored.

The following default parameters are used for plotting and are stored in a json file. (default: Param_plot_curve.json)
They can be changed using a dictionary loaded as **kwargs when necessary:
"figsize"           : a tuple or list with two elements.
"dpi"               : 300 by default.
"title"             : title of the figure. Also used as a part of the filename
"component_name"    : "Ex" or "Ey" or "Ez" or "Hx" or "Hy" or "Hz".
"component_list"    : Should be the subset of ['Abs', 'Re', 'Im'].
"colormap"          : colormap of the image. 'jet' by default.
"grid_linewidth"    : linewidth of the grids. 1 by default
"norm"              : if "zero_in_center", the color map is divided into positive and negative parts with different slopes, making sure 0 is in the middle.
"X_label"           : label for the X axis.
"Y_label"           : label for the Y axis.
"xticks"            : list of the location on X axis where a tick exists.
"yticks"            : list of the location on Y axis where a tick exists.
"xtickslabel"       : list of the labels to replace xticks.
"ytickslabel"       : list of the labels to replace yticks.
"fonttype"          : font used for texts in the figure.
"fontsize"          : fontsize of the texts in the figure.
"cbar_small_ticks"  : the ticks of the color bar are set using the function "Auto_ticks". If True, small_ticks are enabled.
"cbar_label"        : label of the color bar.
"cbar_orientation"  : orientation of the color bar. 'vertical' or 'horizontal'.
"shrink"            : shrink in size of the color bar. value between 0 and 1. default: 0.8.
"pad"               : distance between the color bar and the figure. value between 0 and 1. default: 0.2.
"foldername"        : folder where the figure is saved. "./results/" by default.
'''
def Plot_field_profile(field,field_name,
                    default_param_filename = "Param_plot_field_profile.json",
                    *args, **kwargs):

    with open(default_param_filename,'r', encoding='UTF-8') as f:
        param_dict = json.load(f)

    for key,value in kwargs.items():
        param_dict[key] = value

    component_list = param_dict["component_list"]
    field_list     = []
    if "Abs" in component_list:
        field_list.append(np.abs(field))
    if "Re" in component_list:
        field_list.append(np.real(field))
    if "Im" in component_list:
        field_list.append(np.imag(field))

    fig, ax = plt.subplots(1,len(field_list),
                           figsize=param_dict["figsize"],dpi=param_dict["dpi"])
    plt.subplots_adjust(left=0.05, right=0.95, wspace =0.3, hspace =0.6)

    plt.rcParams["font.family"] = param_dict["fonttype"]
    plt.rcParams.update({'font.size': param_dict["fontsize"]})
    plt.legend()
    if len(component_list) > 1:
        for idx in range(len(component_list)):
            im = ax[idx].imshow(field_list[idx], cmap=param_dict["colormap"])
            ax[idx].set_title(component_list[idx]+'('+field_name+')')
            cbar = fig.colorbar(im, ax=ax[idx], orientation=param_dict["cbar_orientation"],
                                label=param_dict["cbar_label"],
                                shrink=param_dict["shrink"], pad=param_dict["pad"])
            cbar.ax.tick_params(labelsize=param_dict["fontsize"]*0.8)
            cbar.set_label(param_dict["cbar_label"],size=param_dict["fontsize"])
            ax[idx].set_xlabel("X",label)
            ax[idx].set_ylabel("Y")
            ax[idx].invert_yaxis()
            ax[idx].tick_params(axis='both',labelsize=param_dict["fontsize"])
            ax[idx].set_xticks(param_dict["xticks"])
            ax[idx].set_xticklabels(param_dict["xtickslabel"],fontsize=param_dict["fontsize"])
            ax[idx].set_yticks(param_dict["yticks"])
            ax[idx].set_yticklabels(param_dict["ytickslabel"],fontsize=param_dict["fontsize"])
    else:
        im = ax.imshow(field_list[0], cmap=param_dict["colormap"])
        ax.set_title(param_dict["title"]+ ": " +component_list[0]+ '('+field_name+')')
        cbar = fig.colorbar(im, ax=ax, orientation=param_dict["cbar_orientation"],
                            label=param_dict["cbar_label"],
                            shrink=param_dict["shrink"], pad=param_dict["pad"])
        cbar.ax.tick_params(labelsize=param_dict["fontsize"]*0.8)
        cbar.set_label(param_dict["cbar_label"],size=param_dict["fontsize"])
        ax.set_xlabel(param_dict["xlabel"],fontsize= param_dict["fontsize"])
        ax.set_ylabel(param_dict["ylabel"],fontsize= param_dict["fontsize"])
        ax.invert_yaxis()
        ax.tick_params(axis='both',labelsize=param_dict["fontsize"])
        ax.set_xticks(param_dict["xticks"])
        ax.set_xticklabels(param_dict["xtickslabel"],fontsize=param_dict["fontsize"])
        ax.set_yticks(param_dict["yticks"])
        ax.set_yticklabels(param_dict["ytickslabel"],fontsize=param_dict["fontsize"])

    plt.savefig(param_dict["foldername"]+param_dict["title"]+".eps",
                dpi=param_dict["dpi"])
    if param_dict["plot_show"]:
        plt.show()
    plt.close()