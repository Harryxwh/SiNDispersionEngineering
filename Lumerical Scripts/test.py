import os,sys
import numpy as np
sys.path.append("C:\\Program Files\\Lumerical\\v251\\api\\python\\") # locate lumapi files
import lumapi
import matplotlib.pyplot as plt


with lumapi.FDTD() as fdtd:

    # Set up simulation region
    fdtd.addfdtd()
    fdtd.set("x",0)
    fdtd.set("x span",8e-6)
    fdtd.set("y",0)
    fdtd.set("y span",8e-6)
    fdtd.set("z",0.25e-6)
    fdtd.set("z span",0.5e-6)

    # Set up source
    fdtd.addgaussian()
    fdtd.set("injection axis","z")
    fdtd.set("direction","forward")
    fdtd.set("x",0)
    fdtd.set("x span",16e-6)
    fdtd.set("y",0)
    fdtd.set("y span",16e-6)
    fdtd.set("z",0.2e-6)
    fdtd.set("use scalar approximation",1)
    fdtd.set("waist radius w0",2e-6)
    fdtd.set("distance from waist",0)
    fdtd.setglobalsource("wavelength start",1e-6)
    fdtd.setglobalsource("wavelength stop",1e-6)

    # Set up monitor
    fdtd.addpower()
    fdtd.set("monitor type","2D Z-normal")
    fdtd.set("x",0)
    fdtd.set("x span",16e-6)
    fdtd.set("y",0)
    fdtd.set("y span",16e-6)
    fdtd.set("z",0.3e-6)

    # Run simulation
    fdtd.save("fdtd_tutorial.fsp")
    fdtd.run()