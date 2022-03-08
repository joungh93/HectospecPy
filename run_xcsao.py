#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 8 12:32:33 2022
@author: jlee
"""


import numpy as np
from matplotlib import pyplot as plt
import glob, os
from astropy.io import fits


# ----- Directories ----- #
current_dir = os.getcwd()
dir_iraf = "/data/jlee/HSCv6/M81/MMT_2022A/spectra/velocity/"
os.chdir(dir_iraf)


# ----- Importing IRAF tasks ----- #
from pyraf import iraf
from pyraf.iraf import rvsao

iraf.unlearn()
os.chdir(current_dir)
iraf.chdir(current_dir)

try:
    iraf.xcsao
except AttributeError:
    print("Please check again.")


# ----- Listing the spectra files ----- #
dir_root = "/data/jlee/HSCv6/M81/MMT_2022A/process/"
obsID = "2022.0302"
dir_skysub = dir_root + obsID + "/reduction/0100/skysub_target1_1/"
spfile = sorted(glob.glob(dir_skysub+"*.fits"))


# ----- Creating the lists ----- #
if glob.glob("binning") == []:
    os.system("mkdir binning")
else:
    os.system("rm -rfv binning/*")

for i in np.arange(len(spfile)):
    os.system("cp -rpv "+spfile[i]+" binning/"+spfile[i].split('/')[-1].split('.fits')[0]+".bin1.fits")
    
f = open("bin1.list", "w")
for i in np.arange(len(spfile)):
    f.write("binning/"+spfile[i].split('/')[-1].split('.fits')[0]+'.bin1.fits'+"\n")
f.close()

f = open("bin2.list", "w")
for i in np.arange(len(spfile)):
    f.write("binning/"+spfile[i].split('/')[-1].split('.fits')[0]+'.bin2.fits'+"\n")
f.close()

f = open("bin4.list", "w")
for i in np.arange(len(spfile)):
    f.write("binning/"+spfile[i].split('/')[-1].split('.fits')[0]+'.bin4.fits'+"\n")
f.close()


# ----- Templates ----- #
dir_temp = "/home/jlee/anaconda3/envs/geminiconda/iraf_extern/rvsao/templates/"
temp_GC = sorted(glob.glob(dir_temp+"m31*.fits"))
temp_Astar = sorted(glob.glob(dir_temp+"astar*.fits"))
temp_QSO = sorted(glob.glob(dir_temp+"sdss_qso.fits"))
temp_El = sorted(glob.glob(dir_temp+"e*temp.fits"))
temp_Sp = sorted(glob.glob(dir_temp+"sptemp.fits"))
temp_gal = sorted(glob.glob(dir_temp+"h*temp0.fits"))
temps = temp_GC + temp_Astar + temp_QSO + temp_El + temp_Sp + temp_gal

if glob.glob("template_copy") == []:
    os.system("mkdir template_copy")
else:
    os.system("rm -rfv template_copy/*")
    
for i in np.arange(len(temps)):
    os.system("cp -rpv "+temps[i]+" template_copy/")
g = open("templates.list", "w")
for i in np.arange(len(temps)):
    g.write("template_copy/"+temps[i].split("/")[-1]+"\n")
g.close()


# ----- Binning with IRAF/blkavg tasks ----- #
os.system("rm -rfv binning/*.bin2.fits")
iraf.blkavg(input="@bin1.list", output="@bin2.list", b1=2, b2=1, b3=1)

os.system("rm -rfv binning/*.bin4.fits")
iraf.blkavg(input="@bin1.list", output="@bin4.list", b1=4, b2=1, b3=1)


# ----- Running IRAF/xcsao tasks ----- #
iraf.xcsao(spectra="@bin1.list", templates="@templates.list",
	       st_lambda=3800., end_lambda=5400.,
           svel_corr="heliocentric", save_vel="yes")
iraf.xcsao(spectra="@bin2.list", templates="@templates.list",
	       st_lambda=3800., end_lambda=5400.,
           svel_corr="heliocentric", save_vel="yes")
iraf.xcsao(spectra="@bin4.list", templates="@templates.list",
	       st_lambda=3800., end_lambda=5400.,
           svel_corr="heliocentric", save_vel="yes")

