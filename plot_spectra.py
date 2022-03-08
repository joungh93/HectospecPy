#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 6 11:18:20 2022

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os, copy
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
import tqdm


# Directories
dir_root = "/data/jlee/HSCv6/M81/MMT_2022A/process/"
obsID = "2022.0302"
dir_skysub = dir_root + obsID + "/reduction/0100/skysub_target1_1/"
spfile = sorted(glob.glob(dir_skysub+"*.fits"))


# Setting
d, h = fits.getdata(spfile[0], header=True)

knl1 = Gaussian1DKernel(20)
knl2 = Gaussian1DKernel(2*h['CDELT1'])


# Absorption lines
abs_line = [3934.777, 3969.588, 4305.61, 4862.68, 5176.7, 5895.6, 6564.61]
### Ca H&K, G, Hbeta, Na, Mg, Halpha


# Plotting
for i in tqdm.trange(len(spfile)):
    d, h = fits.getdata(spfile[i], header=True)
    wav = h['CRVAL1'] + h['CDELT1']*(np.arange(h['NAXIS1'])+(1-h['CRPIX1']))
    coflx, sqvar = d[0, 0, :], d[3, 0, :]
    fiberID = spfile[i].split('/')[-1].split('.')[0]
    objName = spfile[i].split('/')[-1].split('.')[1]

    xmin, xmax = 3700., 7000.
    if (np.median(coflx[(wav > 6000.) & (wav < 6400.)]) <= 0.):
        med, sig = np.median(coflx[(wav > 6000.) & (wav < 6400.)]), np.std(coflx[(wav > 6000.) & (wav < 6400.)])
        coflx += -med+0.5*sig
        ysc = False
    else:
        ysc = True
    fconv = convolve(coflx, knl1)
    fconv2 = convolve(coflx, knl2)
    nflx = np.median(fconv[(wav > 6000.) & (wav < 6400.)])
    ymin, ymax = -0.09, 1.4*fconv2[(wav > xmin) & (wav < xmax)].max() / nflx
    imin, imax = np.argmin(np.abs(wav-xmin)), np.argmin(np.abs(wav-xmax))
    SN_median = np.median(coflx[imin:imax+1] / sqvar[imin:imax+1])    

    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(wav, fconv2/nflx, '-', color='C0', linewidth=1.5, alpha=0.7)
    for w in abs_line:
        idx_line = np.argmin(np.abs(wav-w))
        ax.plot([w, w],
                [np.median(fconv[idx_line-10:idx_line+10])/nflx+(ymax-ymin)*0.1,
                 np.median(fconv[idx_line-10:idx_line+10])/nflx+(ymax-ymin)*0.2],
                'r-', linewidth=2.5, alpha=0.9)
    ax.set_xlim([xmin, xmax])
    if ysc:
        ax.set_ylim([ymin, ymax])
    ax.set_xlabel(r"Wavelength [${\rm \AA}$] (Observer-frame)", fontsize=15.0)
    ax.set_ylabel("Normalized flux", fontsize=15.0)
    ax.tick_params(axis="both", labelsize=15.0)
    ax.tick_params(width=1.25, length=7.0)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.25)
    ax.text(0.04, 0.95, fiberID+". "+objName, fontsize=16.0, fontweight='bold', color='k',
            ha='left', va='top', transform=ax.transAxes)
    ax.text(0.95, 0.06, f"S/N = {SN_median:.2f}", fontsize=14.0, color='k',
            ha='right', va='bottom', transform=ax.transAxes)
    plt.tight_layout()
    
    plt.savefig("./Coadd_"+fiberID+"."+objName+".png", dpi=300)
#     plt.close()
