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
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle
from astropy import units as u
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.visualization import ZScaleInterval
interval = ZScaleInterval()

import tqdm
import warnings
warnings.filterwarnings('ignore')


# Directories & Files
dir_root = "/data/jlee/DATA/Subaru/HSC/HSCv6/M81/MMT_2022A/process/"
obsID = "2022.0302"
dir_skysub = dir_root + obsID + "/reduction/0100/skysub_target1_1/"
spfile = sorted(glob.glob(dir_skysub+"*.fits"))
msfile = glob.glob(dir_root + obsID + "/reduction/0100/*.ms.fits")[0]
dir_img = "/data/jlee/DATA/Subaru/HSC/HSCv6_data2/M81/Okamoto/F1234/Phot/HSC_sephot/Images/"
imgfile = dir_img+"cbM81_F1.fits"
pixel_scale = 0.168    # arcsec/pixel


# Setting
hdr_sp = fits.getheader(spfile[0], ext=0)
hdr_ms = fits.getheader(msfile, ext=0)
img, hdr = fits.getdata(imgfile, header=True)
wc = WCS(hdr)

dir_imgp = "/data/jlee/DATA/Subaru/HSC/HSCv6/M81/Bell/Red/rerun/object/deepCoadd-results/HSC-G/0/"
img2, hdr2 = fits.getdata(dir_imgp+"6,2/calexp-HSC-G-0-6,2.fits", header=True, ext=1)
wc2 = WCS(hdr2)
img3, hdr3 = fits.getdata(dir_imgp+"7,2/calexp-HSC-G-0-7,2.fits", header=True, ext=1)
wc3 = WCS(hdr3)

knl1 = Gaussian1DKernel(20)
knl2 = Gaussian1DKernel(2*hdr_sp['CDELT1'])


# Absorption lines
abs_line = [3934.777, 3969.588, 4305.61, 4862.68, 5176.7, 5895.6, 6564.61]
### Ca H&K, G, Hbeta, Na, Mg, Halpha


# Plotting
for i in tqdm.trange(len(spfile)):
# i = 0
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

    axins = inset_axes(ax, width="100%", height="100%",
                       bbox_to_anchor=(0.025, 0.625, 0.25, 0.25),
                       bbox_transform=ax.transAxes, borderpad=0)
    axins.tick_params(left=False, right=False, labelleft=False, labelright=False,
                      top=False, bottom=False, labeltop=False, labelbottom=False)
    for axis in ['top','bottom','left','right']:
        axins.spines[axis].set_linewidth(0.8)

    ra0 = Angle(Angle(hdr_ms['APID'+str(int(fiberID))].split(' ')[1], unit=u.hour), unit=u.deg)
    dec0 = Angle(hdr_ms['APID'+str(int(fiberID))].split(' ')[2], unit=u.deg)

    x0, y0 = wc.wcs_world2pix(ra0.value, dec0.value, 1)
    x0, y0 = x0.item(0), y0.item(0)
    x2, y2 = wc2.wcs_world2pix(ra0.value, dec0.value, 1)
    x2, y2 = x2.item(0), y2.item(0)
    x3, y3 = wc3.wcs_world2pix(ra0.value, dec0.value, 1)
    x3, y3 = x3.item(0), y3.item(0)
    rth = 10.0    # arcsec

    if (img[round(y0), round(x0)] > 0.):
        cut_img =  img[round(y0-1-rth/pixel_scale):round(y0-1+rth/pixel_scale)+1,
                       round(x0-1-rth/pixel_scale):round(x0-1+rth/pixel_scale)+1]
    elif ((round(x2-1-rth/pixel_scale) >= 0) & (round(x2-1+rth/pixel_scale)+1 < img2.shape[1]) & 
          (round(y2-1-rth/pixel_scale) >= 0) & (round(y2-1+rth/pixel_scale)+1 < img2.shape[0])):
        cut_img = img2[round(y2-1-rth/pixel_scale):round(y2-1+rth/pixel_scale)+1,
                       round(x2-1-rth/pixel_scale):round(x2-1+rth/pixel_scale)+1]
    elif ((round(x3-1-rth/pixel_scale) >= 0) & (round(x3-1+rth/pixel_scale)+1 < img2.shape[1]) & 
          (round(y3-1-rth/pixel_scale) >= 0) & (round(y3-1+rth/pixel_scale)+1 < img2.shape[0])):
        cut_img = img3[round(y3-1-rth/pixel_scale):round(y3-1+rth/pixel_scale)+1,
                       round(x3-1-rth/pixel_scale):round(x3-1+rth/pixel_scale)+1]
    else:
        cut_img =  img[round(y0-1-rth/pixel_scale):round(y0-1+rth/pixel_scale)+1,
                       round(x0-1-rth/pixel_scale):round(x0-1+rth/pixel_scale)+1]

    vmin, vmax = interval.get_limits(cut_img)
    axins.imshow(cut_img, origin='lower', cmap='gray_r', vmin=vmin, vmax=vmax)
    fib = Circle((round(rth/pixel_scale), round(rth/pixel_scale)), radius=0.75/pixel_scale,
                 linewidth=1.5, edgecolor='magenta', fill=False, alpha=0.9)
    axins.add_artist(fib)

    plt.tight_layout()
    plt.savefig("./Coadd_"+fiberID+"."+objName+".png", dpi=300)
    plt.close()
