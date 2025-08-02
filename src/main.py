import gc
import math
import os
from collections import OrderedDict
from datetime import datetime

import Extensions as ext
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Functions import *
from astropy import units as u, constants as const
from astropy.cosmology import LambdaCDM
from astropy.io import fits
from matplotlib.patches import Ellipse
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredAuxTransformBox
from scipy import integrate
from scipy.optimize import curve_fit

from .Source import Source

data_location = "../data"
Fits_location = "../Fits files"
results_location = "../Results"
Sources_location = results_location+"/Sources"
Images_location = results_location+"/Images"
SEDs_location = results_location+"/SEDs"
Excels_locations = results_location+"/Excels"
analyze_beam_size = False

output_path = ""

def prepare_analysis():
    wl_range = [8, 1000]  # observed wavelengths range in microns
    # Converting the wl range (8-1000 micron) to frequency Range (GHz)
    freq_range = [(ext.wl_to_freq(wl_range[0], "GHz").value), (ext.wl_to_freq(wl_range[1], "GHz").value)]
    # Creating a vector of frequencies in the frequency range
    freq = np.arange(freq_range[1], freq_range[0], step=0.05) * u.GHz
    # Converting the frequency vector to a wavelgnth vector
    wavelengths = ext.freq_to_wl(freq, "GHz").value

    # CE templates
    df_CE = pd.read_excel(data_location+"/Chary_corrected_table_9feb17.xlsx")
    wl_CE = df_CE['Lambda'].to_numpy()[1:]  # CE templates wavelengths
    freq_CE = (const.c * u.s / u.m / np.flip(wl_CE) / 1000)  # converting CE wls to frequencies
    return freq, df_CE,wl_CE,freq_CE

class Plots_Data:
    def __init__(self,flux = None, freq = None, wavelengths = None, wl_CE = None, chary_nuLnu_arr = None, torus = None,
                 torus_wls = None,gb_bf = None,gb_ff = None):
        self.flux = flux
        self.freq = freq
        self.wavelengths = wavelengths
        self.wl_CE = wl_CE
        self.chary_nuLnu_arr = chary_nuLnu_arr
        self.torus = torus
        self.torus_wls = torus_wls
        self.torus_wls = torus_wls
        self.gb_bf = gb_ff



def analysis(sources,saveFile,create_plots,showPlots,create_images,GB_type='BB'):

    freq, df_CE, wl_CE, freq_CE = prepare_analysis()
    
    for source in sources:
        print("Starting source #",source.ind, " (",source.name,")")
        flux = ext.scaled_graybody(freq, source.beta, source.temp*u.K, source.rest_freq, source.UsedFlux)  # GB flux array
        torus, torus_wls = source.getTorus(flux, wavelengths)

        source.ALMA_flux = source.UsedFlux
        source.UsedFlux = source.ALMA_flux*(1-source.torus_fraction)
        source.subtractTorHerschel(torus,torus_wls)

        # FF = free beta and temperature, BF = bound beta and free temperature
        if GB_type in ['FF','BF'] and source.group=='det':
            source.findTBeta()
            source.temp = source.T_ff if GB_type == 'FF' else source.T_bf
            source.beta = source.beta_ff if GB_type == 'FF' else source.beta_bf


        flux = ext.scaled_graybody(freq, source.beta, source.temp*u.K, source.rest_freq, source.UsedFlux)  # re-creating the GB

        source.nuLnu = ext.Flux_to_nuLnu(source.UsedFlux, source.ob_freq, source.dist_cm)  # .value
        source.nuLnu_err = ext.Flux_to_nuLnu(source.Flux_Err, source.ob_freq, source.dist_cm)

        source.ALMA_nuLnu = ext.Flux_to_nuLnu(source.ALMA_flux, source.ob_freq, source.dist_cm)

        closest_wav = np.argmin(np.abs(source.rest_wl - wl_CE))  # finding the right wavelength row
        nuLnu_row_CE = df_CE.iloc[closest_wav + 1, 1:]  # row of nu*L_nu values at selected wavelength
        source.CE_template_num = np.argmin(np.abs(np.log10(source.nuLnu) - nuLnu_row_CE))  # finding the closest value of nu*L_nu
        source.CE_nuLnu = df_CE.iloc[closest_wav + 1, source.CE_template_num + 1]
        chary_nuLnu_arr = df_CE.iloc[1:, source.CE_template_num + 1]  # the distribution which fits the best

        Wise_wls_ob = np.array([3.4, 4.6, 12, 22, 250, 350, 500])
        Wise_wls_rest = ext.wl_obs_to_rest(Wise_wls_ob, source.z)



        if create_plots or create_images:
            path = Sources_location + "/" + source.name
            try:
                os.mkdir(path)
            except FileExistsError:
                pass

        # Creating the SEDs
        if create_plots:
            plots_data_obj = Plots_Data(flux,freq,wavelengths,wl_CE,chary_nuLnu_arr,torus,torus_wls)
            source.findTBeta()
            gb_bf,gb_ff = source.plot_NuLnu_wl_Lsun(flux,freq,wavelengths,wl_CE,chary_nuLnu_arr,torus,torus_wls,showPlots=showPlots)
            plots_data_obj.gb_bf, plots_data_obj.gb_ff = gb_bf,gb_ff
            source.plot_NuLnu_wl_cgs(flux, freq, wavelengths, wl_CE, chary_nuLnu_arr, torus, torus_wls,gb_bf,gb_ff,showPlots=showPlots)
            source.plot_Fnu_freq(flux, freq, wavelengths, wl_CE, chary_nuLnu_arr, torus, torus_wls,gb_bf,gb_ff,showPlots=showPlots)

        source.L_CE = df_CE.iloc[0, source.CE_template_num + 1]  # CE Luminosity (first row in table)

        # integration of total flux:
        source.total_flux = integrate.simps(y=flux / (1 + source.z), x=freq * (10 ** 9))  # integration over obsereved frequencies
        source.total_lum = ext.Flux_to_Lum(source.total_flux, source.dist_cm)  # luminosity in erg/s (L=F*4*pi*(D_l)^2)
        lum_solar = ext.L_to_L_sun(source.total_lum)
        source.SFR_gb = (lum_solar) / (10 ** 10)  # SFR=L/(l_sun*10^10)
        source.SFR_CE = ((10 ** source.L_CE) * u.erg / u.s) / (const.L_sun.to(u.erg / u.s) * (10 ** 10))
        source.gb_nuLnu = ext.Flux_to_nuLnu(ext.scaled_graybody(source.rest_freq * u.GHz, source.beta, source.temp*u.K, source.rest_freq, source.UsedFlux), source.ob_freq,source.dist_cm)

        if create_images:
            source.zoomedInImage()
            source.fullViewImage()

        # calculating the Mass consuming rate of the BH
        source.M_dot_BH = ext.L_AGN_to_M_dot_BH(source.L_AGN)
        source.MBH = ext.L_AGN_to_MBH(source.L_AGN)

        if source.neigh:
            source.companion.find_SFRs(freq,source.ob_freq,source.temp,source.beta,source.dist_cm,df_CE,nuLnu_row_CE,closest_wav)

        if saveFile:
            df_temp = source.saveToFile()
            if source == sources[0]:
                outF = df_temp
            else:
                outF = pd.concat([outF, df_temp])

            if source == sources[-1]:
                File_name = f'{output_path}/Analysis_Results_{datetime.now().strftime("%m%d%Y-%H%M%S")}.xlsx'
                outF.to_excel(File_name, index=False)

        print("Source #", source.ind, " is finished!")


def getPeakInJy(file,sourceList):
    df = pd.read_excel(file)
    arrP = []
    arrErr = []
    for s in range(len(sourceList)):
        arrP.append(df["Peak Value (Jy/beam)"][s]/sourceList[s].beam_pix)
        arrErr.append(df["Peak Error (Jy/beam)"][s]/sourceList[s].beam_pix)
    dfNew = pd.DataFrame({"Peak [Jy]":arrP,"Peak Error [Jy]":arrErr})
    dfNew.to_excel("PeakInJy.xlsx")
    return


def prepare_sources(indices,temp,beta,gb_type,saveFile,create_plots,showPlots,create_images):
    db = pd.read_excel(data_location+"/ALMA x Herschel CASA input - March23.xlsx")
    source_list = []
    for i in indices:
        source_list.append(Source(db,i-1,temp,beta))

    analysis(source_list,saveFile=saveFile,create_plots=create_plots,showPlots=showPlots,create_images=create_images,GB_type=gb_type)


