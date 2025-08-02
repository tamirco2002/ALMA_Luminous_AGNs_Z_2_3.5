from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
from astropy.io import fits
from scipy import integrate
from scipy import interpolate
from scipy.optimize import curve_fit

import Extensions as ext
from src import ALMA_HERSCHEL_INPUT_PATH, CE01_PATH, LANI_TORUS_PATH


# Figure 1 - Separations of quasars and ALMA-detected hosts
def quasars_host_seperation(output_path,sep = False, dates = False):
    db = pd.read_excel(ALMA_HERSCHEL_INPUT_PATH)
    df_CE = pd.read_excel(CE01_PATH)
    df_torus = pd.read_csv(LANI_TORUS_PATH)
    indicators=db['ind'].to_numpy()

    L_AGN = db['log L_AGN'].to_numpy()
    L_SF=db['log L_SF'].to_numpy()
    group=db['group'].to_numpy()

    z = db['z [SDSS/DR16]']  # redshift
    fits_filename = db['Filename']
    Galaxy_name = db['ALMA name']

    RA_SDSS = [str(db.iloc[i]['RA [SDSS/DR16]']) for i in range (len(db['RA [SDSS/DR16]']))]
    Dec_SDSS = [str(db.iloc[i]['Dec [SDSS/DR16]'])for i in range (len(db['Dec [SDSS/DR16]']))]

    RA_CASA = [str(db.iloc[i]['RA [CASA]']) for i in range (len(db['RA [SDSS/DR16]']))]
    Dec_CASA = [str(db.iloc[i]['Dec [CASA]']) for i in range(len(db['Dec [CASA]']))]

    Fits_location = "../data/Fits"
    hdul = [fits.open(Fits_location+"/"+fits_filename[i]) for i in range(len(fits_filename))]
    data = [hdul[i][0].data[0, 0] for i in range(len(hdul))]
    hdr = [hdul[i][0].header for i in range(len(hdul))]

    pixel_to_arc = [ext.deg_to_arc(np.abs(hdr[i]['CDELT1'])) for i in range(len(hdr))]  # 1 pixel = "pixel_to_arc" arc-seconds
    pixel_to_deg = [np.abs(hdr[i]['CDELT1']) for i in range(len(hdr))]  # 1 pixel = "pixel_to_deg" degrees

    coords_sdss = [SkyCoord(RA_SDSS[i],Dec_SDSS[i],unit=(u.hourangle, u.deg)) for i in range(len(RA_SDSS))]
    coords_casa = [SkyCoord(RA_CASA[i],Dec_CASA[i],unit=(u.hourangle, u.deg)) for i in range(len(RA_CASA))]

    x = coords_sdss[0].separation(coords_casa[0])
    seperations = [coords_sdss[i].separation(coords_casa[i]) for i in range(len(coords_sdss))]

    seperations_arc = [ext.deg_to_arc(seperations[i].value) for i in range(len(seperations))]
    seperations_pix = [seperations[i].value / pixel_to_deg[i] for i in range(len(seperations))]

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    arc_kpc = [1 / cosmo.arcsec_per_kpc_proper(z[i]).value for i in range(len(z))]  # 1 arc-sec = "arc_kpc" kpc
    seperations_kpc = [seperations_arc[i]*arc_kpc[i] for i in range(len(seperations_arc))]
    if sep:
        outputDF = pd.DataFrame()
        outputDF["ind"]=db['ind']
        outputDF["group"] = db['group']
        outputDF["ALMA name"] = db['ALMA name']
        outputDF["RA [CASA]"] = db['RA [CASA]']
        outputDF["Dec [CASA]"] = db['Dec [CASA]']
        outputDF["RA [SDSS/DR16]"] = db['RA [SDSS/DR16]']
        outputDF["Dec [SDSS/DR16]"] = db['Dec [SDSS/DR16]']
        outputDF["z [SDSS/DR16]"] = db['z [SDSS/DR16]']
        outputDF["Detected by ALMA"] = db ["Detected by ALMA"]
        outputDF["Seperation [arc-sec]"]=seperations_arc
        outputDF["Seperation [kpc]"] = seperations_kpc
        outputDF["Seperation [pixel]"] = seperations_pix
        outputDF.to_excel(output_path+"\SeperationsFile_"+datetime.now().strftime("%m%d%Y-%H%M%S")+".xlsx")

        df_det = outputDF.loc[outputDF["Detected by ALMA"] == 'det']
        df_undet = outputDF.loc[outputDF["Detected by ALMA"] == 'non-det']

        mean = np.mean(np.array(df_det['Seperation [arc-sec]']))
        std = np.std(np.array(df_det['Seperation [arc-sec]']))
        mean_undet = np.mean(np.array(df_undet['Seperation [arc-sec]']))
        std_undet = np.std(np.array(df_undet['Seperation [arc-sec]']))
        print(np.min(seperations_arc),np.max(seperations_arc))
        bins = np.histogram(seperations_arc, bins=6,range=(0,1.2))[1]
        plt.hist(df_det['Seperation [arc-sec]'], bins=bins,color='green',label="Detected")
        plt.hist(df_undet['Seperation [arc-sec]'], bins=bins,alpha=0.5,color='red',label="Un-detected")

        plt.legend()
        plt.xlabel("Distance between quasar center and galaxy center [arc-sec]")
        plt.ylabel("Number of Sources")
        plt.yticks(np.arange(0, 20, step=2))

        plt.tight_layout()
        plt.savefig(f'{output_path}\Seperation_ArcSec.pdf',format='pdf')
        plt.clf()
        plt.close()

        mean = np.mean(np.array(df_det['Seperation [kpc]']))
        std = np.std(np.array(df_det['Seperation [kpc]']))
        mean_undet = np.mean(np.array(df_undet['Seperation [kpc]']))
        std_undet = np.std(np.array(df_undet['Seperation [kpc]']))

        bins = np.histogram(seperations_kpc,bins=5,range=(0,10))[1]
        plt.hist(df_det['Seperation [kpc]'], bins=bins,color='green',label="Detected")
        plt.hist(df_undet['Seperation [kpc]'], bins=bins,alpha=0.5,color='red',label="Un-detected")

        plt.xlabel("Distance between quasar center and galaxy center [kpc]")
        plt.ylabel("Number of Sources")
        plt.yticks(np.arange(0, 20, step=2))

        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{output_path}\Seperation_KPC.pdf',format='pdf')
        RA_SDSS0 = str(db.iloc[0]['RA [SDSS/DR16]']).split(sep=":")
        Dec_SDSS0 = str(db.iloc[0]['Dec [SDSS/DR16]']).split(sep=":")
        RA_SDSS0 = np.array(RA_SDSS0, dtype=float)
        Dec_SDSS0 = np.array(Dec_SDSS0, dtype=float)
        alpha_SDSS, delta_SDSS = ext.alpha_delta(RA_SDSS0, Dec_SDSS0, 1)

    if dates:
        dateList = [hdr[i]['DATE-OBS'] for i in range(len(hdr))]
        outputDF = pd.DataFrame()
        outputDF["ind"]=db['ind']
        outputDF["group"] = db['group']
        outputDF["ALMA name"] = db['ALMA name']
        outputDF["Observation date"] = dateList
        outputDF.to_excel(f"{output_path}\DatesFile_" + datetime.now().strftime("%m%d%Y-%H%M%S") + ".xlsx")


if __name__ == '__main__':
    out_path = 'Replace with your output path'
    quasars_host_seperation(out_path,sep=True,dates=True)