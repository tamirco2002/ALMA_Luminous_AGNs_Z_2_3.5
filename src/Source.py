

class Source:
    def __init__(self, db, index,temp,beta):
        self.temp = temp
        self.beta = beta

        self.ind = db.iloc[index]['ind']
        self.group = db.iloc[index]['group']
        self.name = db.iloc[index]['ALMA name']
        self.full_name = db.iloc[index]['Name [Netzer+16]']
        self.file = db.iloc[index]['Filename']

        self.RA_CASA = np.array(str(db.iloc[index]['RA [CASA]']).split(sep=":"), dtype=float)
        self.Dec_CASA = np.array(str(db.iloc[index]['Dec [CASA]']).split(sep=":"), dtype=float)
        self.RA_SDSS = np.array(str(db.iloc[index]['RA [SDSS/DR16]']).split(sep=":"), dtype=float)
        self.Dec_SDSS = np.array(str(db.iloc[index]['Dec [SDSS/DR16]']).split(sep=":"), dtype=float)

        # Converting the SDSS and CASA RA/Dec to degrees (alpha/delta)
        self.sign = 1
        if (str(db.iloc[index]['Dec [SDSS/DR16]']).split(sep=":")[0][0] == '-'):
            self.sign = -1
        self.alpha_SDSS, self.delta_SDSS = ext.alpha_delta(self.RA_SDSS, self.Dec_SDSS, self.sign)

        self.sign_CASA = 1
        if (str(db.iloc[index]['Dec [CASA]']).split(sep=":")[0][0] == '-'):
            self.sign_CASA = -1
        self.alpha_CASA, self.delta_CASA = ext.alpha_delta(self.RA_CASA, self.Dec_CASA, self.sign_CASA)

        self.Flux = ext.Jansky_to_erg_s(db.iloc[index]['Flux [Jy, CASA]']).value
        self.Flux_Err = ext.Jansky_to_erg_s(db.iloc[index]['Flux Error [Jy, CASA]']).value
        self.Flux_RMS = ext.Jansky_to_erg_s(db.iloc[index]['Flux RMS [Jy, CASA]']).value
        self.SNR = self.Flux / self.Flux_RMS

        self.z = db.iloc[index]['z [SDSS/DR16]']
        self.ALMA_det = db.iloc[index]['Detected by ALMA']

        if np.abs(self.SNR)<3 or self.ALMA_det == 'non-det':
            self.UsedFlux = np.abs(3*self.Flux_RMS)
        else:
            self.UsedFlux = self.Flux

        self.F_3mic = db.iloc[index]['F 3.4 [mJy]'] * 10 ** -26
        self.F_3mic_Err = db.iloc[index]['dF 3.4'] * 10 ** -26

        self.F_4mic = db.iloc[index]['F 4.6 [mJy]'] * 10 ** -26
        self.F_4mic_Err = db.iloc[index]['dF 4.6'] * 10 ** -26

        self.F_12mic = db.iloc[index]['F 12 [mJy]'] * 10 ** -26
        self.F_12mic_Err = db.iloc[index]['dF 12'] * 10 ** -26

        self.F_22mic = db.iloc[index]['F 22 [mJy]'] * 10 ** -26
        self.F_22mic_Err = db.iloc[index]['dF 22'] * 10 ** -26

        self.F_250mic = db.iloc[index]['F 250 [mJy]'] * 10 ** -26
        self.F_250mic_Err = db.iloc[index]['dF 250'] * 10 ** -26

        self.F_350mic = db.iloc[index]['F 350 [mJy]'] * 10 ** -26
        self.F_350mic_Err = db.iloc[index]['dF 350'] * 10 ** -26

        self.F_500mic = db.iloc[index]['F 500 [mJy]'] * 10 ** -26
        self.F_500mic_Err = db.iloc[index]['dF 500'] * 10 ** -26

        if db.iloc[index]['Neighbour'] == 'det':
            self.neigh = True
            self.companion = Companion(
                F=ext.Jansky_to_erg_s(db.iloc[index]['Flux Neighbour [Jy, CASA]']).value,
                F_err=ext.Jansky_to_erg_s(db.iloc[index]['Flux  Neighbour Error [Jy, CASA]']).value,
                RMS=ext.Jansky_to_erg_s(db.iloc[index]['Flux  Neighbour RMS [Jy, CASA]']).value,
                RA=db.iloc[index]['RA Neighbour [CASA]'],
                Dec=db.iloc[index]['Dec Neighbour [CASA]'],
                main_z=self.z,
            )

            self.neigh_ratio = self.Flux / (self.Flux + self.companion.F)

            self.F_250mic *= self.neigh_ratio
            self.F_250mic_Err *= self.neigh_ratio
            self.F_350mic *= self.neigh_ratio
            self.F_350mic_Err *= self.neigh_ratio
            self.F_500mic *= self.neigh_ratio
            self.F_500mic_Err *= self.neigh_ratio
        else:
            self.neigh = False

        self.L_SF = db.iloc[index]['log L_SF']
        self.L_AGN = db.iloc[index]['log L_AGN']

        self.T_bf = db.iloc[index]['T1']
        self.beta_bf = db.iloc[index]['beta1']

        self.T_ff = db.iloc[index]['T2']
        self.beta_ff = db.iloc[index]['beta2']

        self.get_attributes_from_fits()
        self.get_cosmo_attributes()

    def get_attributes_from_fits(self):
        hdul = fits.open(Fits_location + "/" + self.file)  # opening fits file
        self.data = hdul[0].data[0, 0]  # extracting the data from the file
        hdr = hdul[0].header  # fits header

        self.ob_freq = hdr['RESTFRQ'] / (10 ** 9)  # Observed Frequency
        self.rest_freq = ext.freq_obs_to_rest(self.ob_freq, self.z)  # Rest-Frame frequency (before redshift)
        self.rest_wl = ext.freq_to_wl(self.rest_freq, "GHz").value

        # RA/Dec of observation center location in degrees
        self.alpha_mes, self.delta_mes = hdr['CRVAL1'], hdr['CRVAL2']

        self.center_pix = hdr['CRPIX1'], hdr['CRPIX2']  # Pixel indices of the observation's center pixel
        self.pixel_to_arc = ext.deg_to_arc(np.abs(hdr['CDELT1']))  # 1 pixel = "pixel_to_arc" arc-seconds
        self.pixel_to_deg = np.abs(hdr['CDELT1'])  # 1 pixel = "pixel_to_deg" degrees
        self.CASA_pix = list(reversed(
            ext.equatorial_to_pixel(self.center_pix, self.alpha_CASA, self.delta_CASA, self.alpha_mes, self.delta_mes,
                                    self.pixel_to_deg)[-1]))
        self.Quasar_pix = list(reversed(ext.equatorial_to_pixel(self.center_pix, self.alpha_SDSS, self.delta_SDSS, self.alpha_mes,
                                                               self.delta_mes, self.pixel_to_deg)[-1]))
        # beam size
        BMAJ = hdr['BMAJ']  # beam major-axis (deg)
        BMIN = hdr['BMIN']  # beam minor-axis (deg)
        BMAJ_arc = ext.deg_to_arc(BMAJ)  # major-axis (arc-sec)
        BMIN_arc = ext.deg_to_arc(BMIN)  # minor-axis (arc-sec)

        self.BMAJ_arc = BMAJ_arc
        self.BMIN_arc = BMIN_arc

        self.BMAJ_pix = BMAJ_arc / self.pixel_to_arc  # major-axis (pixels)
        self.BMIN_pix = BMIN_arc / self.pixel_to_arc  # minor-axis (pixels)
        self.BPA = hdr['BPA']  # beam angle (deg)

        ell_area = lambda a, b: np.pi * a * b  # calculate ellipse area
        beam_area = ell_area(BMAJ_arc / 2, BMIN_arc / 2)  # arc-sec^2
        pixel_area = self.pixel_to_arc ** 2  # arc-sec^2
        self.beam_pix = beam_area / pixel_area  # 1 Jy/pixel = 'beam_pix' Jy/beam

    def get_cosmo_attributes(self):
        # Distance of source from Earth and kpc to arc/pixel convertion rates
        cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)  # Cosmological Model
        self.distance = cosmo.luminosity_distance(self.z)  # Distance from source in Mpc
        self.dist_cm = self.distance.to(u.cm).value  # Distance from source in cm
        self.arc_kpc = 1 / cosmo.arcsec_per_kpc_proper(self.z).value  # 1 arc-sec = "arc_kpc" kpc
        self.pix_kpc = self.pixel_to_arc * self.arc_kpc  # 1 pixel = "pix_kpc" kpc

    def Xi_square(self, T_range, beta_range):

        data_freq_rest = []
        data_flux = []
        data_err = []

        if self.F_250mic != 0  and not math.isnan(self.F_250mic):
            data_freq_rest.append(ext.freq_obs_to_rest((ext.wl_to_freq(250, "GHz").value), self.z))
            data_flux.append(self.F_250mic)
            data_err.append(self.F_250mic_Err)

        if self.F_350mic != 0 and not math.isnan(self.F_350mic):
            data_freq_rest.append(ext.freq_obs_to_rest((ext.wl_to_freq(350, "GHz").value), self.z))
            data_flux.append(self.F_350mic)
            data_err.append(self.F_350mic_Err)

        if self.F_500mic != 0 and not math.isnan(self.F_500mic):
            data_freq_rest.append(ext.freq_obs_to_rest((ext.wl_to_freq(500, "GHz").value), self.z))
            data_flux.append(self.F_500mic)
            data_err.append(self.F_500mic_Err)

        if len(data_freq_rest)>0:
            T_min, beta_min, xi_min = ext.Xi_Square(beta_range, T_range, self.rest_freq, self.UsedFlux * (u.erg / u.s),
                                                    data_freq_rest,
                                                    data_flux, data_err, error=self.Flux_Err)
        else:
            T_min, beta_min, xi_min = 0,0,0
        return T_min, beta_min, xi_min

    def getWiseHerschValues(self,wls):

        wls_flux_dict = {
            3.4:(self.F_3mic,self.F_3mic_Err),
            4.6: (self.F_4mic,self.F_4mic_Err),
            12: (self.F_12mic,self.F_12mic_Err),
            22: (self.F_22mic,self.F_22mic_Err),
            250: (self.F_250mic,self.F_250mic_Err),
            350: (self.F_350mic,self.F_350mic_Err),
            500: (self.F_500mic,self.F_500mic_Err),

        }
        data_wl_obs = []
        data_flux = []
        data_err = []
        for w in wls:
            if w in wls_flux_dict.keys():
                if wls_flux_dict[w][0] != 0:
                    data_wl_obs.append(w)
                    data_flux.append(wls_flux_dict[w][0])
                    data_err.append(wls_flux_dict[w][1])

        return data_wl_obs, data_flux, data_err

    def getTorus(self, flux,gb_wls):
        # Torus template
        df_torus = pd.read_csv(data_location+'/SED_table_Lani_2017 - corrected.csv')
        # Stitching the wanted torus template with a different template in the missing edges
        # tor = np.array(list(df_torus["1"][:19]) + list(df_torus["4"][19:-8]))
        # torus_wls = df_torus["wl"][:-8]

        tor = np.array(list(df_torus["1"]))
        torus_wls = df_torus["wl"]
        # Adjusting the torus to go through the 22 micron wise data point
        if self.F_22mic > 0:
            # Convertin Wise Flux to nuLnu
            nuLnu_22 = (ext.Flux_to_nuLnu(self.F_22mic, ext.wl_to_freq(22,"GHz"), self.dist_cm) / const.L_sun.to(
                u.erg / u.s)).value
            nuLnu_err_22 = (ext.Flux_to_nuLnu(self.F_22mic_Err, ext.wl_to_freq(22,"GHz"), self.dist_cm) / const.L_sun.to(
                u.erg / u.s)).value

            torus = ext.fit_torus(nuLnu_22, nuLnu_err_22, ext.wl_obs_to_rest(22,self.z), tor, torus_wls)

        # If 22 micron is missing - adjusting the torus to go through the 12 micron wise data point
        else:
            nuLnu_12 = ext.Flux_to_nuLnu(self.F_12mic, ext.wl_to_freq(12,"GHz"), self.dist_cm) / const.L_sun.to(
                u.erg / u.s).value
            nuLnu_err_12 = ext.Flux_to_nuLnu(self.F_12mic_Err, ext.wl_to_freq(12,"GHz"), self.dist_cm) / const.L_sun.to(
                u.erg / u.s).value

            torus = ext.fit_torus(nuLnu_12, nuLnu_err_12, ext.wl_obs_to_rest(12,self.z), tor, torus_wls)
        torus = np.array(torus)
        hersch_wl_obs, hersch_flux, hersch_f_err = self.getWiseHerschValues([250, 350, 500])
        hersch_nuLnu = (ext.Flux_to_nuLnu(np.array(hersch_flux), ext.wl_to_freq(np.array(hersch_wl_obs),"GHz"),self.dist_cm) / const.L_sun.to(u.erg / u.s)).value

        if self.group == 'det':
            # finding the torus values at the closest rest-frame wavelegnth to Herschel measurements
            indices_wls = [np.argmin(np.abs(ext.wl_obs_to_rest(hersch_wl_obs,self.z)[i] -torus_wls)) for i in range(len(hersch_wl_obs))]
            nuLnu_tor_herschel = np.array([torus[indices_wls[i]] for i in range(len(indices_wls))]) # nuLnu values for torus at herschel wls
            # Subtracting the torus nuLnu's from Herschel measurements
            nuLnu_hersch_minus_tor = np.array([hersch_nuLnu[i] -nuLnu_tor_herschel[i] for i in range(len(indices_wls))])

            # Altering the source saved attributes to not include the torus
            for w in range(len(hersch_wl_obs)):
                f = ext.nuLnu_to_Flux(ext.L_sun_to_L(nuLnu_hersch_minus_tor[w]),
                                      ext.wl_to_freq(hersch_wl_obs[w], freq_units="GHz"), self.dist_cm).value
                if hersch_wl_obs[w] == 250:
                    self.F_250mic = f
                elif hersch_wl_obs[w] == 350:
                    self.F_350mic = f
                elif hersch_wl_obs[w] == 500:
                    self.F_500mic = f


        indices_for_tor = [np.argmin(np.abs(w - gb_wls)) for w in torus_wls]
        fracs = ext.L_sun_to_L(torus) / np.array(ext.Flux_to_nuLnu(flux, ext.wl_to_freq(ext.wl_rest_to_obs(gb_wls, self.z),"GHz"),
                              self.dist_cm).value)[indices_for_tor]
        popt, pcov = curve_fit(ext.powerLaw, torus_wls[111:], fracs[111:])
        self.torus_fraction = ext.powerLaw(self.rest_wl, *popt)
        return torus,torus_wls

    def subtractTorHerschel(self,torus,torus_wls):
        wls = np.array([250,350,500])
        indices = [np.argmin(np.abs(torus_wls-wls[i])) for i in range(len(wls))]

        torus_flux = ext.nuLnu_to_Flux(ext.L_sun_to_L(torus),ext.freq_rest_to_obs(ext.wl_to_freq(torus_wls, freq_units="GHz"),self.z), self.dist_cm)
        torus_vals = torus_flux[indices]
        fracs = np.array(torus_vals)/np.array([self.F_250mic,self.F_350mic,self.F_500mic])

        self.F_250mic = self.F_250mic*(1-fracs[0])
        self.F_250mic_Err = self.F_250mic_Err *(1-fracs[0])
        self.F_350mic = self.F_350mic*(1-fracs[1])
        self.F_350mic_Err = self.F_350mic_Err *(1-fracs[1])
        self.F_500mic = self.F_500mic*(1-fracs[2])
        self.F_500mic_Err = self.F_500mic_Err *(1-fracs[2])

        return

    def plot_NuLnu_wl_Lsun(self,flux,freq,gb_wls,wls_CE,nuLnu_CE,torus,torus_wls,showPlots=False):
        """
        This function plot nu*L_nu [L_sun] vs. Wavelength (8-1000[micron]) for the source
        """
        fig = plt.figure(figsize=(12, 8),dpi=600)
        plt.rc('font', size=16)
        ax = plt.axes()

        N16_wl_obs, N16_flux, N16_err = self.getWiseHerschValues([3.4, 4.6, 12, 22, 250, 350, 500])
        N16_wl_rest = ext.wl_obs_to_rest(N16_wl_obs,self.z)
        N16_nuLnu = (ext.Flux_to_nuLnu(np.array(N16_flux), ext.wl_to_freq(np.array(N16_wl_obs),"GHz"), self.dist_cm) / const.L_sun.to(u.erg / u.s)).value
        N16_nuLnu_err = (ext.Flux_to_nuLnu(np.array(N16_err), ext.wl_to_freq(np.array(N16_wl_obs),"GHz"),self.dist_cm) / const.L_sun.to(u.erg / u.s)).value

        # Netezer Measurements - Torus
        plt.errorbar(N16_wl_rest, N16_nuLnu, yerr=N16_nuLnu_err,label='Netzer+16', fmt='^g', zorder=3)

        # Adding ff and bf Graybodies
        if self.group == 'det':
            if (self.T_bf != 47 or self.beta_bf != 1.6):
                gb_bf = ext.scaled_graybody(freq, self.beta_bf, self.T_bf*u.K, self.rest_freq, self.UsedFlux)
                plt.loglog(gb_wls,
                           ext.Flux_to_nuLnu(gb_bf,ext.wl_to_freq(ext.wl_rest_to_obs(gb_wls * 10 ** 9, self.z)),
                                             self.dist_cm).value / const.L_sun.to(u.erg / u.s).value,':', c="orange",
                           label='Gray-Body T=' + str(np.round(self.T_bf, 1)) + "[K], \u03B2=" + str(
                               np.round(self.beta_bf, 2)), zorder=1)
            if (self.T_ff != 47 or self.beta_ff != 1.6):
                gb_ff = ext.scaled_graybody(freq, self.beta_ff, self.T_ff*u.K, self.rest_freq, self.UsedFlux)
                plt.loglog(gb_wls,
                           ext.Flux_to_nuLnu(gb_ff,ext.wl_to_freq(ext.wl_rest_to_obs(gb_wls * 10 ** 9, self.z)),
                                             self.dist_cm).value / const.L_sun.to(u.erg / u.s).value,':', c="darkviolet",
                           label='Gray-Body T=' + str(np.round(self.T_ff, 1)) + "[K], \u03B2=" + str(np.round(self.beta_ff, 2)), zorder=1)
        else:
            gb_bf = 0
            gb_ff = 0
        # ALMA Measurment

        if (np.abs(self.SNR)<3 or self.ALMA_det == 'non-det'):
            plt.errorbar(self.rest_wl, self.nuLnu / const.L_sun.to(u.erg / u.s).value,yerr=[[(self.nuLnu / const.L_sun.to(u.erg / u.s)).value/3],[0]], marker='s',uplims=True, c="red", ms=6,label='ALMA', zorder=4)
        else:
            plt.errorbar(self.rest_wl, self.nuLnu / const.L_sun.to(u.erg / u.s).value,yerr=self.nuLnu_err/ const.L_sun.to(u.erg / u.s).value, marker='s',color='red', mew=2, ms=6, label='ALMA',zorder=4)

        plt.scatter(self.rest_wl, self.ALMA_nuLnu / const.L_sun.to(u.erg / u.s), marker='.', color='red', s=60,zorder=5)

        plt.plot([self.rest_wl, self.rest_wl],[(self.ALMA_nuLnu / const.L_sun.to(u.erg / u.s)).value, (self.nuLnu / const.L_sun.to(u.erg / u.s)).value], c='red')

        # Cherry & Elbaz
        plt.loglog(wls_CE, 10 ** nuLnu_CE / const.L_sun.to(u.erg / u.s).value, label='Chary & Elbaz 01',zorder=2)

        # Main Graybody - Usually T = 47K, beta = 1.6
        plt.loglog(gb_wls,
                   ext.Flux_to_nuLnu(flux, ext.wl_to_freq(ext.wl_rest_to_obs(gb_wls * 10 ** 9, self.z)),
                                     self.dist_cm).value /
                   const.L_sun.to(u.erg / u.s).value, '--', c='black',
                   label="Gray-Body T=" + str(int(self.temp)) + "[K], \u03B2= " + str(self.beta), zorder=1)

        # Torus
        plt.loglog(torus_wls, torus, label="Torus")

        plt.xlabel('Rest Frame Wavelength, \u03BB [\u03BCm]')  # creating x-label
        plt.ylabel('Luminosity, \u03BD$L_{\u03BD}$ [$L_{\u2609}$]')  # creating y-label
        plt.xlim(0.5, 1000)  # Limitting x-axis
        plt.ylim(10 ** 9, 5 * 10 ** 13)  # Limitting y-axis

        # Changing the order of the legend
        handles, labels = ax.get_legend_handles_labels()

        order = ext.orderSEDsLabels(labels)
        labels = [labels[i] for i in order]
        handles = [handles[i] for i in order]
        ax.legend(handles, labels,loc='best')

        plt.title(self.full_name + " ,z=" + str(np.round(self.z, 3)), loc='left')

        plt.tight_layout()

        # saving/plotting the results
        plt.savefig(SEDs_location +"/"+ self.name + '_SED_nuLnu_w_N16.pdf', format="pdf",dpi=600)
        plt.savefig(Sources_location + '/' + self.name + '/'+self.name+'_SED_nuLnu_w_N16.pdf', format="pdf",dpi=600)
        if showPlots:
            plt.show()

        # closing the last figure
        plt.clf()
        plt.cla()
        plt.close()

        plt.close('all')  # closing all figures used in the function
        gc.collect()  # garbage-collector
        matplotlib.rcParams.update(matplotlib.rcParamsDefault)
        return gb_bf,gb_ff

    def plot_NuLnu_wl_cgs(self,flux,freq,gb_wls,wls_CE,nuLnu_CE,torus,torus_wls,gb_bf,gb_ff,showPlots=False):
        """
        This function plot nu*L_nu [cgs] vs. Wavelength (8-1000[micron]) for the source
        """
        plt.figure(figsize=(12, 8),dpi=600)
        plt.rc('font', size=16)
        ax = plt.axes()

        N16_wl_obs, N16_flux, N16_err = self.getWiseHerschValues([3.4, 4.6, 12, 22, 250, 350, 500])
        N16_wl_rest = ext.wl_obs_to_rest(N16_wl_obs,self.z)
        N16_nuLnu = ext.L_sun_to_L(ext.Flux_to_nuLnu(np.array(N16_flux), ext.wl_to_freq(np.array(N16_wl_obs),"GHz"), self.dist_cm) / const.L_sun.to(u.erg / u.s)).value
        N16_nuLnu_err = ext.L_sun_to_L(ext.Flux_to_nuLnu(np.array(N16_err), ext.wl_to_freq(np.array(N16_wl_obs),"GHz"),self.dist_cm) / const.L_sun.to(u.erg / u.s)).value

        # Netezer Measurements
        plt.errorbar(N16_wl_rest, np.log10(N16_nuLnu),
                     yerr=(N16_nuLnu_err) / N16_nuLnu / 2, label='Netzer+16',
                     fmt='^g', zorder=3)  # remember the division by 2 in the error, maybe its not ok

        # ALMA Measurment
        if (np.abs(self.SNR)<3 or self.ALMA_det == 'non-det'):
            plt.errorbar(self.rest_wl, np.log10(self.nuLnu),yerr=[[0.15],[0]], marker='s',uplims=True, c="red", ms=6, label='ALMA',zorder=4)
        else:
            plt.errorbar(self.rest_wl, np.log10(self.nuLnu),yerr=np.array([[-np.log10(1-self.nuLnu_err/self.nuLnu)],[np.log10(1+self.nuLnu_err/self.nuLnu)]]),marker='s',color='red', mew=2, ms=6, label='ALMA', zorder=4)

        plt.scatter(self.rest_wl, np.log10(self.ALMA_nuLnu), marker='.', color='red', s=60,zorder=5)

        plt.plot([self.rest_wl, self.rest_wl], [np.log10(self.ALMA_nuLnu),
                                                np.log10(self.nuLnu)], c='red')

        # Cherry & Elbaz
        plt.semilogx(wls_CE, nuLnu_CE, label='Chary & Elbaz 01', zorder=2)

        # Main Graybody - Usually T = 47K, beta = 1.6
        plt.semilogx(gb_wls, np.log10(
            ext.Flux_to_nuLnu(flux, ext.wl_to_freq(ext.wl_rest_to_obs(gb_wls * 10 ** 9, self.z)),self.dist_cm).value),
                     '--', c='black',label="Gray-Body T=" + str(int(self.temp)) + "[K], \u03B2= " + str(self.beta), zorder=1)

        # Adding ff and bf Graybodies
        if self.group == 'det':
            if (self.T_bf != 47 * u.K or self.beta_bf != 1.6):
                plt.semilogx(gb_wls, np.log10(
                    ext.Flux_to_nuLnu(gb_bf, ext.wl_to_freq(ext.wl_rest_to_obs(gb_wls * 10 ** 9, self.z)),
                                      self.dist_cm).value), ':', c="orange",
                             label='Gray-Body T=' + str(np.round(self.T_bf, 1)) + "[K], \u03B2=" + str(
                                 np.round(self.beta_bf, 2)), zorder=1)

            if (self.T_ff != 47 or self.beta_ff != 1.6):
                plt.semilogx(gb_wls, np.log10(
                    ext.Flux_to_nuLnu(gb_ff, ext.wl_to_freq(ext.wl_rest_to_obs(gb_wls * 10 ** 9, self.z)),
                                      self.dist_cm).value), ':', c="darkviolet",
                             label='Gray-Body T=' + str(np.round(self.T_ff, 1)) + "[K], \u03B2=" + str(
                                 np.round(self.beta_ff, 2)), zorder=1)

        # Torus
        plt.loglog(torus_wls, np.log10(ext.L_sun_to_L(torus)), label="Torus")

        plt.xlabel('Rest Frame Wavelength, \u03BB [\u03BCm]')  # creating x-label
        plt.ylabel('Log \u03BD$L_{\u03BD}$ [erg/s]')  # creating y-label
        plt.xlim(0.1, 600)  # Limitting x-axis
        plt.ylim(43, 47.5)  # Limitting y-axis
        ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())  # Removing scientific notation from y-axis

        # Changing the order of the legend
        handles, labels = ax.get_legend_handles_labels()
        # if (np.abs(self.SNR)<3 or self.ALMA_det == 'non-det'):
        #     handles = handles[-2:] + [handles[0], handles[-3]] + handles[1:-3]
        #     labels = labels[-2:] + [labels[0], labels[-3]] + labels[1:-3]
        # else:
        #     handles = [handles[0], handles[-1], handles[1], handles[-2]] + handles[2:-2]
        #     labels = [labels[0], labels[-1], labels[1], labels[-2]] + labels[2:-2]
        order = ext.orderSEDsLabels(labels)
        labels = [labels[i] for i in order]
        handles = [handles[i] for i in order]
        ax.legend(handles, labels,loc='lower left')

        plt.title(self.full_name + " ,z=" + str(np.round(self.z, 3)), loc='left')

        plt.tight_layout()

        # saving/plotting the results
        plt.savefig(SEDs_location +"/"+ self.name + '_SED_nuLnu_w_N16_erg_s.pdf', format="pdf",dpi=600)
        plt.savefig(Sources_location + '/' + self.name + '/'+self.name+'_SED_nuLnu_w_N16_erg_s.pdf', format="pdf",dpi=600)
        if showPlots:
            plt.show()

        # closing the last figure
        plt.clf()
        plt.cla()
        plt.close()

        plt.close('all')  # closing all figures used in the function
        gc.collect()  # garbage-collector
        matplotlib.rcParams.update(matplotlib.rcParamsDefault)
        return

    def plot_Fnu_freq(self,flux,freq,gb_wls,wls_CE,nuLnu_CE,torus,torus_wls,gb_bf,gb_ff,showPlots=False):
        """
        This function plots Flux Density [erg/Hz/s/cm^2] vs. Rest-Frame Frequency
        """
        plt.figure(figsize=(12, 8),dpi=600)
        plt.rc('font', size=16)
        ax = plt.axes()

        N16_wl_obs, N16_flux, N16_err = self.getWiseHerschValues([3.4, 4.6, 12, 22, 250, 350, 500])
        N16_wl_rest = np.array(ext.wl_obs_to_rest(N16_wl_obs,self.z))
        N16_freq_rest = ext.wl_to_freq(N16_wl_rest,"GHz").value

        # Netezer Measurements - Torus
        plt.errorbar(N16_freq_rest, N16_flux,yerr=N16_err, label='Netzer+16', fmt='^g', zorder=3)

        # ALMA Measurment
        if (np.abs(self.SNR)<3 or self.ALMA_det == 'non-det'):
            plt.errorbar(self.rest_freq * u.GHz, self.UsedFlux,yerr=[[self.UsedFlux/3],[0]],uplims=True, marker='s', c="red", ms=6, label='ALMA', zorder=4)
        else:
            plt.errorbar(self.rest_freq * u.GHz, self.UsedFlux,yerr=self.Flux_Err,marker='s',color='red', mew=2, ms=6, label='ALMA', zorder=4)

        plt.scatter(self.rest_freq * u.GHz, self.ALMA_flux, marker='.', color='red', s=60, zorder=5)

        plt.plot([self.rest_freq, self.rest_freq], [self.ALMA_flux,
                                                self.UsedFlux], c='red')

        # Cherry & Elbaz
        plt.loglog(ext.wl_to_freq(wls_CE * 10 ** 9).value, ext.nuLnu_to_Flux((10 ** nuLnu_CE),
                    ext.wl_to_freq(ext.wl_rest_to_obs(wls_CE * 10 ** 9,self.z)).value, self.dist_cm),
                   label='Chary & Elbaz 01', zorder=2)

        # Main Graybody - Usually T = 47K, beta = 1.6
        plt.loglog(freq, flux, '--', c='black',
                   label="Gray-Body T=" + str(int(self.temp)) + "[K], \u03B2= " + str(self.beta), zorder=1)


        # Adding ff and bf Graybodies
        if (self.group == 'det'):
            if (self.T_bf != 47 or self.beta_bf != 1.6):
                plt.loglog(freq, gb_bf, ':', c="orange",
                           label='Gray-Body T=' + str(np.round(self.T_bf, 1)) + "[K], \u03B2=" + str(
                               np.round(self.beta_bf, 2)), zorder=1)

            if (self.T_ff != 47 or self.beta_ff != 1.6):
                plt.loglog(freq, gb_ff, ':', c="darkviolet",
                           label='Gray-Body T=' + str(np.round(self.T_ff, 1)) + "[K], \u03B2=" + str(
                               np.round(self.beta_ff, 2)), zorder=1)

        # Torus
        plt.loglog(ext.wl_to_freq(torus_wls, freq_units="GHz"),ext.nuLnu_to_Flux(ext.L_sun_to_L(torus),
                    ext.freq_rest_to_obs(ext.wl_to_freq(torus_wls, freq_units="GHz"),self.z), self.dist_cm), label="Torus")

        plt.xlabel('Rest-Frame Frequency,\u03BD [GHz]')  # creating x-label
        plt.ylabel('Flux Density, $F_{\u03BD}$ [erg/Hz/s/cm\u00b2]')  # creating y-label
        plt.xlim(300, 4 * 10 ** 5)  # Limitting x-axis
        plt.ylim(5 * 10 ** (-28), 3 * 10 ** (-24))  # Limitting y-axis
        # adding small ticks in both axes
        plt.tick_params(axis='both', direction="in", which='both', bottom=True, top=True, left=True, right=True)

        # Changing the order of the legend
        handles, labels = ax.get_legend_handles_labels()
        # if (np.abs(self.SNR)<3 or self.ALMA_det == 'non-det'):
        #     handles = handles[-2:] + [handles[0], handles[-3]] + handles[1:-3]
        #     labels = labels[-2:] + [labels[0], labels[-3]] + labels[1:-3]
        # else:
        #     handles = [handles[0], handles[-1], handles[1], handles[-2]] + handles[2:-2]
        #     labels = [labels[0], labels[-1], labels[1], labels[-2]] + labels[2:-2]
        order = ext.orderSEDsLabels(labels)
        labels = [labels[i] for i in order]
        handles = [handles[i] for i in order]
        ax.legend(handles, labels,loc='best')

        plt.title(self.full_name + " ,z=" + str(np.round(self.z, 3)), loc='left')

        plt.tight_layout()

        # saving/plotting the results
        plt.savefig(SEDs_location +"/"+ self.name + '_SED_Fnu_w_N16.pdf', format="pdf",dpi=600)
        plt.savefig(Sources_location + '/' + self.name + '/'+self.name+'_SED_Fnu_w_N16.pdf', format="pdf",dpi=600)

        if showPlots:
            plt.show()

        # closing the last figure
        plt.clf()
        plt.cla()
        plt.close()

        plt.close('all')  # closing all figures used in the function
        gc.collect()  # garbage-collector
        matplotlib.rcParams.update(matplotlib.rcParamsDefault)
        return

    def zoomedInImage(self):

        fig, ax = plt.subplots(figsize=(12, 12),dpi=600)
        # plt.rc('font', size=24)

        box = AnchoredAuxTransformBox(ax.transData, loc=4)  # creating a box at the corner of the figure
        el = Ellipse((0, 0), width=self.BMAJ_pix, height=self.BMIN_pix, angle=self.BPA + 90)  # creating ellipse with beam size
        box.drawing_area.add_artist(el)  # adding the ellipse to the box
        line = ext.get_scale_line(5, self.pix_kpc)
        ax.add_artist(line)
        ax.add_artist(box)

        RA_src_label, Dec_src_label = ext.equatorial_labels(self.RA_SDSS, self.Dec_SDSS, self.sign)

        ylabels = ["+" + str(-i / 2.5) + "''" for i in range(-3, 0)] + [Dec_src_label] + [str(-i / 2.5) + "''" for i
                                                                                          in range(1, 4)]
        xlabels = [str(np.round(i, 2)) + "s" for i in np.arange(-0.09, 0, 0.03)] + [RA_src_label] + [
            "+" + str(np.round(i, 2)) + "s" for i in np.arange(0.03, 0.12, 0.03)]

        half_image_arc = 1.5

        source_pixel_offset_x = int(self.Quasar_pix[1] - half_image_arc / self.pixel_to_arc) - (
                self.Quasar_pix[1] - half_image_arc / self.pixel_to_arc)
        source_pixel_offset_y = int(self.Quasar_pix[0] - half_image_arc / self.pixel_to_arc) - (
                self.Quasar_pix[0] - half_image_arc / self.pixel_to_arc)

        plt.xticks(
            [(half_image_arc / self.pixel_to_arc + 15 * (i) / self.pixel_to_arc + source_pixel_offset_x + 1 / 2) for i in
             np.arange(-0.09, 0.12, 0.03)],reversed(xlabels),fontsize="small")
        plt.yticks([(half_image_arc / self.pixel_to_arc + (i) / self.pixel_to_arc + source_pixel_offset_y + 1 / 2) for i in
                    np.arange(-1.2, 1.6, 0.4)],reversed(ylabels),fontsize="small")

        sub_image = self.data[int(self.Quasar_pix[1] - half_image_arc / self.pixel_to_arc):int(
            self.Quasar_pix[1] + half_image_arc / self.pixel_to_arc + 2),
                    int(self.Quasar_pix[0] - half_image_arc / self.pixel_to_arc):int(
                        self.Quasar_pix[0] + half_image_arc / self.pixel_to_arc + 2)] / self.beam_pix

        plt.title(self.full_name + "\nz=" + str(np.round(self.z, 3)), loc='left',fontsize=25,y=0.92,pad=-20,weight='bold')
        plt.xlabel("Right Ascension", fontsize=22)
        plt.ylabel("Declination", fontsize=22)

        plt.imshow(sub_image, cmap='binary', vmin=0, vmax=np.nanmax(sub_image))
        cbar = plt.colorbar(fraction=0.046, pad=0.04)
        cbar.ax.set_ylabel("[Jy]", fontsize=22)
        cen_loc = +self.CASA_pix[0] - self.Quasar_pix[0] + (len(sub_image) - 2) / 2 + source_pixel_offset_x, + \
            self.CASA_pix[1] - self.Quasar_pix[1] + (len(sub_image) - 2) / 2 + source_pixel_offset_y
        QSO_loc = (len(sub_image) - 2) / 2 + source_pixel_offset_x, (len(sub_image) - 2) / 2 + source_pixel_offset_y
        plt.errorbar(cen_loc[0], cen_loc[1], xerr=0.1 / self.pixel_to_arc, yerr=0.1 / self.pixel_to_arc, marker='+', color='r',
                     label="CASA")
        plt.errorbar(QSO_loc[0], QSO_loc[1], xerr=0.1 / self.pixel_to_arc, yerr=0.1 / self.pixel_to_arc, marker='+', color='b',
                     label="Quasar")
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        cbar.ax.tick_params(labelsize=13)

        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-5, -5))  # Force exponent notation for values in the range of 10^-5
        cbar.formatter = formatter
        cbar.update_ticks()

        plt.legend()
        ax.invert_yaxis()


        plt.tight_layout()
        plt.savefig(Images_location + '/' + self.name + '_view.pdf', format="pdf", bbox_inches='tight',dpi=600)
        plt.savefig(Sources_location + '/' + self.name +'/'+self.name+ '_view.pdf', format="pdf", bbox_inches='tight',dpi=600)



        plt.clf()
        plt.cla()
        plt.close()

    def fullViewImage(self):
        # finding the expressions of the centers Dec and RA for the figures labels
        RA_center, DEC_center = ext.RA_Dec_deg_to_time(self.alpha_mes, self.delta_mes)
        RA_center_label, Dec_center_label = ext.equatorial_labels(RA_center, DEC_center, self.sign)
        ylabels = ["+" + str(-i) + "''" for i in range(-16, 0, 4)] + [Dec_center_label] + [str(-i) + "''" for i in
                                                                                           range(4, 20, 4)]
        fig, ax = plt.subplots(figsize=(12, 12),dpi=600)
        #fig, ax = plt.subplots(figsize=(12, 12))
        plt.xticks([(self.center_pix[0] + 15 * i / self.pixel_to_arc) for i in np.arange(-1.6, 2.0, 0.4)],
                   reversed([str(np.round(i, 1)) + "s" for i in np.arange(-1.6, 0, 0.4)] + [RA_center_label] + [
                       "+" + str(np.round(i, 1)) + "s" for i in np.arange(0.4, 2.0, 0.4)]), fontsize="small")
        plt.yticks([(self.center_pix[1] + i / self.pixel_to_arc) for i in range(-16, 20, 4)], reversed(ylabels),
                   fontsize="small")

        box = AnchoredAuxTransformBox(ax.transData, loc=4)  # creating a box at the corner of the figure
        # the box is created a seocnd time in order to include it in the entire observations image
        el = Ellipse((0, 0), width=self.BMAJ_pix, height=self.BMIN_pix,angle=self.BPA + 90)  # creating ellipse with beam size #-BPA+90
        box.drawing_area.add_artist(el)  # adding the ellipse to the box
        ax.add_artist(box)
        pownorm = 1.0
        # plt.imshow(self.data / self.beam_pix, cmap='binary', vmin=0, vmax=np.nanmax(self.data / self.beam_pix))
        plt.imshow(self.data / self.beam_pix, cmap='binary', vmin=8*10**-6, vmax=3*10**-5)
        cbar = plt.colorbar(fraction=0.046, pad=0.04)

        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-5, -5))  # Force exponent notation for values in the range of 10^-5
        cbar.formatter = formatter
        cbar.update_ticks()

        circle = plt.Circle(self.Quasar_pix, 30, color='cyan', fill=False)
        ax.add_patch(circle)
        plt.title(self.full_name + "\nz=" + str(np.round(self.z, 3)), loc='left',fontsize=25,y=0.95,pad=-20,weight='bold')
        plt.xlabel("Right Ascension", fontsize=22)
        plt.ylabel("Declination", fontsize=22)
        cbar.ax.set_ylabel("[Jy]",fontsize=22)
        ax.invert_yaxis()

        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        cbar.ax.tick_params(labelsize=13)

        plt.tight_layout()

        plt.savefig(Images_location +"/" + self.name + '_fullview.pdf', format="pdf", bbox_inches='tight',dpi=600)
        plt.savefig(Sources_location + '/' + self.name +'/'+self.name+'_fullview.pdf', format="pdf", bbox_inches='tight',dpi=600)

        plt.clf()  # closing the last figure
        plt.cla()
        plt.close()

    def findTBeta(self):
        T_range = np.arange(20,70,0.5)
        beta_range=np.arange(1,2.5,0.05)
        self.T_bf, self.beta_bf, xi_bf = self.Xi_square(T_range,[1.6])
        self.T_ff, self.beta_ff, xi_ff = self.Xi_square(T_range,beta_range)
        return

    def saveToFile(self):
        file_dict = OrderedDict()
        file_dict = {"ind": [self.ind],
                     "source": [self.name],
                     "group": [self.group],
                     "Alma group": [self.ALMA_det],
                     "redshift(z)": [self.z],
                     "Temperature[K]": [self.temp],
                     "beta": [self.beta],
                     "distance[Mpc]": [self.distance.value],
                     "Observed-frequency[Ghz]": [np.round(self.ob_freq, 2)],
                     "RF-frequency[GHz]": [np.round(self.rest_freq, 2)],
                     "Scale[arcsec/kpc]": [np.round(self.arc_kpc, 2)],
                     "Measured-Flux[erg/s/cm2/Hz]": [self.Flux],
                     "Used-Flux (Measured/3RMS)[erg/s/cm2/Hz]": [self.UsedFlux],
                     "Noise[erg/s/cm2/Hz]": [self.Flux_Err],
                     "RMS[erg/s/cm2/Hz]": [self.Flux_RMS],
                     "SNR": [np.round(self.SNR, 3)],
                     "Torus Fraction": [self.torus_fraction],
                     "Measured-nu*L_nu[erg/s]": [self.nuLnu],
                     "GB-Flux[erg/cm2/s]": [self.total_flux],
                     "CE-Flux[erg/cm2/s]": [(ext.Lum_to_Flux((10 ** self.L_CE),self.dist_cm))],
                     "GB-Luminosity[erg/s]": [self.total_lum],
                     "CE-Luminosity[erg/s]": [10 ** (self.L_CE)],
                     "GB-SFR[M_sun/year]": [np.round(self.SFR_gb, 3)],
                     "CE-SFR[M_sun/year]": [np.round(self.SFR_CE, 3)],
                     "GB-nu*L_nu[erg/s]": [self.gb_nuLnu.value],
                     "CE-nu*L_nu[erg/s]": [10 ** self.CE_nuLnu],
                     "CE-SED#": [self.CE_template_num+1],
                     "log-L_AGN[erg/s]": [self.L_AGN],
                     "log-L_SF[L_sun]": [self.L_SF],
                     "MBH(lower-limit)[M_sun]": [np.round(self.MBH,3)],
                     "M_dot_BH[M_sun/yr]": [self.M_dot_BH],
                     "GB-SFR/M_dot_BH": [np.round(self.SFR_gb / self.M_dot_BH, 3)],
                     "CE-SFR/M_dot_BH": [np.round(self.SFR_CE / self.M_dot_BH, 3)],
                     "T-bf": [self.T_bf],
                     "beta-bf": [self.beta_bf],
                     "T-ff": [self.T_ff],
                     "beta-ff": [self.beta_ff],
                     "Neigh-Flux[erg/s/cm2/Hz]": [""],
                     "neigh-noise[erg/s/cm2/Hz]": [""],
                     "Neigh-SNR": [""],
                     "Source-Neigh-ratio": [""],
                     "SFR-GB-Neigh[M_sun/year]": [""],
                     "SFR-CE-Neigh[M_sun/year]": [""],
                     }
        if self.neigh:
            file_dict["Neigh-Flux[erg/s/cm2/Hz]"] = [self.companion.F]
            file_dict["neigh-noise[erg/s/cm2/Hz]"] = [self.companion.F_err]
            file_dict["Neigh-RMS"] = [self.companion.RMS]
            file_dict["Neigh-SNR"] = [str(np.round(self.companion.SNR, 3))]
            file_dict["Source-Neigh-ratio"] = [str(np.round(self.neigh_ratio, 3))]
            file_dict["SFR-GB-Neigh[M_sun/year]"] = [str(np.round(self.companion.SFR_gb, 3))]
            file_dict["SFR-CE-Neigh[M_sun/year]"] = [str(np.round(self.companion.SFR_CE, 3))]

        df_temp = pd.DataFrame.from_dict(file_dict, orient='columns')
        return df_temp


class Companion():
    def __init__(self,F,F_err,RMS,RA,Dec,main_z):
        self.F = F
        self.F_err = F_err
        self.RMS = RMS
        self.SNR = self.F/self.RMS
        self.RA = RA
        self.Dec = Dec
        self.z = main_z


    def find_SFRs(self,freq,ob_freq,temp,beta,dist_cm,df_CE,nuLnu_row_CE,wl_index):
        flux = ext.scaled_graybody(freq, beta, temp * u.K, ext.freq_obs_to_rest(ob_freq,self.z), self.F)  # GB flux array

        self.nuLnu = ext.Flux_to_nuLnu(self.F, ob_freq, dist_cm)
        self.CE_template_num = np.argmin(np.abs(np.log10(self.nuLnu) - nuLnu_row_CE))  # finding the closest value of nu*L_nu
        self.CE_nuLnu = df_CE.iloc[wl_index + 1, self.CE_template_num + 1]
        chary_nuLnu_arr = df_CE.iloc[1:, self.CE_template_num + 1]  # the distribution which fits the best

        self.L_CE = df_CE.iloc[0, self.CE_template_num + 1]  # CE Luminosity (first row in table)
        self.gb_flux = integrate.simps(y=flux / (1 + self.z),x=freq * (10 ** 9))  # integration over obsereved frequencies
        self.gb_lum = ext.Flux_to_Lum(self.gb_flux, dist_cm)  # luminosity in erg/s (L=F*4*pi*(D_l)^2)
        self.gb_lum_solar = ext.L_to_L_sun(self.gb_lum)
        self.SFR_gb = (self.gb_lum_solar) / (10 ** 10)  # SFR=L/(l_sun*10^10)
        self.SFR_CE = ((10 ** self.L_CE) * u.erg / u.s) / (const.L_sun.to(u.erg / u.s) * (10 ** 10))
        return