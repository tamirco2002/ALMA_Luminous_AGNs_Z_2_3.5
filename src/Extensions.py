import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u, constants as const
from astropy.cosmology import LambdaCDM
from astropy.modeling import models
from matplotlib.lines import Line2D

"""input: blackbody model and frequency in GHz,power law
output : graybody model of this frequency (beta defined at the start)"""
def graybody(blackbody,frequency,beta):
    return blackbody(frequency) * 4 * math.pi * u.sr * (((frequency / u.GHz) * 10 ** 9) ** beta)

def scaled_graybody(frequency, beta, temp, rest_freq, Flux):
    temp_bb = models.BlackBody(temperature=temp, scale=1.0)  # black-body model *without* scaling
    gb = graybody(temp_bb, rest_freq * u.GHz, beta)  # value of original gray-body SED at rest-frame frequency
    gb_scale = float((Flux / gb).value)
    bb = models.BlackBody(temperature=temp, scale=gb_scale)  # black-body model *with* scaling
    return graybody(bb, frequency,beta)

def graybody2(blackbody,frequency,beta):
    return (blackbody(frequency) * 4 * math.pi * u.sr * (((frequency / u.GHz) * 10 ** 9) ** beta)).value

def gb_for_fit(rest_freq, Flux):

    def inner_func(frequency, beta, temp):
        temp_bb = models.BlackBody(temperature=temp*u.K, scale=1.0)  # black-body model *without* scaling
        gb = graybody2(temp_bb, rest_freq * u.GHz, beta)  # value of original gray-body SED at rest-frame frequency
        gb_scale = float((Flux / gb))
        bb = models.BlackBody(temperature=temp*u.K, scale=gb_scale)  # black-body model *with* scaling
        return graybody2(bb, frequency, beta)

    return inner_func


#rest-obs (wl,freq), jansky to erg/s arc to deg,pixel to deg,beam to pix
"""input: frequency [Hz]
"output: wavelength [micron]"""
def freq_to_wl(f,freq_units="Hz"):
    if(freq_units=="GHz"):
        wl = ((const.c.to(u.micron / u.s) / (f*10**9)).value) * u.micron
    if(freq_units=="Hz"):
        wl = ((const.c.to(u.micron / u.s) / f).value) * u.micron
    return wl

"""input: wavelength [micron]
"output: frequency [Hz]"""
def wl_to_freq(wl,freq_units="Hz"):
    if (freq_units == "Hz"):
        f=(const.c.to(u.micron / u.s) / wl)
    if (freq_units == "GHz"):
        f=(const.c.to(u.micron / u.s) / wl)/(10**9)
    return f

"""input: rest frequency and redshift
output: observed frequency"""
def freq_rest_to_obs(f_rest,z):
    return f_rest/(1+z)

"""input: observed frequency and redshift
output: rest frequency"""
def freq_obs_to_rest(f_obs,z):
    return f_obs*(1+z)

"""input: rest wavelength and redshift
output: observed wavelength"""
def wl_rest_to_obs(wl_rest,z):
    return wl_rest*(1+z)

"""input: observed wavelength and redshift
output: rest wavelength"""
def wl_obs_to_rest(wl_obs,z):
    return wl_obs/(1+z)

"""input: size in arc-seconds
output: size in degrees"""
def arc_to_deg(size):
    return size/3600

"""input: size in degrees
output: size in arc-seconds"""
def deg_to_arc(size):
    return size * 3600

#def pixel_conversion(size,ratio):

"""input: Flux in Jansky(unitless parameter)
output: Flux in erg/s (with units)"""
def Jansky_to_erg_s(val):
    return val*10**-23* u.erg / u.cm ** 2 / u.Hz / u.s

"""input: Flux in erg/s (with units)
output: Flux in Jansky(unitless parameter)"""
def erg_s_to_Jansky(val):
    return val/(10**-23* u.erg / u.cm ** 2 / u.Hz / u.s)

"""input: Flux[erg/s/cm**2/Hz], distance[cm]
output: L [erg/s]"""
def Flux_to_Lum(F,d):
    L= F*4*np.pi*d**2
    return L

"""input:  L [erg/s], distance[cm]
output: Flux[erg/s/cm**2/Hz]"""
def Lum_to_Flux(L,d):
    F= L/(4*np.pi*d**2)
    return F

def L_to_L_sun(L):
    new_L= L/const.L_sun.to(u.erg / u.s).value
    return new_L

def L_sun_to_L(L):
    new_L= L*const.L_sun.to(u.erg / u.s).value
    return new_L

"""input:  Flux[erg/s/cm**2/Hz] ,observed frequency[GHz], distance[cm]
output: NuLnu (unitless)"""
def Flux_to_nuLnu(Flux,f_obs,d):
    nuLnu=np.abs(((Flux)*4*np.pi*(d**2)*f_obs*10**9))
    return nuLnu

"""input:  NuLnu (unitless),observed frequency[GHz], distance[cm]
output: Flux[erg/s/cm**2/Hz]"""
def nuLnu_to_Flux(L,f_obs,d):
    F=np.abs(L/(4*np.pi*(d**2)*f_obs*10**9))
    return F

"input: array,size of subarray [pixels],center of image [pixels], pixel to degree ratio, alpha and delta of center and source" \
"output: center of the source\ brightest pixel [pixels] "
def find_center(arr,size,center_sample,pixel_to_deg,alpha_mes,delta_mes,alpha_src,delta_src):
    # The estimated location of the source according to NED, assuming the center of the image is the coordinates from the header
    est_center = equatorial_to_pixel(center_sample,alpha_src,delta_src,alpha_mes,delta_mes,pixel_to_deg)
    # ranges of x and y in which to look for the sources center
    x_check = np.maximum(0, est_center[1] - size), np.minimum(est_center[1] + size, np.shape(arr)[1])
    y_check = np.maximum(0, est_center[0] - size), np.minimum(est_center[0] + size, np.shape(arr)[0])
    slice_arr = arr[int(y_check[0]):int(y_check[1]),int(x_check[0]):int(x_check[1])]  # sub-array to search for the source
    #val = np.nanmax(slice_arr) * 10 ** -23 * u.erg / u.cm ** 2 / u.Hz / u.s  # flux measured at center of the galaxy
    #print(val)
    max_pix = np.unravel_index(np.nanargmax(slice_arr),np.shape(slice_arr))  # indices of the source's center in sub-array
    max_pix = tuple(map(sum, zip(max_pix, (y_check[0], x_check[0]))))  # indices in original array
    return max_pix
    #return est_center

"input: array,radius of search [pixels],center of image [pixels], pixel to degree ratio, alpha and delta of center and source" \
"output: center of the source\ brightest pixel [pixels] "
def find_center_radius(arr,radius,center_sample,pixel_to_deg,alpha_mes,delta_mes,alpha_src,delta_src):
    arr_copy=np.copy(arr)
    est_center = equatorial_to_pixel(center_sample, alpha_src, delta_src, alpha_mes, delta_mes, pixel_to_deg)
    y_grid, x_grid = np.ogrid[0:len(arr_copy), 0:len(arr_copy[0])]
    mask = (x_grid-est_center[1]) ** 2 + (y_grid-est_center[0]) ** 2 > radius ** 2
    arr_min = np.min(arr_copy)
    arr_copy[mask] = arr_min
    max_pix = np.unravel_index(np.nanargmax(arr_copy), arr_copy.shape)
    return max_pix

""" input: area [pixel^2], center pixel (tupple)
output: flux [Jy] inside a square with given area and around the center point."""
def flux_square(data,area,center,beam_pix):
    f=0
    #area_in_pix = area / pix_arc
    r=np.sqrt(area)/(2)
    for i in range(len(data)):
        dist_y = np.abs(center[0] - i)
        if(dist_y<r+0.5):
            for j in range(len(data[i])):
                dist_x=np.abs(center[1]-j)
                if(dist_x<r+0.5):
                    f+=np.nan_to_num(data[i][j])
    return f/beam_pix

""" input: Major-axis size [pixel], Minor-axis size [pixel],center pixel (tupple),ellipse-angle [deg], factor.
output:  flux [Jy] inside an ellipse with given major and minor axis, center point and an angle.
the factor sets the size of the area (factor=1 is the area created by given axis)"""
def flux_ellipse(data,major,minor,center,angle,factor,beam_pix):
    major =major*np.sqrt(factor)
    minor = minor * np.sqrt(factor)
    f=0
    f_arr=[]
    new_center = center[1] * np.sin(angle) + center[0] * np.cos(angle), center[1] * np.cos(angle) - center[0] * np.sin(angle)
    for i in range(len(data)):
        for j in range(len(data[i])):
            i_rot=j*np.sin(angle)+i*np.cos(angle)
            j_rot=j*np.cos(angle)-i*np.sin(angle)
            p = (((i_rot- new_center[0])**2 / minor**2) + ((j_rot- new_center[1])**2 / major**2))
            if (p<=1):
                f+=np.nan_to_num(data[i][j])
                f_arr.append(np.nan_to_num(data[i][j]))

    return f/beam_pix

""" input: Major-axis size [pixel], Minor-axis size [pixel],center pixel (tupple),ellipse-angle [deg], source's 
ellipse factor,seperation ring factor, background ellipse factor.
output:  flux [Jy] inside the source's ellipse,source's flux RMS, flux [Jy] of background ellipse surrounding
the source's ellipse, background's flux RMS.
note: the ellipse's has the given major and minor axis, center point and an angle.
the factor sets the size of the area (factor=1 is the area created by given axis), while factor_src <= factor_bg
(factor_src = factor_bg for no background). """
def aparture(data,major,minor,center,angle,factor_src,factor_sep,factor_bg,beam_pix):
    #source's ellipse major/minor axis
    major_src = major * np.sqrt(factor_src)
    minor_src = minor * np.sqrt(factor_src)
    # gap's ellipse major/minor axis
    major_sep = major * factor_sep
    minor_sep = minor * factor_sep
    # background's ellipse major/minor axis
    major_bg = major * factor_bg
    minor_bg = minor * factor_bg
    #ellipses areas
    src_area=np.pi*major_src*minor_src
    sep_area = np.pi * major_sep * minor_sep
    bg_area=np.pi*major_bg*minor_bg
    #arrays of flux values of pixel inside the ellipses area
    f_arr_src = []
    f_arr_bg =[]
    new_center = center[1] * np.sin(angle) + center[0] * np.cos(angle), center[1] * np.cos(angle) - center[0] * np.sin(angle)
    for i in range(len(data)):
        for j in range(len(data[i])):
            i_rot=j*np.sin(angle)+i*np.cos(angle)
            j_rot=j*np.cos(angle)-i*np.sin(angle)
            p_src = (((i_rot- new_center[0])**2 / minor_src**2) + ((j_rot- new_center[1])**2 / major_src**2))
            p_sep = (((i_rot - new_center[0]) ** 2 / minor_sep ** 2) + ((j_rot - new_center[1]) ** 2 / major_sep ** 2))
            p_bg = (((i_rot - new_center[0]) ** 2 / minor_bg ** 2) + ((j_rot - new_center[1]) ** 2 / major_bg ** 2))
            if (p_src <= 1):
                f_arr_src.append(np.nan_to_num(data[i][j]))
            elif (p_bg <=1 and p_sep>1):
                f_arr_bg.append(np.nan_to_num(data[i][j]))

    #ajusting to Jy/pixel
    f_arr_src=np.array(f_arr_src)/beam_pix
    f_arr_bg=np.array(f_arr_bg)/beam_pix
    #summing for the entire measured flux
    f_src=np.sum(f_arr_src) #source's flux
    f_bg=np.sum(f_arr_bg) #background's flux
    std_bg=np.std(f_arr_bg)
    bg_mean=f_bg/len(f_arr_bg) #average bg flux *per pixel*
    bg_rms=np.sqrt(np.sum(f_arr_bg**2)/len(f_arr_bg)) #BG RMS *per pixel*
    src_rms=np.sqrt(np.sum(f_arr_src**2)/len(f_arr_src)) # Source's RMS *per pixel(
    noise=np.sqrt(len(f_arr_src)*((bg_rms)**2)+(f_src*0.1)**2) #noise = [((#src_pixels)*RMS(BG))^2+((#src_pixels)*RMS(SRC))^2)]^0.5
    noise_without_src = np.sqrt(len(f_arr_src) * ((bg_rms) ** 2))
    SNR=f_src/noise  # SNR = src_Flux/noise
    SNR_without_src=f_src/noise_without_src
    corrected_flux=np.abs(f_src-noise)
    return (f_src, src_rms, f_bg, bg_rms,SNR, corrected_flux,noise,SNR_without_src,noise_without_src)


"""input: range of beta values, range of T values, rest frequency [GHz] (unitless), Flux[erg/s],
add-on measurments frequencies,fluxes and errors
output: Optimal T[K],beta and Xi-square"""
def Xi_Square(beta_range,T_range,rest_freq,ALMA_Flux,freq_mes,flux_mes,err_mes,error):
    T_min = 0
    beta_min = 0
    xi_min = np.inf

    x = np.append(freq_mes, rest_freq)
    y = np.append(flux_mes, ALMA_Flux.value)
    dy = np.append(err_mes, error)

    for t in range(len(T_range)):
        for b in range(len(beta_range)):

            model = scaled_graybody(x*u.GHz, beta_range[b], T_range[t]*u.K, rest_freq, ALMA_Flux).value
            xi = np.sum((((y - model) / dy) ** 2)) / (len(y) + 1)
            if (np.abs(xi) < np.abs(xi_min)):
                xi_min = xi
                T_min = T_range[t]
                beta_min = beta_range[b]
    return T_min,beta_min,xi_min


def fit_torus(wise_nuLnu,wise_err, wise_wl, spec_nuLnu,spec_wl):
    checked_wl = np.argmin(np.abs(spec_wl - wise_wl))
    val=wise_nuLnu/spec_nuLnu[checked_wl]
    return spec_nuLnu*val

# plotting the change in flux for measurments at different beam size
def beam_size_analyzer(lowest,highest,step,data,major,minor,center,angle,beam_pix,path,Galaxy_name):
    beam_range=np.arange(lowest,highest,step)
    beam_flux1=np.array([aparture(data,major,minor,center,angle,i,2*i,beam_pix) for i in beam_range])
    beam_flux=np.array([np.abs(beam_flux1[i][0]-beam_flux1[i][2]) for i in range(len(beam_flux1))])
    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()
    ax1.plot(beam_range,beam_flux,'o')
    ax2.plot(beam_range,beam_flux/beam_flux[0],'o')
    ax2.set_xlabel("beam size")
    ax2.set_ylabel("1 beam flux")
    ax1.set_xlabel("beam size")
    ax1.set_ylabel("flux")
    plt.savefig(path+'/'+Galaxy_name+'flux_grad_test.png')
    return beam_flux

"""input: Right ascention and Declination in the format [RA]=h:m:s , [Dec]=deg:arc-min:arc-sec
output: Right ascention and Declination as alpha and delta in arc-sec"""
def alpha_delta(RA,Dec,sign=1):
    alpha = (RA[0] + RA[1] / 60 + RA[2] / 3600) * 15
    delta=(np.abs(Dec[0])+Dec[1]/(60)+Dec[2]/(3600))*sign
    return alpha,delta

"""input:  Right ascention and Declination as alpha and delta in arc-sec
output:Right ascention and Declination in the format [RA]=h:m:s , [Dec]=deg:arc-min:arc-sec"""
def RA_Dec_deg_to_time(alpha,delta):
    alpha_sec = alpha / 15
    hour_ra = np.floor(alpha_sec)
    min_ra = np.floor((alpha_sec - hour_ra) * 60)
    sec_ra = np.round(((alpha_sec - hour_ra - min_ra / 60) * 3600),4)
    RA = [hour_ra, min_ra, sec_ra]
    deg_dec = np.floor(np.abs(delta))
    amin_dec = np.floor((np.abs(delta) - deg_dec) * 60)
    asec_dec = np.round(np.abs(delta) * 3600 - deg_dec * 3600 - amin_dec * 60,2)

    if (delta < 0):
        deg_dec = deg_dec * (-1)
    Dec=[deg_dec,amin_dec,asec_dec]
    return RA,Dec

"""input:center of image [pixels], alpha and delta of center and source in degree, pixel to degree ratio
output: coordinates of the src in pixels"""
def equatorial_to_pixel(center_sample,alpha_src,delta_src,alpha_cen,delta_cen,pixel_to_deg):
    cos_theta = np.sin(delta_cen*np.pi/180)*np.sin(delta_src*np.pi/180)+np.cos(delta_cen*np.pi/180)*np.cos(delta_src*np.pi/180)*np.cos((alpha_cen-alpha_src)*np.pi/180)
    #print((cos_theta-np.cos(delta_cen*np.pi/180))*(alpha_src - alpha_cen)/pixel_to_deg)
    coordinates = [center_sample[1] + round((delta_src - delta_cen) / pixel_to_deg),
              center_sample[0] - round((alpha_src - alpha_cen)*np.cos(delta_cen*np.pi/180) / pixel_to_deg)]
    unrounded = [center_sample[1] + ((delta_src - delta_cen) / pixel_to_deg),
              center_sample[0] - ((alpha_src - alpha_cen)*np.cos(delta_cen*np.pi/180) / pixel_to_deg)]

    return coordinates,unrounded



    return coordinates,unrounded

"""input:center of image [pixels], alpha and delta of center and source in degree, pixel to degree ratio
output: coordinates of the src in pixels"""
def pixel_to_equatorial(src_pix,center_pix,alpha_cen,delta_cen,pixel_to_deg):
    coordinates=[alpha_cen+pixel_to_deg*(center_pix[0]-src_pix[0]),delta_cen+(src_pix[1]-center_pix[1])*pixel_to_deg]
    return coordinates

"""input: RA and Dec in the format [RA]=h:m:s [Dec]=deg:arc-min:arc-sec,sign of delta
output: labels of the coordinates for the plots"""
def equatorial_labels(RA,Dec,sign):
    Dec_label=str(int(Dec[0]))+u'\N{DEGREE SIGN}'+str(int(Dec[1]))+"'"+str(np.round(Dec[2],2))+"''"
    if(sign==1):
        Dec_label="+"+Dec_label
    elif(Dec[0]==0):
        Dec_label = "-" + Dec_label
    Ra_label = str(int(RA[0])) + "h" + str(int(RA[1])) + "m" + str(np.round(RA[2],2)) + "s"
    return Ra_label,Dec_label


class AnchoredHScaleBar(matplotlib.offsetbox.AnchoredOffsetbox):
    """ size: length of bar in data units
        extent : height of bar ends in axes units """
    def __init__(self, size=1, extent = 0.03, label="", loc=2, ax=None,
                 pad=0.4, borderpad=0.5, ppad = 0, sep=2, prop=None,
                 frameon=True, linekw={}, **kwargs):
        if not ax:
            ax = plt.gca()
        trans = ax.get_yaxis_transform()
        size_bar = matplotlib.offsetbox.AuxTransformBox(trans)
        line = Line2D([0,0],[0,size], **linekw)
        vline1 = Line2D([-extent/2.,extent/2.],[0,0], **linekw)
        vline2 = Line2D([-extent/2.,extent/2.],[size,size], **linekw)
        size_bar.add_artist(line)
        size_bar.add_artist(vline1)
        size_bar.add_artist(vline2)
        txt = matplotlib.offsetbox.TextArea(label, minimumdescent=False)
        self.vpac = matplotlib.offsetbox.VPacker(children=[size_bar,txt],
                                 align="center", pad=ppad, sep=sep)
        matplotlib.offsetbox.AnchoredOffsetbox.__init__(self, loc, pad=pad,
                 borderpad=borderpad, child=self.vpac, prop=prop, frameon=frameon,
                 **kwargs)

def get_scale_line(size_kpc,pix_kpc):
    size_pix=size_kpc/pix_kpc
    ob = AnchoredHScaleBar(size=size_pix, label=str(size_kpc)+" kpc", loc=3, frameon=True,
                       pad=0.6,sep=4, linekw=dict(color="crimson"),)
    return ob


c = const.c.to(u.cm / u.s).value
M_sun = const.M_sun.to(u.g).value
year = 3.154 * 10 ** 7  # seconds in a year


def L_AGN_to_MBH(L_AGN):
    MBH=8*M_sun*(10**-39)*(10 ** L_AGN)/(M_sun)
    return MBH

"""input:log L_AGN [erg/s]
output: M_dot_BH [M_sun/yr]"""
def L_AGN_to_M_dot_BH(L_AGN):
    M_dot_BH = 9 * (10 ** L_AGN) * year / (M_sun * c ** 2)
    return (M_dot_BH)

def M_dot_BH_to_L_AGN(M_dot_BH):
    L_AGN = np.log10(M_dot_BH*(M_sun * c ** 2)/(9*year))
    return (L_AGN)

"""input: SFR [M_sun/yr]
output: log L_Sf [erg/s]"""
def SFR_to_L_SF(SFR):
    L_SF=np.log10(SFR*const.L_sun.to(u.erg/u.s).value*10**10)
    return L_SF

def L_SF_to_SFR(L_SF):
    SFR=10**L_SF/(const.L_sun.to(u.erg/u.s).value*10**10)
    return SFR

def powerLaw(x, a, b):
    return a * x**-b

def orderSEDsLabels(labels):
    finalOrder = []
    for l in range(len(labels)):
        if 'ALMA' in labels[l]:
            finalOrder.append(l)
            break
    for l in range(len(labels)):
        if 'Netzer' in labels[l]:
            finalOrder.append(l)
            break
    for l in range(len(labels)):
        if 'Chary' in labels[l]:
            finalOrder.append(l)
            break
    for l in range(len(labels)):
        if 'Torus' in labels[l]:
            finalOrder.append(l)
            break
    for l in range(len(labels)):
        if '47' in labels[l] and '1.6' in labels[l]:
            finalOrder.append(l)
            break
    for l in range(len(labels)):
        if '47' not in labels[l] and '1.6' in labels[l]:
            finalOrder.append(l)
            break
    for l in range(len(labels)):
        if 'Gray-Body' in labels[l] and '47' not in labels[l] and '1.6' not in labels[l]:
            finalOrder.append(l)
            break
    return finalOrder

def finalMass(M0,t,GR):
    """
    :param M0: Initial mass in solar-masses (log10)
    :param t: Time passes in years (log10)
    :param GR: Growth Rate in solar-masses/yr
    :return: M: Final mass in solar-masses/yr (log10)
    """
    M = np.log10(10**M0+GR*10**t)
    return M



def age_of_universe(redshift):
    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    age = cosmo.age(redshift)
    return age



