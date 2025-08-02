import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import container
from matplotlib.colors import to_rgba

import Extensions as ext


# Figure 6  - AGN luminosities versus star-formation luminosities
def SFR_MBH_dot(output_path,type='GB',GB_type='BB'):
    results_file_path = f'../data/Results/{GB_type}.xlsx'
    db = pd.read_excel(results_file_path)

    SFR = np.log10(np.abs(db[f'{type}-SFR[M_sun/year]'].to_numpy()))

    L_AGN=db['log-L_AGN[erg/s]'].to_numpy()
    hersch=db['group'].to_numpy()
    SNR=db['SNR'].to_numpy()
    Alma_det = db['Alma group']
    L_SFR=db['log-L_SF[L_sun]'].to_numpy()-10

    Errors_ratio_low = -np.log10(1-db['Noise[erg/s/cm2/Hz]']/db['Measured-Flux[erg/s/cm2/Hz]'])
    Errors_ratio_up = np.log10(1+db['Noise[erg/s/cm2/Hz]']/db['Measured-Flux[erg/s/cm2/Hz]'])
    for x in range(len(Errors_ratio_low)):
        if Alma_det[x]!='det':
            Errors_ratio_low[x] = 0.15
            Errors_ratio_up[x] = 0
    labels={}

    font = {'family': 'normal','size': 16}
    matplotlib.rc('font', **font)

    fig,axes=plt.subplots(dpi=600)
    fig.set_figheight(9)
    fig.set_figwidth(12)

    for i in range(len(SFR)):
        lab = ""
        if (hersch[i]=='det'):
            col='r'
            lab+="$\it{Herschel}$ Detected"
        else:
            col='b'
            lab += "$\it{Herschel}$ Non-Detected"

        if(np.abs(SNR[i])>=3 and Alma_det[i]=='det'):
            shape='o'
            fc = col
        else:
            shape='o'
            fc = 'none'
            lab +=" Upper Limit"

        if(lab not in labels):
            if Errors_ratio_up[i] != 0:
                axes.errorbar(L_AGN[i],SFR[i],yerr=[[Errors_ratio_low[i]],[Errors_ratio_up[i]]],marker=shape,color=col,ecolor=to_rgba(col, alpha=0.4),markerfacecolor=fc,ms=8,label=lab,capthick=0.5)
            else:
                axes.errorbar(L_AGN[i], SFR[i], yerr=[[Errors_ratio_low[i]],[Errors_ratio_up[i]]],uplims=True, marker=shape, color=col,ecolor=to_rgba(col, alpha=0.4), markerfacecolor=fc,label=lab, ms=8, capthick=0.5)
            labels[lab]=lab
        else:
            if Errors_ratio_up[i] != 0:
                axes.errorbar(L_AGN[i], SFR[i],yerr=[[Errors_ratio_low[i]],[Errors_ratio_up[i]]], marker=shape, color=col,ecolor=to_rgba(col, alpha=0.4),markerfacecolor=fc, ms=8,capthick=0.5)
            else:
                axes.errorbar(L_AGN[i], SFR[i],yerr=[[Errors_ratio_low[i]],[Errors_ratio_up[i]]],uplims=True, marker=shape, color=col,ecolor=to_rgba(col, alpha=0.4),markerfacecolor=fc, ms=8,capthick=0.5)

        if(np.abs(db['Source-Neigh-ratio'].to_numpy()[i])>0):
            if("has neighbour" not in labels):
                axes.scatter(L_AGN[i],SFR[i],marker='+',c='black',s=45,zorder=3,label="Detected Companion")
                labels["has neighbour"]="has neighbour"
            else:
                axes.scatter(L_AGN[i], SFR[i], marker='+', c='black', s=45, zorder=3)

        if(L_SFR[i]>0):
            if("N16 L_SFR" not in labels):
                axes.scatter(L_AGN[i],L_SFR[i],marker='s',c=col,s=40,label="Netzer+16 ${L}_{SFR}$")
                labels["N16 L_SFR"] = "N16 L_SFR"
            else:
                axes.scatter(L_AGN[i], L_SFR[i], marker='s', c=col, s=40)
            axes.plot([L_AGN[i],L_AGN[i]],[SFR[i],L_SFR[i]],c='black',zorder=-2,ls='--',alpha=0.7)

    plt.ylabel("$log_{10}$ SFR")
    plt.xlabel("$log_{10}$ $L_{AGN}$")

    plt.xlim(45,48)
    plt.ylim(1,4.4)

    ax2 = axes.twinx()
    mn, mx = axes.get_ylim()
    ax2.set_ylim((ext.SFR_to_L_SF(10**mn)), (ext.SFR_to_L_SF(10**mx)))
    ax2.set_ylabel("$log_{10}$ $L_{SF}$")
    ax3 = axes.twiny()
    mn2, mx2 = axes.get_xlim()
    trendx=np.arange(np.log10(ext.L_AGN_to_M_dot_BH(mn2)),np.log10(ext.L_AGN_to_M_dot_BH(mx2)),step=0.01)
    mn3,mx3=ax2.get_ylim()
    trendx2=np.arange(mn2,mx2,step=0.01)

    SFR_range = np.linspace(-1, 5)
    MBH50 = ext.M_dot_BH_to_L_AGN(10 ** SFR_range / 50)
    MBH200 = ext.M_dot_BH_to_L_AGN(10 ** SFR_range / 200)
    MBH500 = ext.M_dot_BH_to_L_AGN(10 ** SFR_range / 500)
    axes.plot(MBH50, SFR_range, ls='--', linewidth=0.5, c="black", zorder=4)
    axes.plot(MBH200, SFR_range, ls='--', linewidth=0.5, c="black", zorder=4)
    axes.plot(MBH500, SFR_range, ls='--', linewidth=0.5, c="black", zorder=4)

    loc = 22
    th500 = axes.text(MBH500[loc+4], SFR_range[loc+4], 'SFR/$\dot{M}_{BH}$=500', fontsize=12,
                      rotation=45, rotation_mode='anchor', transform_rotates_text=True)
    th200 = axes.text(MBH200[loc+2], SFR_range[loc+2], 'SFR/$\dot{M}_{BH}$=200', fontsize=12,
                      rotation=45, rotation_mode='anchor', transform_rotates_text=True)
    th50 = axes.text(MBH50[loc-1], SFR_range[loc-1], 'SFR/$\dot{M}_{BH}$=50', fontsize=12,
                     rotation=45, rotation_mode='anchor', transform_rotates_text=True)

    axes.legend()
    handles, labels = axes.get_legend_handles_labels()
    new_handles = []
    for handle, label in zip(handles, labels):
        if isinstance(handle, container.ErrorbarContainer):
            errorbar_marker = handle[0]
            # Create a new handle with only the marker
            new_handle = plt.Line2D([0], [0], marker=errorbar_marker.get_marker(),
                                    markersize=errorbar_marker.get_markersize(), label=label,markeredgecolor=errorbar_marker.get_markeredgecolor(),
                                    color=errorbar_marker.get_markerfacecolor(),linestyle=None)
            new_handles.append(new_handle)
        else:
            new_handles.append(handle)
    handles = new_handles
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    handles = [handles[-2],handles[2],handles[-1],handles[-3],handles[0],handles[1]]
    labels = [labels[-2],labels[2],labels[-1],labels[-3],labels[0],labels[1]]
    axes.legend(handles,labels)

    ax3.set_xlim(np.log10(ext.L_AGN_to_M_dot_BH(mn2)),np.log10(ext.L_AGN_to_M_dot_BH(mx2)))
    ax3.set_xlabel("$log_{10}$ $\dot{M}_{BH}$")

    plt.tight_layout()
    plt.savefig(f'{output_path}\log_SFR_Log_L_AGN.pdf',format='pdf')


if __name__ == '__main__':
    out_path = 'Replace with your output path'
    out_path = r'C:\Users\Tamir\Downloads'
    SFR_MBH_dot(out_path,type='GB')
