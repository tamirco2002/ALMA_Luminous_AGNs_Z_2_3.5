import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import Extensions as ext


# Figure 10  - SFR/ Mbhdot for the most luminous quasars at $z\lesssim3.5$
def SFR_Mdot_z(output_path,data_path):
    filenames = ['cohen24','netzer16', 'schultz19','Rosario13','schweizer07']
    markers = ['o', 'h' ,'s', 'D','*', 'D']
    colors = ['red','black', 'green','purple', 'lightskyblue', 'orange']
    legend = ['This Work','Netzer+16', 'Schulze+19','Rosario+13','Schweizer+07']

    font = {'family': 'normal',
            'size': 16}

    matplotlib.rc('font', **font)

    fig, ax = plt.subplots(dpi=600)
    fig.set_figheight(9)
    fig.set_figwidth(12)
    all_vals = []
    for i in range(len(filenames)):
        df = pd.read_csv(data_path+filenames[i]+".csv")
        ms = [8 if mark not in ['*','o','D'] else 10 if mark in ['o','D'] else 13 for mark in markers]
        # mfc = [colors[i] if x==1 and markers[i]!='D' else 'none' for x in df['ul']]
        mfc = [colors[i] if markers[i]!='D' else 'none' for x in df['ul']]
        SFR = ext.L_SF_to_SFR(df['logL_ir'])
        M_BH_dot = ext.L_AGN_to_M_dot_BH(df['logL_agn'])
        print(filenames[i]+"\n"+str(np.median(SFR/M_BH_dot)))
        uplim = [False if x==1 else True for x in df['ul']]
        a0 = 0.5
        yerr = [0 if x==1 else a0*np.sinh(np.arcsinh(SFR[j]/M_BH_dot[j]/a0)-0.5/a0) for j,x in enumerate(df['ul'])]
        for j in range(len(df['z'])):
            if j == 0:
                ax.errorbar(df['z'][j],SFR[j]/M_BH_dot[j],yerr=yerr[j],marker=markers[i],ls='',color=colors[i],label=legend[i],
                            zorder=-i,ms=ms[i],uplims=uplim[j],mfc=mfc[j])
            else:
                ax.errorbar(df['z'][j],SFR[j]/M_BH_dot[j],yerr=yerr[j],marker=markers[i],ls='',color=colors[i],
                            zorder=-i,ms=ms[i],uplims=uplim[j],mfc=mfc[j])
            all_vals.append(SFR[j]/M_BH_dot[j])

    all_vals = np.array(all_vals)
    print(np.nanstd(all_vals))
    ax.set_xlabel('z')
    ax.set_ylabel('SFR/$\dot{M}_{BH}$')
    ax.axhline(y = 200,lw=3,ls='--',c='black')
    ax.set_yscale("log")
    from matplotlib.ticker import ScalarFormatter
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_path}\SFR_Mdot_z.pdf',format='pdf')

if __name__ == '__main__':
    out_path = 'Replace with your output path'
    data_folder_path = '../data/SFR_Mdot_Z_Figure/'
    SFR_Mdot_z(out_path,data_folder_path)