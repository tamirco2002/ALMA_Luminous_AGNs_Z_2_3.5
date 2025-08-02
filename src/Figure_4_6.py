# Figure 4.6 (top) in Thesis - Comparison with literature samples
def compapre_LAGN_LSFR(output_path,data_folder_path):
    filenames = ['f6_cohen24','f6_netzer16','f6_schultz19','f6_duras17_newy','f6_scholtz18_update','f6_stanley15','f6_stanley17']
    markers = ['o','o','s','*','h','D','D']
    colors = ['black','black','lightskyblue','green','purple','red','orange']
    legend = ['This Work','Netzer+16','Schulze+19','Duras+17','Scholtz+18','Stanley+15','Stanley+17']

    font = {'family': 'normal','size': 16}

    matplotlib.rc('font', **font)

    fig, axes = plt.subplots(dpi=600)
    fig.set_figheight(9)
    fig.set_figwidth(12)

    ax2 = axes.twinx()

    for i in range(len(filenames)):
        df = pd.read_csv(data_folder_path+filenames[i]+".csv")
        ms = [6 if mark!='*' and mark!='o' else 7.5 if mark=='o' else 9 for mark in markers]


        mfc = [colors[i] if x==1 and markers[i]!='D' else 'none' for x in df['ul']]
        if i ==0:
            mfc = [colors[i] for x in df['ul']]
        if i ==1:
            mfc = ['none' for x in df['ul']]
            ms = [5.5 for mark in markers]

        L_sf = np.log10(ext.L_SF_to_SFR(df['logL_ir']))
        yrr = np.array([df['logL_ir_elow'],df['logL_ir_eup']]).astype(np.float32)
        uplim = [False if x==1 else True for x in df['ul']]
        for j in range(len(L_sf)):
            if uplim[j] == True:
                yrr[0][j] = 0.15
                yrr[1][j] = 0
            if j==0:
                ax2.errorbar(df['logL_agn'][j],df['logL_ir'][j],ls='',yerr=[[yrr[0][j]],[yrr[1][j]]],marker=markers[i],
                             color=colors[i],label=legend[i],mfc=mfc[j],zorder=-i,ms=ms[i],uplims=uplim[j],elinewidth=0.5)
            else:
                ax2.errorbar(df['logL_agn'][j], df['logL_ir'][j], ls='', yerr=[[yrr[0][j]],[yrr[1][j]]], marker=markers[i],
                             color=colors[i], mfc=mfc[j],zorder=-i,ms=ms[i],uplims=uplim[j],elinewidth=0.5,capthick=0.5)

    SFR_range = np.linspace(-1,5)
    MBH50 = ext.M_dot_BH_to_L_AGN(10**SFR_range/50)
    MBH200 = ext.M_dot_BH_to_L_AGN(10 ** SFR_range / 200)
    MBH500 = ext.M_dot_BH_to_L_AGN(10 ** SFR_range / 500)
    axes.plot(MBH50,SFR_range,ls='--',linewidth=0.5,c="black",zorder=4)
    axes.plot(MBH200, SFR_range, ls='--',linewidth=0.5,c="black",zorder=4)
    axes.plot(MBH500, SFR_range, ls='--',linewidth=0.5,c="black",zorder=4)

    xmean = sum(i for i in MBH500) / float(len(MBH500))
    ymean = sum(i for i in SFR_range) / float(len(SFR_range))

    th500 = axes.text(MBH500[10], SFR_range[10], 'SFR/$\dot{M}_{BH}$=500', fontsize=12,
                  rotation=45, rotation_mode='anchor',transform_rotates_text=True)
    th200 = axes.text(MBH200[10], SFR_range[10], 'SFR/$\dot{M}_{BH}$=200', fontsize=12,
                  rotation=45, rotation_mode='anchor',transform_rotates_text=True)
    th50 = axes.text(MBH50[9], SFR_range[9], 'SFR/$\dot{M}_{BH}$=50', fontsize=12,
                  rotation=45, rotation_mode='anchor',transform_rotates_text=True)

    axes.set_ylabel("$log_{10}$ SFR")
    axes.set_xlabel("$log_{10}$ $L_{AGN}$")

    axes.set_xlim(43, 48.2)
    ax2.set_ylim(ext.SFR_to_L_SF(10 ** -0.5), ext.SFR_to_L_SF(10 **4.4))

    mn, mx = ax2.get_ylim()

    axes.set_ylim(np.log10(ext.L_SF_to_SFR(mn)), np.log10(ext.L_SF_to_SFR(mx)))

    ax2.set_ylabel("$log_{10}$ $L_{SF}$")

    ax3 = axes.twiny()

    mn2, mx2 = axes.get_xlim()
    trendx = np.arange(np.log10(ext.L_AGN_to_M_dot_BH(mn2)), np.log10(ext.L_AGN_to_M_dot_BH(mx2)), step=0.01)
    mn3, mx3 = ax2.get_ylim()
    ax3.set_xlim(np.log10(ext.L_AGN_to_M_dot_BH(mn2)),np.log10(ext.L_AGN_to_M_dot_BH(mx2)))
    ax3.set_xlabel("$log_{10}$ $\dot{M}_{BH}$")

    # get handles
    handles, labels = ax2.get_legend_handles_labels()
    # remove the errorbars
    handles = [h[0] for h in handles]
    # use them in the legend
    ax2.legend(handles, labels, numpoints=1,loc=2)
    plt.tight_layout()

    plt.savefig(f'{output_path}\SFR_compare.pdf',format='pdf')

# Figure 4.6 (bottom) in Thesis - Comparison with literature samples
def SFR_Evolution(output_path,data_folder_path):
    filenames = ['f7_cohen24','f7_netzer16','f7_schulze19','f7_duras17','f7_scholtz18','f7_hatziminaoglou18','f7_nguyen20','f7_priddey03']
    markers = ['o','o','s','*','h','^','D','d']
    colors = ['black','black','lightskyblue','green','purple','red','orange','hotpink']
    legend = ['This Work','Netzer+16','Schulze+19','Duras+17','Scholtz+18','Hatziminaoglou+18','Nguyen+20','Priddey+03']

    font = {'family': 'normal',
            'size': 16}

    matplotlib.rc('font', **font)

    fig, axes = plt.subplots(dpi=600)
    fig.set_figheight(9)
    fig.set_figwidth(12)

    for i in range(len(filenames)):
        df = pd.read_csv(data_folder_path+filenames[i]+".csv")
        ms = [6 if mark!='*' and mark!='o' else 7.5 if mark=='o' else 9 for mark in markers]
        mfc = [colors[i] if x==1 else 'none' for x in df['ul']]

        if i ==0:
            mfc = [colors[i] for x in df['ul']]
        if i ==1:
            mfc = ['none' for x in df['ul']]
            ms = [5.5 for mark in markers]

        logSFR = df['logSFR']
        z = df['z']
        try:
            error = np.array(df['logSFR_e']).astype(np.float32)
        except:
            error = np.array([df['logSFR_elow'],df['logSFR_eup']]).T
            new_error = []
            for x in range(len(error)):
                new_error.append([[error[x][0]],[error[x][1]]])
            error = np.array(new_error).astype(np.float32)

        uplim = [False if x == 1 else True for x in df['ul']]

        for j in range(len(logSFR)):
            if uplim[j] == True:
                error[j] = 0.15


            if j==0:
                axes.errorbar(z[j],logSFR[j],ls='',yerr=error[j],marker=markers[i],
                             color=colors[i],label=legend[i],mfc=mfc[j],zorder=-i,ms=ms[i],uplims=uplim[j],elinewidth=0.5)
            else:
                axes.errorbar(z[j],logSFR[j],ls='',yerr=error[j],marker=markers[i],
                             color=colors[i], mfc=mfc[j],zorder=-i,ms=ms[i],uplims=uplim[j],elinewidth=0.5)

    axes.set_ylabel("$log_{10}$ SFR")
    axes.set_xlabel("z")

    # get handles
    handles, labels = axes.get_legend_handles_labels()
    # remove the errorbars
    handles = [h[0] for h in handles]
    # use them in the legend
    axes.legend(handles, labels, numpoints=1,loc=4,bbox_to_anchor=(0.85,0))
    # plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig(f'{output_path}\SFR_Evolution.pdf',format='pdf')


if __name__ == '__main__':
    data_folder_path_f6_top = '../Data/schulze19_external_samples/figure6/'
    data_folder_path_f6_bottom = '../Data/schulze19_external_samples/figure7/'
    out_path == 'Replace with your output path'
    compapre_LAGN_LSFR(output_path,data_folder_path_f6_top)
    SFR_Evolution(output_path,data_folder_path)