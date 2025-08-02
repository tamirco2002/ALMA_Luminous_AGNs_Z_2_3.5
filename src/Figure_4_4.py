
# Figure 4.4 in Thesis - Distributions of SFR/ M_dot_BH
def SFR_Mdot_histogram(output_path,df_path,type='GB'):
    df = pd.read_excel(df_path)
    if type == 'GB':
        all_sources = [df["GB-SFR/M_dot_BH"][s] for s in range(len(df["SNR"]))]
        det_sources = [df["GB-SFR/M_dot_BH"][s] for s in range(len(df["SNR"])) if df["Alma group"][s]=='det' and np.abs(df["SNR"][s])>=3]
        undet_sources = [df["GB-SFR/M_dot_BH"][s] for s in range(len(df["SNR"])) if df["Alma group"][s]!='det' or np.abs(df["SNR"][s])<3]
    else:
        all_sources = [df["CE-SFR/M_dot_BH"][s] for s in range(len(df["SNR"]))]
        det_sources = [df["CE-SFR/M_dot_BH"][s] for s in range(len(df["SNR"])) if df["Alma group"][s]=='det' and np.abs(df["SNR"][s])>=3]
        undet_sources = [df["CE-SFR/M_dot_BH"][s] for s in range(len(df["SNR"])) if df["Alma group"][s]!='det' or np.abs(df["SNR"][s])<3]

    count, bins = np.histogram(det_sources, bins=12,range=(np.min(det_sources),np.max(det_sources)))
    count1, bins1 = np.histogram(undet_sources, bins=10,range=(np.min(undet_sources),np.max(undet_sources)))

    sorted_det= np.sort(det_sources)
    cdf_det = np.arange(1, len(sorted_det) + 1) / len(sorted_det)

    sorted_undet= np.sort(undet_sources)
    cdf_undet = np.arange(1, len(sorted_undet) + 1) / len(sorted_undet)

    font = {'family': 'normal',
            'size': 16}

    matplotlib.rc('font', **font)
    matplotlib.rcParams['xtick.major.size'] = 8
    matplotlib.rcParams['xtick.minor.size'] = 5
    matplotlib.rcParams['xtick.major.width'] = 3
    matplotlib.rcParams['xtick.minor.width'] = 1
    fig, (ax2, ax1) = plt.subplots(nrows=2,dpi=600,sharex=True)
    fig.set_figheight(12)
    fig.set_figwidth(10)

    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    logbins1 = np.logspace(np.log10(bins1[0]), np.log10(bins1[-1]), len(bins1))
    ax1.hist(det_sources, bins=logbins, histtype='step',facecolor='None', edgecolor='green', label="Detected",lw=4)
    undetHist = ax1.hist(undet_sources, bins=logbins1, histtype='step',facecolor='None', edgecolor='red', label="Undetected",lw=4)

    heights = [1,None,None,0.5,None,1.5,None,1.5,None,0.5]
    for bin_edge, h in zip(undetHist[1][1:], heights):
        head_length = (np.log10(bin_edge))**1.9 * 0.2  # Adjust the scaling factor as needed
        length = -(np.log10(bin_edge))**1.9 * 0.7
        if h != None:
            plt.arrow(bin_edge, h, length, 0, head_width=-0.05, head_length=head_length, fc='red', ec='red')

    ax1.set_xlabel("SFR/$\dot{M}_{BH}$")
    ax1.set_ylabel("Number of Sources")
    ax1.set_xscale('log')
    ax1.set_yticks(np.arange(0, 4.5, 1.0))

    ax2.set_xlim(2,)
    ax2.set_ylabel('Cumulative Distribution Function')
    ax2.set_ylim(0, 1)


    # Normalize the CDF by the maximum value of the histogram counts scaled by the bin width
    ax2.plot(np.insert(sorted_det,0,0), np.insert(cdf_det,0,0), label="Detected", ls='-',lw=3, color='green',drawstyle='steps-post')
    ax2.plot(np.insert(sorted_undet,0,0), np.insert(cdf_undet,0,0), label="Undetected", ls='-',lw=3, color='red',drawstyle='steps-post')


    x = np.insert(sorted_det,0,0)
    y = np.insert(cdf_det,0,0)
    index_at_x50 = np.argmin(np.abs(x - 50))
    index_at_x200 = np.argmin(np.abs(x - 200))

    ax2.axvline(x = 50,lw=3,ls='--',c='black',label="SFR/$\dot{M}_{BH}$ = 50")
    ax2.axvline(x = 200,lw=3,ls='-.',c='blue',label="SFR/$\dot{M}_{BH}$ = 200")

    ax1.set_xlim(0)
    ax1.xaxis.set_tick_params(width=3)
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{output_path}\SFR_Mdot_hist.pdf',format='pdf')


if __name__ == '__main__':
    out_path == 'Replace with your output path'
    df_path = pd.read_excel('..data/ALMA x Herschel CASA input - March23.xlsx')
    SFR_histogram(df_path,out_path,type='GB')