# Figure 4.5 in Thesis - Calculated mass growth trajectories
def Growth_Evolution(output_path,AGNs_df_path, rest_df_path,det = True):

    AGNs_df = pd.read_excel(agns_path)
    rest_df = pd.read_excel(rest_path)

    # plt.scatter(AGNs_df["log M_gal"],AGNs_df["log M_BH"],c="brown",marker='.',zorder=-10)
    plt.errorbar(rest_df["log M_gal"], rest_df["log M_BH"],rest_df["e_logM_BH"],ls="",marker='.',color='gray',zorder=-10,alpha=0.5)

    SFR_det = 740
    SFR_undet = 300

    BH_rate_det = 15.8
    BH_rate_undet = 18.3



    if not det:
        SFR_det = SFR_undet
        BH_rate_det = BH_rate_undet
    masses_gal_i = np.array([10,11,10,11],dtype=np.float64)
    masses_BH_i = np.array([9,9,9.5,9.5],dtype=np.float64)
    masses_BH_f = finalMass(M0=masses_BH_i, t=8, GR=BH_rate_det)
    masses_gal_f = finalMass(M0=masses_gal_i, t=8, GR=SFR_det)

    masses_BH_ff = finalMass(M0=masses_BH_i, t=9, GR=BH_rate_det)
    masses_gal_ff = finalMass(M0=masses_gal_i, t=9, GR=SFR_det)

    masses_BH_track = [finalMass(M0=masses_BH_i[j], t=np.linspace(0, 9, 1000), GR=BH_rate_det) for j in range(len(masses_BH_i))]
    masses_gal_track = [finalMass(M0=masses_gal_i[j], t=np.linspace(0,9,1000), GR=SFR_det) for j in range(len(masses_gal_i))]

    masses_BH_track2 = [finalMass(M0=masses_BH_i[j], t=np.linspace(0, 9, 1000), GR=BH_rate_undet) for j in range(len(masses_BH_i))]
    masses_gal_track2 = [finalMass(M0=masses_gal_i[j], t=np.linspace(0,9,1000), GR=SFR_undet) for j in range(len(masses_gal_i))]

    masses_BH_ticks = [finalMass(M0=masses_BH_i[j], t=np.array([np.log10(i*100*10**6) for i in range(2,10)]), GR=BH_rate_det) for j in range(len(masses_BH_i))]
    masses_gal_ticks = [finalMass(M0=masses_gal_i[j], t=np.array([np.log10(i*100*10**6) for i in range(2,10)]), GR=SFR_det) for j in range(len(masses_gal_i))]

    # image_path = get_sample_data('ada.png')

    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 600
    # fig.set_figheight(9)
    # fig.set_figwidth(9)

    plt.scatter(masses_gal_i,masses_BH_i,color='red',marker='*',s=200,label='t = 0')
    plt.scatter(masses_gal_f, masses_BH_f, color='tab:orange',marker='*',s=200,label='t = 100 Myr')
    plt.scatter(masses_gal_ff, masses_BH_ff, color='green', marker='*', s=200,label='t = 1 Gyr')

    colors = ['tab:blue','tab:green','tab:orange','tab:purple']
    for i in range(4):
        plt.plot(masses_gal_track[i],masses_BH_track[i],color=colors[i],zorder=-1)


        plt.plot(masses_gal_ticks[i], masses_BH_ticks[i], '|', markersize=5.5,color=colors[i],zorder=-1)

    # Get the current limits of the plot
    x_limits = plt.xlim()
    y_limits = plt.ylim()

    # Calculate endpoints of the diagonal line based on the current limits
    x_start, x_end = x_limits
    y_start = np.log10(10 ** x_start / 50)

    y_end = np.log10(10 ** x_end / 50)

    # Plot the diagonal line
    plt.plot([x_start, x_end], [y_start, y_end], '--', color='black',lw=1,zorder=-100)
    y_start = np.log10(10 ** x_start / 15)
    y_end = np.log10(10 ** x_end / 15)

    # Plot the diagonal line
    plt.plot([x_start, x_end], [y_start, y_end], '-.', color='black',lw=1,zorder=-100)
    loc = 0

    for val in [200,1000,5000]:
        y_start = np.log10(10 ** x_start / val)

        y_end = np.log10(10 ** x_end / val)
        plt.plot([x_start, x_end], [y_start, y_end], '--', color='black',lw=1,zorder=-100)

        xvals = np.linspace(x_start, x_end, 50)
        vals = np.log10(10 ** xvals / val)

        th = plt.text(xvals[loc + 4], vals[loc + 4], '$M_{*}$/$M_{BH}$='+str(val), fontsize=10,
                        rotation=45, rotation_mode='anchor', transform_rotates_text=True)

    xvals = np.linspace(x_start,x_end,50)
    yvals_50 = np.log10(10 ** xvals / 50)
    yvals_15 = np.log10(10 ** xvals / 15)
    th15 = plt.text(xvals[loc + 4], yvals_15[loc + 4], '$M_{*}$/$M_{BH}$=15', fontsize=10,
                      rotation=45, rotation_mode='anchor', transform_rotates_text=True)
    th50 = plt.text(xvals[loc + 4], yvals_50[loc + 4], '$M_{*}$/$M_{BH}$=50', fontsize=10,
                      rotation=45, rotation_mode='anchor', transform_rotates_text=True)

    plt.xlim(x_limits)
    plt.ylim(y_limits)

    plt.xlabel("$log_{10}$ ($M_{*}$/$M_{\u2609}$)")
    plt.ylabel("$log_{10}$ ($M_{BH}$/$M_{\u2609}$)")
    plt.legend()

    plt.grid()
    plt.tight_layout()
    if det:
        plt.savefig(f'{output_path}\Growth_Evolution_detected.pdf',format='pdf')
    else:
        plt.savefig(f'{output_path}\Growth_Evolution_undetected.pdf',format='pdf')


if __name__ == '__main__':
    out_path == 'Replace with your output path'
    reinesPath = r"D:\אוניברסיטה\אסטרופיזיקה\ALMA-Luminous-AGNs-at-z-2-3.5--\Data\Reines"
    AGNs_df_path = reinesPath + "\\table1.xlsx"
    rest_df_path = reinesPath + "\\table3.xlsx"
    SFR_Mdot_z(output_path,AGNs_df_path, rest_df_path)
