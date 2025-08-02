import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src import ALMA_HERSCHEL_INPUT_PATH


# Figure 5 - Distribution of SFRs
def SFR_histogram(output_path,type='GB',GB_type='BB'):
    results_file_path = f'../data/Results/{GB_type}.xlsx'
    df = pd.read_excel(results_file_path)
    if type == 'GB':
        all_sources = [np.log10(df["GB-SFR[M_sun/year]"][s]) for s in range(len(df["SNR"]))]
        det_sources = [np.log10(df["GB-SFR[M_sun/year]"][s]) for s in range(len(df["SNR"])) if df["Alma group"][s]=='det' and np.abs(df["SNR"][s])>=3]
        undet_sources = [np.log10(df["GB-SFR[M_sun/year]"][s]) for s in range(len(df["SNR"])) if df["Alma group"][s]!='det' or np.abs(df["SNR"][s])<3]
    else:
        all_sources = [np.log10(df["CE-SFR[M_sun/year]"][s]) for s in range(len(df["SNR"]))]
        det_sources = [np.log10(df["CE-SFR[M_sun/year]"][s]) for s in range(len(df["SNR"])) if df["Alma group"][s]=='det' and np.abs(df["SNR"][s])>=3]
        undet_sources = [np.log10(df["CE-SFR[M_sun/year]"][s]) for s in range(len(df["SNR"])) if df["Alma group"][s]!='det' or np.abs(df["SNR"][s])<3]

    undet_sources = undet_sources[:-2]
    det_sources = det_sources[:-2]

    count, bins = np.histogram(all_sources, bins=[1.7,2,2.3,2.6,2.9,3.2,3.5,3.7])

    sorted_det= np.sort(det_sources)
    cdf_det = np.arange(1, len(sorted_det) + 1) / len(sorted_det)

    sorted_undet= np.sort(undet_sources)
    cdf_undet = np.arange(1, len(sorted_undet) + 1) / len(sorted_undet)

    fig, ax1 = plt.subplots(dpi=600)

    ax1.hist(det_sources, bins=bins, color='limegreen', label="Detected")
    ax1.hist(undet_sources, bins=bins, histtype='step', edgecolor='red', label="Undetected")
    ax1.set_xlabel("log SFR [$M_{\u2609}$/yr]")
    ax1.set_ylabel("Number of Sources")

    ax2 = ax1.twinx()
    ax2.set_ylabel('Cumulative Distribution Function')
    ax2.set_ylim(0, 1)


    # Normalize the CDF by the maximum value of the histogram counts scaled by the bin width
    ax2.plot(np.insert(sorted_det,0,1.7), np.insert(cdf_det,0,0), label="CDF", ls='--',lw=2.5, color='green',drawstyle='steps-post')
    ax2.plot(np.insert(sorted_undet,0,1.7), np.insert(cdf_undet,0,0), label="CDF", ls='--',lw=2.5 , color='red',drawstyle='steps-post')

    ax1.legend(loc='upper left')
    ax1.set_xlim(1.6,3.8)
    plt.tight_layout()
    plt.savefig(f'{output_path}\SFR_distribution.pdf',format='pdf')

if __name__ == '__main__':
    out_path = 'Replace with your output path'
    SFR_histogram(out_path,type='GB')
