import io
import math
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
from scipy.stats import norm
from sklearn.metrics import confusion_matrix

matplotlib.use('pdf')


def plot_efficiency_significance(mode, tsd, significance, efficiency, data_range_array):
    fig, ax1 = plt.subplots()

    color = 'tab:blue'

    ax1.set_xlabel('BDT Score')
    ax1.set_ylabel('Significance', color=color)
    ax1.plot(tsd, significance, color=color)
    ax1.tick_params(axis='y', labelcolor=color, direction='in')

    ax2 = ax1.twinx()

    color = 'tab:red'
    # we already handled the x-label with ax1
    ax2.set_ylabel('BDT efficiency', color=color)
    ax2.plot(tsd, efficiency, color=color)
    ax2.tick_params(axis='y', labelcolor=color, direction='in')

    fig.tight_layout()

    fig_eff_path = os.environ['HYPERML_FIGURES_{}'.format(
        mode)]+'/Significance'
    if not os.path.exists(fig_eff_path):
        os.makedirs(fig_eff_path)

    fig_name = '/sign_eff_ct{}{}_pT{}{}_cen{}{}.pdf'.format(
        data_range_array[0],
        data_range_array[1],
        data_range_array[2],
        data_range_array[3],
        data_range_array[4],
        data_range_array[5])
    plt.savefig(fig_eff_path + fig_name)
    plt.close()


def plot_significance_scan(
        max_index, significance, significance_error, expected_signal, bkg_df, score_list, data_range_array,
        n_ev, mode, split='', mass_bins=40, hist_range = [2.96,3.04]):

    label = 'Significance x Efficiency'

    raw_yield = expected_signal[max_index]
    max_score = score_list[max_index]

    selected_bkg = bkg_df.query('score>@max_score')

    bkg_counts, bins = np.histogram(
        selected_bkg['m'], bins=mass_bins, range=hist_range)

    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    signal_counts_norm = norm.pdf(bin_centers, loc=2.992, scale=0.0025)
    signal_counts = raw_yield * signal_counts_norm / sum(signal_counts_norm)

    side_map = (bin_centers < 2.98) + (bin_centers > 3.005)
    bins_side = bin_centers[side_map]
    mass_map = np.logical_not(side_map)

    bkg_side_counts = bkg_counts[side_map]

    bkg_roi_shape = np.polyfit(bins_side, bkg_side_counts, 2)
    bkg_roi_counts = np.polyval(bkg_roi_shape, bin_centers)

    tot_counts = (bkg_roi_counts + signal_counts)[mass_map]

    fig, axs = plt.subplots(1, 2, figsize=(12, 4))

    axs[0].set_xlabel('Score')
    axs[0].set_ylabel(label)
    axs[0].tick_params(axis='x', direction='in')
    axs[0].tick_params(axis='y', direction='in')
    axs[0].plot(score_list, significance, 'b',
                label='Expected {}'.format(label))

    significance = np.asarray(significance)
    significance_error = np.asarray(significance_error)

    low_limit = significance - significance_error
    up_limit = significance + significance_error

    axs[0].fill_between(score_list, low_limit, up_limit,
                        facecolor='deepskyblue', label=r'$ \pm 1\sigma$')
    axs[0].grid()
    axs[0].legend(loc='upper left')

    bkg_side_error = np.sqrt(bkg_side_counts)
    tot_counts_error = np.sqrt(np.absolute(tot_counts))

    bins_mass = bin_centers[mass_map]

    axs[1].errorbar(bins_side, bkg_side_counts, yerr=bkg_side_error,
                    fmt='.', ecolor='k', color='b', elinewidth=1., label='Data')
    axs[1].errorbar(bins_mass, tot_counts, yerr=tot_counts_error,
                    fmt='.', ecolor='k', color='r', elinewidth=1., label='Pseudodata')
    axs[1].plot(bin_centers, bkg_roi_counts, 'g-', label='Background fit')

    x = np.linspace(2.9923 - 3 * 0.0025, 2.9923 + 3 * 0.0025, 1000)
    gauss_signal_counts = norm.pdf(x, loc=2.992, scale=0.0025)
    gauss_signal_counts = (raw_yield / sum(signal_counts_norm)) * \
        gauss_signal_counts + np.polyval(bkg_roi_shape, x)

    axs[1].plot(x, gauss_signal_counts, 'y', color='orange',
                label='Signal model (Gauss)')
    axs[1].set_xlabel(r'$m_{\ ^{3}He+\pi^{-}}$')
    axs[1].set_ylabel(r'Events /  ${}\ \rm{{MeV}}/c^{{2}}$'.format(mass_bins))
    axs[1].tick_params(axis='x', direction='in')
    axs[1].tick_params(axis='y', direction='in')
    axs[1].legend(loc='best', frameon=False)
    plt.ylim(bottom=0)

    s = sum(tot_counts) - sum(bkg_roi_counts[mass_map])
    b = sum(bkg_roi_counts[mass_map])

    sign_score = s / np.sqrt(s + b)

    plt.suptitle(r'%1.f$\leq ct \leq$%1.f %1.f$\leq \rm{p}_{T} \leq$%1.f  Cut Score=%0.2f  Significance=%0.2f  Raw yield=%0.2f' % (
        data_range_array[0], data_range_array[1], data_range_array[2], data_range_array[3], max_score,  sign_score, raw_yield))

    # text = '\n'.join(
    #     r'%1.f GeV/c $ \leq \rm{p}_{T} < $ %1.f GeV/c ' % (data_range_array[0], data_range_array[1]),
    #     r' Significance/Sqrt(Events) = %0.4f$x10^{-4}$' % (max_significance / np.sqrt(n_ev) * 1e4))

    # props = dict(boxstyle='round', facecolor='white', alpha=0)

    # axs[1].text(0.37, 0.95, text, transform=axs[1].transAxes, verticalalignment='top', bbox=props)

    fig_name = 'Significance_ct{}{}_pT{}{}_cen{}{}{}.pdf'.format(
        data_range_array[0],
        data_range_array[1],
        data_range_array[2],
        data_range_array[3],
        data_range_array[4],
        data_range_array[5],
        split)

    fig_sig_path = os.environ['HYPERML_FIGURES_{}'.format(
        mode)]+'/Significance'
    if not os.path.exists(fig_sig_path):
        os.makedirs(fig_sig_path)

    plt.savefig(fig_sig_path + '/' + fig_name)
    plt.close()


def plot_confusion_matrix(y_true, df, mode, score,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues, fig_name='confusion.pdf'):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # if the score is closer to max then to min it's recognised as signal
    y_pred = [1 if i > score else 0 for i in df['score']]

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    classes = ['Background', 'Signal']
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.

    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")

    fig.tight_layout()

    fig_sig_path = os.environ['HYPERML_FIGURES_{}'.format(mode)]+'/Confusion'
    if not os.path.exists(fig_sig_path):
        os.makedirs(fig_sig_path)

    plt.savefig(fig_sig_path + '/' + fig_name)
    plt.close()

    return ax
