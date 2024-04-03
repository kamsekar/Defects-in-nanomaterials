# ana.rebeka.kamsek@ki.si, 2023

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib


def bar_plots(mass, specific, ecsa, labels=None, colors=None):
    """Plots bar plots for mass and specific activities and ECSA, including error bars.

    :param mass: mass activities in mA/g_Pt (no_of_samples, no_of_measurements)
    :param specific: specific activities in mA/cm^2 (no_of_samples, no_of_measurements)
    :param ecsa: electrochemically active surface area in m^2/g_Pt (no_of_samples, no_of_measurements)
    :param labels: optional list with labels for all samples
    :param colors: optional list of colors for all samples
    """

    no_of_samples = mass.shape[0]

    if labels is None:
        labels = ["sample no." + str(i + 1) for i in range(no_of_samples)]
    if colors is None:
        color_options = ["lightsalmon", "indianred", "silver", "maroon"]
        colors = []
        for i in range(no_of_samples):
            colors.append(color_options[i % len(color_options)])

    mass_avg = np.average(mass, axis=1)
    specific_avg = np.average(specific, axis=1)
    ecsa_avg = np.average(ecsa, axis=1)

    mass_std = np.std(mass, axis=1)
    specific_std = np.std(specific, axis=1)
    ecsa_std = np.std(ecsa, axis=1)

    plt.figure(figsize=(12, 4))
    plt.rcParams.update({'font.size': 14})

    plt.subplot(131)
    plt.bar(np.arange(1, 5), mass_avg, yerr=mass_std, color=colors, alpha=0.8)
    plt.ylabel(r"${j_m}$ / mA g$^{-1}_{Pt}$")
    plt.xticks((1, 2, 3, 4), labels=labels)
    plt.gca().tick_params(axis='both', which='major', labelsize=11)

    plt.subplot(132)
    plt.bar(np.arange(1, 5), specific_avg, yerr=specific_std, color=colors, alpha=0.8)
    plt.ylabel(r"${j_k}$ / mA cm$^{-2}}$")
    plt.xticks((1, 2, 3, 4), labels=labels)
    plt.gca().tick_params(axis='both', which='major', labelsize=11)

    plt.subplot(133)
    plt.bar(np.arange(1, 5), ecsa_avg, yerr=ecsa_std, color=colors, alpha=0.8)
    plt.ylabel(r"${ECSA_{CO}}$ / m$^2$ g$^{-1}_{Pt}}$")
    plt.xticks((1, 2, 3, 4), labels=labels)
    plt.gca().tick_params(axis='both', which='major', labelsize=11)

    plt.tight_layout()
    plt.show()


def scatter_plot(activity, disorder, colors=None):
    """A scatter plot for specific activities and degrees of disorder including error bars.

    :param activity: array with specific activities in mA/cm^2 and their standard deviations (no_of_samples, 2)
    :param disorder: array with fractions of the disordered phase and their st. devs. (no_of_samples, 2)
    :param colors: optional list of colors for all samples
    """

    if colors is None:
        color_options = ["lightsalmon", "indianred", "silver", "maroon"]
        no_of_samples = activity.shape[0]
        colors = []
        for i in range(no_of_samples):
            colors.append(color_options[i % len(color_options)])

    plt.figure(figsize=(5, 4))
    plt.rcParams.update({'font.size': 14})

    plt.subplot(111)
    plt.scatter(disorder[:, 0], activity[:, 0], c=colors, alpha=0.8, s=15)
    plt.errorbar(disorder[:, 0], activity[:, 0], xerr=disorder[:, 1], yerr=activity[:, 1],
                 fmt='none', ecolor=colors, alpha=0.8)

    plt.ylabel(r"${j_k}$ / mA cm$^{-2}}$")
    plt.xlabel(r"fraction of disordered phase")
    plt.ylim(bottom=0)
    plt.xlim(-0.05, 1.05)

    plt.gca().tick_params(axis='both', which='major', labelsize=11)
    plt.tight_layout()
    plt.show()


def rde_measurement_results(orr, tafel, co, labels=None, colors=None):
    """Plots ORR and CO stripping measurements as well as a calculated Tafel plot.

    Assumes data points in potential-current measurements are spaced 0.001 V apart.
    :param orr: currents in mA from >=two ORR cycles for all measurements (no_of_measurements, data_points, 2)
    :param tafel: calculated Tafel data for all measurements (no_of_measurements, tafel_length, 2)
    :param co: currents in mA from two CO stripping cycles for all measurements (no_of_measurements, data_points, 2)
    :param labels: optional list with labels for all samples
    :param colors: optional list of colors for all samples
    """

    no_of_measurements = orr.shape[0]
    if labels is None:
        labels = ["sample no." + str(i + 1) for i in range(no_of_measurements)]
    if colors is None:
        color_options = ["lightsalmon", "indianred", "silver", "maroon"]
        colors = []
        for i in range(no_of_measurements):
            colors.append(color_options[i % len(color_options)])

    # divide the ORR currents by electrode surface area in cm^2 to get current densities
    orr[:, :, 1] /= 0.196

    # take only one half of the second cycle for ORR measurements and subtract the background values
    lower_potential_limit = 0.05  # V
    upper_potential_limit = 1  # V
    half_cycle = int((upper_potential_limit - lower_potential_limit) * 1000)
    orr_bg = orr[:, 2 * half_cycle:3 * half_cycle, :] - co[:, 2 * half_cycle:3 * half_cycle, :]

    plt.figure(figsize=(12, 4))

    plt.subplot(131)
    for i in range(no_of_measurements):
        plt.plot(orr_bg[i, :, 0], orr_bg[i, :, 1], alpha=0.8, linewidth=1, color=colors[i], label=labels[i])

    plt.rcParams.update({'font.size': 14})
    plt.xlabel("E vs RHE [V]")
    plt.ylabel(r"j [$\mathrm{mA/cm^2_{geo}}$]")
    plt.legend(fontsize=11)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.subplot(132)
    for i in range(no_of_measurements):
        plt.plot(tafel[i, :, 0], tafel[i, :, 1], alpha=0.8, linewidth=1, color=colors[i], label=labels[i])

    plt.ylim(0.9, 1.0)
    plt.xscale('log', base=10)
    plt.xlim(0.03, 6)

    plt.rcParams.update({'font.size': 14})
    plt.xlabel(r"log ($\mathrm{j_k}$) [$\mathrm{mA/cm^2}$]")
    plt.ylabel("E vs RHE [V]")
    plt.yticks((0.90, 0.92, 0.94, 0.96, 0.98, 1.00), fontsize=11)
    plt.xticks(fontsize=12)

    plt.subplot(133)
    for i in range(no_of_measurements):
        plt.plot(co[i, :, 0], co[i, :, 1], alpha=0.8, linewidth=1, color=colors[i], label=labels[i])

    plt.rcParams.update({'font.size': 14})
    plt.xlabel("E vs RHE [V]")
    plt.ylabel("i [mA]")
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)

    plt.tight_layout()
    plt.show()


def plot_diffraction_patterns(patterns, labels=None, offset=0, empty_legend=False):
    """Plots a number of diffraction patterns with optional labels and y-offsetting.

    :param patterns: data array with diffraction patterns (no_of_patterns, no_of_data_points, 2)
    :param labels: list with labels corresponding to patterns
    :param offset: optional margin between diffraction patterns on the plot
    :param empty_legend: option to display a legend even when no specific labels are passed
    """

    matplotlib.rcParams.update({'font.size': 14})
    colormap = cm.get_cmap('hot', 10)

    if labels is None:
        labels = ["structure no." + str(i + 1) for i in range(patterns.shape[0])]

    plt.figure()
    for i in range(patterns.shape[0]):
        plt.plot(patterns[i, :, 0], patterns[i, :, 1] + offset * i, color=colormap(i * 2), label=labels[i])

    plt.xlim(20, 90)
    plt.ylim(0, np.amax(patterns))

    plt.xlabel(r'$2\theta\ [°]$')
    plt.ylabel("intensity (a.u.)")
    if labels is not None or empty_legend:
        plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()


def exp_fit_difference(exp, fit, offset=5000):
    """Plots an experimental diffraction pattern and its numerical fit, and their difference.

    Assumes patterns with 2theta values spanning at least over (20°, 140°),
    requires xlim and xticks corrections otherwise.
    :param exp: experimental diffraction pattern (no_of_data_points, 2)
    :param fit: fit to the experimental pattern (no_of_data_points, 2)
    :param offset: margin between diffraction patterns and their difference on the plot
    """

    matplotlib.rcParams.update({'font.size': 12})

    plt.figure()
    plt.plot(exp[:, 0], exp[:, 1], color='steelblue', linewidth=1.2, alpha=0.7)
    plt.plot(exp[:, 0], fit[:, 1], color='darkred', linewidth=0.6, alpha=0.8, ls='dotted')
    plt.plot(exp[:, 0], exp[:, 1] - fit[:, 1] - offset, color='gray', alpha=0.8, linewidth=0.8)

    plt.xlim(19, 141)
    plt.xticks((20, 40, 60, 80, 100, 120, 140))
    plt.yticks([])

    plt.xlabel(r'$2\theta\ [°]$')
    plt.ylabel("intensity (a.u.)")
    plt.tight_layout()
    plt.show()


def temp_profile(profile):
    """Plots the temperature-time profile for a HT-XRD measurement.

    :param profile: array (data, 3) with times, temperatures, and measurement indicators (1 for recording, 0 otherwise)
    """

    # extract the points from the profile where diffraction patterns are recorded, indicated with ones
    points = np.where(profile[:, 2] == 1)

    plt.figure()
    plt.plot(profile[:, 1], profile[:, 0], c='black', alpha=0.5, label='temperature profile')
    plt.scatter(profile[points, 1] + 0.083, profile[points, 0], s=15, c='black', label='XRD measurement')
    plt.xlabel("time [h]")
    plt.ylabel("temperature [°C]")
    plt.legend()
    plt.show()


def heatmap_plotting(patterns, temperatures, smoothing=False):
    """Plots sqrt values of in-situ HT-XRD measurements as a heatmap with corresponding temperature info.

    Assumes patterns with 2theta values spanning at least over (20°, 60°),
    requires xlim and xticks corrections otherwise.
    :param patterns: diffraction data (no_of_scans, no_of_data_points, 2)
    :param temperatures: a list with temperatures in °C for each diffraction pattern
    :param smoothing: optional convolve filter on the data to make it look more smooth
    """

    patterns = np.abs(patterns)

    no_of_scans = patterns.shape[0]
    x = patterns[0, :, 0]
    diffs = patterns[:, :, 1]

    temperatures = np.array(temperatures)
    # flip the order of given temperatures because a heatmap is constructed from the ground up
    temperatures = np.flip(temperatures)

    matplotlib.rcParams.update({'font.size': 13})

    plt.figure()
    plt.subplot(111)
    if smoothing:
        N = 3
        end_shape = (diffs.shape[0], diffs.shape[1] - N + 1)
        smooth_diffs = np.zeros(end_shape)

        for i in range(no_of_scans):
            smooth_diffs[i] = np.convolve(diffs[i], np.ones(N) / N, mode='valid')

        plt.imshow(np.sqrt(smooth_diffs), cmap='plasma', aspect='auto', extent=[min(x), max(x), 0, no_of_scans])
    else:
        plt.imshow(np.sqrt(diffs), cmap='plasma', aspect='auto', extent=[min(x), max(x), 0, no_of_scans])

    plt.ylabel("temperature [°C]")
    plt.yticks(ticks=np.arange(0.5, 0.5 + no_of_scans), labels=temperatures)
    plt.xlabel(r'$2\theta\ [°]$')
    plt.xticks((20, 40, 60))

    # clip the color values for better visualization
    plt.colorbar()
    plt.clim(np.amin(np.sqrt(diffs)), 0.7 * np.amax(np.sqrt(diffs)))
    plt.show()
