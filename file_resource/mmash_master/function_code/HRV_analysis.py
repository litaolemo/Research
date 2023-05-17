import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from matplotlib.patches import Ellipse
from scipy import interpolate

def _create_timestamp_list(nn_intervals):
    """
    Creates corresponding time interval for all nn_intervals
    Parameters
    ---------
    nn_intervals : list
        List of Normal to Normal Interval.
    Returns
    ---------
    nni_tmstp : list
        list of time intervals between first NN-interval and final NN-interval.
    """
    # Convert in seconds
    nni_tmstp = np.cumsum(nn_intervals) / 1000

    # Force to start at 0
    return nni_tmstp - nni_tmstp[0]


def _create_interpolated_timestamp_list(nn_intervals, sampling_frequency = 7):
    """
    Creates the interpolation time used for Fourier transform's method
    Parameters
    ---------
    nn_intervals : list
        List of Normal to Normal Interval.
    sampling_frequency : int
        Frequency at which the signal is sampled.
    Returns
    ---------
    nni_interpolation_tmstp : list
        Timestamp for interpolation.
    """
    time_nni = _create_timestamp_list(nn_intervals)
    # Create timestamp for interpolation
    nni_interpolation_tmstp = np.arange(0, time_nni[-1], 1 / float(sampling_frequency))
    return nni_interpolation_tmstp

from collections import namedtuple

WELCH_METHOD = "welch"
VlfBand = namedtuple("Vlf_band", ["low", "high"])
LfBand = namedtuple("Lf_band", ["low", "high"])
HfBand = namedtuple("Hf_band", ["low", "high"])


def _get_freq_psd_from_nn_intervals(nn_intervals, method = WELCH_METHOD,
                                    sampling_frequency = 4,
                                    interpolation_method = "linear",
                                    vlf_band = VlfBand(0.003, 0.04),
                                    hf_band = HfBand(0.15, 0.40)):
    """
    Returns the frequency and power of the signal.
    Parameters
    ---------
    nn_intervals : list
        list of Normal to Normal Interval
    method : str
        Method used to calculate the psd. Choice are Welch's FFT or Lomb method.
    sampling_frequency : int
        Frequency at which the signal is sampled. Common value range from 1 Hz to 10 Hz,
        by default set to 7 Hz. No need to specify if Lomb method is used.
    interpolation_method : str
        Kind of interpolation as a string, by default "linear". No need to specify if Lomb
        method is used.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.
    Returns
    ---------
    freq : list
        Frequency of the corresponding psd points.
    psd : list
        Power Spectral Density of the signal.
    """

    timestamp_list = _create_timestamp_list(nn_intervals)

    if method == WELCH_METHOD:
        # ---------- Interpolation of signal ---------- #
        funct = interpolate.interp1d(x=timestamp_list, y=nn_intervals, kind=interpolation_method)

        timestamps_interpolation = _create_interpolated_timestamp_list(nn_intervals, sampling_frequency)
        nni_interpolation = funct(timestamps_interpolation)

        # ---------- Remove DC Component ---------- #
        nni_normalized = nni_interpolation - np.mean(nni_interpolation)

        #  --------- Compute Power Spectral Density  --------- #
        freq, psd = signal.welch(x=nni_normalized, fs=sampling_frequency, window='hann',
                                 nfft=4096)

    elif method == LOMB_METHOD:
        freq, psd = LombScargle(timestamp_list, nn_intervals,
                                normalization='psd').autopower(minimum_frequency=vlf_band[0],
                                                               maximum_frequency=hf_band[1])
    else:
        raise ValueError("Not a valid method. Choose between 'lomb' and 'welch'")

    return freq, psd


def plot_psd(nn_intervals, method = "welch", sampling_frequency = 7,
             interpolation_method = "linear", vlf_band: namedtuple = VlfBand(0.003, 0.04),
             lf_band: namedtuple = LfBand(0.04, 0.15), hf_band: namedtuple = HfBand(0.15, 0.40),
             save_fig = False):
    """
    Function plotting the power spectral density of the NN Intervals.
    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval.
    method : str
        Method used to calculate the psd. Choice are Welch's FFT (welch) or Lomb method (lomb).
    sampling_frequency : int
        frequence at which the signal is sampled. Common value range from 1 Hz to 10 Hz, by default
        set to 7 Hz. No need to specify if Lomb method is used.
    interpolation_method : str
        kind of interpolation as a string, by default "linear". No need to specify if lomb method is
        used.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    lf_band : tuple
        Low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.
    """

    freq, psd = _get_freq_psd_from_nn_intervals(nn_intervals=nn_intervals, method=method,
                                                sampling_frequency=sampling_frequency,
                                                interpolation_method=interpolation_method)

    # Calcul of indices between desired frequency bands
    vlf_indexes = np.logical_and(freq >= vlf_band[0], freq < vlf_band[1])
    lf_indexes = np.logical_and(freq >= lf_band[0], freq < lf_band[1])
    hf_indexes = np.logical_and(freq >= hf_band[0], freq < hf_band[1])

    frequency_band_index = [vlf_indexes, lf_indexes, hf_indexes]
    label_list = ["VLF component", "LF component", "HF component"]

    # Plot parameters
    # sns.set_style("ticks")
    fig,ax = plt.subplots(figsize=(6,4))
    plt.xlabel("Frequency (Hz)", fontsize=15)
    plt.ylabel("PSD (s2/ Hz)", fontsize=15)

    if method == "lomb":
        plt.title("Lomb's periodogram", fontsize=20)
        for band_index, label in zip(frequency_band_index, label_list):
            plt.fill_between(freq[band_index], 0, psd[band_index] / (1000 * len(psd[band_index])), label=label)
        plt.legend(prop={"size": 15}, loc="best")

    elif method == "welch":
        plt.title("FFT Spectrum : Welch's periodogram", fontsize=20)
        for band_index, label in zip(frequency_band_index, label_list):
            plt.fill_between(freq[band_index], 0, psd[band_index] / (1000 * len(psd[band_index])), label=label)
        plt.legend(prop={"size": 15}, loc="best")
        plt.xlim(0, hf_band[1])
    else:
        raise ValueError("Not a valid method. Choose between 'lomb' and 'welch'")

    fig.tight_layout()

    if save_fig == True:
        fig.savefig('psd.pdf')

    plt.show()
    

def plot_poincare(nn_intervals, plot_sd_features = True, save_fig = False):
    """
    Pointcare / Lorentz Plot of the NN Intervals
    Arguments
    ---------
    nn_intervals : list
        list of NN intervals
    plot_sd_features : bool
        Option to show or not SD1 and SD2 features on plot. By default, set to True.
    Notes
    ---------
    The transverse axis (T) reflects beat-to-beat variation
    the longitudinal axis (L) reflects the overall fluctuation
    """

    # For Lorentz / poincaré Plot
    ax1 = nn_intervals[:-1]
    ax2 = nn_intervals[1:]

    # compute features for ellipse's height, width and center
    dict_sd1_sd2 = get_poincare_plot_features(nn_intervals)
    sd1 = dict_sd1_sd2["sd1"]
    sd2 = dict_sd1_sd2["sd2"]
    mean_nni = np.mean(nn_intervals)

    # Plot options and settings
    # sns.set_style("ticks")
    fig,ax = plt.subplots(figsize=(6,6))
    plt.title("Poincaré Plot", fontsize=20)
    plt.xlabel('NN_n (ms)', fontsize=15)
    plt.ylabel('NN_n+1 (ms)', fontsize=15)
    plt.xlim(min(nn_intervals) - 10, max(nn_intervals) + 10)
    plt.ylim(min(nn_intervals) - 10, max(nn_intervals) + 10)
    plt.plot([0,1500],[0,1500],'--',c='k')

    # Poincaré Plot
    ax.scatter(ax1, ax2, c='b', s=20, alpha=0.4)

    if plot_sd_features:
        # Ellipse plot settings
        ells = Ellipse(xy=(mean_nni, mean_nni), width=2 * sd2 + 1,
                       height=2 * sd1 + 1, angle=45, linewidth=3, edgecolor = 'k',
                       fill=False)
        ax.add_patch(ells)

        ells = Ellipse(xy=(mean_nni, mean_nni), width=2 * sd2,
                       height=2 * sd1, angle=45)
        ells.set_alpha(0.05)
        ells.set_facecolor("blue")
        ax.add_patch(ells)

        # Arrow plot settings
        sd1_arrow = ax.arrow(mean_nni, mean_nni, -sd1 * np.sqrt(2) / 2, sd1 * np.sqrt(2) / 2,
                             linewidth=4, ec='r', fc="r", label="SD1")
        sd2_arrow = ax.arrow(mean_nni, mean_nni, sd2 * np.sqrt(2) / 2, sd2 * np.sqrt(2) / 2,
                             linewidth=4, ec='orange', fc="orange", label="SD2")

        plt.legend(handles=[sd1_arrow, sd2_arrow], fontsize=12, loc="best")

    fig.tight_layout()

    if save_fig == True:
        fig.savefig('poincare.pdf')

    plt.show()
    
    
def get_poincare_plot_features(nn_intervals):
    """
    Function returning a dictionary containing 3 features from non linear domain
    for HRV analyses.
    Known practise is to use this function on short term recordings, from 5 minutes window.
    Parameters
    ---------
    nn_intervals : list
        Normal to Normal Interval
    Returns
    ---------
    poincare_plot_features : dict
        Dictionary containing non linear domain features for hrV analyses. There
        are details about each features are given below.
    Notes
    ---------
    - **sd1** : The standard deviation of projection of the Poincaré plot on the line \
    perpendicular to the line of identity.
    - **sd2** : SD2 is defined as the standard deviation of the projection of the Poincaré \
    plot on the line of identity (y=x).
    - **ratio_sd2_sd1** : Ratio between SD2 and SD1.
    References
    ----------
    .. [4] Pre-ictal heart rate variability assessment of epileptic seizures by means of linear \
    and non- linear analyses, Soroor Behbahani, Nader Jafarnia Dabanloo et al - 2013
    """
    diff_nn_intervals = np.diff(nn_intervals)
    # measures the width of poincare cloud
    sd1 = np.sqrt(np.std(diff_nn_intervals, ddof=1) ** 2 * 0.5)
    # measures the length of the poincare cloud
    sd2 = np.sqrt(2 * np.std(nn_intervals, ddof=1) ** 2 - 0.5 * np.std(diff_nn_intervals, ddof=1) ** 2)
    ratio_sd2_sd1 = sd2 / sd1

    poincare_plot_features = {
        'sd1': sd1,
        'sd2': sd2,
        'ratio_sd2_sd1': ratio_sd2_sd1
    }

    return poincare_plot_features


def get_time_domain_features(nn_intervals):
    """
    Returns a dictionary containing time domain features for HRV analysis.
    Mostly used on long term recordings (24h) but some studies use some of those features on
    short term recordings, from 1 to 5 minutes window.
    Parameters
    ----------
    nn_intervals : list
        list of Normal to Normal Interval
    Returns
    -------
    time_domain_features : dict
        dictionary containing time domain features for HRV analyses. There are details
        about each features below.
    Notes
    -----
    Here are some details about feature engineering...
    - **mean_nni**: The mean of RR-intervals.
    - **sdnn** : The standard deviation of the time interval between successive normal heart beats \
    (i.e. the RR-intervals).
    - **sdsd**: The standard deviation of differences between adjacent RR-intervals
    - **rmssd**: The square root of the mean of the sum of the squares of differences between \
    adjacent NN-intervals. Reflects high frequency (fast or parasympathetic) influences on hrV \
    (*i.e.*, those influencing larger changes from one beat to the next).
    - **median_nni**: Median Absolute values of the successive differences between the RR-intervals.
    - **nni_50**: Number of interval differences of successive RR-intervals greater than 50 ms.
    - **pnni_50**: The proportion derived by dividing nni_50 (The number of interval differences \
    of successive RR-intervals greater than 50 ms) by the total number of RR-intervals.
    - **nni_20**: Number of interval differences of successive RR-intervals greater than 20 ms.
    - **pnni_20**: The proportion derived by dividing nni_20 (The number of interval differences \
    of successive RR-intervals greater than 20 ms) by the total number of RR-intervals.
    - **range_nni**: difference between the maximum and minimum nn_interval.
    - **cvsd**: Coefficient of variation of successive differences equal to the rmssd divided by \
    mean_nni.
    - **cvnni**: Coefficient of variation equal to the ratio of sdnn divided by mean_nni.
    - **mean_hr**: The mean Heart Rate.
    - **max_hr**: Max heart rate.
    - **min_hr**: Min heart rate.
    - **std_hr**: Standard deviation of heart rate.
    References
    ----------
    .. [1] Heart rate variability - Standards of measurement, physiological interpretation, and \
    clinical use, Task Force of The European Society of Cardiology and The North American Society \
    of Pacing and Electrophysiology, 1996
    """

    diff_nni = np.diff(nn_intervals)
    length_int = len(nn_intervals)

    # Basic statistics
    mean_nni = np.mean(nn_intervals)
    median_nni = np.median(nn_intervals)
    range_nni = max(nn_intervals) - min(nn_intervals)

    sdsd = np.std(diff_nni)
    rmssd = np.sqrt(np.mean(diff_nni ** 2))

    nni_50 = sum(np.abs(diff_nni) > 50)
    pnni_50 = 100 * nni_50 / length_int


    # Features only for long term recordings
    sdnn = np.std(nn_intervals, ddof=1)  # ddof = 1 : unbiased estimator => divide std by n-1

    # Heart Rate equivalent features
    heart_rate_list = np.divide(60000, nn_intervals)
    mean_hr = np.mean(heart_rate_list)
    min_hr = min(heart_rate_list)
    max_hr = max(heart_rate_list)
    std_hr = np.std(heart_rate_list)

    time_domain_features = {
        'mean_nni': mean_nni,
        'sdnn': sdnn,
        'pnni_50': pnni_50,
        'rmssd': rmssd,
        'median_nni': median_nni,
        'range_nni': range_nni,
        'mean_hr': mean_hr,
        "max_hr": max_hr,
        "min_hr": min_hr,
        "std_hr": std_hr,
    }

    return time_domain_features


def get_frequency_domain_features(nn_intervals, method = WELCH_METHOD,
                                  sampling_frequency = 4, interpolation_method = "linear",
                                  vlf_band = VlfBand(0.003, 0.04),
                                  lf_band = LfBand(0.04, 0.15),
                                  hf_band = HfBand(0.15, 0.40)):
    """
    Returns a dictionary containing frequency domain features for HRV analyses.
    To our knowledge, you might use this function on short term recordings, from 2 to 5 minutes  \
    window.
    Parameters
    ---------
    nn_intervals : list
        list of Normal to Normal Interval
    method : str
        Method used to calculate the psd. Choice are Welch's FFT or Lomb method.
    sampling_frequency : int
        Frequency at which the signal is sampled. Common value range from 1 Hz to 10 Hz,
        by default set to 4 Hz. No need to specify if Lomb method is used.
    interpolation_method : str
        kind of interpolation as a string, by default "linear". No need to specify if Lomb
        method is used.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    lf_band : tuple
        Low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.
    Returns
    ---------
    frequency_domain_features : dict
        Dictionary containing frequency domain features for HRV analyses. There are details
        about each features below.
    Notes
    ---------
    Details about feature engineering...
    - **total_power** : Total power density spectral
    - **vlf** : variance ( = power ) in HRV in the Very low Frequency (.003 to .04 Hz by default). \
    Reflect an intrinsic rhythm produced by the heart which is modulated primarily by sympathetic \
    activity.
    - **lf** : variance ( = power ) in HRV in the low Frequency (.04 to .15 Hz). Reflects a \
    mixture of sympathetic and parasympathetic activity, but in long-term recordings, it reflects \
    sympathetic activity and can be reduced by the beta-adrenergic antagonist propanolol.
    - **hf**: variance ( = power ) in HRV in the High Frequency (.15 to .40 Hz by default). \
    Reflects fast changes in beat-to-beat variability due to parasympathetic (vagal) activity. \
    Sometimes called the respiratory band because it corresponds to HRV changes related to the \
    respiratory cycle and can be increased by slow, deep breathing (about 6 or 7 breaths per \
    minute) and decreased by anticholinergic drugs or vagal blockade.
    - **lf_hf_ratio** : lf/hf ratio is sometimes used by some investigators as a quantitative \
    mirror of the sympatho/vagal balance.
    - **lfnu** : normalized lf power.
    - **hfnu** : normalized hf power.
    References
    ----------
    .. [1] Heart rate variability - Standards of measurement, physiological interpretation, and \
    clinical use, Task Force of The European Society of Cardiology and The North American Society \
    of Pacing and Electrophysiology, 1996
    .. [2] Signal Processing Methods for Heart Rate Variability - Gari D. Clifford, 2002
    """

    # ----------  Compute frequency & Power spectral density of signal  ---------- #
    freq, psd = _get_freq_psd_from_nn_intervals(nn_intervals=nn_intervals, method=method,
                                                sampling_frequency=sampling_frequency,
                                                interpolation_method=interpolation_method,
                                                vlf_band=vlf_band, hf_band=hf_band)

    # ---------- Features calculation ---------- #
    freqency_domain_features = _get_features_from_psd(freq=freq, psd=psd,
                                                      vlf_band=vlf_band,
                                                      lf_band=lf_band,
                                                      hf_band=hf_band)

    return freqency_domain_features


def _get_features_from_psd(freq, psd, 
                           vlf_band = VlfBand(0.003, 0.04),
                           lf_band = LfBand(0.04, 0.15),
                           hf_band = HfBand(0.15, 0.40)):
    """
    Computes frequency domain features from the power spectral decomposition.
    Parameters
    ---------
    freq : array
        Array of sample frequencies.
    psd : list
        Power spectral density or power spectrum.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    lf_band : tuple
        Low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.
    Returns
    ---------
    freqency_domain_features : dict
        Dictionary containing frequency domain features for HRV analyses. There are details
        about each features given below.
    """

    # Calcul of indices between desired frequency bands
    vlf_indexes = np.logical_and(freq >= vlf_band[0], freq < vlf_band[1])
    lf_indexes = np.logical_and(freq >= lf_band[0], freq < lf_band[1])
    hf_indexes = np.logical_and(freq >= hf_band[0], freq < hf_band[1])

    # Integrate using the composite trapezoidal rule
    lf = np.trapz(y=psd[lf_indexes], x=freq[lf_indexes])
    hf = np.trapz(y=psd[hf_indexes], x=freq[hf_indexes])

    # total power & vlf : Feature often used for  "long term recordings" analysis
    vlf = np.trapz(y=psd[vlf_indexes], x=freq[vlf_indexes])
    total_power = vlf + lf + hf

    lf_hf_ratio = lf / hf
    lfnu = (lf / (lf + hf)) * 100
    hfnu = (hf / (lf + hf)) * 100

    freqency_domain_features = {
        'lf': lf,
        'hf': hf,
        'lf_hf_ratio': lf_hf_ratio,
        'lfnu': lfnu,
        'hfnu': hfnu,
        'total_power': total_power,
        'vlf': vlf
    }

    return freqency_domain_features

def plot_HRV(df_window):
    
    """
    Returns Poincarrè plot and spectral analysis plot
    """
    # DRAW PLOTS
    fig, (ax1,ax2) = plt.subplots(ncols=2, nrows=1,figsize=(10,7))

    # POINCARè PLOT
    nn_intervals = list(1000*df_window['ibi_s'].dropna().values)
    val1 = nn_intervals[:-1]
    val2 = nn_intervals[1:]

    # compute features for ellipse's height, width and center
    dict_sd1_sd2 = get_poincare_plot_features(nn_intervals)
    sd1 = dict_sd1_sd2["sd1"]
    sd2 = dict_sd1_sd2["sd2"]
    mean_nni = np.mean(nn_intervals)

    ax1.scatter(val1, val2, c='b', s=20, alpha=0.4)

    ells = Ellipse(xy=(mean_nni, mean_nni), width=2 * sd2 + 1,
            height=2 * sd1 + 1, angle=45, linewidth=3, edgecolor = 'k',
            fill=False)
    ax1.add_patch(ells)

    ells = Ellipse(xy=(mean_nni, mean_nni), width=2 * sd2,
                    height=2 * sd1, angle=45)
    ells.set_alpha(0.05)
    ells.set_facecolor("blue")
    ax1.add_patch(ells)

    # Arrow plot settings
    sd1_arrow = ax1.arrow(mean_nni, mean_nni, -sd1 * np.sqrt(2) / 2, sd1 * np.sqrt(2) / 2,
                          linewidth=4, ec='r', fc="r", label="SD1")
    sd2_arrow = ax1.arrow(mean_nni, mean_nni, sd2 * np.sqrt(2) / 2, sd2 * np.sqrt(2) / 2,
                          linewidth=4, ec='orange', fc="orange", label="SD2")

    ax1.plot([np.min(nn_intervals)-20,np.max(nn_intervals)+20],[np.min(nn_intervals)-20,np.max(nn_intervals)+20],'k--')
    ax1.legend(handles=[sd1_arrow, sd2_arrow], fontsize=12, loc="best")

    ax1.set_title("Poincaré Plot", fontsize=20)
    ax1.set_xlabel('NN_n (ms)', fontsize=15)
    ax1.set_ylabel('NN_n+1 (ms)', fontsize=15)


    # PSD PLOT
    VlfBand = namedtuple("Vlf_band", ["low", "high"])
    LfBand = namedtuple("Lf_band", ["low", "high"])
    HfBand = namedtuple("Hf_band", ["low", "high"])

    vlf_band = VlfBand(0.003, 0.04)
    lf_band = LfBand(0.04, 0.15)
    hf_band = HfBand(0.15, 0.40)

    freq, psd = _get_freq_psd_from_nn_intervals(nn_intervals=nn_intervals, method = "welch", sampling_frequency = 7,
                                                    interpolation_method = "linear")

    # Calcul of indices between desired frequency bands
    vlf_indexes = np.logical_and(freq >= vlf_band[0], freq < vlf_band[1])
    lf_indexes = np.logical_and(freq >= lf_band[0], freq < lf_band[1])
    hf_indexes = np.logical_and(freq >= hf_band[0], freq < hf_band[1])

    frequency_band_index = [vlf_indexes, lf_indexes, hf_indexes]
    label_list = ["VLF component", "LF component", "HF component"]

    # Plot parameters
    # sns.set_style("ticks")
    ax2.set_xlabel("Frequency (Hz)", fontsize=15)
    ax2.set_ylabel("PSD (s2/ Hz)", fontsize=15)

    ax2.set_title("FFT Spectrum", fontsize=20) # Welch's periodogram
    for band_index, label in zip(frequency_band_index, label_list):
        ax2.fill_between(freq[band_index], 0, psd[band_index] / (1000 * len(psd[band_index])), label=label)
    ax2.legend(prop={"size": 12}, loc="best")
    ax2.set_xlim(0, hf_band[1])

    fig.tight_layout()
    plt.show()