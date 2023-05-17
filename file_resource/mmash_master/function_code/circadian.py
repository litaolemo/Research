def fit_sin(tt, yy,plot=False):
  
    '''
    Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period"
    Parameters
    ---------
    tt : list
        List of timestamp expressed in seconds
    yy: list
        List of values to fit circadian rhythm (e.g., Heart Rate, intra-beats intervals).
    Returns
    ---------
    res : dictionaire
        dictionaires cotaining fitting parmeters.
    '''
    
    import scipy
    import numpy
    from sklearn.metrics import r2_score
    import pandas as pd
    import matplotlib.pyplot as plt

    tt = numpy.array(tt)
    yy = numpy.array(yy)
    freq = 24*3600 # dayly frequency (24h in seconds)
    guess_amp = numpy.std(yy) * 2.**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.*numpy.pi*freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return  A * numpy.sin((t/(freq))*2*numpy.pi + p) + c 

    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*numpy.pi)

    fitfunc = lambda t: A * numpy.sin((t/freq)*2*numpy.pi + p) + c

    coefficient_of_dermination = r2_score(yy, fitfunc(tt))
    acrophase_time = pd.DataFrame([[x,y] for x,y in zip(tt,fitfunc(tt))],columns=['t','val']).set_index('t')
    acrophase_timestamp = acrophase_time[acrophase_time.val == acrophase_time.max().values[0]].index[0]

    res = {"amp": abs(A), "phase": p, "APhase": acrophase_timestamp, "offset": c, "r2":coefficient_of_dermination,"tt":tt,'ff':fitfunc(tt),"yy":yy}
    
    if plot==True:
        fig,ax = plt.subplots(figsize=(6,5))
        ax.plot(tt, yy, "-k", linewidth=1, alpha=0.3)
        ax.plot(tt, fitfunc(tt), "r-", label="circadian rhythm", linewidth=1)
        plt.xticks(numpy.arange(130000)[::20000][2:],['11AM day 1','5:30PM day 1','10PM day 1', '3:30AM day 2', '9AM day 2'],rotation=90)
        plt.legend()
        fig.tight_layout()
        plt.show()

    return res