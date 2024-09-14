import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt


# Cover function for leastsq more suited for physics labs.
# notice the main feature is an adjust list.
def curfit(x, y, sigy, a0, function_name=None, adj=None):
    """
      chisq,a,sigmaa,nresids = curfit(x,y,sigy,a0,function_name=None,adj=None):

    input:
      x,y are the data.
      sigy are the errors in y and must be given.
      a0 is the initial guess, which must be given.
      function_name is the fitting function of the form y=f(x,p) where
          p are the parameters (like a0). By default it is named funct(x,p).
      adj is an array of indices into a0 that will be varied or fit to. The
          rest are held constant. Default is all parameters.
    results:
      chisq is the reduced chi-square of the fit.
      a are  the fitted result.
      sigmaa are the errors in the parameters with same meaning as in sigy.
          That means if sigy is 2 sigma so is sigmaa.
      nresids are the normalized residuals (y-f(x,a))/sigy.
    notes:
      the covariance matrix (cv) is available in the code and could be
          returned if needed by a simple change in code.
    """
    a = a0.copy()
    if adj is None:
        adj = np.arange(len(a), dtype="int32")
    if function_name is None:
        function_name = funct
    afit, cv, idt, m, ie = leastsq(
        _residuals,
        a[adj],
        args=(x, y, sigy, a, adj, function_name),
        full_output=True,
        ftol=0.001,
        epsfcn=0.00025**2,
    )
    # ftol=.001,epsfcn=.000025**2)
    # coarser numerical derivatives usually work better.
    # print("ie=%d"%ie,m)
    a[adj] = afit
    realcv = np.identity(a0.size)
    nresids = idt["fvec"]
    # yfit=y-nresids*sigy
    chisq = np.sum(nresids**2) / (len(y) - len(adj))
    if cv is None:
        raise RuntimeError("Problem with fit: is number of parameters correct?")
    realcv[np.ix_(adj, adj)] = cv
    sigmaa = np.zeros(len(a0))
    sigmaa[adj] = np.sqrt(np.diag(cv))
    return (chisq, a, sigmaa, nresids)


# hidden residuals function for leastsq
def _residuals(p, x, y, sigy, pall, adj, fun):
    pall[adj] = p
    return (fun(x, pall) - y) / sigy


def fitpr(chisq, a, sigmaa, title=None, lbl=None):
    """str = fitpr(chisq,a,sigmaa,title=None,lbl=None)
    Generate string to nicely print out results of a fit.
    """
    # get fitted results.
    if lbl is None:
        lbl = []
        for i in xrange(a.size):
            lbl.append("A%(#)02d" % {"#": i})
    # print resuls of a fit.
    str = ""
    if title is not None:
        str = title + " "
    str += "chisq=%(c).4f" % {"c": chisq}
    for i in range(a.size):
        str += "\n     %(lbl)8s =%(m)10.4f +/- %(s).4f" % {
            "lbl": lbl[i],
            "m": a[i],
            "s": sigmaa[i],
        }
    return str


# easy plot for fit
def fitplot(x, y, sigy, nresids, pl=plt):
    """Plot the data and the function."""
    yfit = y - sigy * nresids
    plt.plot(x, yfit)
    plt.errorbar(x, y, fmt="o", yerr=sigy)


# define a default function, a straight line.
def funct(x, p):
    return p[0] * x + p[1]
