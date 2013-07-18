"""

   Common statistics library

"""


# python libs
from math import *
from itertools import izip
import cmath
import random
import os

# rasmus libs
from rasmus import util
from rasmus import tablelib




def logprod(lst):
    """Computes the product of a list of numbers"""
    return sum(log(i) for i in lst)


def prod(lst):
    """Computes the product of a list of numbers"""
    p = 1.0
    for i in lst:
        p *= i
    return p


def mean(vals):
    """Computes the mean of a list of numbers"""
    n = 0
    s = 0.0
    for i in vals:
        s += i
        n += 1
    return s / float(n)

def median(vals):
    """Computes the median of a list of numbers"""
    lenvals = len(vals)
    sortvals = sorted(vals)
    
    if lenvals % 2 == 0:
        return (sortvals[lenvals / 2] + sortvals[lenvals / 2 - 1]) / 2.0
    else:
        return sortvals[lenvals / 2]

def mode(vals):
    """Computes the mode of a list of numbers"""
    top = 0
    topkey = None
    for key, val in util.hist_dict(vals).iteritems():
        if val > top:
            top = val
            topkey = key
    return topkey


def msqerr(vals1, vals2):
    """Mean squared error"""
    
    assert len(vals1) == len(vals2), "lists are not the same length"
    
    
    return mean([(vals1[i] - vals2[i]) ** 2 
                 for i in xrange(len(vals1))])
    
    

def variance(vals):
    """Variance"""
    u = mean(vals)
    return sum((x - u)**2 for x in vals) / float(len(vals)-1)

def sdev(vals):
    """Standard deviation"""
    return sqrt(variance(vals))

def serror(vals):
    """Stanadrd error"""
    return sdev(vals) / sqrt(len(vals))

def covariance(lst1, lst2):
    """Covariance"""
    m1 = mean(lst1)
    m2 = mean(lst2)
    tot = 0.0
    for i in xrange(len(lst1)):
        tot += (lst1[i] - m1) * (lst2[i] - m2)    
    return tot / (len(lst1)-1)


def covmatrix(mat):
    """Covariance Matrix"""
    size = len(mat)
    return [[cov(mat[i], mat[j]) for j in range(size)]
            for i in range(size)]

def corrmatrix(mat):
    """Correlation Matrix"""
    size = len(mat)
    return [[corr(mat[i], mat[j]) for j in range(size)]
            for i in range(size)]


def corr(lst1, lst2):
    """Pearson's Correlation"""
    num = covariance(lst1, lst2)
    denom = float(sdev(lst1) * sdev(lst2))
    if denom != 0:
        return num / denom
    else:
        return 1e1000


def corr_pvalue(r, n):
    """Returns the signficance of correlation > r with n samples"""

    t = r / sqrt((1 - r*r) / float(n - 2))
    return rpy.r.pt(-t, n-2)
    

def qqnorm(data, plot=None):
    """Quantile-quantile plot"""
    
    from rasmus import gnuplot

    data2 = sorted(data)
    norm = [random.normalvariate(0, 1) for x in range(len(data2))]
    norm.sort()
    
    if plot == None:
        return gnuplot.plot(data2, norm)
    else:
        plot.plot(data2, norm)
        return plot


def entropy(probs, base=2):
    """Shannon's entropy"""

    return - sum(p * log(p, base) for p in probs if p > 0.0)

def cross_entropy(p, q, base=2):
    try:
        return - sum(i * log(j, base) for i,j in izip(p, q) if i > 0.0)
    except OverflowError:
        return util.INF

def kl_div(p, q):
    """Compute the KL divergence for two discrete distributions"""
    return cross_entropy(p, q) - entropy(p)

def akaike_ic(lnl, k):
    """Akaike information criterion"""
    return 2 * k - 2 * lnl

def akaike_icc(lnl, n, k):
    """Akaike information criterion with second order correction
       Good for small sample sizes
    """
    return akaike_ic(lnl, k) + 2*k*(k+1) / (n - k - 1)


def bayesian_ic(lnl, n, k):
    """Bayesian information criterion

       lnl -- ln(L)
       n   -- number of data points
       k   -- number of parameters
    """
    return -2 * lnl + k * log(n)



def fitLine(xlist, ylist):
    """2D regression"""
    
    xysum = 0
    xxsum = 0
    n = len(xlist)        
    for i in range(n):
        xysum += xlist[i] * ylist[i]
        xxsum += xlist[i] * xlist[i]
    avgx = mean(xlist)
    avgy = mean(ylist)

    if (xxsum - n*avgx*avgx) == 0:
        slope = 1e10
    else:
        slope = (xysum - n*avgx*avgy) / float(xxsum - n*avgx*avgx)

    inter = (avgy*xxsum - avgx*xysum) / float(xxsum - n*avgx*avgx)

    return (slope, inter)


def fitLineError(xlist, ylist, slope, inter):
    """Returns the Mean Square Error of the data fit"""
    error = 0
    n = len(xlist)
    
    for i in range(n):
        error += ((xlist[i]*slope + inter) - ylist[i]) ** 2
    return error / n


def pearsonsRegression(observed, expected):
    """Pearson's coefficient of regression"""
    
    # error sum of squares
    ess = sum((a - b)**2 for a, b in izip(observed, expected))
    
    # total sum of squares
    u = mean(observed)
    tss = sum((a - u)**2 for a in observed)
    
    r2 = 1 - ess / tss
    return r2

    
def pearsonsRegressionLine(x, y, m, b):
    observed = y
    expected = [m*i + b for i in x]
    return pearsonsRegression(observed, expected)


def rank(vals, x, norm=False, sort=True):
    """
    Returns the rank of x in list vals
    rank(x) = i if vals[i-1] <= x < vals[i]

    x    -- value to rank within values
    vals -- list of values to compute the rank of
    sort -- if True, vals will be sorted first
    norm -- if True, return normalized ranks (i.e. percentiles)
    """

    if sort:
        vals = sorted(vals)
    n = len(vals)

    for r, v in enumerate(vals):
        if v > x:
            break
    else:
        r = n

    if norm:
        r /= float(n + 1)

    return r


def percentile(vals, perc, rounding=-1, sort=True):
    """Give the value at a percentile 'perc'

       vals     -- list of values
       perc     -- perctile
       rounding -- round down if -1 or round up for 1
       sort     -- if True, sort vals first
    """
    
    if sort:
        vals2 = sorted(vals)
    else:
        vals2 = vals
    n = len(vals2)
    if rounding == -1:
        return vals2[util.clamp(int(perc * n), 0, n-1)]
    elif rounding == 1:
        return vals2[util.clamp(int(ceil(perc * n)), 0, n-1)]
    else:
        raise Exception("rounding must be 1 or -1")


def dither(vals, radius):
    return [x + random.uniform(-radius, radius) for x in vals]


def logadd(lna, lnb):
    """Adding numbers in log-space"""
    
    diff = lna - lnb
    if diff < 500:
        return log(exp(diff) + 1.0) + lnb
    else:
        return lna


def logsum(vals):
    SUM_LOG_THRESHOLD = -15
    maxval = vals[0]
    maxi = 0

    # find maxval
    for i in range(1, len(vals)):
        if vals[i] > maxval:
            maxval = vals[i]
            maxi = i

    expsum = 1.0
    for i in xrange(len(vals)):
        if i != maxi and vals[i] - maxval > SUM_LOG_THRESHOLD:
            expsum += exp(vals[i] - maxval)
  
    return maxval + log(expsum)


def logsub(lna, lnb):
    """
    subtracting numbers in log-space

    must have lna > lnb
     """

    diff = lna - lnb
    if diff < 500:
        diff2 = exp(diff) - 1.0
        if diff2 == 0.0:
            return -util.INF
        else:
            return log(diff2) + lnb
    else:
        return lna
    

def logadd_sign(sa, lna, sb, lnb):
    """Adding numbers in log-space"""

    if sa > 0 and sb > 0:
        return 1, logadd(lna, lnb)

    elif sa == 0:
        return sb, lnb

    elif sb == 0:
        return sa, lna

    elif sa < 0 and sb < 0:
        return -1, logadd(lna, lnb)

    elif sa > 0 and sb < 0:
        if lna > lnb:
            return 1, logsub(lna, lnb)
        elif lna == lnb:
            return 0, -util.INF
        else:
            return -1, logsub(lnb, lna)

    elif sa < 0 and sb > 0:
        if lna > lnb:
            return -1, logsub(lna, lnb)
        elif lna == lnb:
            return 0, -util.INF
        else:
            return 1, logsub(lnb, lna)

    else:
        raise Exception("unhandled case")


def smooth(vals, radius):
    """
    return an averaging of vals using a radius
    
    Note: not implemented as fast as possible
    runtime: O(len(vals) * radius)
    """
    
    vals2 = []
    vlen = len(vals)
    
    for i in xrange(vlen):
        radius2 = min(i, vlen - i - 1, radius)
        vals2.append(mean(vals[i-radius2:i+radius2+1]))
    
    return vals2




def iter_window_index(x, xdist, esp=None):
    """
    iterates a sliding window over x with width 'xdist'

    returns an iterator over list of indices in x that represent windows
    
    x must be sorted least to greatest
    """

    vlen = len(x)
    #if esp is None:
    #    esp = min(x[i+1] - x[i] for i in range(vlen-1)
    #              if x[i+1] - x[i] > 0) / 2.0
    
    # simple case
    if vlen == 0:
        return
    
    start = x[0]
    end = x[-1]
    window = [0]
    
    low = start
    high = start + xdist
    lowi = 0 # inclusive
    highi = 0 # inclusive

    # move up high boundary
    while highi+1 < vlen and x[highi+1] < high:
        highi += 1

    yield (lowi, highi, low, high)
    
    while highi+1 < vlen:
        low_step = x[lowi] - low    # dist until expell
        high_step = x[highi+1] - high # dist until include

        # advance though duplicates
        if low_step == 0:
            lowi += 1
            continue
        
        if high_step == 0:
            highi += 1
            continue

        # determine new low high boundary
        if low_step <= high_step:
            low = x[lowi] #+ min(esp, (high_step - low_step) / 2.0)
            high = low + xdist            
            lowi += 1
            
        if high_step <= low_step:
            highi += 1
            if highi >= vlen: break
            high = x[highi] #+ min(esp, (low_step - high_step) / 2.0)
            low = high - xdist

        assert abs((high - low) - xdist) < .001, (low, high)
        
        yield (lowi, highi, low, high)


def iter_window_index_step(x, size, step, minsize=0):

    vlen = len(x)
    start = x[0]
    end = x[-1]

    low = start
    high = start + size
    i = 1

    lowi = 0
    highi = 0
    
    # move up high boundary
    while highi+1 < vlen and x[highi+1] < high:
        highi += 1

    while highi < vlen and high < end:
        if highi - lowi >= minsize:
            yield lowi, highi, low, high
        low = start + i * step
        high = low + size
        i += 1

        # move up low boundary
        while lowi < vlen and x[lowi] < low:
            lowi += 1

        # move up high boundary
        while highi+1 < vlen and x[highi+1] < high:
            highi += 1
        
    

def iter_window(x, xdist, func=lambda win: win, minsize=0, key=lambda x: x):
    """
    iterates a sliding window over x with radius xradius
    
    x must be sorted least to greatest
    """

    for lowi, highi, low, high in iter_window_index(map(key, x), xdist):
        if highi - lowi >= minsize:
            yield (high + low)/2.0, func(x[lowi:highi])


def iter_window_step(x, width, step, func=lambda win: win, minsize=0):
    """
    iterates a sliding window over x with width 'width'
    
    x must be sorted least to greatest

    return an iterator with (midx, func(x[lowi:highi]))
    """
    
    for lowi, highi, low, high in iter_window_index_step(x, width, step, minsize):
        yield (high + low) / 2.0, func(x[lowi:highi])






def smooth2(x, y, xradius, minsize=0, sort=False):
    """
    return an averaging of x and y using xradius
    
    x must be sorted least to greatest
    """

    vlen = len(x)
    assert vlen == len(y)
    
    # simple case
    if vlen == 0:
        return [], []
    
    if sort:
        x, y = util.sort_many(x, y)
    
    x2 = []
    y2 = []
    
    start = min(x)
    end = max(x)
    xtot = x[0]
    ytot = y[0]
    
    low = 0
    high = 0
    
    for i in xrange(vlen):
        xi = x[i]
    
        xradius2 = min(xi - start, end - xi, xradius)
    
        # move window
        while x[low] < xi - xradius2:
            xtot -= x[low]
            ytot -= y[low]
            low += 1
        while x[high] < xi + xradius2:
            high += 1
            xtot += x[high]
            ytot += y[high]
        
        denom = float(high - low + 1)
        if denom >= minsize:
            x2.append(xtot / denom)
            y2.append(ytot / denom)
    
    return x2, y2


def factorial(x, k=1):
    """Simple implementation of factorial"""
    
    n = 1
    for i in xrange(int(k)+1, int(x)+1):
        n *= i
    return n


def logfactorial(x, k=1):
    """returns the log(factorial(x) / factorial(k)"""
    
    n = 0
    for i in xrange(int(k)+1, int(x)+1):
        n += log(i)
    return n


def choose(n, k):
    if n == 0 and k == 0:
        return 1
        
    if n < 0 or k < 0 or k > n:
        return 0
    
    # optimization for speed
    if k > n/2:
        k = n - k
    
    t = 1.0
    n2 = n + 1.0
    for i in xrange(1, k+1):
        t *= (n2 - i) / i
    return int(t + 0.5)
    #return factorial(n, n - k) / factorial(k)

def fchoose(n, k):
    if n == 0 and k == 0:
        return 1
        
    if n < 0 or k < 0 or k > n:
        return 0
    
    # optimization for speed
    if k > n/2:
        k = n - k
    
    t = 1.0
    n2 = n + 1.0
    for i in xrange(1, k+1):
        t *= (n2 - i) / i
    return t

def logchoose(n, k):
    if n == 0 and k == 0:
        return 0.0
        
    if n < 0 or k < 0 or k > n:
        return -util.INF
    
    # optimization for speed
    if k > n/2:
        k = n - k
    
    t = 0.0
    n2 = n + 1.0
    for i in xrange(1, k+1):
        t += log((n2 - i) / i)
    return t


def multinomial(vals):
    n = sum(vals)

    res = logfactorial(n)
    for v in vals:
        res -= logfactorial(v)
    return int(exp(res) + .05)


def logmultinomial(vals):
    n = sum(vals)

    res = logfactorial(n)
    for v in vals:
        res -= logfactorial(v)
    return res


def sample(weights):
    """
    Randomly choose an int between 0 and len(probs)-1 using
    the weights stored in list probs.
    
    item i will be chosen with probability weights[i]/sum(weights)
    """

    total = sum(weights)
    pick = random.random() * total
    x = 0
    for i in xrange(len(weights)):
        x += weights[i]
        if x >= pick:
            return i
    return len(weights) - 1
    

def rhyper(m, n, M, N, report=0):
    '''
    calculates cumulative probability based on
    hypergeometric distribution
    over/under/both (report = 0/1/2)
    (uses R through RPy)
    
    N = total balls in urn
    M = total white balls in urn
    n = drawn balls from urn
    m = drawn white balls from urn
    
    '''

    from rpy import r

    
    assert( (type(m) == type(n) == type(M) == type(N) == int)
            and m <= n and m <= M and n <= N)
    
    
    
    if report == 0:
        #p-val for over-repr.
        return r.phyper(m-1, M, N-M, n, lower_tail=False)
    elif report == 1:
        #p-val for under-repr.
        return r.phyper(m, M, N-M, n)
    elif report == 2:
        #tuple (over, under)
        return r.phyper(m-1, M, N-M, n, lower_tail=False), r.phyper(m, M, N-M, n)
    else:
        raise "unknown option"


def cdf(vals):
    """Computes the CDF of a list of values"""
    
    vals = sorted(vals)
    tot = float(len(vals))
    x = []
    y = []
    
    for i, x2 in enumerate(vals):
        x.append(x2)
        y.append(i / tot)
        
    return x, y
    
    
def enrichItems(in_items, out_items, M=None, N=None, useq=True, extra=False):
    """Calculates enrichment for items within an in-set vs and out-set.
       Returns a sorted table.
    """
    # DEPRECATED
    # TODO: remove this function

    
    # count items
    counts = util.Dict(default=[0, 0])
    for item in in_items:
        counts[item][0] += 1
    for item in out_items:
        counts[item][1] += 1
    
    if N is None:
        N = len(in_items) + len(out_items)
    if M is None:
        M = len(in_items)
    
    tab = tablelib.Table(headers=["item", "in_count", "out_count", 
                                  "pval", "pval_under"])
    
    # do hypergeometric
    for item, (a, b) in counts.iteritems():
        tab.add(item=item,
                in_count=a,
                out_count=b,
                pval=rhyper(a, a+b, M, N),
                pval_under=rhyper(a, a+b, M, N, 1))
    
    # add qvalues
    if useq:
        qval = qvalues(tab.cget("pval"))
        qval_under = qvalues(tab.cget("pval_under"))
        
        tab.add_col("qval", data=qval)
        tab.add_col("qval_under", data=qval_under)
    
    if extra:
        tab.add_col("in_size", data=[M]*len(tab))
        tab.add_col("out_size", data=[N-M]*len(tab))
        tab.add_col("item_ratio", data=[
            row["in_count"] / float(row["in_count"] + row["out_count"])
            for row in tab])
        tab.add_col("size_ratio", data=[
            M / float(N) for row in tab])
        tab.add_col("fold", data=[row["item_ratio"] / row["size_ratio"]
                                 for row in tab])
    
    tab.sort(col='pval')
    return tab


def qvalues(pvals):
    import rpy
    ret = rpy.r.p_adjust(pvals, "fdr")
    return ret



#=============================================================================
# Distributions
#

def uniformPdf(x, params):
    a, b = params
    if x < a or x > b:
        return 0.0
    else:
        return 1.0 / (b - a)


def binomialPdf(k, params):
    p, n = params
    return choose(n, k) * (p ** k) * ((1.0-p) ** (n - k))

def gaussianPdf(x, params):
    return 1/sqrt(2*pi) * exp(- x**2 / 2.0)

def normalPdf(x, params):
    mu, sigma = params
    # sqrt(2*pi) = 2.5066282746310002
    return exp(- (x - mu)**2 / (2.0 * sigma**2)) / (sigma * 2.5066282746310002)

def normalCdf(x, params):
    mu, sigma = params
    return (1 + erf((x - mu)/(sigma * sqrt(2)))) / 2.0

def logNormalPdf(x, params):
    """mu and sigma are the mean and standard deviation of the 
       variable's logarithm"""
    
    mu, sigma = params
    return 1/(x * sigma * sqrt(2*pi)) * \
           exp(- (log(x) - mu)**2 / (2.0 * sigma**2))

def logNormalCdf(x, params):
    """mu and sigma are the mean and standard deviation of the 
       variable's logarithm"""
    
    mu, sigma = params
    return (1 + erf((log(x) - mu)/(sigma * sqrt(2)))) / 2.0


def poissonPdf(x, params):
    lambd = params[0]
    
    if x < 0 or lambd <= 0:
        return 0.0
    
    a = 0
    for i in xrange(1, int(x)+1):
        a += log(lambd / float(i))
    return exp(-lambd + a)


def poissonCdf(x, params):
    """Cumulative distribution function of the Poisson distribution"""
    # NOTE: not implemented accurately for large x or lambd
    lambd = params[0]
    
    if x < 0:
        return 0
    else:
        return (gamma(floor(x+1)) - gammainc(floor(x + 1), lambd)) / \
               factorial(floor(x))


def poissonvariate(lambd):
    """Sample from a Poisson distribution"""
    l = -lambd
    k = 0
    p = 0.0

    while 1:
        k += 1
        p += log(random.random())
        if p < l:
            return k - 1


def exponentialPdf(x, params):
    lambd = params[0]
    
    if x < 0 or lambd < 0:
        return 0.0
    else:
        return lambd * exp(-lambd * x)


def exponentialCdf(x, params):
    lambd = params[0]
    
    if x < 0 or lambd < 0:
        return 0.0
    else:
        return 1.0 - exp(-lambd * x)


def exponentialvariate(lambd):
    return -log(random.random()) / lambd

def gammaPdf(x, params):
    alpha, beta = params
    if x <= 0 or alpha <= 0 or beta <= 0:
        return 0.0
    else:
        return (exp(-x * beta) * (x ** (alpha - 1)) * (beta ** alpha)) / \
           gamma(alpha)


def loggammaPdf(x, params):
    alpha, beta = params
    if x <= 0.0 or alpha <= 0.0 or beta <= 0.0:
        return -util.INF
    else:
        return -x*beta + (alpha - 1)*log(x) + alpha*log(beta) - gammaln(alpha)

def gammaPdf2(x, params):
    alpha, beta = params
    if x <= 0 or alpha <= 0 or beta <= 0:
        return 0.0
    else:
        return exp(loggammaPdf(x, params))


def gammaCdf(x, params):
    alpha, beta = params
    if x <= 0:
        return 0
    else:
        return gammainc(alpha, x * beta) / gamma(alpha)

def invgammaPdf(x, params):
    a, b = params

    if x <=0 or a <= 0 or b <= 0:
        return 0.0
    else:
        return (b**a) / gamma(a) * (1.0/x)**(a + 1) * exp(-b/x)

def loginvgammaPdf(x, params):
    a, b = params
    if x < 0 or a < 0 or b < 0:
        return -util.INF
    else:
        return a*log(b) - gammaln(a) + (a+1)*log(1.0/x) -b/x




def betaPdf2(x, params):
    """A simpler implementation of beta distribution but will overflow
       for values of alpha and beta near 100
    """
    
    alpha, beta = params
    if 0 < x < 1 and alpha > 0 and beta > 0:
        return gamma(alpha + beta) / (gamma(alpha)*gamma(beta)) * \
               x ** (alpha-1) * (1-x)**(beta-1)
    else:
        return 0.0

def betaPdf(x, params):
    alpha, beta = params
    
    if 0 < x < 1 and alpha > 0 and beta > 0:
        return exp(gammaln(alpha + beta) - (gammaln(alpha) + gammaln(beta)) + \
                   (alpha-1) * log(x) +  (beta-1) * log(1-x))
    else:
        return 0.0
    


def betaPdf3(x, params):
    alpha, beta = map(int, params)
    if 0 < x < 1 and alpha > 0 and beta > 0:
        n = min(alpha-1, beta-1)
        m = max(alpha-1, beta-1)
        
        prod1 = 1
        for i in range(1,n+1):
            prod1 *= ((n+i)*x*(1-x))/i
        
        prod2 = 1
        if alpha > beta:
            for i in range(n+1, m+1):
                prod2 *= ((n+i)*x)/i
        else:
            for i in range(n+1, m+1):
                prod2 *= ((n+i)*(1-x))/i
        
        return prod1 * prod2 * (alpha + beta - 1)
    else:
        return 0.0


def negbinomPdf(k, r, p):
    return exp(gammaln(r+k) - gammaln(k+1) - gammaln(r) +
               r*log(p) + k * log(1-p))



def gamma(x):
    """
    Lanczos approximation to the gamma function. 
    
    found on http://www.rskey.org/gamma.htm   
    """
    
    ret = 1.000000000190015 + \
          76.18009172947146 / (x + 1) + \
          -86.50532032941677 / (x + 2) + \
          24.01409824083091 / (x + 3) + \
          -1.231739572450155 / (x + 4) + \
          1.208650973866179e-3 / (x + 5) + \
          -5.395239384953e-6 / (x + 6)
    
    return ret * sqrt(2*pi)/x * (x + 5.5)**(x+.5) * exp(-x-5.5)



def gammaln(xx):
    """
    From numerical alogrithms in C

    float gammln(float xx)
    Returns the value ln[(xx)] for xx > 0.
    {
        Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
        accuracy is good enough.
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
             24.01409824083091,-1.231739572450155,
             0.1208650973866179e-2,-0.5395239384953e-5};
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
    }
    """

    cof = [76.18009172947146,-86.50532032941677,
         24.01409824083091,-1.231739572450155,
         0.1208650973866179e-2,-0.5395239384953e-5]
    
    y = x = xx
    tmp = x + 5.5
    tmp -= (x + 0.5) * log(tmp)
    ser = 1.000000000190015
    
    for j in range(6):
        y += 1
        ser += cof[j] / y
    
    return - tmp + log(2.5066282746310005 * ser / x)
    



GAMMA_INCOMP_ACCURACY = 1000
def gammainc(a, x):
    """Lower incomplete gamma function"""
    # found on http://www.rskey.org/gamma.htm
    
    ret = 0
    term = 1.0/x
    for n in xrange(GAMMA_INCOMP_ACCURACY):
        term *= x/(a+n)
        ret += term
        if term < .0001:
            break
    return x**a * exp(-x) * ret


def erf(x):
    # http://www.theorie.physik.uni-muenchen.de/~serge/erf-approx.pdf
    
    a = 8/(3*pi) * (pi - 3)/(4 - pi)
    axx = a * x * x
    
    if x >= 0:
        return sqrt(1 - exp(-x*x * (4.0/pi + axx)/(1 + axx)))
    else:
        return - sqrt(1 - exp(-x*x * (4.0/pi + axx)/(1 + axx)))
    


def chiSquare(rows, expected=None, nparams=0):
    # ex: rows = [[1,2,3],[1,4,5]]
    assert util.equal(map(len, rows))

    if 0 in map(sum,rows): return 0,1.0
    cols = zip(* rows)
    if 0 in map(sum,cols): return 0,1.0

    if not expected:
        expected = make_expected(rows)

    chisq = 0
    for obss,exps in zip(rows,expected):
        for obs, exp in zip(obss, exps):
            chisq += ((obs-exp)**2)/exp

    df = max(len(rows)-1, 1)*max(len(rows[0])-1, 1) - nparams

    p = chi_square_lookup(chisq,df)

    return chisq,p


def make_expected(rows):
    rowtotals = map(sum, rows)
    coltotals = map(sum, zip(* rows))
    grandtotal = float(sum(rowtotals))

    expected = []
    for row,rowtotal in zip(rows,rowtotals):
        expected_row = []
        for obs, coltotal in zip(row, coltotals):
            exp = rowtotal * coltotal / grandtotal
            expected_row.append(exp)
        expected.append(expected_row)
    return expected


def chiSquareFit(xbins, ybins, func, nsamples, nparams, minsamples=5):
    sizes = [xbins[i+1] - xbins[i] for i in xrange(len(xbins)-1)]
    sizes.append(sizes[-1]) # NOTE: assumes bins are of equal size
    
    # only focus on bins that are large enough
    counts = [ybins[i] * sizes[i] * nsamples for i in xrange(len(xbins)-1)]
    
    expected = []
    for i in xrange(len(xbins)-1):
        expected.append((func(xbins[i]) + func(xbins[i+1]))/2.0 * 
                         sizes[i] * nsamples)
        
    # ensure we have enough expected samples in each bin
    ind = util.find(util.gefunc(minsamples), expected)
    counts = util.mget(counts, ind)
    expected = util.mget(expected, ind)
    
    if len(counts) == 0:
        return [0, 1], counts, expected
    else:
        return chiSquare([counts], [expected], nparams), counts, expected


chi_square_table = {
    1: [1.64, 2.71, 3.84, 5.02, 6.64, 10.83],
    2: [3.22, 4.61, 5.99, 7.38, 9.21, 13.82],
    3: [4.64, 6.25, 7.82, 9.35, 11.34, 16.27],
    4: [5.99, 7.78, 9.49, 11.14, 13.28, 18.47],
    5: [7.29, 9.24, 11.07, 12.83, 15.09, 20.52],
    6: [8.56, 10.64, 12.59, 14.45, 16.81, 22.46],
    7: [9.80, 12.02, 14.07, 16.01, 18.48, 24.32],
    8: [11.03, 13.36, 15.51, 17.53, 20.09, 26.12],
    9: [12.24, 14.68, 16.92, 19.02, 21.67, 27.88],
    10: [13.44, 15.99, 18.31, 20.48, 23.21, 29.59],
    11: [14.63, 17.28, 19.68, 21.92, 24.72, 31.26],
    12: [15.81, 18.55, 21.03, 23.34, 26.22, 32.91],
    13: [16.98, 19.81, 22.36, 24.74, 27.69, 34.53],
    14: [18.15, 21.06, 23.68, 26.12, 29.14, 36.12],
    15: [19.31, 22.31, 25.00, 27.49, 30.58, 37.70],
    16: [20.47, 23.54, 26.30, 28.85, 32.00, 39.25],
    17: [21.61, 24.77, 27.59, 30.19, 33.41, 40.79],
    18: [22.76, 25.99, 28.87, 31.53, 34.81, 42.31],
    19: [23.90, 27.20, 30.14, 32.85, 36.19, 43.82],
    20: [25.04, 28.41, 31.41, 34.17, 37.57, 45.31],
    21: [26.17, 29.62, 32.67, 35.48, 38.93, 46.80],
    22: [27.30, 30.81, 33.92, 36.78, 40.29, 48.27],
    23: [28.43, 32.01, 35.17, 38.08, 41.64, 49.73],
    24: [29.55, 33.20, 36.42, 39.36, 42.98, 51.18],
    25: [30.68, 34.38, 37.65, 40.65, 44.31, 52.62],
    26: [31.79, 35.56, 38.89, 41.92, 45.64, 54.05],
    27: [32.91, 36.74, 40.11, 43.19, 46.96, 55.48],
    28: [34.03, 37.92, 41.34, 44.46, 48.28, 56.89],
    29: [35.14, 39.09, 42.56, 45.72, 49.59, 58.30],
    30: [36.25, 40.26, 43.77, 46.98, 50.89, 59.70]
}


def chi_square_lookup(value, df):
    
    ps = [0.20, 0.10, 0.05, 0.025, 0.01, 0.001]
    
    if df <= 0:
        return 1.0    
    
    row = chi_square_table[min(df, 30)]

    for i in range(0,len(row)):
        if row[i] >= value:
            i = i-1
            break
    
    if i == -1: return 1
    else: return ps[i]



def spearman(vec1, vec2):
    """Spearman's rank test"""
    
    assert len(vec1) == len(vec2), "vec1 and vec2 are not the same length"
    
    n = len(vec1)
    rank1 = util.sortranks(vec1)
    rank2 = util.sortranks(vec2)
    
    R = sum((vec1[i] - vec2[i])**2 for i in xrange(n))
    
    Z = (6*R - n*(n*n - 1)) / (n*(n + 1) * sqrt(n - 1))
    
    return Z
    


# input:
#   xdata, ydata  - data to fit
#   func          - a function of the form f(x, params)
#
def fitCurve(xdata, ydata, func, paramsInit):   
    import scipy
    import scipy.optimize

    y = scipy.array(ydata)
    p0 = scipy.array(paramsInit)
    
    def error(params):
        y2 = scipy.array(map(lambda x: func(x, params), xdata))
        return y - y2

    params, msg = scipy.optimize.leastsq(error, p0)
    
    resid = error(params)
    
    return list(params), sum(resid*resid)

    
def fitDistrib(func, paramsInit, data, start, end, step, perc=1.0):
    xdata, ydata = util.distrib(data, low=start, width=step)
    ydata = [i / perc for i in ydata]
    xdata = util.histbins(xdata)
    params, resid = fitCurve(xdata, ydata, func, paramsInit)
    return params, resid
    

def plotfuncFit(func, paramsInit, xdata, ydata, start, end, step, plot = None,
                **options):
    from rasmus import gnuplot

    if not plot:
        plot = gnuplot.Gnuplot()
    
    options.setdefault('style', 'boxes')
    
    params, resid = fitCurve(xdata, ydata, func, paramsInit)
    plot.plot(util.histbins(xdata), ydata, **options)
    plot.plotfunc(lambda x: func(x, params), start, end, step)
    
    return plot, params, resid
    

def plotdistribFit(func, paramsInit, data, start, end, step, plot = None,
                   **options):
    xdata, ydata = util.distrib(data, low=start, width=step)
    return plotfuncFit(func, paramsInit, xdata, ydata, start, end, step/10, plot,
                       **options)




def chi_square_fit(cdf, params, data, ndivs=20, minsamples=5, plot=False,
                   start=-util.INF, end=util.INF):

    from rasmus import gnuplot
    import scipy
    import scipy.stats

    # determine ndiv and binsize
    binsize = len(data) / ndivs
    if binsize < minsamples:
        ndivs = len(data) / minsamples
        binsize = len(data) / ndivs

    data = sorted(data)
    bins = [data[i:i+binsize] for i in xrange(0, len(data), binsize)]
    obs = scipy.array(map(len, bins))
    ind = util.find(lambda x: x[-1] >= start and x[0] <= end, bins)
    obs = util.mget(obs, ind)
    
    x = [bin[0] for bin in bins]
    expected = [len(data) * cdf(x[1], params)]
    expected.extend([len(data) *
                     (cdf(x[i+1], params) - cdf(x[i], params))
                     for i in range(1, len(x)-1)])
    expected.append(len(data) * (1.0 - cdf(x[-1], params)))
    expected = scipy.array(util.mget(expected, ind))
    
    chi2, pval = scipy.stats.chisquare(obs, expected)

    if plot:        
        p = gnuplot.plot(util.mget(x, ind), obs)
        p.plot(util.mget(x, ind), expected)
    
    return chi2, pval


def fit_distrib(cdf, params_init, data, ndivs=20, minsamples=5,
                start=-util.INF, end=util.INF):

    import scipy
    import scipy.optimize
    import scipy.stats

    # determine ndiv and binsize
    binsize = len(data) / ndivs
    if binsize < minsamples:
        ndivs = len(data) / minsamples
        binsize = len(data) / ndivs

    data = sorted(data)
    bins = [data[i:i+binsize] for i in xrange(0, len(data), binsize)]
    obs = scipy.array(map(len, bins))
    ind = util.find(lambda x: x[-1] >= start and x[0] <= end, bins)
    obs = util.mget(obs, ind)
    
    def optfunc(params):
        x = [bin[0] for bin in bins]
        expected = [len(data) * cdf(x[1], params)]
        expected.extend([len(data) *
                         (cdf(x[i+1], params) - cdf(x[i], params))
                         for i in range(1, len(x)-1)])
        expected.append(len(data) * (1.0 - cdf(x[-1], params)))
        expected = scipy.array(util.mget(expected, ind))
        
        chi2, pval = scipy.stats.chisquare(obs, expected)
        return chi2

    params = scipy.optimize.fmin(optfunc, params_init, disp=False)
    chi2, pval = chi_square_fit(cdf, params, data, ndivs, minsamples)

    return list(params), pval

    



def solveCubic(a, b, c, real=True):
    """solves x^3 + ax^2 + bx + c = 0 for x"""
    
    p = b - a*a / 3.0
    q = c + (2*a*a*a - 9*a*b) / 27.0
    
    # special case: avoids division by zero later on
    if p == q == 0:
        return [- a / 3.0]
    
    # 
    # u = (q/2 +- sqrt(q^2/4 + p^3/27))^(1/3)
    #
    
    # complex math is used to find complex roots
    sqrteqn = cmath.sqrt(q*q/4.0 + p*p*p/27.0)
    
    # find fist cube root
    u1 = (q/2.0 + sqrteqn)**(1/3.0)
    
    # special case: avoids division by zero later on
    if u1 == 0:
        u1 = (q/2.0 - sqrteqn)**(1/3.0)
    
    # find other two cube roots
    u2 = u1 * complex(-.5, -sqrt(3)/2)
    u3 = u1 * complex(-.5, sqrt(3)/2)
    
    # finds roots of cubic polynomial
    root1 = p / (3*u1) - u1 - a / 3.0
    root2 = p / (3*u2) - u2 - a / 3.0
    root3 = p / (3*u3) - u3 - a / 3.0
    
    if real:
        return [x.real 
                for x in [root1, root2, root3]
                if abs(x.imag) < 1e-10]
    else:
        return [root1, root2, root3]


def _solveCubic_test(n=100):

    def test(a, b, c):
        xs = solveCubic(a, b, c)

        for x in xs:
            y = x**3 + a*x*x + b*x + c
            assert abs(y) < 1e-4, y

    test(0, 0, 0)
    test(0, 1, 1)
    test(0, 0, 1)

    for i in xrange(n):
        
        a = random.normalvariate(10, 5)
        b = random.normalvariate(10, 5)
        c = random.normalvariate(10, 5)

        test(a, b, c)


def bisect_root(f, x0, x1, err=1e-7):
    """Find a root of a function func(x) using the bisection method"""
    f0 = f(x0)
    f1 = f(x1)
    
    while (x1 - x0) / 2.0 > err:
        x2 = (x0 + x1) / 2.0
        f2 = f(x2)
        
        if f0 * f2 > 0:
            x0 = x2
            f0 = f2
        else:
            x1 = x2
            f1 = f2

    return (x0 + x1) / 2.0





#=============================================================================
# testing
    
if __name__ == "__main__":
    

    # iter_window
    from rasmus import util
    from rasmus import gnuplot

    vals = sorted([random.random() * 20 for x in range(600)])

    vals += sorted([40 + random.random() * 20 for x in range(600)])

    '''    
    win = filter(lambda x: len(x) > 0,
                 list(iter_window_index(vals, 5)))

    p = util.plot(util.cget(win, 2))#, style="lines")
    p.enableOutput(False)
    p.plot(util.cget(win, 3)) #, style="lines")    

    for i, y in enumerate(vals):
        p.plot([i, len(vals)], [y, y], style="lines")
    p.enableOutput(True)
    p.replot()
    '''

    def mean2(v):
        if len(v) == 0:
            return 0.0
        else:
            return mean(v)

    x, y = zip(* iter_window_step(vals, 5, 1, len))
    gnuplot.plot(x, y)
    



#=============================================================================
# OLD CODE

'''
def smooth_old(x, radius):
    """
    return an averaging of vals using a radius
    
    Note: not implemented as fast as possible
    runtime: O(len(vals) * radius)
    """
    
    vlen = len(x)
    
    # simple case
    if vlen == 0:
        return []
    
    x2 = []
    
    tot = x[0]
    
    low = 0
    high = 0
    
    for i in range(vlen):
        xi = x[i]
    
        xradius2 = min(i, vlen - i - 1, xradius)
    
        # move window
        while x[low] < xi - xradius2:
            xtot -= x[low]
            ytot -= y[low]
            low += 1
        while x[high] < xi + xradius2:
            high += 1
            xtot += x[high]
            ytot += y[high]
        
        denom = float(high - low + 1)
        x2.append(xtot / denom)
        y2.append(ytot / denom)
    
    return x2, y2
'''
