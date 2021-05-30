import numpy as np
from scipy.sparse import diags, spdiags, vstack
from scipy.sparse.linalg import spsolve


def beads(y, d, fc, r, Nit, lam0, lam1, lam2, pen, conv=None):
    """

    Baseline estimation and denoising using sparsity (BEADS)
    The notation follows the rule in the original MATLAB code.

    INPUT
        y: Noisy observation.
        d: Filter order (d = 1 or 2).
        fc: Filter cut-off frequency (cycles/sample) (0 < fc < 0.5).
        r: Asymmetry ratio for penalty function.
        Nit: Number of iteration (usually 10 to 30 is enough).
        lam0, lam1, lam2: Regularization parameters.
        pen: Penalty function, 'L1_v1' or 'L1_v2'.

        Original feature in this translation.
        conv: Smoothing factor for differential matrix D. This stabilizes
              outputs with slightly diffrent values for lam1 and lam2.
              Must be integer. 3 to 5 is recommended.
              

    OUTPUT
        x: Estimated sparse-derivative signal.
        f: Estimated baseline.
        cost: Cost function history

    Reference:
        Chromatogram baseline estimation and denoising using sparsity (BEADS)
        Xiaoran Ning, Ivan W. Selesnick, Laurent Duval
        Chemometrics and Intelligent Laboratory Systems (2014)
        doi: 10.1016/j.chemolab.2014.09.014
        Available online 30 September 2014

    Maintainer:
        Kotaro Saito
    
    Initial translator:
        Hisao Chun-Yi


    """
    # The following parameter may be altered.
    EPS0 = 1e-6  # cost smoothing parameter for x (small positive value)
    EPS1 = 1e-6  # cost smoothing parameter for derivatives(small positive value)

    if pen is 'L1_v1':
        phi = lambda xx: np.sqrt(np.power(abs(xx), 2) + EPS1)
        wfun = lambda xx: 1. / np.sqrt(np.power(abs(xx), 2) + EPS1)
    elif pen is 'L1_v2':
        phi = lambda xx: abs(xx) - EPS1 * np.log(abs(xx) + EPS1)
        wfun = lambda xx: 1. / (abs(xx) + EPS1)
    else:
        ValueError('penalty must be L1_v1, L1_v2')
    
    #  equation (25)
    theta = lambda xx: sum(xx[(xx > EPS0)]) - r * sum(xx[(xx < -EPS0)]) \
                       + sum((1+r)/(4*EPS0) * xx[abs(xx) <= EPS0] ** 2 \
                       + (1-r)/2 * xx[abs(xx) <= EPS0] + EPS0*(1+r)/4)

    y = np.reshape(a=y, newshape=(len(y), 1))
    x = y
    cost = []
    N = len(y)
    A, B = BAfilt(d, fc, N)
    H = lambda xx: B.dot(linv(A, xx))
    D1, D2 = make_diff_matrices(N)
    D = vstack([D1, D2])  # scipy.sparse.vstack, not np.vstack
    BTB = B.transpose().dot(B)

    w = np.vstack(([lam1 * np.ones((N-1, 1)), lam2 * np.ones((N-2, 1))]))
    b = (1-r) / 2 * np.ones((N, 1))
    d = BTB.dot(linv(A, y)) - lam0 * A.transpose().dot(b)

    gamma = np.ones((N, 1))

    for i in range(1, Nit+1):
        # print('step: ', i)
        if type(conv) is int:
            diff = np.convolve(D.dot(x.squeeze()), np.ones(conv)/conv, mode='same')[:, np.newaxis]
        else:
            diff = D.dot(x)
        wf = w * wfun(diff)
        Lmda = spdiags(wf.transpose(), 0, 2 * N - 3, 2 * N - 3)

        k = np.array(abs(x) > EPS0)  # return index 1d
        gamma[~k] = ((1 + r) / 4) / abs(EPS0)
        gamma[k] = ((1 + r) / 4) / abs(x[k])
        Gamma = spdiags(gamma.transpose(), 0, N, N)

        M = 2 * lam0 * Gamma + (D.transpose().dot(Lmda)).dot(D).transpose()
        x = A.dot(linv(BTB + A.transpose().dot(M.dot(A)), d))

        a = y - x
        cost.append(
            0.5 * sum(abs(H(a)) ** 2)
            + lam0 * theta(x)
            + lam1 * sum(phi(np.diff(x.squeeze())))
            + lam2 * sum(phi(np.diff(x.squeeze(), 2))))
        pass

    f = y - x - H(y - x)

    return x.squeeze(), f.squeeze(), cost

def BAfilt(d, fc, N):
    """
     --- local function ----

    function [A, B] = BAfilt(d, fc, N)
     [A, B] = BAfilt(d, fc, N)

     Banded matrices for zero-phase high-pass filter.
     The matrices are 'sparse' data type in MATLAB.
     In this code, output matrices are scipy.sparse.dia_matrix.

     INPUT
       d  : degree of filter is 2d (use d = 1 or 2)
       fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
       N  : length of signal
    """

    b1 = [1, -1]
    for i in range(1, d):
        b1 = np.convolve(a=b1, v=[-1, 2, -1])

    b = np.convolve(a=b1, v=[-1, 1])

    omc = 2 * np.pi * fc
    t = np.power(((1 - np.cos(omc)) / (1 + np.cos(omc))), d)

    a = 1
    for i in range(1, d+1):  # for i = 1:d
        a = np.convolve(a=a, v=[1, 2, 1])
    
    a = b + t * a
    xa, xb = (a*np.ones((N, 1))).transpose(), (b*np.ones((N, 1))).transpose()
    dr = np.arange(-d, d+1)
    A = spdiags(xa, dr, N, N)  # A: Symmetric banded matrix
    B = spdiags(xb, dr, N, N)  # B: banded matrix

    return A, B


# left inverse
def linv(a, b):
    '''
    a: sparse matrix
    b: vector

    '''
    return spsolve(a.tocsc(), b).reshape(len(b), 1)

def make_diff_matrices(N):
    D1 = diags([-1, 1], [0, 1], shape=(N - 1, N))
    D2 = diags([1, -2, 1], [0, 1, 2], shape=(N - 2, N))

    return D1, D2
