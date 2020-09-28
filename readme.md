# BEADS: Baseline Estimation And Denoising using Sparsity 
This is a Python translation of an awesome baseline estimation algorithm "BEADS" [originally written in MATLAB](https://jp.mathworks.com/matlabcentral/fileexchange/49974-beads-baseline-estimation-and-denoising-with-sparsity). The initial translation had been done by @hsiaocy. 

# Installation
Simply run `pip install pybeads` in your terminal.

# Usage

```
import pybeads as be

y = (your data containing signal, noise and background)

# try changing fc and lam0-2, amp if these dont' fit your data
fc = 0.006
d = 1
r = 6
amp = 0.8
lam0 = 0.5 * amp
lam1 = 5 * amp
lam2 = 4 * amp
Nit = 15
pen = 'L1_v2'

signal_est, bg_est, cost = be.beads(y, d, fc, r, Nit, lam0, lam1, lam2, pen, conv=None)
```

```
# Quick visualization code
signal_est, bg_est, cost = be.beads(y, d, fc, r, Nit, lam0, lam1, lam2, pen, conv=None)
fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
fig.subplots_adjust(hspace=0)
axes[0].plot(y, c='k', label='original data')
axes[0].plot(bg_est, c='r', label='BG estimated by BEADS')
axes[0].legend()
axes[1].plot(signal_est, label='signal estimated by BEADS')
axes[1].legend()
axes[2].plot(y-signal_est-bg_est, label='noise estimated by BEADS')
axes[2].legend()
```

Parameter notation and usage are same as the original MATLAB code's. Please see "Example" tab of [this page](https://jp.mathworks.com/matlabcentral/fileexchange/49974-beads-baseline-estimation-and-denoising-with-sparsity). If you want to know the details of the parameters, please refer to the original paper in Resources.

## Sample data
Real chromatograms with eight different background levels (probably including the one tested in the original paper shown in Figure 12) and additional computed white noise are available in MATLAB format at https://jp.mathworks.com/matlabcentral/fileexchange/49974-beads-baseline-estimation-and-denoising-with-sparsity, but MathWork account is required to download the zip file. For those who don't want to create a MathWork account just for this, I converted them into CSV format and included in this Python repo (`sample_data/chromatograms_and_noise.csv`). You can play with this test data with Jupyter notebook `pybeads_demo.ipynb`.

_Redistributed with permission by courtesy of Laurent Duval_. 

## A new parameter, `conv`
The main function `beads` in this package has a parameter called `conv` which does not exist in the MATLAB code. This parameter is used for a smoothing feature for derivatives (moving average using convolution, to be correct). I noticed that the MATLAB implementation sometimes gives completely different results for regularization parameters with a slight difference, for example `lam1=0.001` and `lam1=0.0011`. This unpredictable behavior is suppressed when derivatives are smoothed. When you face such instability, I reccomend you to set `conv` to 3 or 5. Default is `None` which means "no smoothing".

## A tip for real data
When you apply BEADS to your experimetal data, probably you'll find something is wrong at the both ends of your data which does not happen to the example data. The trick in the exmaple is that the both ends are smoothly decreasing to zero. To solve the issue, all you have to do is to do the same trick to your data by extending the both ends smoothly to zero. I reccomend to use a sigmoid function with an appropriate scaling factor for x, e.g. 30 to 100 or even bigger depending on the number of points of your data. Below is an exmaple.

```
import numpy as np

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# y is your data

xscale_l, xscale_r = 100, 100
dx = 1
y_l = y[0]*sigmoid(1/xscale_l*np.arange(-5*xscale_l, 5*xscale_l, dx))
y_r = y[-1]*sigmoid(-1/xscale_r*np.arange(-5*xscale_r, 5*xscale_r, dx))
y_ext = np.hstack([y_l, y, y_r])
len_l, len_o, len_r = len(y_l), len(y), len(y_r)
plt.plot(range(len_l, len_l+len_o), y)
plt.plot(y_l, 'C1')
plt.plot(range(len_l+len_o, len_l+len_o+len_r), y_r, 'C1')
```

# Resources
- [Original paper (2014)](https://doi.org/10.1016/j.chemolab.2014.09.014)
- Preprint on [laurent-duval.eu](http://www.laurent-duval.eu/Articles/Ning_X_2014_j-chemometr-intell-lab-syst_chromatogram_bedusbeads-preprint.pdf)
- [Project website](http://www.laurent-duval.eu/siva-beads-baseline-background-removal-filtering-sparsity.html)
- [Original MATLAB code](https://jp.mathworks.com/matlabcentral/fileexchange/49974-beads-baseline-estimation-and-denoising-with-sparsity)
