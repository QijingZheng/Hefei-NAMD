#!/usr/bin/env python

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def gaussian(x,c):
    return np.exp(-x**2/(2*c**2))

def dephase(Et, dt=1.0):
    '''
    Calculate the autocorrelation function (ACF), dephasing function, and FT of
    ACF.

    The dephasing function was calculated according to the following formula

    G(t) = (1 / hbar**2) \int_0^{t_1} dt_1 \int_0^{t_2} dt_2 <E(t_2)E(0)>
    D(t) = exp(-G(t))

    where Et is the difference of two KS energies in unit of eV, <...> is the
    ACF of the energy difference and the brackets denote canonical averaging.

    Fourier Transform (FT) of the normalized ACF gives the phonon influence
    spectrum, also known as the phonon spectral density.

    I(\omega) \propto | FT(Ct / Ct[0]) |**2

    Jaeger, Heather M., Sean Fischer, and Oleg V. Prezhdo. "Decoherence-induced surface hopping." JCP 137.22 (2012): 22A545.
    '''

    from scipy.integrate import cumtrapz
    from scipy.fftpack import fft

    hbar = 0.6582119513926019       # eV fs

    Et = np.asarray(Et)
    Et -= np.average(Et)

    # Autocorrelation Function (ACF) of Et
    Ct = np.correlate(Et, Et, 'full')[Et.size:] / Et.size
    
    # Cumulative integration of the ACF
    Gt = cumtrapz(cumtrapz(Ct, dx=dt, initial=0), dx=dt, initial=0)
    Gt /= hbar**2
    # Dephasing function
    Dt = np.exp(-Gt)

    # FT of normalized ACF
    Iw = np.abs(fft(Ct / Ct[0]))**2

    return Ct, Dt, Iw

energy = np.loadtxt('energy.dat')
nbasis = energy.shape[1]
matrix = np.zeros((nbasis, nbasis), dtype=float)

dt = 1.0 # fs

for ii in range(nbasis):
    for jj in range(ii):
        Et = energy[:, ii] - energy[:, jj]
        T = np.arange(Et.size-1) * dt

        Ct, Dt, Iw = dephase(Et)
        popt, pcov = curve_fit(gaussian, T, Dt)
        Dt_fit = gaussian(T, *popt)
        matrix[ii,jj] =  popt[0]
        matrix[jj,ii] =  matrix[ii,jj]

np.savetxt('DEPHTIME', matrix, fmt='%10.4f')
# inp = np.loadtxt('kaka.dat')
# Et = (inp[:,1] - inp[:,0])
# dat = Et / Et[0] * 0.4
# from scipy.interpolate import interp1d
# Tnew = np.arange(1993)
#
# func = interp1d(np.arange(dat.size) * 2000 / dat.size, Et)
# Et = func(Tnew) / 100
#
# # Et = np.loadtxt('delt_E.dat')
#
# Ct, Gt, Dt = dephase(Et)
#
# np.savetxt('dep3.dat', Dt)
#
# dt = 1.0 # fs
# T = np.arange(Et.size) * dt
#
# ############################################################
# import matplotlib.pyplot as plt
# import matplotlib as mpl
#
# mpl.rcParams['axes.unicode_minus'] = False
#
# fig, axes = plt.subplots(nrows=3, ncols=1,
#                          sharex=False,
#                          sharey=False)
# plt.subplots_adjust(left=0.15, right=0.95,
#                     bottom=0.10, top=0.95,
#                     wspace=0.10, hspace=0.20)
# fig.set_size_inches((4.5, 6))
#
# inp = [Et, Ct, Dt]
# ylabs = [r'$\Delta$E (eV)', 'ACF', 'Dephasing Function']
# # xlims = [(0, T.max()), (0, T.max()), (0, 50)]
# clrs = ['k', 'r', 'b']
#
# for ii in range(3):
#     ax = axes.flatten()[ii]
#     dat = inp[ii]
#
#     ax.minorticks_on()
#
#     N = min(T.size, dat.size)
#     ax.plot(T[:N], dat[:N], ls='-', lw=1.0,
#             # marker='o', markevery=3, mew=0.0, ms=3,
#             color=clrs[ii])
#
#     if ii == 2:
#         ax.set_xlabel('Time [fs]', labelpad=10, fontsize='small')
#     ax.set_ylabel(ylabs[ii], labelpad=5, fontsize='small') 
#
#     if ii == 2:
#         ax.set_xlim(0, 100)
#
#     ax.tick_params(which='both', labelsize='x-small')
#
# plt.savefig('kaka.png', dpi=300)
# # plt.show()
