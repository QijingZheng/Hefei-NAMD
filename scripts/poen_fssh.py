#!/usr/bin/env python
############################################################
import os, re
import numpy as np
from glob import glob

############################################################
def EnergyFromPro(infile='PROCAR'):
    """
    Extract band energies from PROCAR.
    """
    print infile 
    assert os.path.isfile(infile), '%s cannot be found!' % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # when the band number is too large, there will be no space between ";" and
    # the actual band number. A bug found by Homlee Guo.
    # Here, #kpts, #bands and #ions are all integers
    nkpts, nbands, nions = [int(xx) for xx in re.sub('[^0-9]', ' ', FileContents[1]).split()]

    energies = np.asarray([line.split()[-4] for line in FileContents
                            if 'occ.' in line], dtype=float)
    nspin = energies.shape[0] / (nkpts * nbands)
    energies.resize(nspin, nkpts, nbands)

    return energies

def parallel_energy(runDirs, nproc=None):
    '''
    calculate localization of some designated in parallel.
    '''
    import multiprocessing
    nproc = multiprocessing.cpu_count() if nproc is None else nproc
    pool = multiprocessing.Pool(processes=nproc)

    results = []
    for rd in runDirs:
        res = pool.apply_async(EnergyFromPro, (rd + '/PROCAR',))
        results.append(res)

    en = []
    for ii in range(len(results)):
        tmp_en = results[ii].get()
        en.append(tmp_en)

    return np.array(en)

############################################################
# calculate spatial localization
############################################################
nsw     = 7000
nproc   = 8
prefix  = 'NAMD/run/'
runDirs = [prefix + '{:04d}'.format(ii + 1) for ii in range(nsw)]

if os.path.isfile('all_en.npy'):
    energies = np.load('all_en.npy')
else:
    # for gamma point version, no-spin
    energies = parallel_energy(runDirs, nproc=nproc)
    print energies.shape
    energies = energies[:, 0,0, :]

    np.save('all_en.npy', energies)

############################################################
# load FSSH result files
############################################################

########################################
bmin     = 152
bmax     = 236
namdTime = 6000
potim    = 1.0
inpFiles = glob('NAMD/SHPROP.*')
########################################

if not os.path.isfile('po.npy'):

    iniTimes = [int(F.split('.')[-1]) for F in inpFiles]

    dat = np.array([np.loadtxt(F) for F in inpFiles])
    dat = np.average(dat, axis=0)

    Ci_t = dat[:,2:]
    En_t = np.zeros_like(Ci_t)
    Time = np.zeros_like(Ci_t)
    #
    energies = energies[:, bmin:bmax]
    EVBM = np.average(energies[:,-1])

    for start in iniTimes:
        end = start + namdTime
        En_t += energies[start:end,:] - EVBM
    else:
        En_t /= len(iniTimes)

    for i in range(bmax - bmin):
        Time[:,i] = np.arange(namdTime)

    dat[:,1] -= EVBM
    x, y = dat[:,:2].T

    np.save('po.npy', (Time, En_t, Ci_t))
    np.savetxt('average_energy.dat', dat[:,:2], fmt='%8.4f')

else:
    inp = np.load('po.npy')
    Time = inp[0,:]
    En_t = inp[1,:]
    Ci_t = inp[2,:]
    x, y = np.loadtxt('average_energy.dat', unpack=True)

############################################################
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['axes.unicode_minus'] = False

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure()
fig.set_size_inches(4.8, 3.0)

ax  = plt.subplot()
divider = make_axes_locatable(ax)
ax_cbar = divider.append_axes('right', size="5%", pad=0.02)

############################################################
line, = ax.plot(x, y,
                ls='--',
                color='blue',
                lw=1.5,
                alpha=0.6)

kmap = ax.scatter(Time, En_t, c=Ci_t,
                  cmap='hot_r',
                  vmin=0,
                  vmax=1,
                  s=15, alpha=0.8, lw=0.0)
cbar = plt.colorbar(kmap, cax=ax_cbar,
                    orientation='vertical',
                    ticks=np.linspace(0, 1, 6, endpoint=True))

ax.legend([line,], ['Average Hole Energy', ],
          fancybox=True,
          loc='lower right',
          framealpha=0.7,
          fontsize=9)

############################################################
ax.set_xlim(0, namdTime)
# ax.set_ylim(-1.0, 1.0)

ax.set_xlabel('Time [fs]',   fontsize='small', labelpad=5)
ax.set_ylabel('Energy [eV]', fontsize='small', labelpad=5)

########################################
plt.tight_layout(pad=0.2)
plt.savefig('kpoen.png', dpi=360)
