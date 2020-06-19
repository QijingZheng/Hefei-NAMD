#!/usr/bin/env python

from subprocess import call
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle


def get_molecule_positions(theta=0, shrink=0.25, r=0.6):

    L = 1.0
    bond_ratio = 0.55
    phi = 103.99987509868838 / 180 * np.pi
    theta = theta / 180 * np.pi
    M = np.array([
        [np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]
    ])
    pos = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [np.cos(phi), np.sin(phi)],
    ])
    pos = np.dot(M, pos.T).T * shrink
    r1 = L / (bond_ratio + 1)
    r2 = r1 * bond_ratio

    return np.array([r1, r2, r2]) * shrink, pos


figure = plt.figure(
    figsize=(3.0, 4.0)
)
ax = plt.subplot()
ax.set_aspect('equal')
ax.axis('off')

############################################################
R0 = 1.0
C0 = Circle(
    (0.0, 0.0), radius=R0,
    facecolor='k',
    alpha=0.7,
)
C1 = Circle(
    (0.0, 0.0), radius=0.55 * R0,
    facecolor='w'
)

ax.add_patch(C0)
ax.add_patch(C1)

############################################################

MOLCOLORS = ['r', 'b', 'b']
R, P = get_molecule_positions(theta=-10)
for ii in range(3):
    cc = Circle(
        P[ii], radius=R[ii],
        facecolor=MOLCOLORS[ii],
        alpha=0.6
    )
    ax.add_patch(cc)

R, P = get_molecule_positions(theta=-120, shrink=0.15)
P += np.array([-0.2, -0.2])
for ii in range(3):
    cc = Circle(
        P[ii], radius=R[ii],
        facecolor=MOLCOLORS[ii],
        alpha=0.5
    )
    ax.add_patch(cc)

R, P = get_molecule_positions(theta=120, shrink=0.1)
P += np.array([0.2, -0.2])
for ii in range(3):
    cc = Circle(
        P[ii], radius=R[ii],
        facecolor=MOLCOLORS[ii],
        alpha=0.4
    )
    ax.add_patch(cc)
############################################################

NAMD = [x.upper() for x in 'namd']
start_angle = 60.0 / 180. * np.pi
for ii in range(4):
    tc = np.exp(1j * (start_angle + ii * np.pi / 2))
    cc = Circle(
        (tc.real, tc.imag), radius=0.4,
        facecolor='white',
        # alpha=0.9
    )
    ax.add_patch(cc)
    cc = Circle(
        (tc.real, tc.imag), radius=0.32,
        facecolor='red',
        alpha=0.5
    )
    ax.add_patch(cc)
    ax.text(tc.real, tc.imag - 0.05, NAMD[ii],
            ha="center",
            va="center",
            fontsize=52,
            # family='monospace',
            color='white',
            fontweight='bold',
            transform=ax.transData,
            fontname='Cooper Black',
            # bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            )

############################################################

ax.plot([-1.4, 1.4], [-1.4, -1.4], 'k-', lw=5.0)

############################################################
TEXT = 'HEFEI'
TPOS = np.linspace(-1.1, 1.1, 5, endpoint=True)
TROT = [-19.05462068,   1.69593734, -11.93570325,  -0.39501574,  22.31639076]
# TROT = [-16.59700324, -7.07637502, 29.74676562,  7.06633301, 22.90883817]
# TROT = np.random.uniform(-1, 1, 5) * 30
# print(TROT)

TCLR = plt.rcParams['axes.prop_cycle'].by_key()['color']

for ii in range(5):
    h = TPOS[ii]
    v = -2.05 + np.sin(ii * np.pi / 2) * 0.2
    # v = -1.95
    T = plt.text(h, v, TEXT[ii],
        ha="center",
        va="center",
        fontsize=52,
        # family='monospace',
        color='white', alpha=0.95,
        rotation=TROT[ii],
        fontweight='bold',
        transform=ax.transData,
        fontname='Cooper Black',
        # bbox=dict(pad=0.1, lw=0.0, boxstyle='circle', facecolor='green', alpha=0.6)
    )
    # ax.plot([h,], [v], marker='o', ms=5)

    # box = T.get_window_extent().inverse_transformed(ax.transData)
    # box_w = box.x1 - box.x0
    # box_h = box.y1 - box.y0
    # RR = Rectangle(
    #     (box.x0, box.y0 + 0.1), box_w, box_h,
    #     facecolor='green'
    # )
    # ax.add_patch(RR)

    cc = Circle(
        (h, v + 0.05), radius=0.3,
        # facecolor='green',
        facecolor=TCLR[ii],
        alpha=0.6
    )
    ax.add_patch(cc)


ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-2.5, 1.4)

plt.tight_layout(pad=0.1)
plt.savefig('logo.png', dpi=300)
# plt.show()

# from subprocess import call
# call('feh -xdF logo.png'.split())
