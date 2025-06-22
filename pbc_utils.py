import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
import numpy as np
import pickle
import glob
import cmath
import argparse
from scipy.fft import fft, ifft

def make_movie(pars, movie, x, snap):
    #Create a figure with a side panel using gridspec
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])  # 3:1 ratio for plot and side panel
    ax_plot = plt.subplot(gs[0])  # Main plot area
    ax_panel = plt.subplot(gs[1])  # Side panel for parameters
    #List of parameters
    display_params = {
        r'$a$': pars.a,
        r'$b$': pars.b,
        r'$\chi^\prime$':pars.chiprime,
        r'$K$': pars.K,
        r'$\chi$': pars.chi,
        r'$\alpha$': pars.alpha,
        # r'$s_0$': pars.s0,
        r'$s_1$': pars.s1,
        r'$s_2$': pars.s2,
        r'$\zeta$': pars.zeta,
        r'$\sigma$': pars.sigma,
        r'$\kappa$': pars.kappa,
    }

    #Draw the side panel once
    ax_panel.axis("off")  # Turn off axes
    # ax_panel.set_title("Parameters", fontsize=16, pad=10)
    for i, (key, value) in enumerate(display_params.items()):
        ax_panel.text(0.05, 0.95 - i * 0.1, f"{key}: {value:.3f}", fontsize=16, transform=ax_panel.transAxes)

    #Initialize plot lines for h, psi1, and psi2
    line_h, = ax_plot.plot([], [], 'k-', label=r'$h(x)$', linewidth=3)  #h in black
    line_psi1, = ax_plot.plot([], [], 'g-', label=r'$\psi_1(x)$', linewidth=3)  #psi1 in green
    line_psi2, = ax_plot.plot([], [], 'b-', label=r'$\psi_2(x)$', linewidth=3)  #psi2 in blue

    ax_plot.set_xlim(x[0], x[-1])
    ax_plot.set_ylim(-3.0,3.0) 
    ax_plot.set_xlabel(r'$x$')
    ax_plot.set_ylabel(r'Field values')
    ax_plot.legend()
    ax_plot.grid()

    N = pars.N
    #Update function for animation
    def update(frame):
        print(f'alpha = {pars.alpha}, L = {pars.L:.2f}, frame {frame} out of {len(snap.t)}')
        # Update the data for each plot line
        line_h.set_data(x, snap.h[frame])
        line_psi1.set_data(x, snap.psi1[frame])
        line_psi2.set_data(x, snap.psi2[frame])
        # Update the title with the current time and averages of protein densitites
        psi1avg = snap.psi1avg[frame]
        psi2avg = snap.psi2avg[frame]
        ax_plot.set_title(f"Time: {snap.t[frame]:.2f}, $avg(\psi_1) = {psi1avg:.2e}$, $avg(\psi_2) = {psi2avg:.2e}$", fontsize=14)

        return line_h, line_psi1, line_psi2

    ani = FuncAnimation(fig, update, frames=len(snap.t), blit=True, repeat=True)
    ani.save(movie, writer="ffmpeg", fps=10)

def make_kymograph(pars, figname, x, snap):
    # Unpack parameters
    chi, alpha, s0, s1, s2, zeta = pars.chi, pars.alpha, pars.s0, pars.s1, pars.s2, pars.zeta
    a, b, chiprime, K, sigma, kappa = pars.a, pars.b, pars.chiprime, pars.K, pars.sigma, pars.kappa
    # Make panels
    fig = plt.figure(figsize=(9, 4))
    gs = gridspec.GridSpec(3, 2, width_ratios=[20, 1], height_ratios=[1, 1, 1], wspace=0.05, hspace=0.2)

    # Main axes
    ax0 = fig.add_subplot(gs[0, 0], )
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0, sharey=ax0)
    ax2 = fig.add_subplot(gs[2, 0], sharex=ax0, sharey=ax0)

    # Colorbar axes
    cax0 = fig.add_subplot(gs[0, 1])  # colorbar for h
    cax1 = fig.add_subplot(gs[1:, 1])  # shared colorbar for psi1 and psi2

    fig.suptitle(
        fr'$\chi$ = {chi}, $\alpha$ = {alpha}, $s0$ = {s0}, $s_1$ = {s1}, $s_2$ = {s2}, $\zeta$ = {zeta}, '
        fr'$a$ = {a}, $b$ = {b}, $\chi^\prime$ = {chiprime}, $K$ = {K}, $\sigma$ = {sigma}, $\kappa$ = {kappa}',
        fontsize=10
    )
    # Compute global vmin/vmax for consistent color scaling
    vmin = min(np.min(snap.psi1), np.min(snap.psi2))
    vmax = max(np.max(snap.psi1), np.max(snap.psi2))
    # vmin = -1.
    # vmax = 1.
    # Plot the fields
    kym_h = ax0.pcolormesh(snap.t, x, np.transpose(snap.h), shading='nearest', rasterized=True, vmin=np.min(snap.h), vmax=np.max(snap.h),cmap='Spectral')
    kym_psi1 = ax1.pcolormesh(snap.t, x, np.transpose(snap.psi1), shading='nearest', rasterized=True, vmin=vmin, vmax=vmax,cmap='Blues')
    kym_psi2 = ax2.pcolormesh(snap.t, x, np.transpose(snap.psi2), shading='nearest', rasterized=True, vmin=vmin, vmax=vmax,cmap='Blues')

    ax0.tick_params(labelbottom=False)
    ax1.tick_params(labelbottom=False)

    ax0.text(0.05, 1.02, r"$h(x,t)$", transform=ax0.transAxes, fontsize=9, va='bottom', ha='left')
    ax1.text(0.05, 1.02, r"$\psi_1(x,t)$", transform=ax1.transAxes, fontsize=9, va='bottom', ha='left')
    ax2.text(0.05, 1.02, r"$\psi_2(x,t)$", transform=ax2.transAxes, fontsize=9, va='bottom', ha='left')

    # Colorbars
    cb0 = fig.colorbar(kym_h, cax=cax0)
    cb1 = fig.colorbar(kym_psi2, cax=cax1)  # can use either psi1 or psi2 (same scale)

    # Axis labels
    ax2.set_xlabel(r'$t$')
    for ax in [ax0, ax1, ax2]:
        ax.set_ylabel(r'$x$')

    plt.savefig(figname)
    plt.close(fig) 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def make_kymograph_amp(pars, figname, x, snap):
    # Unpack parameters
    chi, alpha, s0, s1, s2, zeta = pars.chi, pars.alpha, pars.s0, pars.s1, pars.s2, pars.zeta
    a, b, chiprime, K, sigma, kappa = pars.a, pars.b, pars.chiprime, pars.K, pars.sigma, pars.kappa

    # Compute modulus squared of complex field
    density = np.array(snap.psi1)**2 + np.array(snap.psi2)**2

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 3))
    c = ax.pcolormesh(snap.t, x, density.T, shading='nearest', cmap='viridis', rasterized=True)

    fig.suptitle(
    fr'$\chi$ = {chi}, $\alpha$ = {alpha}, $s0$ = {s0}, $s_1$ = {s1}, $s_2$ = {s2}, $\zeta$ = {zeta}, '
    fr'$a$ = {a}, $b$ = {b}, $\chi^\prime$ = {chiprime}, $K$ = {K}, $\sigma$ = {sigma}, $\kappa$ = {kappa}',
    fontsize=10
    )

    # Labels
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$x$')
    ax.set_title(r'$\psi_1^2 + \psi_2^2$')

    # Colorbar
    fig.colorbar(c, ax=ax, label=r'$|\varphi|^2$')

    # Save and close
    plt.tight_layout()
    plt.savefig(figname)
    plt.close(fig)

def make_structure_factor(pars, figname, x, snap, t1, t2):
    # Unpack parameters
    chi, alpha, s0, s1, s2, zeta = pars.chi, pars.alpha, pars.s0, pars.s1, pars.s2, pars.zeta
    a, b, chiprime, K, sigma, kappa = pars.a, pars.b, pars.chiprime, pars.K, pars.sigma, pars.kappa

    # Time window indices
    t = np.array(snap.t)
    t1_idx = np.searchsorted(t, t1)
    t2_idx = np.searchsorted(t, t2)

    # Setup
    dx = x[1] - x[0]
    N = len(x)
    q = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    q = np.fft.fftshift(q)
    dq = q[1] - q[0]
    S = np.zeros(N)

    # Accumulate power spectrum over time
    for i in range(t1_idx, t2_idx):
        ftpsi1 = fft(snap.psi1[i])
        ftpsi2 = fft(snap.psi2[i])
        S += np.abs(ftpsi1)**2 + np.abs(ftpsi2)**2

    S /= (t2_idx - t1_idx)
    S = np.fft.fftshift(S)

    # Restrict to positive q ≤ 1
    q_mask = (q > 0) & (q <= 1)
    # q_mask = (q > 0) 
    q_plot = q[q_mask]
    S_plot = S[q_mask]

    # Plot
    plt.figure(figsize=(6.5, 4.5))
    plt.loglog(q_plot, S_plot, color='darkblue', linewidth=1.5)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$S(q)$")
    plt.grid(True, which="both", linestyle=":", alpha=0.5)
    plt.title(
        fr"$\chi$ = {chi}, $\alpha$ = {alpha}, $s_0$ = {s0}, $s_1$ = {s1}, $s_2$ = {s2}, $\zeta$ = {zeta}, "
        fr"$a$ = {a}, $b$ = {b}, $\chi^\prime$ = {chiprime}, $K$ = {K}, $\sigma$ = {sigma}, $\kappa$ = {kappa}",
        fontsize=10
    )
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

def twave(kappa,sigma,alpha,zeta,s,dirname,q,x):
    #Generates travelling wave profile e.g. to be used for initial condition 

    files = sorted(glob.glob(dirname + '/frame*.p'))
    nsteps = len(files)

    data = pickle.load(open(f'{dirname}/frame{nsteps-10}.p','rb'))
    ftpsi1 = fft(data['psi1'])
    idx = np.argmax(ftpsi1)

    qselect = q[idx]
    qselect = 0.5

    # phi0 = np.sqrt(1 - qselect**2) #can be plus or minus
    # #Helper functions:
    # cbrt2 = 2**(1/3)
    # fsq = (-qselect**2*kappa - sigma + 2*qselect**2*phi0**2 + 0.5*phi0 - qselect**2*phi0**2 - 0.25*phi0**4)**2
    # D = -18*alpha*fsq*zeta**4 + 9*alpha*kappa*2*qselect**2*s**2*zeta**5 - 2*alpha**3*zeta**6
    # N = -alpha**2*zeta**4 + 3*zeta**2*(fsq + kappa**2*qselect**2*s**2*zeta)

    # omega = -alpha*qselect**4/3 - qselect**2*cbrt2*N/(3*zeta**2*(D + np.sqrt(D**2 + 4*N**3))**(1/3)) + qselect**2*(D + np.sqrt(D**2 + 4*N**3))**(1/3)/(3*cbrt2*zeta**2)
    # h0 = -phi0*(omega + alpha*qselect**2)/(np.sqrt(2)*kappa*s*qselect**4)

    B = zeta - qselect**2 + 0.5
    C = zeta*(alpha - 1 + qselect**2) + qselect**2*kappa + sigma
    sqrt_term = cmath.sqrt(B**2 - 4*C)
    phi0sq_plus = 0.5 * (-B + sqrt_term)
    phi0 = cmath.sqrt(phi0sq_plus)

    g = -qselect**4*kappa - qselect**2*sigma + 2*qselect**4*phi0**2 + qselect**2*0.5*phi0**2 - qselect**4*phi0**2 - qselect**2*phi0**4
    omega = -g/zeta
    h0 = -phi0*(omega + alpha*qselect**2)/(np.sqrt(2)*kappa*s*qselect**4)

    psi1 = np.zeros(len(x))
    psi2 = np.zeros(len(x))
    h = np.zeros(len(x))
    t = 0
    for i in range(0,len(x)):
        x_single = x[i]  
        psi1[i] = phi0*np.cos(qselect*x_single - omega*t)
        psi2[i] = phi0*np.sin(qselect*x_single - omega*t)
        h[i] = h0*np.sin(qselect*x_single - omega*t)

    return psi1, psi2, h

def nonlin_terms(psi1, psi2, ftpsi1, ftpsi2, fth, q, q2, a, b, K, chi, chiprime, s0, s1, s2, kappa):

    #Compute real space first and second derivatives
    hx = ifft(1j*q*fth).real
    hxx = ifft(-q2*fth).real
    psi1x = ifft(1j*q*ftpsi1).real
    psi2x = ifft(1j*q*ftpsi2).real

    # Nonlinear terms in Fourier space
    Nh = 0. #Code will be updated to include the nonlinear h terms once the paper is available on arXiv
    Npsi1 = -q2*b*fft(psi1**3) - 2*q2*chiprime*fft(psi1*psi2**2)
    Npsi2 = -q2*b*fft(psi2**3) - 2*q2*chiprime*fft(psi2*psi1**2)

    return Nh, Npsi1, Npsi2

def get_args():
    parser = argparse.ArgumentParser()

    #Floats
    parser.add_argument('-a', type=float, default=-1.0, help='Parameter a')
    parser.add_argument('-b', type=float, default=1.0, help='Parameter b')
    parser.add_argument('-cp', '--chiprime', type=float, default=0.0, help='Fourth order interaction')
    parser.add_argument('-K', type=float, default=1.0, help='Parameter K')
    parser.add_argument('-s', '--sigma', type=float, default=0.01, help='Surface tension')
    parser.add_argument('-k', '--kappa', type=float, default=0.1, help='Curvature')
    parser.add_argument('-z', '--zeta', type=float, default=1.0, help='Friction')
    parser.add_argument('-chi', type=float, default=0.2, help='Reciprocal interaction')
    parser.add_argument('-alpha', type=float, default=0., help='Non-reciprocal interaction')
    parser.add_argument('-s1', type=float, default=0.1, help='Curvature coupling, species 1')
    parser.add_argument('-s2', type=float, default=0.1, help='Curvature coupling, species 2')
    parser.add_argument('-s0', type=float, default=0., help='Equilibrium spontaneous curvature')
    parser.add_argument('-psi1bar', type=float, default=0., help='average density of species 1')
    parser.add_argument('-psi2bar', type=float, default=0., help='avergae density of species 2')
    parser.add_argument('-syslength', type=float, default=40, help='system length')
    parser.add_argument('-initnoise', type=float, default=40, help='inital condition noise level')

    #Switches
    parser.add_argument('-vc', type=int, default=1, help='Volume conserve or no')
    parser.add_argument('-ds', type=int, default=1, help='Save data or no')
    parser.add_argument('-ad', type=int, default=1, help='Adapt time step or no')
    parser.add_argument('-vn', type=int, default=0, help='Vary initial noise step or no')
    parser.add_argument('-mov', type=int, default=0, help='Make movie')
    parser.add_argument('-kym', type=int, default=1, help='Make kymograph')
    parser.add_argument('-amp', type=int, default=1, help='Make kymograph')
    parser.add_argument('-struct', type=int, default=0, help='Plot structure factor')

    #Other options
    parser.add_argument('-t1', type=float, default=0, help='Start of time averaging')
    parser.add_argument('-t2', type=float, default=0, help='End of time averaging')
    parser.add_argument('-tread', type=float, default=0, help='When to start reading data from')
    parser.add_argument('-res', type=float, default=6.4, help='Number of grid points = resolution*syslength')

    return parser.parse_args()