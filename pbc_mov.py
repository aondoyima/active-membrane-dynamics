import pickle
from types import SimpleNamespace as nt
import pbc_utils
import numpy as np
import argparse
import glob

#Get parameters from argparser
args = pbc_utils.get_args()

#Global options
vol_conserve = bool(args.vc)
vary_noise = bool(args.vn)

#Plotting options
mov = bool(args.mov)
kym = bool(args.kym)
struct = bool(args.struct)
kym_amp = bool(args.amp)

#Extract arguments
a, b, chiprime, K, sigma, kappa = args.a, args.b, args.chiprime, args.K, args.sigma, args.kappa
chi, alpha, s0, s1, s2, zeta = args.chi, args.alpha, args.s0, args.s1, args.s2, args.zeta
psi1bar, psi2bar = args.psi1bar, args.psi2bar

L = args.syslength*np.pi

print(f'chi = {chi}, alpha = {alpha}')

#Directory structure
if vol_conserve:
    dirname = f'./results_volconsv_L={L:.2f}_zeta={zeta}/s0={s0}_s1={s1}_s2={s2}/chi={chi}_alpha={alpha}/a={a}_b={b}_K={K}_kappa={kappa}_chiprime={chiprime}_sigma={sigma}/'
else:
    dirname = f'./results_freevol_L={L:.2f}_zeta={zeta}/s0={s0}_s1={s1}_s2={s2}/chi={chi}_alpha={alpha}/a={a}_b={b}_K={K}_kappa={kappa}_chiprime={chiprime}_sigma={sigma}/'

if vary_noise:
    dirname = f'./results_noise_L={L:.2f}_zeta={zeta}/s0={s0}_s1={s1}_s2={s2}/chi={chi}_alpha={alpha}/a={a}_b={b}_K={K}_kappa={kappa}_chiprime={chiprime}_sigma={sigma}/noise={args.initnoise}'

#Count the number of frames saved
files = sorted(glob.glob(dirname + '/frame*.p'))
nsteps = len(files) #Counting starts from zero so the index of the final file is nsteps - 1.

#Choose what time to start reading data from
last_step = pickle.load(open(f'{dirname}/frame{nsteps-1}.p','rb'))
t_max = last_step['t']
t_interval = t_max/(nsteps-1)
i = int(args.tread/t_interval)

#Create a namesapace compatible with pbc_utils.make_movie.
snap = nt(t = [], h = [], psi1 = [], psi2 = [], psi1avg = [], psi2avg = [])
while i < nsteps:
    data = pickle.load(open(f'{dirname}/frame{i}.p','rb'))
    snap.t.append(data['t'])
    snap.h.append(data['h'])
    snap.psi1.append(data['psi1'])
    snap.psi2.append(data['psi2'])
    snap.psi1avg.append(data['psi1avg'])
    snap.psi2avg.append(data['psi2avg'])
    i += 1

params = pickle.load(open(f'{dirname}/params.p','rb'))
x = params['x']

pars = nt(
    N = np.size(x), a = a, b = b, K = K, s0 = s0, s1 = s1, s2 = s2, zeta = zeta, L = L,
    sigma = sigma, kappa = kappa, chi = chi, alpha = alpha, chiprime = chiprime, psi1bar = psi1bar, psi2bar = psi2bar
    )

if kym:
    figname = f'{dirname}/kym.pdf'
    pbc_utils.make_kymograph(pars,figname,x,snap)

if kym_amp:
    figname = f'{dirname}/kym_amp.pdf'
    pbc_utils.make_kymograph_amp(pars,figname,x,snap)

if struct:
    figname = f'{dirname}/sf.pdf'
    pbc_utils.make_structure_factor(pars, figname, x, snap, args.t1, args.t2)

if mov:
    movname = f'{dirname}/mov.mp4'
    pbc_utils.make_movie(pars,movname,x,snap)



