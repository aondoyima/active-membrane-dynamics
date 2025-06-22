#This code uses an implicit Crank-Nicholson-Adams-Bashford algorithm to solve the coupled equations for two interacting enzymes on a 1D membrane. 

import numpy as np
import os
import pickle
import pbc_advance
from scipy import integrate
from types import SimpleNamespace as nt
from scipy.fft import fft, ifft
import pbc_utils

#Get parameters from argparser
args = pbc_utils.get_args()

#Global options
vol_conserve = bool(args.vc)
data_save = bool(args.ds)
adapt_step = bool(args.ad)

#Algorithm parameters
L, N, dt_min, dt_max, tmax, tol, factor_min, factor_max = args.syslength*np.pi, 256, 1e-12, 1e-1, 20000, 1e-5, 0.1, 10.

#Set resolution for varying system length
if L > 40*np.pi:
    N = int(args.res*args.syslength)

#Chemical parameters
a, b, K, chi, chiprime, alpha = args.a, args.b, args.K, args.chi, args.chiprime, args.alpha
#Surface parameters
s0, s1, s2, kappa, sigma, zeta = args.s0, args.s1, args.s2, args.kappa, args.sigma, args.zeta
#Initial condition parameters
psi1bar, psi2bar = args.psi1bar, args.psi2bar

#Real space 
dx = L / N
x = np.linspace(-L/2, L/2, N, endpoint=False) #Should not include end points - see Navdeep's notes

#Fourier space 
q = np.fft.fftfreq(N, d=dx) * 2 * np.pi 
q2 = q**2
q4 = q**4
q6 = q**6
dt = dt_min

qmax = q.max()*2/3
dealias = np.array(np.abs(q) < qmax)
# dealias = 1

#Initial noise level
initnoise = args.initnoise

#Directory structure
if vol_conserve:
    dirname = f'./results_volconsv_L={L:.2f}_zeta={zeta}/s0={s0}_s1={s1}_s2={s2}/chi={chi}_alpha={alpha}/a={a}_b={b}_K={K}_kappa={kappa}_chiprime={chiprime}_sigma={sigma}/'
else:
    dirname = f'./results_freevol_L={L:.2f}_zeta={zeta}/s0={s0}_s1={s1}_s2={s2}/chi={chi}_alpha={alpha}/a={a}_b={b}_K={K}_kappa={kappa}_chiprime={chiprime}_sigma={sigma}/'

os.makedirs(dirname, exist_ok=True)

#Save parameters
if data_save:
    params = {
            'x': x, 'a': a, 'chiprime': chiprime, 'K': K, 
            'sigma': sigma, 'kappa': kappa, 'chi': chi,
            'alpha': alpha, 's1': s1, 's2': s2, 'b': b,
            'zeta': zeta, 's0': s0, 
        }
    pickle.dump(params, open(f'{dirname}/params.p', 'wb'))

#Initial conditions
if vol_conserve:
    h = np.ones(N) + initnoise*(np.random.rand(N) - 0.5) 
else:
    h = np.zeros(N) + initnoise*(np.random.rand(N) - 0.5) 

psi1 = initnoise*(np.random.rand(N) - 0.5) + psi1bar
psi2 = initnoise*(np.random.rand(N) - 0.5) + psi2bar

# dirname_look = f'./results_volconsv_L={L:.2f}_zeta={zeta}/s0={s0}_s1={s1}_s2={s2}/chi={chi}_alpha={alpha}/a={a}_b={b}_K={K}_kappa={kappa}_chiprime={chiprime}_sigma={sigma}/'
# psi1, psi2, h = pbc_utils.twave(kappa,sigma,alpha,zeta,s1,dirname_look,q,x)

#Fourier transform fields
ftpsi1 = fft(psi1)
ftpsi2 = fft(psi2)
fth = fft(h)

#Volume conserving lagrange multiplier
lambda_v = 0.

#Store fields
#Create namesapce to store fields
sys_old = nt(psi1 = psi1, psi2 = psi2, h = h, ftpsi1 = ftpsi1, ftpsi2 = ftpsi2, fth = fth)
sys_now = nt(psi1 = [], psi2 = [], h = [], ftpsi1 = [], ftpsi2 = [], fth = [])
sys_new = nt(psi1 = [], psi2 = [], h = [], ftpsi1 = [], ftpsi2 = [], fth = [])

#Saving
framecounter = 0
next_save_time = 0. 
save_interval = 2.0 #make sure this is always much bigger than dt, or your saving will mess up.

#For adaptive step
counter_factor = 0
counter_max = 10

#Start time
t = 0. 

#Integrate initally with a one step method
sys_now = pbc_advance.stepCN(
                             dt, sys_old, q, q2, q4, q6, dealias,
                             a, b, K, chi, chiprime, alpha,
                             s0, s1, s2, kappa, sigma, zeta,
                            )
#t+=dt

while t < tmax:
    #Choose whether to use adaptive time stepping
    if not adapt_step:
        #Advance by one timestep
        sys_new = pbc_advance.stepCNAB(
                                       dt, sys_now, sys_old, q, q2, q4, q6, dealias,
                                       a, b, K, chi, chiprime, alpha,
                                       s0, s1, s2, kappa, sigma, zeta,
                                      )

        #Update system and time
        sys_old = nt(
        psi1 = sys_now.psi1.copy(),
        psi2 = sys_now.psi2.copy(),
        h = sys_now.h.copy(),
        ftpsi1 = sys_now.ftpsi1.copy(),
        ftpsi2 = sys_now.ftpsi2.copy(),
        fth = sys_now.fth.copy()
        )

        sys_now = nt(
        psi1 = sys_new.psi1.copy(),
        psi2 = sys_new.psi2.copy(),
        h = sys_new.h.copy(),
        ftpsi1 = sys_new.ftpsi1.copy(),
        ftpsi2 = sys_new.ftpsi2.copy(),
        fth = sys_new.fth.copy()
        )
        t += dt

    elif adapt_step:
        #Advance once by two timesteps
        sys_star = pbc_advance.stepCNAB(
                                        2*dt, sys_now, sys_old, q, q2, q4, q6, lambda_v, dealias,
                                        a, b, K, chi, chiprime, alpha,
                                        s0, s1, s2, kappa, sigma, zeta,
                                       )

        #Adance twice by one timestep 
        sys_temp = pbc_advance.stepCNAB(
                                        dt, sys_now, sys_old, q, q2, q4, q6, lambda_v, dealias,
                                        a, b, K, chi, chiprime, alpha,
                                        s0, s1, s2, kappa, sigma, zeta,
                                       )
        sys_new = pbc_advance.stepCNAB(
                                       dt, sys_temp, sys_now, q, q2, q4, q6, lambda_v, dealias,
                                       a, b, K, chi, chiprime, alpha,
                                       s0, s1, s2, kappa, sigma, zeta,
                                      )
        
        #Extract fields
        psi1_star, psi2_star, h_star = sys_star.psi1, sys_star.psi2, sys_star.h 
        psi1_new, psi2_new, h_new = sys_new.psi1, sys_new.psi2, sys_new.h

        #Calcualte error for adaptive step
        rel_err_psi1 = (1/3)*np.max(np.abs(psi1_new - psi1_star))/np.max(np.abs(psi1_new))
        rel_err_psi2 = (1/3)*np.max(np.abs(psi2_new - psi2_star))/np.max(np.abs(psi2_new))
        rel_err_h = (1/3)*np.max(np.abs(h_new - h_star))/np.max(np.abs(h_new))
        rel_err = max(rel_err_psi1, rel_err_psi2, rel_err_h)
        factor = (tol/rel_err)**(1/3)

        #Reject and try again with smaller timestep
        if rel_err > tol:
            # print(f'Trying with smaller time step dt = {dt:.2e}')
            if factor > 0.5:
                dt *= 0.5
            else:
                dt *= max(factor,factor_min)
                if dt < dt_min:
                    dt = dt_min
            continue
        
        #Accept and update system
        sys_old = nt(
        psi1 = sys_now.psi1.copy(),
        psi2 = sys_now.psi2.copy(),
        h = sys_now.h.copy(),
        ftpsi1 = sys_now.ftpsi1.copy(),
        ftpsi2 = sys_now.ftpsi2.copy(),
        fth = sys_now.fth.copy()
        )

        sys_now = nt(
        psi1 = sys_new.psi1.copy(),
        psi2 = sys_new.psi2.copy(),
        h = sys_new.h.copy(),
        ftpsi1 = sys_new.ftpsi1.copy(),
        ftpsi2 = sys_new.ftpsi2.copy(),
        fth = sys_new.fth.copy()
        )

        #Update volume conserving lagrange multiplier. It is fixed at zero if the volume conserve option is off.
        if vol_conserve:
            hxx = ifft(-q2*sys_now.fth).real
            psi1x = ifft(1j*q*sys_now.ftpsi1).real
            psi2x = ifft(1j*q*sys_now.ftpsi2).real
            lambda_v = 2*K*integrate.simpson(hxx*(psi1x**2 + psi2x**2),x=x)*fft(np.ones(N))/L

        #Update time and time step
        t += 2*dt
        dt *= min(factor,factor_max)
        if dt > dt_max:
            dt = dt_max
        if dt < dt_min:
            dt = dt_min

    if t >= next_save_time:  #this will save at the first instant that t exceeds next_save_time (next_save_time is updated immediately after saving)

        psi1, psi2, h = sys_now.psi1, sys_now.psi2, sys_now.h
        ftpsi1, ftpsi2, fth = sys_now.ftpsi1, sys_now.ftpsi2, sys_now.fth

        #Calculate average psi and dominant mode
        psi1avg = (1/L)*integrate.simpson(psi1,x=x)
        psi2avg = (1/L)*integrate.simpson(psi2,x=x)
        idx = np.argmax(sys_now.ftpsi1)

        #Print useful variables
        print(f't = {t:.2f}, dt = {dt:.2e}')
        print(f'chi = {chi:.2f}, chiprime = {chiprime:.2f} alpha = {alpha:.2f}, s0 = {s0}, s1 = {s1:.2f}, s2 = {s2:.2f}, a = {a:.2f}, b = {b:.2f}, kappa = {kappa:.2f}, zeta = {zeta:.2f}')
        print(f'rel_err = {rel_err:.2e}, factor = {factor:.3f}, psi1avg = {psi1avg:.2e}, pis2avg = {psi2avg:.2e}, Lselect = {2*np.pi/np.abs(q[idx]):.2f}')
        print(f'delta_psi1 = {np.linalg.norm(sys_now.psi1 - sys_old.psi1):.2e}, max_ftpsi1 = {np.max(np.abs(sys_now.ftpsi1)):.2e}, maxpsi1 = {np.max(psi1):.2f} \n')

        #Save the data 
        if data_save:
            data = {'t': t, 'x': x, 'psi1': psi1, 'psi2': psi2, 'h': h, 'ftpsi1': ftpsi1, 'ftpsi2': ftpsi2, 'fth': fth, 'psi1avg': psi1avg, 'psi2avg': psi2avg, 'dt': dt}
            pickle.dump(data, open(f'{dirname}/frame{framecounter}.p', 'wb'))
        
        framecounter += 1
        next_save_time += save_interval