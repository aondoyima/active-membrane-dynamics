import numpy as np
from scipy.fft import fft, ifft
from types import SimpleNamespace as nt
import pbc_utils

def stepCN(
          dt, sys, q, q2, q4, q6, dealias,
          a, b, K, chi, chiprime, alpha,
          s0, s1, s2, kappa, sigma, zeta,
          ):

    inv_zeta = 1./zeta

    A = 1 + 0.5*dt*(a*q2 + K*q4)
    D = 1 + inv_zeta*0.5*dt*(sigma*q2 + kappa*q4)

    #Unpack
    psi1, psi2 = sys.psi1, sys.psi2
    ftpsi1, ftpsi2, fth = sys.ftpsi1, sys.ftpsi2, sys.fth

    #Nonlinear terms
    Nh, Npsi1, Npsi2 = pbc_utils.nonlin_terms(psi1, psi2, ftpsi1, ftpsi2, fth, q, q2, a, b, K, chi, chiprime, s0, s1, s2, kappa)

    #Remove high wave number modes
    Nh *= dealias
    Npsi1 *= dealias
    Npsi2 *= dealias

    #Linear terms
    Rh = fth + inv_zeta*dt*(-0.5*(kappa*q4 + sigma*q2)*fth + 0.5*kappa*q2*(s1*ftpsi1 + s2*ftpsi2) + Nh)
    Rpsi1 = ftpsi1 + dt*(-0.5*(a*q2 + K*q4)*ftpsi1 - 0.5*(chi + alpha)*q2*ftpsi2 + 0.5*kappa*s1*q2*fth + Npsi1)
    Rpsi2 = ftpsi2 + dt*(-0.5*(a*q2 + K*q4)*ftpsi2 - 0.5*(chi - alpha)*q2*ftpsi1 + 0.5*kappa*s2*q2*fth + Npsi2)

    #Solution matrix
    M11 = A*D - 0.25*dt**2*kappa**2*s1**2*q6*inv_zeta
    M22 = A*D - 0.25*dt**2*kappa**2*s2**2*q6*inv_zeta
    M12 = 0.5*dt*(chi + alpha)*q2*D - 0.25*dt**2*kappa**2*s1*s2*q6*inv_zeta
    M21 = 0.5*dt*(chi - alpha)*q2*D - 0.25*dt**2*kappa**2*s1*s2*q6*inv_zeta

    Jpsi1 = Rpsi1*D + 0.5*dt*kappa*s1*q4*Rh
    Jpsi2 = Rpsi2*D + 0.5*dt*kappa*s2*q4*Rh

    denominator = M11*M22 - M12*M21
    numerator_psi1 = M22*Jpsi1 - M12*Jpsi2
    numerator_psi2 = M11*Jpsi2 - M21*Jpsi1

    ftpsi1_new = numerator_psi1/denominator
    ftpsi2_new = numerator_psi2/denominator
    fth_new = (Rh + 0.5*dt*kappa*q2*(s1*ftpsi1_new + s2*ftpsi2_new))

    #Transform back to real space
    h_new = ifft(fth_new).real
    psi1_new = ifft(ftpsi1_new).real
    psi2_new = ifft(ftpsi2_new).real

    sys_new = nt(psi1 = psi1_new, psi2 = psi2_new, h = h_new, ftpsi1 = ftpsi1_new, ftpsi2 = ftpsi2_new, fth = fth_new)
    
    return sys_new

def stepCNAB(
            dt, sys_now, sys_old, q, q2, q4, q6, lambda_v, dealias,
            a, b, K, chi, chiprime, alpha,
            s0, s1, s2, kappa, sigma, zeta,
            ):
    
    inv_zeta = 1./zeta

    A = 1 + 0.5*dt*(a*q2 + K*q4)
    D = 1 + inv_zeta*0.5*dt*(sigma*q2 + kappa*q4)
    
    #Unpack
    psi1, psi2 = sys_now.psi1, sys_now.psi2
    ftpsi1, ftpsi2, fth = sys_now.ftpsi1, sys_now.ftpsi2, sys_now.fth

    psi1_old, psi2_old = sys_old.psi1, sys_old.psi2
    ftpsi1_old, ftpsi2_old, fth_old = sys_old.ftpsi1, sys_old.ftpsi2, sys_old.fth

    #Nonlinear terms
    Nh_now, Npsi1_now, Npsi2_now = pbc_utils.nonlin_terms(psi1, psi2, ftpsi1, ftpsi2, fth, q, q2, a, b, K, chi, chiprime, s0, s1, s2, kappa)
    Nh_old, Npsi1_old, Npsi2_old = pbc_utils.nonlin_terms(psi1_old, psi2_old, ftpsi1_old, ftpsi2_old, fth_old, q, q2, a, b, K, chi, chiprime, s0, s1, s2, kappa)

    #Remove high wavenumber modes
    Nh_now *= dealias
    Npsi1_now *= dealias
    Npsi2_now *= dealias
    Nh_old *= dealias
    Npsi1_old *= dealias
    Npsi2_old *= dealias

    #Linear terms
    Rh = fth + inv_zeta*dt*(-0.5*(kappa*q4 + sigma*q2)*fth + 0.5*kappa*q2*(s1*ftpsi1 + s2*ftpsi2) + 0.5*(3*Nh_now - Nh_old) + lambda_v)
    Rpsi1 = ftpsi1 + dt*(-0.5*(a*q2 + K*q4)*ftpsi1 - 0.5*(chi + alpha)*q2*ftpsi2 + 0.5*kappa*s1*q2*fth + 0.5*(3*Npsi1_now - Npsi1_old))
    Rpsi2 = ftpsi2 + dt*(-0.5*(a*q2 + K*q4)*ftpsi2 - 0.5*(chi - alpha)*q2*ftpsi1 + 0.5*kappa*s2*q2*fth + 0.5*(3*Npsi2_now - Npsi2_old))

    #Solution matrix
    M11 = A*D - 0.25*dt**2*kappa**2*s1**2*q6*inv_zeta
    M22 = A*D - 0.25*dt**2*kappa**2*s2**2*q6*inv_zeta
    M12 = 0.5*dt*(chi + alpha)*q2*D - 0.25*dt**2*kappa**2*s1*s2*q6*inv_zeta
    M21 = 0.5*dt*(chi - alpha)*q2*D - 0.25*dt**2*kappa**2*s1*s2*q6*inv_zeta

    Jpsi1 = Rpsi1*D + 0.5*dt*kappa*s1*q4*Rh
    Jpsi2 = Rpsi2*D + 0.5*dt*kappa*s2*q4*Rh

    denominator = M11*M22 - M12*M21
    numerator_psi1 = M22*Jpsi1 - M12*Jpsi2
    numerator_psi2 = M11*Jpsi2 - M21*Jpsi1

    ftpsi1_new = numerator_psi1/denominator
    ftpsi2_new = numerator_psi2/denominator
    fth_new = (Rh + 0.5*dt*kappa*q2*(s1*ftpsi1_new + s2*ftpsi2_new))/D

    #Transform back to real space
    h_new = ifft(fth_new).real
    psi1_new = ifft(ftpsi1_new).real
    psi2_new = ifft(ftpsi2_new).real

    sys_new = nt(psi1 = psi1_new, psi2 = psi2_new, h = h_new, ftpsi1 = ftpsi1_new, ftpsi2 = ftpsi2_new, fth = fth_new)
    
    return sys_new

    