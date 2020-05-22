import numpy as np
import AAKwrapper as AAK

def P_OMS(f):
    return (1.5e-11)**2*(1+(2e-3/f)**4)

def P_acc(f):
    return (3e-15)**2*(1+(0.4e-3/f)**2)*(1+(f/8e-3)**4)

def S_n(f, L=2.5e9, f_star=19.09e-3):
    return 10/(3*L**2)*(P_OMS(f)+4*P_acc(f)/(2*np.pi*f)**4)*(1+6/10*(f/f_star)**2)

fit_pars = {
    0.5:{
        'alpha': 0.133,
        'beta': 243,
        'kappa': 482,
        'gamma': 917,
        'f_k': 2.58e-3
        },
    1:{
        'alpha': 0.171,
        'beta': 292,
        'kappa': 1020,
        'gamma': 1680,
        'f_k': 2.15e-3
        },
    2:{
        'alpha': 0.165,
        'beta': 299,
        'kappa': 611,
        'gamma': 1340,
        'f_k': 1.73e-3
        },
    4:{
        'alpha': 0.138,
        'beta': -221,
        'kappa': 521,
        'gamma': 1680,
        'f_k': 1.13e-3
        }
    }

def S_c(f, A=9e-45, dur=4):
    if dur not in [0.5,1,2,4]:
        raise ValueError("dur needs to be 0.5, 1, 2, or 4 years")
    alpha = fit_pars[dur]['alpha']
    beta = fit_pars[dur]['beta']
    kappa = fit_pars[dur]['kappa']
    gamma = fit_pars[dur]['gamma']
    f_k = fit_pars[dur]['f_k']
    return A*f**(-7/3)*np.exp(-f**alpha+beta*f*np.sin(kappa*f))*(1+np.tanh(gamma*(f_k-f)))

def LISA_Noise(f, conf=True, L=2.5e9, f_star=19.09e-3, A=9e-45, dur=4):
    if conf:
        return S_n(f, L=L, f_star=f_star)+S_c(f, A=A, dur=dur)
    else:
        return S_n(f, L=L, f_star=f_star)

def AAK_Wave(mu, M, a, p0, e0, i0, D):
	pars = {
		'backint': False,
        'LISA': False,
        'length': 3000000,
        'dt': 10,
        'p': p0,
        'T': 1,
        'f': 2e-3,
        'T_fit': 1,
        'mu': mu,
        'M': M,
        's': a,
        'e': e0,
        'iota': i0,
        'gamma': 0,
        'psi': 0,
        'theta_S': np.pi/3,
        'phi_S': 0,
        'theta_K': np.pi/4,
        'phi_K': np.pi/4,
        'alpha': 0,
        'D': D
        }
	t, hI, hII, timing = AAK.wave(pars)
	return np.vstack([t,hI,hII])