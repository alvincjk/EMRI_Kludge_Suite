"""
Converter to convert parameters in a sampler to physical parameters
needed for the likelihood calculation.
"""
import numpy as np
from scipy import constants as ct
from scipy import stats
import warnings


class Converter:
    def __init__(self, key_order, **kwargs):

        self.conversions = []
        if "ln_mT" in key_order:
            # replace ln_mT and mr with m1 and m2 respectively
            self.ind_ln_mT = key_order.index("ln_mT")
            self.conversions.append(self.ln_mT)

        if "ln_mu" in key_order:
            # replace ln_mT and mr with m1 and m2 respectively
            self.ind_ln_mT = key_order.index("ln_mu")
            self.conversions.append(self.ln_mu)

        if "cos_iota" in key_order:
            self.ind_iota = key_order.index("cos_iota")
            self.conversions.append(self.cos_iota)

        if "sin_theta_S" in key_order:
            self.ind_theta_S = key_order.index("sin_theta_S")
            self.conversions.append(self.sin_theta_S)

        if "sin_theta_K" in key_order:
            self.ind_theta_K = key_order.index("sin_theta_K")
            self.conversions.append(self.sin_theta_K)

        if "ln_distance" in key_order:
            self.ind_ln_distance = key_order.index("ln_distance")
            self.conversions.append(self.ln_distance)

    def ln_mT(self, x):
        x[self.ind_ln_mT] = np.exp(x[self.ind_ln_mT])
        return x

    def ln_mu(self, x):
        x[self.ind_ln_mu] = np.exp(x[self.ind_ln_mu])
        return x

    def ln_distance(self, x):
        x[self.ind_ln_distance] = (
            np.exp(x[self.ind_ln_distance]) * 1e9 * ct.parsec
        )  # Gpc
        return x

    def cos_iota(self, x):
        x[self.ind_iota] = np.arccos(x[self.ind_iota])
        return x

    def sin_theta_S(self, x):
        x[self.ind_theta_S] = np.arcsin(x[self.ind_theta_S])
        return x

    def sin_theta_K(self, x):
        x[self.ind_theta_K] = np.arcsin(x[self.ind_theta_K])
        return x


class Recycler:
    def __init__(self, key_order, **kwargs):

        # Setup of recycler
        self.recycles = []

        if "phi_S" in key_order:
            # assumes beta is also there
            self.ind_phi_S = key_order.index("phi_S")
            self.recycles.append(self.phi_S)

        if "phi_K" in key_order:
            # assumes beta is also there
            self.ind_phi_K = key_order.index("phi_K")
            self.recycles.append(self.phi_K)

        if "gamma" in key_order:
            self.ind_gamma = key_order.index("gamma")
            self.recycles.append(self.gamma)

        if "psi" in key_order:
            self.ind_psi = key_order.index("psi")
            self.recycles.append(self.psi)

        if "alph" in key_order:
            self.ind_alph = key_order.index("alph")
            self.recycles.append(self.alph)

    def phi_S(self, x):
        x[self.ind_phi_S] = x[self.ind_phi_S] % (2 * np.pi)
        return x

    def phi_K(self, x):
        x[self.ind_phi_K] = x[self.ind_phi_K] % (2 * np.pi)
        return x

    def phiRef(self, x):
        x[self.ind_phiRef] = x[self.ind_phiRef] % (2 * np.pi)
        return x

    def gamma(self, x):
        x[self.ind_gamma] = x[self.ind_gamma] % (2 * np.pi)
        return x

    def psi(self, x):
        x[self.ind_psi] = x[self.ind_psi] % (2 * np.pi)
        return x

    def alph(self, x):
        x[self.ind_alph] = x[self.ind_alph] % (2 * np.pi)
        return x

    def recycle(self, x):
        for func in self.recycles:
            x = func(x)
        return x
