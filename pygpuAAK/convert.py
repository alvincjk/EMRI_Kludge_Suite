"""
Converter to convert parameters in a sampler to physical parameters
needed for the likelihood calculation.
"""
import numpy as np
from scipy import constants as ct
from scipy import stats
from astropy import units
import warnings


def modpi(phase):
    # from sylvan
    return phase - np.floor(phase / np.pi) * np.pi


def mod2pi(phase):
    # from sylvan
    return phase - np.floor(phase / (2 * np.pi)) * 2 * np.pi


source_recycle_guide = {
    "emri": {
        "phi_S": 2 * np.pi,
        "phi_K": 2 * np.pi,
        "phiRef": 2 * np.pi,
        "gamma": 2 * np.pi,
        "psi": 2 * np.pi,
        "alpha": 2 * np.pi,
    },
}


class Converter:
    def __init__(
        self,
        source_type,
        key_order,
        dist_unit="Gpc",
        t0=None,
        transform_frame=None,
        **kwargs
    ):

        variables_to_recycle = source_recycle_guide[source_type].keys()
        recycle_guide = source_recycle_guide[source_type]

        if t0 is not None:
            self.t0 = t0 * YRSID_SI

        self.dist_unit = getattr(units, dist_unit)
        self.dist_conversion = self.dist_unit.to(units.m)

        self.recycles = []
        self.conversions = []
        self.inds_list = []
        self.keys = []
        for i, key_in in enumerate(key_order):

            quant = key_in.split("_")[-1]
            setattr(self, "ind_" + quant, i)

            self.inds_list.append(getattr(self, "ind_" + quant))
            self.keys.append(quant)
            if key_in.split("_")[0] in ["ln", "cos", "sin"]:
                self.conversions.append(
                    [getattr(self, key_in.split("_")[0] + "_convert"), np.array([i])]
                )

            if quant == "distance":
                self.conversions.append([self.convert_distance, np.array([i])])

            if quant in variables_to_recycle:
                wrap_value = recycle_guide[quant]
                if wrap_value == 2 * np.pi:
                    self.recycles.append([self.wrap_2pi, np.array([i])])

                elif wrap_value == np.pi:
                    self.recycles.append([self.wrap_pi, np.array([i])])

                else:
                    raise ValueError("wrap_value must be preprogrammed available")

    def ln_convert(self, x):
        return np.exp(x)

    def cos_convert(self, x):
        return np.arccos(x)

    def sin_convert(self, x):
        return np.arcsin(x)

    def convert_distance(self, x):
        return x * self.dist_conversion

    def convert(self, x):
        for func, inds in self.conversions:
            x[inds] = np.asarray(func(*x[inds]))
        return x

    def wrap_2pi(self, x):
        return x % (2 * np.pi)

    def wrap_pi(self, x):
        return x % (np.pi)

        """if x[self.ind_beta] < -np.pi/2 or x[self.ind_beta] > np.pi/2:
            # assumes beta = 0 at ecliptic plane [-pi/2, pi/2]
            x_trans = np.cos(x[self.ind_beta])*np.cos(x[self.ind_lam])
            y_trans = np.cos(x[self.ind_beta])*np.sin(x[self.ind_lam])
            z_trans = np.sin(x[self.ind_beta])

            x[self.ind_lam] = np.arctan2(y_trans, x_trans)
            x[self.ind_beta] = np.arcsin(z_trans/np.sqrt(x_trans**2 + y_trans**2 + z_trans**2))  # check this with eccliptic coordinates
        """

    def recycle(self, x):
        for func, inds in self.recycles:
            x[inds] = np.asarray(func(*x[inds]))
        return x
