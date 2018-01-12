import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pp

import AAKwrapper

import time

print("Loading AAK from {}.".format(AAKwrapper.__path__))

pars = {'length': 1000000,
        'dt': 5.,
        'p': 5.,
        'T': 1.,
        'f': 2.e-3,
        'T_fit': 1.,
        'mu': 1.e1,
        'M': 1.e6,
        's': 0.5,
        'e': 0.1,
        'iota': 0.524,
        'gamma': 0.,
        'psi': 0.,
        'theta_S': 0.785,
        'phi_S': 0.785,
        'theta_K': 1.05,
        'phi_K': 1.05,
        'alpha': 0.,
        'D': 1.}

print("Computing...")

start = time.time()

t, hI, hII = AAKwrapper.AAK(pars)

end = time.time()
print("Time taken: {}".format(end - start))

print("Plotting...")

pp.ioff()
fig = pp.figure()
pp.plot(t[:10000], hI[:10000])
pp.plot(t[:10000], hII[:10000])
pp.savefig("AAKdemo.pdf")
pp.close(fig)

print("All done. Please see {}.".format('AAKdemo.pdf'))