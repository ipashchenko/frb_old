import os
import sys
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
import numpy as np
from frb.frames import SimFrame, DataFrame
import george
from george.kernels import Matern32Kernel
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


if __name__ == '__name__':
    # Generating fake FRB
    frame = SimFrame(128, 1000, 1676., 0., 16. / 128., 0.001)
    frame.add_pulse(0.3, 0.5, 0.0015, dm=2400.)
    frame.add_noise(0.5)
    frame.plot()

    # Fitting noise with GP
    fname = '/home/ilya/code/frb/data/data.txt'
    frame = DataFrame(fname, 1684, 0, 0.5, 0.001)
    # Set up GP
    kernel = Matern32Kernel([10., 0.05], ndim=2)
    gp = george.GP(kernel)
    # Coordinates
    x = list()
    for nu in frame.nu:
        for t in frame.t:
            x.append([nu, t])
    x = np.array(x)

    gp.compute(x[::100], 0.3)
    mu, cov = gp.predict(frame.values.flatten()[::100], x[::50])