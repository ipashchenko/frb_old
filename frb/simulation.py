import george
import numpy as np
import matplotlib.pyplot as plt
from george import kernels


class Frame(object):
    """
    Class that represents a set of regulary spaced frequency channels with
    regulary measured values.
    """
    def __init__(self, n_nu, n_t, nu_0, t_0, dnu, dt):
        self.n_nu = n_nu
        self.n_t = n_t
        self.nu_0 = nu_0
        self.t_0 = t_0
        self.values = np.zeros(n_nu * n_t, dtype=float).reshape((n_nu, n_t,))
        nu = np.arange(n_nu)
        t = np.arange(n_t)
        self.nu = nu_0 - nu * dnu
        self.t = t_0 + t * dt

    def slice(self, channels, times):
        """
        Slice frame using specified channels and/or times.
        """
        pass

    def add_pulse(self, t_0, amp, width, dm=0.):
        """
        Add pulse to data.

        :param t_0:
            Arrival time of pulse at highest frequency channel.
        :param amp:
            Amplitude of pulse.
        :param width:
            Width of pulse (in time domain).
        :param dm: (optional)
            Dispersion measure of pulse. (Default: ``0.``)

        """
        # Calculate arrival times for all channels
        t0_all = (t_0 * np.ones(self.n_nu)[:, np.newaxis] +
                  dm * (1. / self.nu ** 2. -
                        1. / self.nu_0 ** 2.))[0]
        pulse = amp * np.exp(-0.5 * (self.t -
                                     t0_all[:, np.newaxis]) ** 2 / width)
        self.values += pulse

    def add_noise(self, std, mean=0., scale=None, amp=None):
        """
        Add noise to frame using specified gaussian process or simple gaussian
        noise.
        :param std:
        :param scale:
        :param amp:
        :return:
        """
        pass

    def de_disperse(self, dm):
        pass

    def average_in_time(self, times=None):
        pass

    def average_in_freq(self, freqs=None):
        pass

    def plot(self, freqs=None, times=None):
        pass

if __name__== '__main__':

    frame = Frame(512, 300, 1.6 * 10 ** 3, 0., 15625.0 * 10 ** (-6), 0.01)
    # Plot first (highest frequency) channel
    plt.plot(frame.t, frame.values[0])
    frame.add_pulse(1.0, 3., 0.1)
    plt.plot(frame.t, frame.values[100])
    plt.plot(frame.t, frame.values[300])




