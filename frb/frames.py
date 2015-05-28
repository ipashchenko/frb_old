import numpy as np
try:
    import george
    from george import kernels
except ImportError:
    george = None
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

vround = np.vectorize(round)
vint = np.vectorize(int)

class Frame(object):
    """
    Basic class that represents a set of regulary spaced frequency channels with
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
        self.nu = (nu_0 - nu * dnu)[::-1]
        self.t = t_0 + t * dt
        self.dt = dt

    def slice(self, channels, times):
        """
        Slice frame using specified channels and/or times.
        """
        raise NotImplementedError

    def de_disperse(self, dm):
        """
        De-disperse frame using specified value of DM.

        :param dm:
            Dispersion measure to use in de-dispersion [cm^3 / pc].

        """
        # MHz ** 2 * cm ** 3 * s / pc
        k = 1. / (2.410331 * 10 ** (-4))

        # Calculate shift of time caused by de-dispersion for all channels
        dt_all = k * dm * (1. / self.nu ** 2. - 1. / self.nu_0 ** 2.)
        # Find what number of time bins corresponds to this shifts
        nt_all = vint(vround(dt_all / self.dt))
        # Roll each axis (freq. channel) to each own number of time steps.
        for i in range(self.n_nu):
            self.values[i] = np.roll(self.values[i], -nt_all[i])

    def average_in_time(self, plot=False):
        """
        Average frame in time.

        :return:
            Numpy array with length equals number of frequency channels.
        """
        result = np.mean(self.values, axis=1)
        if plt is not None and plot:
            plt.plot(np.arange(self.n_nu), result, '.k')
            plt.xlabel("frequency channel #")
        return result

    def average_in_freq(self, plot=False):
        """
        Average frame in frequency.

        :return:
            Numpy array with length equals number of time steps.
        """
        result = np.mean(self.values, axis=0)
        if plt is not None and plot:
            plt.plot(np.arange(self.n_t), result, '.k')
            plt.xlabel("time steps")
        return result

    # TODO: if one choose what channels to plot - use ``extent`` kwarg.
    def plot(self, freqs=None, times=None, plot_indexes=True):
        if plt is not None:
            plt.imshow(self.values, interpolation='none', aspect='auto',
                       cmap=plt.cm.Reds)
            plt.colorbar()
            if not plot_indexes:
                raise NotImplementedError("Ticks haven't implemented yet")
                # plt.xticks(np.linspace(0, 999, 10, dtype=int),
                # frame.t[np.linspace(0, 999, 10, dtype=int)])
            plt.xlabel("time")
            plt.ylabel("frequency ch. #")


class DataFrame(Frame):
    """
    Class that represents the frame of real data.

    :param fname:
        Name of txt-file with rows representing frequency channels and columns -
        1d-time series of data for each frequency channel.

    """
    def __init__(self, fname, nu_0, t_0, dnu, dt):
        values = np.loadtxt(fname, unpack=True)
        n_nu, n_t = np.shape(values)
        super(DataFrame, self).__init__(n_nu, n_t, nu_0, t_0, dnu, dt)
        self.values += values


class SimFrame(Frame):
    """
    Class that represents the simulation of data.

    """
    def add_noise(self, amp, std, mean=0., kamp=None, kscale=None, kmean=None):
        """
        Add noise to frame using specified gaussian process or simple gaussian
        noise.

        :param amp:
            Amplitude of gaussian uncorrelated noise.
        :param std:
            Std of gaussian uncorrelated noise.
        :param mean: (optional)
            Mean of gaussian uncorrelated noise. (default: ``0.``)
        :param kamp: (optional)
            Amplitude of GP kernel. (default: ``1.0``)
        :param kscale:
            Scale of GP kernel. If ``None`` then don't add correlated noise.
            (default: ``None``)
        :param kmean: (optional)
            Mean of GP kernel. (default: ``0.0``)

        """
        noise = np.random.normal(amp, std, size=(self.n_t * self.n_nu)).reshape(np.shape(self.values))
        gp = george.GP(kamp * kernels.ExpSquaredKernel(kscale))

    def add_pulse(self, t_0, amp, width, dm=0.):
        """
        Add pulse to frame.

        :param t_0:
            Arrival time of pulse at highest frequency channel [s].
        :param amp:
            Amplitude of pulse.
        :param width:
            Width of gaussian pulse [s] (in time domain).
        :param dm: (optional)
            Dispersion measure of pulse [cm^3 / pc]. (Default: ``0.``)

        """
        # MHz ** 2 * cm ** 3 * s / pc
        k = 1. / (2.410331 * 10 ** (-4))

        # Calculate arrival times for all channels
        t0_all = (t_0 * np.ones(self.n_nu)[:, np.newaxis] +
                  k * dm * (1. / self.nu ** 2. -
                            1. / self.nu_0 ** 2.))[0]
        pulse = amp * np.exp(-0.5 * (self.t -
                                     t0_all[:, np.newaxis]) ** 2 / width ** 2.)
        self.values += pulse


if __name__ == '__main__':

    frame = SimFrame(128, 1000, 1600., 0., 16. / 128., 0.001)
    # Plot first (highest frequency) channel
    frame.add_pulse(0.5, 3., 0.003, dm=1000.)
    frame.plot()
    plt.plot(frame.t, frame.values[10])
    plt.plot(frame.t, frame.values[30])
    frame.de_disperse(dm=1000.)
    frame.plot()
