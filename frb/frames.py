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


# TODO: add masking edges
class Frame(object):
    """
    Basic class that represents a set of regulary spaced frequency channels with
    regulary measured values.

    :param n_nu:
        Number of spectral channels.
    :param n_t:
        Number of time steps.
    :param nu_0:
        Frequency of highest frequency channel [MHz].
    :param t0:
        Time of first measurement.
    :param dnu:
        Width of spectral channel [MHz].
    :param dt:
        Time step [s].

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
        self.dnu = dnu

    def slice(self, channels, times):
        """
        Slice frame using specified channels and/or times.
        """
        raise NotImplementedError

    def de_disperse(self, dm, replace=True):
        """
        De-disperse frame using specified value of DM.

        :param dm:
            Dispersion measure to use in de-dispersion [cm^3 / pc].
        :param replace: (optional)
            Replace instance's frame values with de-dispersed ones? (default:
            ``True``)

        """
        # MHz ** 2 * cm ** 3 * s / pc
        k = 1. / (2.410331 * 10 ** (-4))

        # Calculate shift of time caused by de-dispersion for all channels
        dt_all = k * dm * (1. / self.nu ** 2. - 1. / self.nu_0 ** 2.)
        # Find what number of time bins corresponds to this shifts
        nt_all = vint(vround(dt_all / self.dt))
        # Roll each axis (freq. channel) to each own number of time steps.
        # TODO: vectorize - check ``np.roll`` docs
        values = list()
        for i in range(self.n_nu):
            values.append(np.roll(self.values[i], -nt_all[i]))
        values = np.vstack(values)
        if replace:
            self.values = values[:, :]
        return values

    def average_in_time(self, values=None, plot=False):
        """
        Average frame in time.

        :param values: ``(n_t, n_nu)`` (optional)
            Numpy array of Frame values to average. If ``None`` then use current
            instance's values. (default: ``None``)
        :param plot: (optional)
            Plot figure? If ``False`` then only return array. (default:
            ``False``)

        :return:
            Numpy array with length equals number of frequency channels.
        """
        if values is None:
            values = self.values
        result = np.mean(values, axis=1)
        if plt is not None and plot:
            plt.plot(np.arange(self.n_nu), result, '.k')
            plt.xlabel("frequency channel #")
        return result

    def average_in_freq(self, values=None, plot=False):
        """
        Average frame in frequency.

        :param values: ``(n_t, n_nu)`` (optional)
            Numpy array of Frame values to average. If ``None`` then use current
            instance's values. (default: ``None``)
        :param plot: (optional)
            Plot figure? If ``False`` then only return array. (default:
            ``False``)

        :return:
            Numpy array with length equals number of time steps.
        """
        if values is None:
            values = self.values
        result = np.mean(values, axis=0)
        if plt is not None and plot:
            plt.plot(np.arange(self.n_t), result, '.k')
            plt.xlabel("time steps")
        return result

    # TODO: if one choose what channels to plot - use ``extent`` kwarg.
    def plot(self, plot_indexes=True):
        if plt is not None:
            plt.imshow(self.values, interpolation='none', aspect='auto')
            plt.colorbar()
            if not plot_indexes:
                raise NotImplementedError("Ticks haven't implemented yet")
                # plt.xticks(np.linspace(0, 999, 10, dtype=int),
                # frame.t[np.linspace(0, 999, 10, dtype=int)])
            plt.xlabel("time steps")
            plt.ylabel("frequency ch. #")
            plt.show()

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

    def save_to_txt(self, fname):
        np.savetxt(fname, self.values.T)

    def add_noise(self, std, kamp=None, kscale=None, kmean=0.0):
        """
        Add noise to frame using specified gaussian process or simple
        rayleigh-distributed noise.

        Correlated noise is correlated along the frequency axis.

        :param std:
            Std of rayleigh-distributed uncorrelated noise.
        :param kamp: (optional)
            Amplitude of GP kernel. If ``None`` then don't add correlated noise.
            (default: ``None``)
        :param kscale:
            Scale of GP kernel [MHz]. If ``None`` then don't add correlated
            noise.  (default: ``None``)
        :param kmean: (optional)
            Mean of GP kernel. (default: ``0.0``)

        """
        noise = np.random.rayleigh(std,
                                   size=(self.n_t *
                                         self.n_nu)).reshape(np.shape(self.values))
        self.values += noise
        if kscale is not None and kamp is not None:
            gp1 = george.GP(kamp * kernels.ExpSquaredKernel(kscale))
            gp2 = george.GP(kamp * kernels.ExpSquaredKernel(kscale))
            for i in xrange(self.n_t):
                gp_samples = np.sqrt(gp1.sample(self.nu) ** 2. +
                                     gp2.sample(self.nu) ** 2.)
                self.values[:, i] += gp_samples


class DataFrame(Frame):
    """
    Class that represents the frame of real data.

    :param fname:
        Name of txt-file with rows representing frequency channels and columns -
        1d-time series of data for each frequency channel.

    """
    def __init__(self, fname, nu_0, t_0, dnu, dt, n_nu_discard=0):
        values = np.loadtxt(fname, unpack=True)
        n_nu, n_t = np.shape(values)
        super(DataFrame, self).__init__(n_nu - n_nu_discard, n_t,
                                        nu_0 - n_nu_discard * dnu / 2., t_0,
                                        dnu, dt)
        if n_nu_discard:
            self.values += values[n_nu_discard / 2 : -n_nu_discard / 2, :]
        else:
            self.values += values


class SimFrame(Frame):
    """
    Class that represents the simulation of data.

    """
    pass
