import math
import numpy as np
from utils import delta_dm_max
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
        self.dnu = dnu

    def slice(self, channels, times):
        """
        Slice frame using specified channels and/or times.
        """
        raise NotImplementedError

    def de_disperse(self, dm, in_place=True):
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
        # TODO: vectorize - check ``np.roll`` docs
        values = list()
        for i in range(self.n_nu):
            values.append(np.roll(self.values[i], -nt_all[i]))
        values = np.vstack(values)
        if in_place:
            self.values = values[:,:]
        return values

    def average_in_time(self, values=None, plot=False):
        """
        Average frame in time.

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
    def plot(self, freqs=None, times=None, plot_indexes=True):
        if plt is not None:
            plt.imshow(self.values, interpolation='none', aspect='auto',
                       cmap=plt.cm.Reds)
            plt.colorbar()
            if not plot_indexes:
                raise NotImplementedError("Ticks haven't implemented yet")
                # plt.xticks(np.linspace(0, 999, 10, dtype=int),
                # frame.t[np.linspace(0, 999, 10, dtype=int)])
            plt.xlabel("time steps")
            plt.ylabel("frequency ch. #")
            plt.show()


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
        noise = np.random.normal(amp, std,
                                 size=(self.n_t *
                                       self.n_nu)).reshape(np.shape(self.values))
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

    # frame = SimFrame(128, 1000, 1600., 0., 16. / 128., 0.001)
    # # Plot first (highest frequency) channel
    # frame.add_pulse(0.5, 3., 0.003, dm=1000.)
    # frame.plot()
    # plt.plot(frame.t, frame.values[10])
    # plt.plot(frame.t, frame.values[30])
    # frame.de_disperse(dm=1000.)
    # frame.plot()
    # B1641-45 test drive:)
    frame = DataFrame('data.txt', 1676., 0., 16. / 32, 0.001)
    frame.plot()
    dm_min = 0.
    dm_max = 1000.
    dm_delta = delta_dm_max(frame.nu_0, frame.nu_0 - frame.n_nu * frame.dnu,
                            frame.dt)

    n_dm = (dm_max - dm_min) / dm_delta
    n_dm = int(math.ceil(n_dm))
    max_snrs = list()
    dm_used = list()
    for i in xrange(n_dm):
        dm = dm_min + (i + 1) * dm_delta
        print "Searching DM = ", dm, " cm^3 / pc"
        _frame = frame.de_disperse(dm=dm, in_place=False)
        _frame_t = frame.average_in_freq(_frame, plot=False)
        max_snr = max(abs((_frame_t - np.mean(_frame_t)) / np.std(_frame_t)))
        print "Max SNR = ", max_snr
        max_snrs.append(max_snr)
        dm_used.append(dm)

    dm_result = np.mean(np.array(dm_used)[np.where(np.array(max_snrs) ==
                                                   max(max_snrs))])
    print "DM = ",  dm_result

    # Plot and print results
    plt.close()
    plt.figure()
    plt.subplot(211)
    plt.plot(dm_used, max_snrs, '.k')
    plt.xlabel("DM-correction")
    plt.ylabel("max SNR")
    plt.axvline(dm_result, color='r', lw=3, label=str(round(dm_result, 1)) +
                " cm^3 / pc")
    plt.legend(loc='best')
    plt.subplot(212)
    values = frame.de_disperse(dm=dm_result)
    plt.imshow(values, interpolation='none', aspect='auto', cmap=plt.cm.Reds)
    plt.colorbar()
    plt.xlabel("time steps")
    plt.ylabel("frequency ch. #")
    plt.show()
