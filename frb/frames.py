import gc
import numpy as np
import pyfits as pf
import pickle_method
from multiprocessing import Pool
from collections import OrderedDict
from utils import vint, vround, delta_dm_max

try:
    import george
    from george import kernels
except ImportError:
    george = None
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


class Frame(object):
    """
    Basic class that represents a set of regulary spaced frequency channels with
    regulary measured values (time sequence of autospectra).

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

    def add_values(self, array):
        """
        Add dyn. spectra in form of numpy array (#ch, #t,) to instance.

        :param array:
            Array-like of dynamical spectra (#ch, #t,).
        """
        array = np.atleast_2d(array)
        assert self.values.shape == array.shape
        self.values += array

    def slice(self, channels, times):
        """
        Slice frame using specified channels and/or times.
        """
        raise NotImplementedError

    # FIXME: at small ``dt`` it uses too small DM-step for my laptop RAM:)
    def de_disperse(self, dm, replace=False):
        """
        De-disperse frame using specified value of DM.

        :param dm:
            Dispersion measure to use in de-dispersion [cm^3 / pc].
        :param replace: (optional)
            Replace instance's frame values with de-dispersed ones? (default:
            ``False``)

        """
        # MHz ** 2 * cm ** 3 * s / pc
        k = 1. / (2.410331 * 10 ** (-4))

        # Calculate shift of time caused by de-dispersion for all channels
        dt_all = k * dm * (1. / self.nu ** 2. - 1. / self.nu_0 ** 2.)
        # Find what number of time bins corresponds to this shifts
        nt_all = vint(vround(dt_all / self.dt))
        # Roll each axis (freq. channel) to each own number of time steps.
        values = list()
        for i in range(self.n_nu):
            values.append(np.roll(self.values[i], -nt_all[i]))
        values = np.vstack(values)
        #values = roll2d(self.values, -nt_all, axis=1)

        if replace:
            self.values = values[:, :]
        return values

    # FIXME: at small ``dt`` it uses too small DM-step for my laptop RAM:)
    def _de_disperse_freq_average(self, dm):
        """
        De-disperse frame using specified value of DM and average in frequency.

        :param dm:
            Dispersion measure to use in de-dispersion [cm^3 / pc].

        :notes:
            This method avoids creating ``(n_nu, n_t)`` arrays and must be
            faster for data with big sizes. But it returns already frequency
            averaged de-dispersed dyn. spectra.

        """
        # MHz ** 2 * cm ** 3 * s / pc
        k = 1. / (2.410331 * 10 ** (-4))

        # Calculate shift of time caused by de-dispersion for all channels
        dt_all = k * dm * (1. / self.nu ** 2. - 1. / self.nu_0 ** 2.)
        # Find what number of time bins corresponds to this shifts
        nt_all = vint(vround(dt_all / self.dt))
        # Container for summing de-dispersed frequency channels
        values = np.zeros(self.n_t)
        # Roll each axis (freq. channel) to each own number of time steps.
        for i in range(self.n_nu):
            values += np.roll(self.values[i], -nt_all[i])

        return values / self.n_nu

    def average_in_time(self, values=None, plot=False):
        """
        Average frame in time.

        :param values: ``(n_nu, n_t)`` (optional)
            Numpy array of Frame values to average. If ``None`` then use current
            instance's values. (default: ``None``)
        :param plot: (optional)
            Plot figure? If ``False`` then only return array. (default:
            ``False``)

        :return:
            Numpy array with length equals the number of frequency channels.
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
    def plot(self, plot_indexes=True, savefig=None):
        if plt is not None:
            plt.figure()
            plt.matshow(self.values, aspect='auto')
            plt.colorbar()
            if not plot_indexes:
                raise NotImplementedError("Ticks haven't implemented yet")
                # plt.xticks(np.linspace(0, 999, 10, dtype=int),
                # frame.t[np.linspace(0, 999, 10, dtype=int)])
            plt.xlabel("time steps")
            plt.ylabel("frequency ch. #")
            plt.title('Dynamical spectra')
            if savefig is not None:
                plt.savefig(savefig, bbox_inches='tight')
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
            if not george:
                raise Exception("Install george for correlated noise option.")
            gp1 = george.GP(kamp * kernels.ExpSquaredKernel(kscale))
            gp2 = george.GP(kamp * kernels.ExpSquaredKernel(kscale))
            for i in xrange(self.n_t):
                gp_samples = np.sqrt(gp1.sample(self.nu) ** 2. +
                                     gp2.sample(self.nu) ** 2.)
                self.values[:, i] += gp_samples

    def _step_dedisperse(self, dm):

        """
        Method that de-disperses frame using specified value of DM and frequency
        averages the result.

        :param dm:
        :return:
        """
        values = self.de_disperse(dm)
        return self.average_in_freq(values)

    def create_dm_grid(self, dm_min, dm_max, dm_delta=None):
        """
        Method that create DM-grid for current frame.
        :param dm_min:
        :param dm_max:
        :param dm_delta:
        :return:
        """
        if dm_delta is None:
            # Find step for DM grid
            # Seems that ``5`` is good choice (1/200 of DM range)
            dm_delta = 5 * delta_dm_max(self.nu_0,
                                        self.nu_0 - self.n_nu * self.dnu,
                                        self.dt)

        # Create grid of searched DM-values
        return np.arange(dm_min, dm_max, dm_delta)

    def grid_dedisperse(self, dm_grid, savefig=None, threads=1):
        """
        Method that de-disperse ``Frame`` instance with range values of
        dispersion measures and average them in frequency to obtain image in
        (t, DM)-plane.

        :param dm_grid:
            Array-like of value of DM on which to de-disperse [cm^3/pc].
        :param savefig: (optional)
            File to save picture.
        :param threads: (optional)
            Number of threads used for parallelization with ``mupltiprocessing``
            module. If > 1 then it isn't used. (default: 1)

        """
        pool = None
        if threads > 1:
            pool = Pool(threads, maxtasksperchild=100)

        if pool:
            m = pool.map
        else:
            m = map

        # Accumulator of de-dispersed frequency averaged frames
        frames = list(m(self._de_disperse_freq_average, dm_grid.tolist()))
        frames = np.array(frames)

        if pool:
            # Close pool
            pool.close()
            pool.join()

        # Plot results
        if savefig is not None:
            plt.imshow(frames, interpolation='none', aspect='auto')
            plt.xlabel('De-dispersed by DM freq.averaged dyn.spectr')
            plt.ylabel('DM correction')
            plt.yticks(np.linspace(0, len(dm_grid) - 10, 5, dtype=int),
                       vint(dm_grid[np.linspace(0, len(dm_grid) - 10, 5,
                                                dtype=int)]))
            plt.colorbar()
            plt.savefig(savefig, bbox_inches='tight')
            plt.show()
            plt.close()

        return frames


class CompositeFrame(object):
    """
    Class that represents Frame partitioned in some number of frequency bands.

    :param n_bands:
        Number of subbands.
    :param n_nu_band:
        Number of frequency channels in each subband.
    :param n_t:
        Number of time steps.
    :param nu_0_bands:
        Frequency of the highest frequency channel or iterable of frequencies of
        the highest frequency channels in each subband. If former then for each
        of ``n_nu_band`` subband the highest frequency is determined from simple
        partition of full bandwidth (that is ``dnu * n_nu_band``) in
        ``n_nu_band`` parts. [MHz].
    :param t0:
        Time of first measurement.
    :param dnu:
        Width of single spectral channel [MHz].
    :param dt:
        Time step [s].

    """
    def __init__(self, n_bands, n_nu_band, n_t, nu_0_bands, t_0, dnu, dt):
        try:
            self.nu_0_bands = sorted(nu_0_bands)
        except TypeError:
            self.nu_0_bands = [nu_0_bands - i * n_nu_band for i in
                               range(n_bands)]
        self.n_nu_band = n_nu_band
        self.n_bands = n_bands
        self.n_t = n_t
        self.t_0 = t_0
        self.dnu = dnu
        self.dt = dt
        self.frames = OrderedDict()
        for nu_0 in self.nu_0_bands:
            self.frames.update({nu_0: Frame(n_nu_band, n_t, nu_0, t_0, dnu,
                                            dt)})

    def create_dm_grid(self, dm_min, dm_max, dm_delta=None):
        """
        Method that create DM-grid for current frame.
        :param dm_min:
        :param dm_max:
        :param dm_delta:
        :return:
        """
        if dm_delta is None:
            # Find step for DM grid
            # Seems that ``5`` is good choice (1/200 of DM range)
            nu_max = max(self.nu_0_bands)
            nu_min = min(self.nu_0_bands) - self.n_nu_band * self.dnu
            dm_delta = 5 * delta_dm_max(nu_max, nu_min, self.dt)

        # Create grid of searched DM-values
        return np.arange(dm_min, dm_max, dm_delta)

    def grid_dedisperse(self, dm_min, dm_max, dm_delta=None, threads=1):
        # Find grid of DM-value for combined frame
        dm_grid = self.create_dm_grid(dm_min, dm_max, dm_delta=dm_delta)

        # Prepare container of (t, DM)-plane values
        frame_t_dedm = np.zeros((len(dm_grid), self.n_t,), dtype=float)
        # Create grid of searched DM-values
        for nu_0, frame in self.frames.items():
            print "De-dispersing frame with nu_0 = ", nu_0
            frame_t_dedm += frame.grid_dedisperse(dm_grid, savefig=None,
                                                  threads=threads)
            gc.collect()
        return frame_t_dedm / self.n_bands


# TODO: should i use just one class ``Frame`` but different io-methods?
# TODO: Create subclass for FITS input.
class DataFrame(Frame):
    """
    Class that represents the frame of real data.

    :param fname:
        Name of txt-file with rows representing frequency channels and columns -
        1d-time series of data for each frequency channel.

    """
    def __init__(self, fname, nu_0, t_0, dnu, dt, n_nu_discard=0):
        # Assert even number of channels to discard
        assert not int(n_nu_discard) % 2

        values = np.loadtxt(fname, unpack=True)
        n_nu, n_t = np.shape(values)
        super(DataFrame, self).__init__(n_nu - n_nu_discard, n_t,
                                        nu_0 - n_nu_discard * dnu / 2., t_0,
                                        dnu, dt)
        if n_nu_discard:
            self.add_values(values[n_nu_discard / 2 : -n_nu_discard / 2, :])
            # self.values += values[n_nu_discard / 2 : -n_nu_discard / 2, :]
        else:
            self.add_values(values)
            # self.values += values


class CompositeDataFrame(CompositeFrame):
    def __init__(self, fname, nu_0, t_0, dnu, dt, n_bands=1):
        values = np.loadtxt(fname, unpack=True)
        n_nu, n_t = np.shape(values)
        assert n_nu % n_bands == 0
        n_nu_bands = n_nu / n_bands
        super(CompositeDataFrame, self).__init__(n_bands, n_nu_bands, n_t,
                                                 nu_0, t_0, dnu, dt)
        for i, nu_0 in enumerate(self.nu_0_bands):
            self.frames[nu_0].add_values(values[i * n_nu_bands:
                                                (i + 1) * n_nu_bands, :])


if __name__ == '__main__':
    import time
    fname = '/home/ilya/code/frb/data/out_crab_full_64x1'
    # cdframe = CompositeDataFrame(fname, 1684., 0., 16. / 64., 0.001, n_bands=4)
    # t0 = time.time()
    # frames_t_dedm1 = cdframe.grid_dedisperse(0, 1000., threads=4)
    # t1 = time.time()
    # print " Composite de-disperse took ", t1 - t0

    frame = DataFrame(fname, 1684., 0., 16. / 64., 0.001)
    dm_grid = frame.create_dm_grid(0, 1000.)
    t0 = time.time()
    frames_t_dedm2 = frame.grid_dedisperse(dm_grid, threads=4)
    t1 = time.time()
    print " Ordinary de-disperse took ", t1 - t0
