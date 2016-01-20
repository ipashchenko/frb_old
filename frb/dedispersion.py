import numpy as np
from objects import BatchedTDMIO

vint = np.vectorize(int)
vround = np.vectorize(round)
k = 1. / (2.410331 * 10 ** (-4))


def de_disperse(image, nu_0, d_nu, d_t, dm_values):
    """
    De-disperse dynamical spestra with grid of user specifies values of DM.

    :param image:
        2D numpy array of dynamical spectra (#freq, #t).
    :param nu_0:
        Frequency of highest frequency channel [MHz].
    :param d_nu:
        Width of spectral channel [MHz].
    :param d_t:
        Time step [s].
    :param dm_values:
        Array-like of DM values to dedisperse [cm^3 /pc].

    :return:
        2D numpy array (#DM, #t)

    """
    dm_values = np.array(dm_values)
    n_nu, n_t = image.shape
    nu = np.arange(n_nu, dtype=float)
    nu = (nu_0 - nu * d_nu)[::-1]
    print "Pre-calculating cumulative sums..."
    cumsums = np.cumsum(image[::-1, :], axis=0)
    print "Pre-calculating their diffs..."
    dcumsums = np.roll(cumsums, 1, axis=1) - cumsums
    # dcumsums = np.hstack((np.zeros(n_nu, dtype=float).reshape(n_nu, 1),
    #                       dcumsums))


    ## MHz ** 2 * cm ** 3 * s / pc
    k = 1. / (2.410331 * 10 ** (-4))
    # Calculate shift of time caused by de-dispersion for all channels and all
    # values of DM
    # (#DM, #nu)
    dt_all = k * dm_values[:, np.newaxis] * (1. / nu ** 2. - 1. / nu_0 ** 2.)
    # Find what number of time bins corresponds to this shifts
    nt_all = vint(vround(dt_all / d_t))[:, ::-1]

    # Create array for TDM
    values = np.zeros((len(dm_values), n_t), dtype=float)
    # Fill DM=0 row
    values[0] = cumsums[-1]

    # Cycle over DM values and fill TDM array for others DM values
    for i, nt in enumerate(nt_all[1:]):
        # print "Current shifts are : ", nt
        # Find at which frequency channels time shifts have occurred
        indx = np.array(np.where(nt[1:] - nt[:-1] == 1)[0].tolist() +
                        [n_nu - 1])
        # indx = np.where(nt[1:] - nt[:-1] == 1)[0]
        # print "Current indexes of change time shifts : ", indx
        # print "First add shifted on ", -nt[-1], "cumsum[-1]"
        result = np.roll(cumsums[-1], -nt[-1])
        for ix, j in enumerate(indx[:-1]):
            # print "Shift number ", ix, " on channel ", j, " with value ", -nt[j]
            result += np.roll(dcumsums[j], -nt[j])
        values[i + 1] = result

    return values


def split_de_disperse(image, nu_min, nu_max, d_t, dm_values, n_split):

    dm_values = np.array(dm_values)
    # NUmber of spectral channels/ time bins
    n_nu, n_t = image.shape
    # Width of single spectral channel
    d_nu = (nu_max - nu_min) / n_nu
    # Number of spectral channels in one split
    n_nu_split = n_nu / n_split
    # Width of one split
    d_nu_split = d_nu * n_nu_split

    # Calculate t-DM plane image for first split of highest frequencies first
    print "De-disperse first split"
    tdm = de_disperse(image[: n_nu_split, :], nu_max, d_nu, d_t, dm_values)

    # Calculate t-DM images for other splits and shift them and add to original
    for i in range(n_split)[1:]:
        print "De-disperse {} split of {}".format(i + 1, n_split)
        nu_max_current = nu_max - i * d_nu_split
        tdm_ = de_disperse(image[i * n_nu_split: (i + 1) * n_nu_split, :],
                           nu_max_current, d_nu, d_t, dm_values)
        # Shift each DM value along t direction on value
        # 4.15 DM * (1/v_max_current ** 2 - 1/v_max**2)
        dt_all = k * dm_values * (1. / nu_max_current ** 2. - 1. / nu_max ** 2.)
        # Find what number of time bins corresponds to this shifts
        nt_all = vint(vround(dt_all / d_t))[::-1]
        for j, nt in enumerate(nt_all):
            tdm[j] += np.roll(tdm_[j], -nt)

    return tdm




if __name__ == '__main__':
    import time
    from frames import DataFrame
    from frames import Frame
    # from objects import BatchedTDMIO
    # fname = '/home/ilya/code/frb/data/630_sec_wb_raes08a_128ch.npy'
    # fname = '/home/ilya/code/frb/data/crab_600sec_64ch_1ms.npy'
    # frame = DataFrame(fname, 1684., 0., 16. / 64., 0.001)
    # frame = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    # frame.add_pulse(100., 0.7, 0.003, 300.)
    # frame.add_pulse(200., 0.7, 0.001, 400.)
    # frame.add_pulse(300., 0.7, 0.006, 200.)
    # frame.add_pulse(400., 2.7, 0.002, 500.)
    # frame.add_pulse(500., 0.7, 0.003, 600.)
    # frame.add_pulse(600., 0.7, 0.004, 300.)
    # frame.add_pulse(300., 1.0, 0.003, dm=500.)
    # dm_values = frame.create_dm_grid(0., 300., 10.)
    # t0 = time.time()
    print "Creating frame"
    frame = Frame(512, 600000, 1684., 0., 16./512, 1./1000)
    print "Adding pulse"
    # frame.add_pulse(100., 0.09, 0.003, 100.)
    # frame.add_pulse(200., 0.09, 0.003, 200.)
    frame.add_pulse(300., 0.09, 0.003, 300.)
    # frame.add_pulse(400., 0.09, 0.003, 500.)
    # frame.add_pulse(500., 0.09, 0.003, 700.)
    print "Adding noise"
    frame.add_noise(0.5)
    dm_values = frame.create_dm_grid(0., 1000., 10)
    # t0 = time.time()
    # result = de_disperse(frame.values, 1684, 16./512, 0.001, dm_values)
    # t1 = time.time()
    # print "Dedispersion took", t1 - t0
    t0 = time.time()
    result = split_de_disperse(frame.values, 1668., 1684., 0.001, dm_values, 4)
    t1 = time.time()
    print "Dedispersion took", t1 - t0
    # t0 = time.time()
    # result = frame.grid_dedisperse(dm_values)
    # t1 = time.time()
    # print "Grid-dedispersion with 1 thread took", t1 - t0
    # t0 = time.time()
    # result = frame.grid_dedisperse(dm_values, threads=4)
    # t1 = time.time()
    # print "Grid-dedispersion with 4 threads took", t1 - t0
    # btdmi = BatchedTDMIO(result, frame.t, dm_values, 99.85, d_dm=150, dt=0.003)
    # t0 = time.time()
    # xy = btdmi.run(batch_size=100000)
    # t1 = time.time()
    # print xy
    # print "Search of pulses took", t1 - t0
