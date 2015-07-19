import numpy as np

vint = np.vectorize(int)
vround = np.vectorize(round)


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
        # Find at wich frequency channels time shifts have occurred
        indx = np.where(nt[1:] - nt[:-1] == 1)[0]
        result = cumsums[indx[0]] + np.roll(cumsums[-1] - cumsums[indx[-1]],
                                            -nt[-1])
        diff = [np.roll(cumsums[indx[j + 1]] - cumsums[indx[j]],
                        -nt[indx[j + 1]]) for j in range(len(indx) - 1)]
        result += np.sum(diff, axis=0)
        values[i + 1] = result

    return values


if __name__ == '__main__':
    import time
    from frames import DataFrame
    # from objects import BatchedTDMIO
    # fname = '/home/ilya/code/frb/data/630_sec_wb_raes08a_128ch.npy'
    fname = '/home/ilya/code/frb/data/out_crab_full_64x1.npy'
    frame = DataFrame(fname, 1684., 0., 16. / 64., 0.001)
    # frame = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    # frame.add_pulse(100., 0.7, 0.003, 300.)
    # frame.add_pulse(200., 0.7, 0.001, 400.)
    # frame.add_pulse(300., 0.7, 0.006, 200.)
    # frame.add_pulse(400., 0.7, 0.002, 500.)
    # frame.add_pulse(500., 0.7, 0.003, 600.)
    # frame.add_pulse(600., 0.7, 0.004, 300.)
    # frame.add_pulse(300., 1.0, 0.003, dm=500.)
    dm_values = frame.create_dm_grid(0., 1000., 30.)
    t0 = time.time()
    result = de_disperse(frame.values, 1684, 16./64, 0.001, dm_values)
    # result = de_disperse(frame.values, 1684, 16./128, 0.001, dm_values)
    t1 = time.time()
    print t1 - t0
    # t0 = time.time()
    # result = frame.grid_dedisperse(dm_values, threads=4)
    # t1 = time.time()
    # print t1 - t0
    # btdmi = BatchedTDMIO(result, frame.t, dm_values, 99.85, d_dm=350, dt=0.003)
    # t0 = time.time()
    # xy = btdmi.run(batch_size=100000)
    # t1 = time.time()
    # print xy
    # print t1 - t0
