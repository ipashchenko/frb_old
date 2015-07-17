import ctypes
import multiprocessing
import numpy as np

vint = np.vectorize(int)
vround = np.vectorize(round)

def de_disperse(image, nu_0, d_nu, d_t, dm_range):
    """
    De-disperse dynamical spestra with grid of user specifies values.
    """
    n_nu, n_t = image.shape
    nu = np.arange(n_nu, dtype=float)
    nu = (nu_0 - nu * d_nu)[::-1]
    print "Pre-calculating cumulative sums..."
    cumsums = np.cumsum(image[::-1, :], axis=0)
    # print "Pre-calculating diffs of cumulative sums..."
    # dcumsums = cumsums[:, 1:] - cumsums[:, :-1]
    # zero_column = np.zeros(n_nu, dtype=float).reshape((n_nu, 1))
    # dcumsums = np.hstack((zero_column, dcumsums))

    ## MHz ** 2 * cm ** 3 * s / pc
    k = 1. / (2.410331 * 10 ** (-4))
    #dm = d_t / (k * (1. / (nu[0]) ** 2.) - 1. / nu[-1] ** 2.)
    #print "Minimum DM = ", dm

    #dm_array = np.array([i * dm for i in range(dt_max)])
    #print "DM grid = ", dm_array

    # Calculate shift of time caused by de-dispersion for all channels and all
    # values of DM
    # (#DM, #nu)
    dt_all = k * dm_range[:, np.newaxis] * (1. / nu ** 2. - 1. / nu_0 ** 2.)
    # Find what number of time bins corresponds to this shifts
    nt_all = vint(vround(dt_all / d_t))[:, ::-1]

    # Create array for TDM
    n_dm = len(dm_range)
    values = np.zeros((n_dm, n_t), dtype=float)
    values[0] = cumsums[-1]

    # shared_array_base = multiprocessing.Array(ctypes.c_float, n_nu * n_t)
    # _values = np.ctypeslib.as_array(shared_array_base.get_obj()).reshape((n_dm,
    #                                                                       n_t,))
    # Cycle over DM values and fill TDM array
    for i, nt in enumerate(nt_all[1:]):
        # print i
        # print nt
        # Find at witch frequency channels time shifts have occured
        indx = np.where(nt[1:] - nt[:-1] == 1)[0]
        # print indx
        # print "================"
        # if not channels.size:
        #     values[i + 1] = values[i]
        #     continue
        # shifts = np.arange(len(channels))
        # Sum differences of shifted cumsums
        # result = np.sum(np.array([np.roll(cumsums[ch + 1] - cumsums[ch],
        #                                          nt[ch]) for ch in
        #                                  channels[1:-1]]), axis=0)
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
    fname = '/home/ilya/code/frb/data/630_sec_wb_raes08a_128ch.npy'
    frame = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    frame.add_pulse(300., 2.5, 0.003, dm=800.)
    dm_grid = frame.create_dm_grid(0., 1000., 35.)
    t0 = time.time()
    result = de_disperse(frame.values, 1684, 16./128, 0.001, dm_grid)
    t1 = time.time()
    print t1 - t0
