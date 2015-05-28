def delta_dm_max(nu_max, nu_min, dt):
    """
    Return difference in DM that corresponds to arrival time shift between
    highest and lowest frequency channels equals to time resolution.

    :param nu_max:
        Frequency of highest frequency channel [MHz].
    :param nu_min:
        Frequency of lowest frequency channel [MHz].
    :param dt:
        Time interval between measurements [s].
    :return:
        Maximum difference of DM that we should use for FRB grid search.

    >>> delta_dm_max(1600., 1600. - 16., 3. / 1000)
    1.8147253416989504

    """
    # MHz ** 2 * cm ** 3 * s / pc
    k = 1. / (2.410331 * 10 ** (-4))
    return dt / (k * (1. / (nu_min) ** 2.) - 1. / nu_max ** 2.)


def time_interval(size_gb, n_nu, nu_max, nu_min, dt, dm_max):
    """
    Return time interval [s] that could be spanned by 3D-numpy float array with
    ``n_nu`` frequency channels with given max. and min. frequencies and max.
    DM used in grid search.

    :param size_gb:
        Size of array [GB]
    :param n_nu:
        Number of freq. channels.
    :param nu_max:
        Frequency of highest channel [MHz].
    :param nu_min:
        Frequency of lowest channel [MHz].
    :param dt:
        Time resolution [s].
    :param dm_max:
        Maximum DM to search using grid of DM values.
    :return:
        Time interval that could be spanned.

    >>> time_interval(1., 128, 1600., 1600. - 16., 3./1000, 1000)
    10.63407441016334

    """
    # Maximum delta for grid search
    d_dm = delta_dm_max(nu_max, nu_min, dt)
    # Length of DM grid.
    n_dm = int(dm_max / d_dm)
    # 4 bytes in ``float``
    nt = size_gb * 10 ** 9 / (4 * n_nu * n_dm)
    return dt * nt