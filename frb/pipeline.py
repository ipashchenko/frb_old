import glob
import pyfits as pf
from frames import DataFrame
from search_utils import grid_dedisperse_frame
from objects import Objects


# TODO: To compare pulse candidates from different antennas we need time of
# pulse at some fixed frequency. Or use the same set ups for each antenna (that
# is not flexible enough) or bring all times to common frequency (that could be
# not precise)
def search_antenna(antenna, experiment, dm_min=0, dm_max=1000, path=None):
    """
    Function that search pulse candidates in FITS-files for given experiment
    name and antenna name.

    :param antenna:
        Antenna name.
    :param experiment:
        Experiment name.
    :param dm_min: (optional)
        Minimum DM value to search [cm^3/pc]. (default: ``0``)
    :param dm_max: (optional)
        Maximum DM value to search [cm^3/pc]. (default: ``1000``)
    :param path: (optional)
        Path to files. If ``None`` then search current directory for FITS-files.
        (default: ``None``)
    :return:
    """
    if path:
        path = path + ".FITS"
    fnames = glob.glob(path)
    antenna_fnames = list()
    # Accumulate file names with given antenna/experiment
    for fname in fnames:
        hdu = pf.open(fname)
        # Some artificial FITS-keywords:)
        if (hdu.header['EXPER'] == experiment and
                    hdu.header['TELESC'] == antenna):
            antenna_fnames.append(fname)

    candidates = None
    # For each FITS-file with given antenna name search candidates and save
    for fname in antenna_fnames:
        pars = get_pars(fname)
        frame = DataFrame(fname, *pars)
        dm_grid, dm_t_frame = grid_dedisperse_frame(frame, dm_min=dm_min,
                                                    dm_max=dm_max)
        new_candidates = Objects(dm_t_frame, dm_grid, frame.t)
        if candidates is None:
            candidates = new_candidates
        else:
            candidates += new_candidates

    candidates.save_txt(experiment + '_' + antenna + '.txt')


def get_pars(fname):
    """
    Function that gets parameters needed for processing from FITS headers.
    :param fname:
        FITS-file name.
    :return:
    """
    pass


def search_experiment_antennnas(experiment, antennas, d_t=1., d_dm=150.,
                                path=None):
    """
    Function that search in txt-files for given experiment name and antenna
    names and find close in (t, DM)-space pulse candidates.

    :param experiment:
        Experiment name.
    :param antenna:
        Iterable of antenna names.
    :param d_t: (optional)
        Difference in time of pulse at highest frequency.
    :param path: (optional)
        Path to files. If ``None`` then search current directory for FITS-files.
        (default: ``None``)
    :return:
        Creates txt-file with successful pulse candidates.
    """
    pass