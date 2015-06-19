import glob
import os
import pyfits as pf
import numpy as np
from frames import DataFrame
from objects import TDMImageObjects
from search_close import find_close
import json


# TODO: To compare pulse candidates from different antennas we need time of
# pulse at some fixed frequency. Or use the same set ups for each antenna (that
# is not flexible enough) or bring all times to common frequency (that could be
# not precise)
def search_antenna(antenna, experiment, dm_min, dm_max, d_t, d_dm,
                   dm_delta=None, path=None, outpath=None):
    """
    Function that search pulse candidates in FITS-files for given experiment
    name and antenna name.

    :param antenna:
        Antenna name.
    :param experiment:
        Experiment name.
    :param dm_min:
        Minimum DM value to search [cm^3/pc].
    :param dm_max:
        Maximum DM value to search [cm^3/pc].
    :param path: (optional)
        Path to files. If ``None`` then search current directory for FITS-files.
        (default: ``None``)
    :param outpath: (optional)
        Path to place where keep results. If None then use cwd. (default: None)
    :return:
    """
    path = path or os.getcwd()
    if not path.endswith('/'):
        path += '/'
    path = path + "*.FITS"
    fnames = glob.glob(path)
    if not fnames:
        raise Exception("No FITS-files found in ", path)
    antenna_fnames = list()
    # Accumulate file names with given antenna/experiment
    for fname in fnames:
        hdu = pf.open(fname)
        # Some artificial FITS-keywords:)
        if (hdu.header['EXPER'] == experiment and
                    hdu.header['TELESC'] == antenna):
            antenna_fnames.append(fname)
    if not antenna_fnames:
        raise Exception("No FITS-files for ", antenna, " in ", experiment,
                        " in ", path)

    candidates = None
    # For each FITS-file with given antenna name search candidates and save
    for fname in antenna_fnames:
        pars = get_pars(fname)
        frame = DataFrame(fname, *pars)
        dm_grid, dm_t_frame = frame.grid_dedisperse(dm_min=dm_min,
                                                    dm_max=dm_max,
                                                    dm_delta=dm_delta)
        new_candidates = TDMImageObjects(dm_t_frame, frame.t, dm_grid, d_t,
                                         d_dm)
        if candidates is None:
            candidates = new_candidates
        else:
            candidates += new_candidates

    outpath = outpath or os.getcwd()
    if not outpath.endswith('/'):
        outpath += '/'
    candidates.save_txt(outpath + experiment + '_' + antenna + '.txt', 'x', 'y')


def get_pars(fname):
    """
    Function that gets parameters needed for processing from FITS headers.
    :param fname:
        FITS-file name.
    :return:
    """
    hdulist = pf.open(fname)
    # We don't access ``data`` attribute - just fetch some header info.
    header = hdulist[0].header


def search_experiment_antennnas(experiment, antennas, d_t, d_dm, dm_min, dm_max,
                                dm_delta, path=None, outpath=None):
    """
    Function that search in txt-files for given experiment name and antenna
    names and find close in (t, DM)-space pulse candidates.

    :Note:
        Times ``t`` in files must be for the same frequency

    :param experiment:
        Experiment name.
    :param antennas:
        Iterable of antenna names.
    :param d_t:
        Difference in time of pulse at highest frequency.
    :param path: (optional)
        Path to files. If ``None`` then search current directory for FITS-files.
        (default: ``None``)
    :param outpath: (optional)
        Path to place where keep results. If None then use cwd. (default: None)
    :return:
        Creates txt-file with successful pulse candidates.
    """
    for antenna in antennas:
        search_antenna(antenna, experiment, dm_min, dm_max, d_t, d_dm,
                       dm_delta=dm_delta, path=path, outpath=outpath)


def compare_antennas_candidates(experiment, d_t, d_dm, antennas=None,
                                path=None, outpath=None):
    """
    Function that uses text files with candidates from several antennas to find
    coinciding up to user-specified tolerance levels for ``t`` & ``DM``.

    :param experiment:
        Experiment name.
    :param d_t:
        Tolerance in time [s].
    :param d_dm:
        Tolerance in dispersion measure [cm^3/pc].
    :param antennas: (optional)
        Iterable of antenna names to consider. If None then use all available
        files for given experiment. (default: None)
    :param path: (optional)
        Path to text files. If None then search in cwd. (default: None)
    :param outpath: (optional)
        Path to place where to keep results. If ``None`` then use cwd.
        (default: None)
    :return:
        Text file with antennas (at least 2) and (t, DM)-coordinates of
        confirmed candidates in ``outpath`` directory.
    """
    path = path or os.getcwd()
    if not path.endswith('/'):
        path += '/'
    path = path + experiment + "_" + "*.txt"
    fnames = glob.glob(path)
    if not fnames:
        raise Exception("No files for " + experiment + " found in ", path)

    antennas_results = list()
    for fname in fnames:
        antenna_name = fname.split('_')[-1].split('.')[0]
        values = np.loadtxt(fname)
        antennas_results.append({antenna_name: values})

    results = find_close(antennas_results, {0: d_t, 1: d_dm})


