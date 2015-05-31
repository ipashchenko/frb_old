import os
import sys
import math
import argparse
import numpy as np
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
from frb.frames import DataFrame
from frb.search_utils import search_frame
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


if __name__ == '__main__':
    parser =\
        argparse.ArgumentParser(description='Search dispersed signals in'
                                'dynamic spectra')

    parser.add_argument('fname', type=str,
                        metavar='txt-file with dynamical spectra',
                        help='- path to file with dynamical spectra')
    parser.add_argument('-nu_max', action='store', dest='nu_max', type=float,
                        help='- frequency of highest freq. channel')
    parser.add_argument('-t0', action='store', dest='t0', default=0.0,
                        type=float, help='- beginning time')
    parser.add_argument('-dnu', action='store', dest='dnu', type=float,
                        help='- frequency width of each channel')
    parser.add_argument('-dt', action='store', dest='dt', type=float,
                        help='- time resolution')
    parser.add_argument('-dm_min', action='store', dest='dm_min', type=float,
                        help='- mininum DM to search')
    parser.add_argument('-dm_max', action='store', dest='dm_max', type=float,
                        help='- maxmum DM to search')
    parser.add_argument('-savefig', action='store', nargs='?', default=None,
                        type=str, metavar='path to file', help='- file to save'
                                                               ' figure')

    args = parser.parse_args()

    frame = DataFrame(args.fname, args.nu_max, args.t0, args.dnu, args.dt)
    frame.plot()
    dm_min = args.dm_min
    dm_max = args.dm_max

    search_frame(frame, dm_min, dm_max, savefig=args.savefig)
