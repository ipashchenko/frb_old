import os
import sys
import math
import argparse
import numpy as np
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
from frb.frames import DataFrame
from frb.utils import delta_dm_max
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

    # Print and (optionally) plot results
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
    if args.savefig:
        plt.savefig(args.savefig, bbox_inches='tight')
        plt.close()
