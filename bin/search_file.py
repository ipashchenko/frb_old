import os
import sys
import argparse
import numpy as np
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
from frb.frames import DataFrame
from frb.objects import BatchedTDMIO
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


if __name__ == '__main__':
    parser =\
        argparse.ArgumentParser(description='Search for dispersed signals in'
                                'dynamic spectra')

    parser.add_argument('fname', type=str,
                        metavar='txt-file with dynamical spectra',
                        help='- path to file with dynamical spectra')
    parser.add_argument('-nu_max', action='store', dest='nu_max', type=float,
                        help='- frequency of highest freq. channel [MHz]')
    parser.add_argument('-t0', action='store', dest='t0', default=0.0,
                        type=float, help='- beginning time [s].'
                                         'Default: 0.0')
    parser.add_argument('-dnu', action='store', dest='dnu', type=float,
                        help='- frequency width of each channel [MHz]')
    parser.add_argument('-dt', action='store', dest='dt', type=float,
                        help='- time resolution [s]')
    parser.add_argument('-dm_min', action='store', dest='dm_min', type=float,
                        help='- mininum DM to search [cm^3/pc]')
    parser.add_argument('-dm_max', action='store', dest='dm_max', type=float,
                        help='- maxmum DM to search [cm^3/pc]')
    parser.add_argument('-perc', action='store', dest='perc', default=99.95,
                        type=float, help='- percentile of image values that is'
                                         ' used to blank image before'
                                         ' searching for objects.'
                                         ' Default: 99.5')
    parser.add_argument('-d_dm', action='store', dest='d_dm', default=100.,
                        type=float, help='- width of feature [cm^3/pc] in'
                                         '(t, DM)-space along DM-axis to treat'
                                         ' it as candidate.'
                                         ' Default: 100.')
    parser.add_argument('-d_t', action='store', dest='d_t', default=0.003,
                        type=float, help='- width of feature [s] in'
                                         '(t, DM)-space along DM-axis to treat'
                                         ' it as candidate.'
                                         ' Default: 0.005')
    parser.add_argument('-savefig_dyn', action='store', nargs='?',
                        default=None, type=str, metavar='path to file',
                        help='- file to save dyn.spectr plot')
    parser.add_argument('-savefig_dedm', action='store', nargs='?',
                        default=None, type=str, metavar='path to file',
                        help='- file to save dedispersioned freq. averaged dyn.'
                             ' spectra')
    parser.add_argument('-save_result', action='store', nargs='?',
                        default=None, type=str, metavar='path to file',
                        help='- file to save (t, DM)-coordinates of candidates')
    parser.add_argument('-threads', action='store', dest='threads', default=1,
                        type=int, help='- number of threads used for '
                                         'parallelization'
                                         ' Default: 1')
    parser.add_argument('-batchsize', action='store', dest='batchsize',
                        default=100000, type=int, help='- size of batch size '
                                                       'in t-direction.'
                                                       'Default: 100000')

    args = parser.parse_args()

    frame = DataFrame(args.fname, args.nu_max, args.t0, args.dnu, args.dt)
    dm_min = args.dm_min
    dm_max = args.dm_max

    dm_grid, frames_t_dedm = frame.grid_dedisperse(dm_min, dm_max,
                                                   savefig=args.savefig_dedm,
                                                   threads=args.threads)
    btdmio = BatchedTDMIO(frames_t_dedm, frame.t, dm_grid, perc=args.perc,
                          d_dm=args.d_dm, dt=args.d_t)
    candidates = btdmio.run(batch_size=args.batchsize)
    print "Candidates: t [s], DM [cm^3/pc]:"
    print "==============================="
    print candidates
    print "==============================="
    if args.save_result:
        np.savetxt(args.save_result, candidates)
