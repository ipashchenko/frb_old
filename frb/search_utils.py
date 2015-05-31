import os
import sys
import math
import numpy as np
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
from frb.utils import delta_dm_max
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def search_frame(frame, dm_min, dm_max, savefig=None):
    dm_delta = delta_dm_max(frame.nu_0, frame.nu_0 - frame.n_nu * frame.dnu,
                            frame.dt)

    n_dm = (dm_max - dm_min) / dm_delta
    n_dm = int(math.ceil(n_dm))
    max_snrs = list()
    dm_used = list()
    for i in xrange(n_dm):
        dm = dm_min + (i + 1) * dm_delta
        print "Searching DM = ", dm, " cm^3 / pc"
        _frame = frame.de_disperse(dm=dm, replace=False)
        _frame_t = frame.average_in_freq(_frame)
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
    plt.show()
    if savefig is not None:
        plt.savefig(savefig, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    from frb.frames import SimFrame
    frame = SimFrame(128, 1000, 1676., 0., 16. / 128., 0.001)
    frame.add_pulse(0.3, 0.5, 0.0015, dm=2400.)
    frame.save_to_txt('fake_clean.txt')
    frame.plot()
    frame.add_noise(3.00)
    frame.save_to_txt('fake_dirty.txt')
    search_frame(frame, 0, 3000.)
