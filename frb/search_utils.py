import os
import sys
import math
import numpy as np
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
from frb.utils import delta_dm_max
from frb.frames import vint
from scipy.ndimage.morphology import generate_binary_structure
from scipy.ndimage.measurements import find_objects, label, maximum_position
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


# TODO: add ranking of output pulses based on square of it's label
def search_frame(frame, dm_min, dm_max, savefig=None):
    # Find step for DM grid
    # Seems that ``5`` is good choice (1/200 of DM range)
    dm_delta = 5 * delta_dm_max(frame.nu_0, frame.nu_0 - frame.n_nu * frame.dnu,
                                frame.dt)

    # Find length of DM grid
    # Actually, to build grid i can use
    # np.arange(dm_min, dm_max, dm_delta)
    n_dm = (dm_max - dm_min) / dm_delta
    n_dm = int(math.ceil(n_dm))
    dm_used = list()
    frames_dedm = list()
    for i in xrange(n_dm):
        dm = dm_min + (i + 1) * dm_delta
        print "Searching DM = ", dm, " cm^3 / pc"
        _frame = frame.de_disperse(dm=dm, replace=False)
        _frame_t = frame.average_in_freq(_frame)
        frames_dedm.append(_frame_t)
        dm_used.append(dm)

    dm_used = np.array(dm_used)
    frames_dedm = np.array(frames_dedm)

    return dm_used, frames_dedm

    #  Print and (optionally) plot results
    # if savefig is not None:
    #     plt.savefig(savefig, bbox_inches='tight')
    #     plt.close()


if __name__ == '__main__':
    from frb.frames import DataFrame
    fname = '/home/ilya/code/frb/data/data.txt'
    # pulses at t=0.2 and 0.65 with DM=480
    frame = DataFrame(fname, 1684., 0., 16. / 32., 0.001)
    # frame = SimFrame(128, 1000, 1676., 0., 16. / 128., 0.001)
    frame.add_pulse(0.3, 0.05, 0.0015, dm=200.)
    frame.add_pulse(0.9, 0.1, 0.003, dm=700.)
    # plt.close()
    # plt.imshow(frame.values, interpolation='none', aspect='auto')
    # plt.xlabel('Time')
    # plt.ylabel('Freq. channel')
    # plt.colorbar()
    # plt.savefig('pulse_clean1.png')
    frame.add_noise(0.1)
    # plt.close()
    # plt.imshow(frame.values, interpolation='none', aspect='auto')
    # plt.xlabel('Time')
    # plt.ylabel('Freq. channel')
    # plt.colorbar()
    # plt.savefig('pulse_dirty1.png')
    dm_used, frames_dedm = search_frame(frame, 0, 1000.)
    # plt.close()
    # plt.imshow(frames_dedm, interpolation='none', aspect='auto')
    # plt.xlabel('De-dispersed by DM freq.averaged frame')
    # plt.ylabel('DM correction')
    # plt.yticks(np.linspace(0, len(dm_used) - 10, 5, dtype=int),
    #            vint(dm_used[np.linspace(0, len(dm_used) - 10, 5, dtype=int)]))
    # plt.colorbar()
    # plt.savefig('stripes1.png')

    # Now find stripes in "DM-correction vs. freq.averaged" plane
    threshold = np.percentile(frames_dedm.flatten(), 98)
    a = frames_dedm.copy()
    # Keep only tail of distribution with signal (and correlated noise:)
    a[a < threshold] = 0
    s = generate_binary_structure(2, 2)
    # Label image
    labeled_array, num_features = label(a, structure=s)

    objects = find_objects(labeled_array)
    # Find only objects that are elongated along DM axis
    dm_objects = list()
    for obj in objects:
        freq_slice = obj[0]
        if int(freq_slice.stop - freq_slice.start) > 50:
            dm_objects.append(obj)

    # Find what labels number corresponds to such objects
    n_labels = list()
    for obj in dm_objects:
        # Counts of all labels but zero
        counts = np.bincount(labeled_array[obj].flatten())[1:]
        n_label = int(np.where(counts == max(counts))[0] + 1)
        n_labels.append(n_label)

    # center_of_mass(frames_dedm, labeled_array, n_labels)
    pos = maximum_position(frames_dedm, labels=labeled_array, index=n_labels)
    print 'Found ', len(n_labels), ' pulses'
    for i in range(len(n_labels)):
        print "at time ", pos[i][1] / 1000., ' with DM = ', dm_used[pos[i][0]]
