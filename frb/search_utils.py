import os
import sys
import numpy as np
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
from frb.utils import delta_dm_max, vint
from scipy.ndimage.morphology import generate_binary_structure
from scipy.ndimage.measurements import find_objects, label, maximum_position
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def grid_dedisperse_frame(frame, dm_min, dm_max, savefig=None):
    # Find step for DM grid
    # Seems that ``5`` is good choice (1/200 of DM range)
    dm_delta = 5 * delta_dm_max(frame.nu_0, frame.nu_0 - frame.n_nu * frame.dnu,
                                frame.dt)

    # Create grid of searched DM-values
    dm_grid = np.arange(dm_min, dm_max, dm_delta)
    # Accumulator of de-dispersed frequency averaged frames
    frames_t_dedm = list()
    for dm in dm_grid:
        _frame = frame.de_disperse(dm=dm, replace=False)
        _frame_t = frame.average_in_freq(_frame)
        frames_t_dedm.append(_frame_t)

    frames_t_dedm = np.array(frames_t_dedm)

    # Plot results
    if savefig is not None:
        plt.imshow(frames_t_dedm, interpolation='none', aspect='auto')
        plt.xlabel('De-dispersed by DM freq.averaged dyn.spectr')
        plt.ylabel('DM correction')
        plt.yticks(np.linspace(0, len(dm_grid) - 10, 5, dtype=int),
                   vint(dm_grid[np.linspace(0, len(dm_grid) - 10, 5,
                                            dtype=int)]))
        plt.colorbar()
        plt.savefig(savefig, bbox_inches='tight')
        plt.show()
        plt.close()

    return dm_grid, frames_t_dedm


# TODO: add ranking of output pulses based on square of it's size
# TODO: ``dm_range`` here in pixels.
def find_pulses(dm_grid, frames_t_dedm, perc=99.95, dm_range=50):
    """
    Search frequency averaged de-dispersed dynamical spectra dependency on DM
    to find broadband pulses.

    :param dm_grid:
        Array-like of used DM values to produce ``frames_t_dedm``.
    :param frames_t_dedm:
        Sequence of frequency averaged de-dispersed dynamical spectra.
    :return:
        Pulse objects.
    """
    # TODO: fit Rayleigh to density of ``frames_t_dedm`` values to find
    # ``threshold``
    # Now find stripes in "DM-correction vs. freq.averaged" plane
    # TODO: I should find no more then ``nobjects`` objects
    threshold = np.percentile(frames_t_dedm.ravel(), perc)
    a = frames_t_dedm.copy()
    # Keep only tail of distribution with signal (and correlated noise:)
    a[a < threshold] = 0
    s = generate_binary_structure(2, 2)
    # Label image
    labeled_array, num_features = label(a, structure=s)
    print "Checking ", num_features, " objects..."

    objects = find_objects(labeled_array)
    # Find only objects that are elongated along DM axis
    dm_objects = list()
    for obj in objects:
        dm_slice = obj[0]
        # TODO: convert DM to pixel numbers
        # Classification step
        if int(dm_slice.stop - dm_slice.start) > dm_range:
            dm_objects.append(obj)

    # Find what labels number corresponds to such objects
    n_labels = list()
    for obj in dm_objects:
        # Counts of all labels but zero
        counts = np.bincount(labeled_array[obj].ravel())[1:]
        n_label = int(np.where(counts == max(counts))[0] + 1)
        n_labels.append(n_label)

    # center_of_mass(frames_t_dedm, labeled_array, n_labels) is poor choice
    pos = maximum_position(frames_t_dedm, labels=labeled_array, index=n_labels)
    print 'Found ', len(n_labels), ' pulses'
    for i in range(len(n_labels)):
        print "at time ", pos[i][1] / 1000., ' with DM = ', dm_grid[pos[i][0]]


def n_objects(image, threshold):
    """
    Find number of objects in image using specified threshold value.
    :param image:
        2D-Numpy array of image
    :param threshold:
        Threshold value. Image regions with values less then ``threshold`` will
        be zero'd.
    :return:
        Number of objects found.

    """
    a = image.copy()
    a[a < threshold] = 0
    s = generate_binary_structure(2, 2)
    # Label image
    labeled_array, num_features = label(a, structure=s)
    return num_features


if __name__ == '__main__':
    from frb.frames import DataFrame
    fname = '/home/ilya/code/frb/data/90_sec_wb_raes08a_128ch'
    frame_orig = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    frame = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    # frame.add_pulse(30., 0.3, 0.0015, dm=200.)
    # frame.add_pulse(60., 0.2, 0.003, dm=700.)
    # frame.add_pulse(10., 0.3, 0.003, dm=500.)
    # frame.add_pulse(40., 0.2, 0.003, dm=400.)
    # for i in range(100):
    #     frame.add_pulse(np.random.uniform(0, 90), np.random.uniform(0.05, 0.3),
    #                     np.random.uniform(0.0005, 0.004),
    #                     dm=np.random.uniform(0, 1000))
    dm_grid, frames_t_dedm = grid_dedisperse_frame(frame, 0, 1000.)
    # find_pulses(dm_grid, frames_t_dedm)
