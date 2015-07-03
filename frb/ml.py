import numpy as np
from scipy.ndimage.measurements import (maximum_position, label, find_objects,
                                        mean, minimum, sum, variance, maximum,
                                        median, center_of_mass)
from scipy.ndimage.morphology import generate_binary_structure


def ml(image, t, dm, perc):
    threshold = np.percentile(image.ravel(), perc)
    a = image.copy()
    # Keep only tail of image values distribution with signal
    a[a < threshold] = 0
    s = generate_binary_structure(2, 2)
    # Label image
    labeled_array, num_features = label(a, structure=s)
    # Find objects
    objects = find_objects(labeled_array)
    labels = np.arange(num_features) + 1
    # Calculate features of objects
    max_pos = maximum_position(image, labels=labeled_array, index=labels)
    means = mean(image, labels=labeled_array, index=labels)
    mins = minimum(image, labels=labeled_array, index=labels)
    maxs = maximum(image, labels=labeled_array, index=labels)
    meds = median(image, labels=labeled_array, index=labels)
    cmass = center_of_mass(image, labels=labeled_array, index=labels)
    sums = sum(image, labels=labeled_array, index=labels)
    vars = variance(image, labels=labeled_array, index=labels)
    dx = [int(obj[1].stop - obj[1].start) for obj in objects]
    dy = [int(obj[0].stop - obj[0].start) for obj in objects]
    return num_features, objects, labeled_array, dx, dy, means, mins, maxs,\
           meds, sums, vars, cmass, max_pos


if __name__ == '__main__':
    import time
    from frames import DataFrame
    fname = '/home/ilya/code/frb/data/crab_600sec_64ch_1ms.npy'
    frame = DataFrame(fname, 1684., 0., 16. / 64., 0.001)
    # frame.add_pulse(10., 9.0, 0.001, dm=1000.)
    # frame.add_pulse(20., 2.0, 0.003, dm=500.)
    # frame.add_pulse(30., 1.0, 0.003, dm=500.)
    # frame.add_pulse(40., 0.5, 0.003, dm=500.)
    # frame.add_pulse(50., 0.25, 0.003, dm=500.)
    # frame.add_pulse(60., 0.125, 0.003, dm=500.)
    # frame.add_pulse(70., 0.0625, 0.003, dm=500.)
    # frame.add_pulse(80., 0.03125, 0.003, dm=500.)
    t0 = time.time()
    dm_grid = frame.create_dm_grid(0., 1000.)
    frames_t_dedm = frame.grid_dedisperse(dm_grid, threads=4)
    num_features, objects, labeled_array, dx, dy, means, mins, maxs, \
           meds, sums, vars, cmass, max_pos = ml(frames_t_dedm, frame.t,
                                                 dm_grid, 99.75)
