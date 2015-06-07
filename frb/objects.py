import numpy as np
from scipy.ndimage.measurements import maximum_position, label, find_objects,\
    mean, maximum, minimum
from scipy.ndimage.morphology import generate_binary_structure
from sklearn.ensemble import GradientBoostingClassifier
from frames import DataFrame
from search_utils import grid_dedisperse_frame


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


class ObjectExplorer(object):
    def __init__(self, objects):
        self.objects = objects


class Objects(object):
    """
    Class that describes collections of labeled regions in image.
    """
    def __init__(self, image, dm_grid, perc=99.95):
        self.dm_grid = dm_grid
        self.dm_step = dm_grid[1] - dm_grid[0]
        threshold = np.percentile(image.ravel(), perc)
        a = image.copy()
        # Keep only tail of distribution with signal (and correlated noise:)
        a[a < threshold] = 0
        s = generate_binary_structure(2, 2)
        # Label image
        labeled_array, num_features = label(a, structure=s)
        # Find objects
        # TODO: Check that objects in list are sorted (thus, 1 goes first, ...)
        objects = find_objects(labeled_array)
        # Container of object's properties
        _objects = np.empty(num_features, dtype=[('label', 'int'),
                                                 ('dx', '<f8'),
                                                 ('dy', '<f8'),
                                                 ('min', '<f8'),
                                                 ('mean', '<f8'),
                                                 ('max', '<f8'),
                                                 ('max_pos', 'int',
                                                  (2,))])

        labels = np.arange(num_features) + 1
        max_pos = maximum_position(image, labels=labeled_array,
                                        index=labels)
        _mean = mean(image, labeled_array, index=labels)
        _max = maximum(image, labeled_array, index=labels)
        _min = minimum(image, labeled_array, index=labels)
        dx = [int(obj[1].stop - obj[1].start) for obj in objects]
        dy = [int(obj[0].stop - obj[0].start) for obj in objects]

        # Filling objects structured array
        _objects['label'] = labels
        _objects['dx'] = dx
        _objects['dy'] = dy
        _objects['min'] = _min
        _objects['mean'] = _mean
        _objects['max'] = _max
        _objects['max_pos'] = max_pos
        self.objects = _objects

    @property
    def dx(self):
        return self.objects['dx']

    @property
    def dy(self):
        return self.objects['dy']

    @property
    def d_t(self):
        return self.objects['dx']

    @property
    def d_dm(self):
        return self.objects['dy'] * (self.dm_step)

    @property
    def labels(self):
        return self.objects['labels']

    @property
    def max_pos(self):
        return self.objects['max_pos']

    @property
    def max(self):
        return self.objects['max']

    @property
    def min(self):
        return self.objects['min']

    @property
    def mean(self):
        return self.objects['mean']

    @property
    def dm(self):
        return self.dm_grid[self.objects['max_pos'][:, 0]]

    @property
    def t(self):
        return self.objects['max_pos'][:, 1]


if __name__ == '__main__':
    fname = '/home/ilya/code/frb/data/90_sec_wb_raes08a_128ch'
    frame = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    frame.add_pulse(10., 0.3, 0.003, dm=500.)
    frame.add_pulse(20., 0.275, 0.003, dm=500.)
    frame.add_pulse(30., 0.25, 0.003, dm=500.)
    frame.add_pulse(40., 0.225, 0.003, dm=500.)
    frame.add_pulse(50., 0.2, 0.003, dm=500.)
    frame.add_pulse(60., 0.175, 0.003, dm=500.)
    frame.add_pulse(70., 0.15, 0.003, dm=500.)
    frame.add_pulse(80., 0.125, 0.003, dm=500.)
    dm_grid, frames_t_dedm = grid_dedisperse_frame(frame, 0, 1000.)
    objects = Objects(frames_t_dedm, dm_grid)
