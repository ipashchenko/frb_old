import numpy as np
from scipy.ndimage.measurements import maximum_position, label, find_objects
from scipy.ndimage.morphology import generate_binary_structure
from frames import DataFrame
from search_utils import grid_dedisperse_frame
from itertools import combinations


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


class ImageObjects(object):
    """
    Abstract class for finding and handling image objects.

    :param image:
        2D numpy array with image.
    :param perc:
        Percent of image values to blank while labeling image with objects.

    """
    def __init__(self, image, perc):
        threshold = np.percentile(image.ravel(), perc)
        a = image.copy()
        # Keep only tail of distribution with signal (and correlated noise:)
        a[a < threshold] = 0
        s = generate_binary_structure(2, 2)
        # Label image
        labeled_array, num_features = label(a, structure=s)
        # Find objects
        objects = find_objects(labeled_array)
        # Container of object's properties
        _objects = np.empty(num_features, dtype=[('label', 'int'),
                                                 ('dx', '<f8'),
                                                 ('dy', '<f8'),
                                                 ('max_pos', 'int',
                                                  (2,))])

        labels = np.arange(num_features) + 1
        dx = [int(obj[1].stop - obj[1].start) for obj in objects]
        dy = [int(obj[0].stop - obj[0].start) for obj in objects]

        # Filling objects structured array
        _objects['label'] = labels
        _objects['dx'] = dx
        _objects['dy'] = dy
        self.objects = _objects
        self._classify()
        self._sort()
        self.max_pos = self._find_positions(image, labeled_array)

    def _find_positions(self, image, labeled_array):
        return maximum_position(image, labels=labeled_array, index=self.label)

    def _sort(self):
        raise NotImplementedError

    def _classify(self):
        """
        Method that select only image objects with desirable properties.
        """
        raise NotImplementedError

    def plot(self, image, labels=None):
        """
        Overplot image with found labelled objects,
        """
        pass

    def save_txt(self, fname):
        """
        Save found object to text file.
        """
        np.savetxt(fname, self.objects)

    @property
    def dx(self):
        return self.objects['dx']

    @property
    def dy(self):
        return self.objects['dy']

    @property
    def label(self):
        return self.objects['label']

    @property
    def max_pos(self):
        return self.objects['max_pos']

    @max_pos.setter
    def max_pos(self, max_pos):
        self.objects['max_pos'] = max_pos

    def __len__(self):
        return len(self.objects)


class Objects(object):
    """
    Class that describes collections of labeled regions in image.

    :param image:
        Numpy 2D array with image values.
    :param dm_grid:
        Array-like of y-coordinates.
    :param t_grid:
        Array-like of x-coordinates.

    :notes:
        Instances of this class doesn't need original image so it could be
        deleted right after instantiation.
    """
    # TODO: Add parameters of classification - should it be function that will
    # be ``map``-ed to sequence of candidates or ``filter`` them?
    # TODO: Abstract class shouldn't use ``dm_grid`` & ``t_grid``.
    def __init__(self, image, dm_grid, t_grid, perc=99.95):
        self.t_grid = t_grid
        self.dm_grid = dm_grid
        self.dm_step = dm_grid[1] - dm_grid[0]
        self.t_step = t_grid[1] - t_grid[0]
        threshold = np.percentile(image.ravel(), perc)
        a = image.copy()
        # Keep only tail of distribution with signal (and correlated noise:)
        a[a < threshold] = 0
        s = generate_binary_structure(2, 2)
        # Label image
        labeled_array, num_features = label(a, structure=s)
        # Find objects
        objects = find_objects(labeled_array)
        # Container of object's properties
        _objects = np.empty(num_features, dtype=[('label', 'int'),
                                                 ('dx', '<f8'),
                                                 ('dy', '<f8'),
                                                 # ('min', '<f8'),
                                                 # ('mean', '<f8'),
                                                 # ('max', '<f8'),
                                                 ('max_pos', 'int',
                                                  (2,))])

        labels = np.arange(num_features) + 1
        dx = [int(obj[1].stop - obj[1].start) for obj in objects]
        dy = [int(obj[0].stop - obj[0].start) for obj in objects]

        # Filling objects structured array
        _objects['label'] = labels
        _objects['dx'] = dx
        _objects['dy'] = dy
        self.objects = _objects
        self._classify()
        self._sort()
        self.max_pos = self._find_positions(image, labeled_array)

    def _find_positions(self, image, labeled_array):
        return maximum_position(image, labels=labeled_array, index=self.label)

    # TODO: Sort using square of size (dx * dy)
    def _sort(self):
        self.objects = self.objects[np.lexsort((self.dx, self.dy))[::-1]]

    def _classify(self, d_dm=150, dt=0.005):
        """
        Method that select only candidates which have dimensions > ``d_dm``
        [cm*3/pc] and > ``dt`` [s]
        :param d_dm:
        :param dt:
        :return:
        """
        self.objects = self.objects[np.logical_and(self.d_dm > d_dm,
                                                   self.d_t > dt)]
        self._sort()

    def plot(self, image, labels=None):
        """
        Overplot image with found labelled objects,
        """
        pass

    def save_txt(self, fname):
        np.savetxt(fname, self.objects)

    # TODO: I need it only in subclasses of ABC (when ``t_grid`` is available).
    def __add__(self, other):
        values = other.objects.copy()
        # Keep each own's numbering to show it later.
        # values['label'] += len(self.objects)
        self.objects = np.concatenate((self.objects, values))
        self._sort()
        return self

    @property
    def dx(self):
        return self.objects['dx']

    @property
    def dy(self):
        return self.objects['dy']

    # TODO: Should be available if grids are given
    @property
    def d_t(self):
        return self.objects['dx'] * self.t_step

    # TODO: Should be available if grids are given
    @property
    def d_dm(self):
        return self.objects['dy'] * self.dm_step

    @property
    def label(self):
        return self.objects['label']

    @property
    def max_pos(self):
        return self.objects['max_pos']

    @max_pos.setter
    def max_pos(self, max_pos):
        self.objects['max_pos'] = max_pos

    # TODO: Should be available if grids are given
    @property
    def dm(self):
        return self.dm_grid[self.max_pos[:, 0]]

    # TODO: Should be available if grids are given
    @property
    def t(self):
        return self.t_grid[self.max_pos[:, 1]]

    # TODO: Should be available if grids are given
    @property
    def tdm(self):
        return np.vstack((self.t, self.dm)).T

    def __len__(self):
        return len(self.objects)


# FIXME: Currently it saves only one of two close candidates. Should save 2.
# TODO: Create abstract function for case of 2 arrays.
def find_close_2(obj1, obj2, dt=0.1, dm=200):
    """
    Function that finds close objects in (t, DM)-space for 2 ``Objects``
    instances.

    :param dt:
        Threshold for time difference to consider objects close.
    :param dm:
        Threshold for DM difference to consider objects close.
    :return:
    """
    successful_candidates = list()
    for candidate_tdm in obj1.tdm:
        diff = abs(obj2.tdm - candidate_tdm)
        indxs = np.logical_and(diff[:, 0] < dt, diff[:, 1] < dm)
        if indxs.any():
            successful_candidates.append(obj2.tdm[indxs])

    return np.vstack(successful_candidates)


def find_close_many(objects, dt=0.1, dm=200):
    """
    Function that finds close objects in (t, DM)-space.

    :param dt:
        Threshold for time difference to consider objects close.
    :param dm:
        Threshold for DM difference to consider objects close.
    :return:
    """
    # Remove empty objects (with no candidates)
    for object in objects:
        if not object:
            objects.remove(object)
    successful_candidates = list()
    for pair in combinations(objects, 2):
        # Reduce the length of for-loop in ``find_close_2``
        if len(pair[0]) > len(pair[1]):
            pair = pair[::-1]
        successful_candidates.append(find_close_2(pair[0], pair[1], dt=dt,
                                                  dm=dm))
    return unique_rows(np.vstack(successful_candidates))


def unique_rows(a):
    """
    Find unique rows in 2D numpy array.
    see Stackexchange.
    """
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1)
    return a[ui]


if __name__ == '__main__':
    fname = '/home/ilya/code/frb/data/90_sec_wb_raes08a_128ch'
    # fname = '/home/ilya/code/frb/data/BIG.txt'
    frame1 = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    frame1.add_pulse(10., 0.3, 0.003, dm=500.)
    frame1.add_pulse(20., 0.275, 0.003, dm=500.)
    frame1.add_pulse(30., 0.25, 0.003, dm=500.)
    frame1.add_pulse(40., 0.225, 0.003, dm=500.)
    frame1.add_pulse(50., 0.2, 0.003, dm=500.)
    frame1.add_pulse(60., 0.175, 0.003, dm=500.)
    frame1.add_pulse(70., 0.15, 0.003, dm=500.)
    frame1.add_pulse(80., 0.125, 0.003, dm=500.)
    #frame2 = DataFrame(fname, 1684., 1., 16. / 128., 0.001)
    #frame2.add_pulse(15., 0.3, 0.003, dm=500.)
    #frame2.add_pulse(20., 0.275, 0.003, dm=500.)
    #frame2.add_pulse(35., 0.25, 0.003, dm=500.)
    #frame2.add_pulse(40., 0.225, 0.003, dm=500.)
    #frame2.add_pulse(55., 0.2, 0.003, dm=500.)
    #frame2.add_pulse(60., 0.175, 0.003, dm=500.)
    #frame2.add_pulse(75., 0.15, 0.003, dm=500.)
    #frame2.add_pulse(80., 0.125, 0.003, dm=500.)
    #frame3 = DataFrame(fname, 1684., 1., 16. / 128., 0.001)
    #frame3.add_pulse(15., 0.3, 0.003, dm=500.)
    #frame3.add_pulse(20., 0.275, 0.003, dm=500.)
    #frame3.add_pulse(39., 0.25, 0.003, dm=500.)
    #frame3.add_pulse(40., 0.225, 0.003, dm=500.)
    #frame3.add_pulse(51., 0.2, 0.003, dm=500.)
    #frame3.add_pulse(65., 0.175, 0.003, dm=500.)
    #frame3.add_pulse(75., 0.15, 0.003, dm=500.)
    #frame3.add_pulse(87., 0.125, 0.003, dm=500.)
    dm_grid1, frames_t_dedm1 = frame1.grid_dedisperse(0, 1000.)
    #dm_grid2, frames_t_dedm2 = frame2.grid_dedisperse(0, 1000.)
    #dm_grid3, frames_t_dedm3 = frame3.grid_dedisperse(0, 1000.)
    objects1 = Objects(frames_t_dedm1, dm_grid1, frame1.t)
    #objects2 = Objects(frames_t_dedm2, dm_grid2, frame2.t)
    #objects3 = Objects(frames_t_dedm3, dm_grid3, frame3.t)
    #result = find_close_many([objects1, objects2, objects3], dt=1., dm=150.)