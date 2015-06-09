import numpy as np
from scipy.ndimage.measurements import maximum_position, label, find_objects
from scipy.ndimage.morphology import generate_binary_structure
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


# TODO: Sometimes we need more features for classification. Currently it uses
# only ``dx``, ``dy`` for that. Should be an option to include others (``max``,
# ``mean``, ...) but in that case calculations of them for all labelled regions
# will take a lot of time (there are thousands of objects in my use case).
# Current implementation allows this option by means of ``_classify`` method. It
# should calculate necessary features for objects in ``objects`` array using
# original image and labeled array passed as arguments.
class BasicImageObjects(object):
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
        # Keep only tail of image values distribution with signal
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
        self._classify(image, labeled_array)
        self._sort()
        # Fetch positions of objects only on classified ones
        self.max_pos = self._find_positions(image, labeled_array)

    def _find_positions(self, image, labeled_array):
        return maximum_position(image, labels=labeled_array, index=self.label)

    def _sort(self):
        """
        Method that sorts image objects somehow.
        """
        raise NotImplementedError

    def _classify(self, *args, **kwargs):
        """
        Method that selects only image objects with desirable properties.
        You can use any features of ``objects`` attribute to do classification.
        Or you can fetch any other features (like ``max`` or ``variance``)
        using passed as arguments to method image and labeled_array.
        """
        raise NotImplementedError

    def plot(self, image, labels=None):
        """
        Overplot image with found labelled objects,
        """
        pass

    def save_txt(self, fname):
        """
        Save image object's parameters to text file.
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


class ImageObjects(BasicImageObjects):
    """
    Class that handles image objects for images with specified x,y -
    coordinates.

    """
    def __init__(self, image, x_grid, y_grid, perc):
        super(ImageObjects, self).__init__(image, perc)
        self.x_grid = x_grid
        self.y_grid = y_grid
        self.x_step = x_grid[1] - x_grid[0]
        self.y_step = y_grid[1] - y_grid[0]

    def __add__(self, other):
        values = other.objects.copy()
        # Keep each own's numbering to show it later.
        # values['label'] += len(self.objects)
        self.objects = np.concatenate((self.objects, values))
        self._sort()
        return self

    @property
    def d_x(self):
        return self.objects['dx'] * self.x_step

    @property
    def d_y(self):
        return self.objects['dy'] * self.y_step

    @property
    def y(self):
        return self.y_grid[self.max_pos[:, 0]]

    @property
    def x(self):
        return self.x_grid[self.max_pos[:, 1]]

    @property
    def xy(self):
        return np.vstack((self.x, self.y)).T


class TDMImageObjects(ImageObjects):

    def _sort(self):
        self.objects = self.objects[np.lexsort((self.dx, self.dy))[::-1]]

    def _classify(self, d_dm=150, dt=0.005):
        """
        Method that select only candidates which have dimensions > ``d_dm``
        [cm*3/pc] and > ``dt`` [s]
        :param d_dm:
            Value of DM spanned by object to count it as candidate for pulse
            [cm^3/pc].
        :param dt:
            Value of t spanned by object to count it as candidate for pulse
            [s].
        """
        self.objects = self.objects[np.logical_and(self.d_y > d_dm,
                                                   self.d_x > dt)]


if __name__ == '__main__':
    fname = '/home/ilya/code/frb/data/90_sec_wb_raes08a_128ch'
    # fname = '/home/ilya/code/frb/data/BIG.txt'
    frame = DataFrame(fname, 1684., 0., 16. / 128., 0.001)
    frame.add_pulse(10., 0.3, 0.003, dm=500.)
    frame.add_pulse(20., 0.275, 0.003, dm=500.)
    frame.add_pulse(30., 0.25, 0.003, dm=500.)
    frame.add_pulse(40., 0.225, 0.003, dm=500.)
    frame.add_pulse(50., 0.2, 0.003, dm=500.)
    frame.add_pulse(60., 0.175, 0.003, dm=500.)
    frame.add_pulse(70., 0.15, 0.003, dm=500.)
    frame.add_pulse(80., 0.125, 0.003, dm=500.)
    dm_grid, frames_t_dedm = frame.grid_dedisperse(0, 1000.)
    objects = TDMImageObjects(frames_t_dedm, dm_grid, frame.t)
