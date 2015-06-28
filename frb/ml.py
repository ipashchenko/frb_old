import numpy as np
from scipy.ndimage.measurements import (maximum_position, label, find_objects,
                                        mean, minimum, sum, variance)
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
    max_pos = maximum_position(image, labels=labeled_array, index=label)

    dx = [int(obj[1].stop - obj[1].start) for obj in objects]
    dy = [int(obj[0].stop - obj[0].start) for obj in objects]
