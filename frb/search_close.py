import numpy as np
from itertools import combinations


def find_close_2(arr1, arr2, ind_dict):
    """
    :return:
    """
    successful_candidates = list()
    for el in arr1:
        diff = arr2 - el
        indxs = list()
        for key, value in ind_dict.items():
            indxs.append(abs(diff[:, key]) < value)
        indxs = np.vstack(indxs).T
        for i, indx in enumerate(indxs):
            if np.alltrue(indx):
                successful_candidates.append((el, arr2[i],))

    return successful_candidates


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
