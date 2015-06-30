FRB
=====

## Searching FRB in auto-spectral data

frb - tool used to search Fast Radio Bursts in auto-spectral data

## Requirements:

numpy, scipy, matplotlib.pyplot (for plots)

## Searching for pulses in txt-format data file:

``user@host:~$ python search_file.py data.dat -nu_max NU_MAX -dnu DNU -dt DT
-dm_min DM_MIN -dm_max DM_MAX [-batchsize BATCHSIZE] [-d_t D_T] [-d_dm D_DM] [-perc PERC] [-savefig_dyn
fig.png] [-savefig_dedm fig.png] [-threads THREADS]``

Parameters:

- ``data.dat`` - binary or text file with time sequence of dynamical spectra.
    If text file => # of columns = # of frequency channels and # of rows = # of
    time measurements. If binary then ``np.shape(np.load('data.dat'))`` = (# of
    freq. channels, # of time measurements,)

- ``-nu_max`` - frequency of highest frequency channel [MHz].

- ``-dnu`` - frequency width of single frequency channel [MHz].

- ``-dt`` - time step (resolution) [s].

- ``dm_min`` - minimal value of DM window to search [cm^3 / pc].

- ``dm_max`` - maximum value of DM window to search [cm^3 / pc].

- ``batchsize`` - size of image in t-direction, that will be searched for
    candidates in batches. Default: 100000

- ``d_t`` - width of feature [s] in (t, DM)-space along t-axis to treat it as
    candidate. Default: 0.003

- ``d_dm`` - width of feature [cm^3/pc] in (t, DM)-space along DM-axis to treat
    it as candidate. Default: 100.

- ``perc`` - percentile of image values that is used to blank image before
    searching for objects. Default: 99.5

- ``savefig_dyn`` - file name for saving picture of dynamical spectra.

- ``savefig_dedm`` - file name for saving picture of de-dispersed frequency
    averaged dynamical spectra.

- ``save_result`` - file name to save (t, DM)-coordinates of found candidates.
    [s, cm^3/pc]
    
- ``threads`` - number of threads used for parallelization of grid
    de-dispersion. Default is ``1`` (don't use parallelization).

## Notes

Algorithm searches for extended regions in image of de-dispersed frequency
averaged dynamical spectra (that is (t, DM)-plane). Currently there are 3
tunable parameters:

- ``perc``. Current experience suggests values ``99.9 - 99.95``. Low value could
    bring many false features in (t, DM)-space, but as long as we can compare
    results using different telescopes it doesn't seem to be an issue. Value
    that is too high can split characteristic x-shaped dispersed signal
    features.

- ``d_t`` & ``d_dm``. Some experience with fake FRB injected in real data have
    shown that 2 most informative features in classification of FRB candidates
    on (t, DM)-plane are their widths in t- and DM-directions. Currently it is
    the only method of classification that has been implemented. Nonetheless one
    can use any other features and their own algorithms of classification by
    overriding/extending ``TDMImageObjects._classify`` method that gets
    ``image`` and ``labelled array`` as first two positional arguments.

Current implementation allows parallelization of grid de-dispersion step using
``multiprocessing`` module.

Using text files requires additional RAM that can be a problem for large data
sets.


License
-------

Copyright 2015 Ilya Pashchenko.

frb is free software made available under the MIT License. For details see the
LICENSE file.
