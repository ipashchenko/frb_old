FRB
=====

**Searching FRB in RA data**

frb - tool used to search FRB in Radioastron data

Documentation
-------------

Requirements:
^^^^^^^^^^^^^
numpy, scipy, matplotlib.pyplot (for plots)

Searching for pulses in txt-format data file:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``user@host:~$ python search_file.txt data.txt -nu_max NU_MAX -dnu DNU -dt DT
-dm_min DM_MIN -dm_max DM_MAX -savefig_dyn fig.png -savefig_dedm fig.png``

Parameters:

- ``-nu_max`` - frequency of highest frequency channel [MHz].

- ``-dnu`` - frequency with of single frequency channel [MHz].

- ``-dt`` - time step (resolution) [s].

- ``dm_min`` - minimal value of DM window to search [cm^3 / pc].

- ``dm_max`` - maximum value of DM window to search [cm^3 / pc].

- ``d_t`` - width of feature [s] in (t, DM)-space along DM-axis to treat it as
  candidate. Default: 0.005

- ``d_dm`` - width of feature [cm^3/pc] in (t, DM)-space along DM-axis to treat
  it as candidate. Default: 100.

- ``perc`` - percentile of image values that is used to blank image before
 searching for objects. Default: 99.5

- ``savefig_dyn`` - file name for saving picture of dynamical spectra.

- ``savefig_dedm`` - file name for saving picture of de-dispersed frequency
  averaged dynamical spectra.

- ``save_result`` - file name to save (t, DM)-coordinates of found candidates.
  [s, cm^3/pc]

Notes
^^^^^

Algorithm search for extended regions in image of de-dispersed frequency
averaged dynamical spectra (that is (t, DM)-plane). There are 3 tunable
parameters:

- ``perc``. Current experience suggest values ``99.9 - 99.95``. Low value could
    bring many false features in (t, DM)-space, but as long as we can compare
    results using different telescopes it doesn't seem to be an issue. Value
    that is too high can split characteristic x-shaped FRB-features.

- ``d_t`` & ``d_dm``. Some experience with fake FRB injected in real data have
    showen that 2 most informative features in classification of FRB candidates
    on (t, DM)-plane are their widths in t- and DM-directions. Currently it is
    the only method of classification that is implemented. Nonetheless one can
    use any other features and their own algorithms by overriding/extending
    ``TDMImageObjects._classify`` method that gets ``image`` and ``labelled
    array`` as first two positional arguments.



License
-------

Copyright 2015 Ilya Pashchenko.

frb is free software made available under the MIT License. For details see the
LICENSE file.
