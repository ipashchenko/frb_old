uvmod
=====

**Searching FRB in RA data**

frb - tool used to search FRB in Radioastron data

Documentation
-------------

Requirements:
^^^^^^^^^^^^^
numpy, pylab (for plots)

Searching for dispersioned pulse in txt-format data file:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``user@host:~$ python search_file.txt data.txt -nu_max NU_MAX -dnu DNU -dt DT
-dm_min DM_MIN -dm_max DM_MAX -savefig fig.png``

Parameters:

- ``-nu_max`` - frequency of highest frequency channel [MHz].

- ``-dnu`` - frequency with of single frequency channel [MHz].

- ``-dt`` - time step (resolution) [s].

- ``dm_min`` - minimal value of DM window to search [cm^3 / pc].

- ``dm_max`` - maximum value of DM window to search [cm^3 / pc].

- ``savefig`` - file name for saving picture of result.

- Currently only grid searches specified range of D.

Notes:

- Currently, ``search_file.py`` only finds DM value that gives maximum SNR in
  frequency average data (after de-dispersion on that value).

License
-------

Copyright 2015 Ilya Pashchenko.

frb is free software made available under the MIT License. For details see the
LICENSE file.
