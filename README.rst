FRB
=====

**Searching FRB in RA data**

frb - tool used to search FRB in Radioastron data

Documentation
-------------

Requirements:
^^^^^^^^^^^^^
numpy, matplotlib.pyplot (for plots)

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

- ``savefig_dyn`` - file name for saving picture of dynamical spectra.

- ``savefig_dedm`` - file name for saving picture of de-dispersed frequency
  averaged dynamical spectra.

License
-------

Copyright 2015 Ilya Pashchenko.

frb is free software made available under the MIT License. For details see the
LICENSE file.
