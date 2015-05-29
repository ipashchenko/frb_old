import os
import sys
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
from frb.frames import SimFrame
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


frame = SimFrame(128, 1000, 1676., 0., 16. / 128., 0.001)
frame.add_pulse(0.3, 1.3, 0.0015, dm=700.)
frame.add_noise(0.5, )
