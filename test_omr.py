from pylab import *

from omr import OmR

o = OmR()
figure()
plot(o.freq, o.left.output.spectrum_sym(pi/2), label='reflected phase')
legend()
show()
