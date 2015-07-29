from pylab import *

from optomeca.cavity import Cavity

cav = Cavity(transmission_input=30e-6,
             transmission_output=30e-6,
             losses=0.1e-6,
             delta_hz=0,
             freq=linspace(0,3e6, 1000))

cav.left.input.squeezing = 10
plot(cav.freq, cav.right.output.spectrum_sym(0), label='transmitted')
plot(cav.freq, cav.left.output.spectrum_sym(0), label='reflected')
legend()
show()
