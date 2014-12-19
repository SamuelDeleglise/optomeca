from pylab import *
import os.path as osp

from omr import OmR
import curve
from pandas import Series

current_dir = osp.join(osp.dirname(__file__))
if current_dir=='':
    current_dir = '.'
### Arcizet 2006
arcizet = OmR()
arcizet.set_params(m=190e-9,
        f_mech=814e3,
        q=1.e4,
        lambda_nm = 1064.,
        losses=1.e-6,
        transmission_input=2*pi/3e4,
        transmission_output=1e-6,
        length=0.0024,
        i_incident=2e-4,
        temp=300,    
        delta_hz=1.)
arcizet.freq=linspace(813.e3, 816.e3, 1000)

kappa_hz = arcizet.kappa_hz

name = "Arcizet2006"
figure(name)
for power in array((1.4,2.5,4.5,6.5,9.5)):
    shifts = []
    cooling = []
    phis = linspace(-4,4,100)
    arcizet.delta_hz = 0
    arcizet.i_intra = power
    for index, ph in enumerate(phis): #(0.03, 0.06, 0.09, 0.11, 0.13)):#
        arcizet.delta_hz = ph*kappa_hz/2
        y = arcizet.right.output.spectrum_sym(pi/2)#dummy_spec()
        c = curve.Curve()
        c.set_data(Series(abs(y), index=arcizet.freq))
        f,a = c.fit('lorentz')
        print f.params
        shifts.append(f.params['x0'])#arcizet.freq[y.argmax()])
        cooling.append(2*f.params['bandwidth'])#y.max())
    #gca().set_yscale('log')
    subplot(211)
    plot(phis, array(shifts) - arcizet.omega_m/(2*pi),label=str(power) + ' W')
    xlabel('phi')
    ylabel('freq shift (Hz)')
    subplot(212)
    plot(phis, array(cooling), label=str(power) + ' W')
xlabel('phi')
ylabel('Total effective damping (Hz)')
legend()
show()

savefig(current_dir + '/figures/' + name + '.pdf')
savefig(current_dir + '/figures/' + name + '.png')
