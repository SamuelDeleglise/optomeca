from pylab import *

from input_output import InOutSpec

kappa = 1e6*2*pi
delta = 10*2*pi
gamma_m = 10*2*pi
omega_m = 1e6*2*pi
n_th = 30.5
g = 2*pi*2000

class DummySyst(InOutSpec):
    def get_evolution_mat(self):
        return array([[-kappa/2.,0,delta,0],
           [0,-gamma_m/2.,0,omega_m],
           [-delta,-2*g, -kappa/2.,0],
           [-2*g, -omega_m,0,-gamma_m/2.]])
    def init_specific(self):
        self.add_in_out_port('left', coupling_hz=kappa/3, index_x=0, index_y=2)
        self.add_in_out_port('right', coupling_hz=kappa/3, index_x=0, index_y=2)
        self.add_in_out_port('losses', coupling_hz=kappa/3, index_x=0, index_y=2)
        self.add_in_out_port('mech', coupling_hz=gamma_m, index_x=1, index_y=3)
        self.mech.input.set_thermal_state(n_th)

sp = DummySyst()
sp.freq = linspace(0.9999e6, 1.0001e6, 400)


sp.add_intra_field('mech_intra', index_x=1, index_y=3)


plot(sp.freq, sp.mech_intra.spectrum_sym(pi/3), label='mechanics')
plot(sp.freq, sp.left.output.spectrum_sym(pi/2), label='reflected phase')
plot(sp.freq, sp.left.output.spectrum_sym(0), label='reflected intensity')

legend()
show()
