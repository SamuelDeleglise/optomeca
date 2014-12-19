from pylab import *
import os.path as osp

from omr import OmR


current_dir = osp.join(osp.dirname(__file__))
if current_dir=='':
    current_dir = '.'
N_TEMPERATURES = 10

o = OmR()
o.set_params(m=6.75e-12,
        f_mech=1e6,
        q=1000.,
        lambda_nm = 1064.,
        losses=1.e-6,
        transmission_input=300e-6,
        length=0.00764,
        i_incident=1e-3,    
        delta_hz=1.)
tau = 1./(2*pi*1.7e6)
o.tau_in = tau/0.4
o.tau_ex = tau/0.6
o.n_th = 0
o.freq = linspace(0.99e6, 1.01e6, 400)
o.delta_hz = 0





figure('SB asymmetry on resonance', figsize=(10,13))

subplot(211)
o.n_th = 0
o.cooperativity = 0.1
title('Cooperativity = 0.1, n_th = 0')
plot(o.freq, o.right.output.spectrum_ss('upper'), '-r', label='Upper sideband')
plot(o.freq, o.right.output.spectrum_ss('lower'), '-b', label='Lower sideband')
xlabel('Frequency (Hz)')
ylabel('Single sided spectrum (Vacua)')
legend(loc='best')


subplot(212)
ax=gca()
ax.set_color_cycle([cm.bwr(k) for k in linspace(0,1,N_TEMPERATURES)])

coops = logspace(-2, 2, 100)
o.freq = array([1e6])
for n_th in [0.01*2**n for n in range(N_TEMPERATURES)]:
    asyms = []
    for cooperativity in coops:
        o.n_th = n_th
        o.cooperativity = cooperativity
        asym = ((o.right.output.spectrum_ss('upper').max() - 1) - (o.right.output.spectrum_ss('lower').max() - 1))/\
        ((o.right.output.spectrum_ss('upper').max() - 1) + (o.right.output.spectrum_ss('lower').max() - 1))
        asyms.append(asym)
    plot(coops, asyms, label='n_th=' + str(n_th))
gca().set_xscale('log')
xlabel('Cooperativity')
ylabel('Sideband asymmetry')
legend(loc='best')
show()

savefig(current_dir + '/figures/SB asymmetry on resonance.png')
savefig(current_dir + '/figures/SB asymmetry on resonance.pdf')


figure('SB asymmetry Painter style',figsize=(10,13))

subplot(211)
o.rsf = 10
o.n_th = 0
o.cooperativity = 0.01
title('Cooperativity = 0.01, n_th = 0')
o.freq = linspace(o.f_mech - 3000, o.f_mech+3000, 200)
o.delta_hz = o.f_mech
plot(o.freq, o.right.output.spectrum_sym(pi/2), '-b', label='drive on blue side (lower sb. resonant)')
o.delta_hz = -o.f_mech
plot(o.freq, o.right.output.spectrum_sym(pi/2), '-r', label='drive on red side (upper sb. resonant)')
xlabel('Frequency (Hz)')
ylabel('Phase quadrature spectrum (Vacua)')

legend(loc='best')

subplot(212)
ax=gca()
ax.set_color_cycle([cm.bwr(k) for k in linspace(0,1,N_TEMPERATURES)])

o.freq = array([o.f_mech])
o.rsf = 10


n_ths = [0.01*2**n for n in range(N_TEMPERATURES)]
coops = logspace(-3, 0, 100)
for n_th in n_ths:
    asyms = []
    for cooperativity in coops:
        o.n_th = n_th
        o.cooperativity = cooperativity
        o.delta_hz = o.f_mech
        red = o.right.output.spectrum_sym(pi/2).max() - 1.
        o.delta_hz = -o.f_mech
        blue = o.right.output.spectrum_sym(pi/2).max() - 1.
        asyms.append((blue - red)/(red + blue))
    plot(coops, asyms, label='n_th=' + str(n_th))
gca().set_xscale('log')
xlabel('Cooperativity')
ylabel('Sideband asymmetry')
legend(loc='best')
show()
savefig(current_dir +  '/figures/SB_asymmetry_Painter_style.png')
savefig(current_dir +  '/figures/SB_asymmetry_Painter_style.pdf')
