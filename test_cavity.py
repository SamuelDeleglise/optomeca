from pylab import *
import os.path as osp

from cavity import Cavity

current_dir = osp.join(osp.dirname(__file__))
if current_dir=='':
    current_dir = '.'

fig = figure('Squeezing in cavity')
cav = Cavity()
cav.set_params(transmission_input=30e-6,
                transmission_output=20e-6,
                losses=0.1e-6,
                delta_hz=1e-6,
                freq=linspace(0,5e5, 1000))
#cav.transmission_input = 0.00021776668116086753
#cav.transmission_output = 0.0003266500217413013
#cav.losses = 1e-6
cav.left.input.set_squeezed_state(10.)
cav.update()
plot(cav.freq, cav.right.output.spectrum_sym(pi/2), label='transmitted')
plot(cav.freq, cav.left.output.spectrum_sym(pi/2), label='reflected')
legend()
xlabel('Frequency(Hz)')
ylabel('Phase noise (Vacua)')
title('10dB input squeezing')
show()


fig.savefig(current_dir +  '/figures/cavity_squeezing_input.png')
fig.savefig(current_dir +  '/figures/cavity_squeezing_input.pdf')


fig = figure("Phase of fields cavity")

reflected = []
transmitted = []
deltas = linspace(-5*cav.bandwidth_hz, 5*cav.bandwidth_hz, 300)
for delta_hz in deltas:
    cav.delta_hz = delta_hz
    reflected.append(cav.left.output.mean_field)#/cav.intra.mean_field))
    transmitted.append(cav.right.output.mean_field)#/cav.intra.mean_field))
reflected = array(reflected)
transmitted = array(transmitted)
plot(deltas, angle(reflected), label='reflected_phase')
plot(deltas, angle(transmitted), label='transmitted_phase')
xlabel('Frequency (Hz)')
ylabel('Phase(rad.)')
       
legend()
show()

fig = figure('mean_fields', figsize=(12.9, 10.9))
last=None
for index, ((input_trans, out_trans, loss), title_) in enumerate((((1e-6,5e-6, 0), 'undercoupled'),
                                              ((5e-6, 5e-6,0), 'critically coupled'),
                                              ((5e-6, 0.5e-6,0), 'over coupled'),
                                              ((1e-6, 1e-6, 5e-6), 'loss dominated'))):
    last = subplot(2, 2 ,index + 1, sharex=last)
    cav.transmission_input = input_trans
    cav.transmission_output = out_trans
    cav.losses = loss
    title(title_)
    reflected = []
    transmitted = []
    deltas = linspace(-4*cav.bandwidth_hz, 4*cav.bandwidth_hz, 64)
    for delta_hz in deltas:
        cav.delta_hz = delta_hz
        reflected.append(cav.left.output.mean_field)#/cav.intra.mean_field))
        transmitted.append(cav.right.output.mean_field)#/cav.intra.mean_field))
    reflected = array(reflected)
    transmitted = array(transmitted)

    paths = scatter(real(reflected), imag(reflected),marker='s', s=40,linewidths=0, label='reflected', c=[cm.bwr(x) for x in (deltas - deltas.min())/(deltas.max() - deltas.min())])
    paths = scatter(real(transmitted), imag(transmitted),marker='^',s=40,linewidths=0,label='transmitted',c=[cm.bwr(x) for x in (deltas - deltas.min())/(deltas.max() - deltas.min())])
    xlabel('Real')
    ylabel('Imag')
    paths.set_cmap(cm.bwr)
    paths.set_array(deltas/cav.bandwidth_hz)
    cb = colorbar(paths)
    cb.set_label('$\Delta/$bandwidth')

    legend()
    gca().set_aspect('equal')
show()

fig.savefig(current_dir +  '/figures/cavity_mean_field.png')
fig.savefig(current_dir +  '/figures/cavity_mean_field.pdf')
