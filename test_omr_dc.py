from pylab import *
import os.path as osp

from omr import OmR
current_dir = osp.join(osp.dirname(__file__))
if current_dir=='':
    current_dir = '.'

cav = OmR() #\epsilon_p*\epsilon_d
cav.set_params(m=6.75e-12,
        f_mech=1.5243e6,
        q=100.,
        lambda_nm = 1064.,
        losses=1.e-6,
        transmission_input=30e-6,
        transmission_output=10.e-6,
        length=0.00764,
        i_incident=1e-3,
        temp=3.8e-4,    
        delta_hz=1.,
        detection_losses=0.304)

cav.gamma_m_hz = 2560.
tau = 1./(2*pi*1.7e6)
cav.tau_input = tau/0.4
cav.tau_output = tau/0.6
cav.n_photons = 1.1e8
cav.freq = linspace(1.5e6, 1.56e6, 1000)
figure('mean_fields_omr', figsize=(12.9, 10.9))
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
    deltas = linspace(-25*cav.bandwidth_hz, 25*cav.bandwidth_hz, 300)
    for delta_hz in deltas:
        cav.delta_hz = delta_hz
        reflected.append(cav.left.output.mean_field)#/cav.intra.mean_field))
        transmitted.append(cav.right.output.mean_field)#/cav.intra.mean_field))
    reflected = array(reflected)
    transmitted = array(transmitted)

    plot(real(reflected), imag(reflected),'o', label='reflected')
    plot(real(transmitted), imag(transmitted),'o', label='transmitted')
    xlabel('Real')
    ylabel('Imag')
    legend()
    gca().set_aspect('equal')
    show()




#savefig(current_dir +  '/figures/regal_fig2.png')
