from pylab import *
import os.path as osp

from omr import OmR

current_dir = osp.join(osp.dirname(__file__))
if current_dir=='':
    current_dir = '.'


#### PRX Regal
o = OmR() #\epsilon_p*\epsilon_d
o.set_params(m=6.75e-12,
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

o.gamma_m_hz = 2560.
tau = 1./(2*pi*1.7e6)
o.tau_input = tau/0.4
o.tau_output = tau/0.6
o.n_photons = 1.1e8
o.freq = linspace(1.5e6, 1.56e6, 1000)
#o.right.output._get_losses = lambda:



if True:
    figname = 'PRX Regal fig1'
    close(figname)
    figure(figname, figsize=(15,10))
    dets = (3e3, 6e3, 13e3, 20e3)


    for index, delta_hz in enumerate(dets):
        subplot(len(dets),2,2*index+1)
        o.delta_hz = delta_hz
        o.n_photons = 1.1e8
        title("$\Delta/2\pi$ = " + str(delta_hz) + ' hz')
        ylabel('Intensity noise (vacua)')
        plot(o.freq, abs(o.right.output.spectrum_sym(0)))
        ylim(0.6,1.2)


    xlabel('Frequency (Hz)')

    subplot(122)
    dets = linspace(0, 3000e3, 10)
    opts = []
    opts_all_quad = []
    for det in dets:
        o.delta_hz = det
        opts.append(abs(o.right.output.spectrum_sym(0.)).min())
        opts_all_quad.append(o.right.output.squeezing.min())
    plot(dets, opts,'o', label='Intensity squeezing')
    plot(dets, opts_all_quad,'o', label='Any quadrature')
    ylabel('Min. Intensity noise (vacua)')
    xlabel('$\Delta/2\pi$ (Hz)')
    legend(loc='best')
    show()
    savefig(current_dir +  '/figures/regal_fig1.png')
    savefig(current_dir +  '/figures/regal_fig1.pdf')

if True:
    o.freq=linspace(1.45e6, 1.65e6, 300)

    for label in ('trans', 'refl'):
        title(label)
    for index, delta_hz in enumerate((10, 42e3, 100e3)):
        
        o.delta_hz = delta_hz
        o.n_photons = 1.1e8

        sq = array([o.right.output.spectrum_sym(angle) for angle in linspace(-135*pi/180, 45*pi/180,100)])
        sqr = array([o.left.output.spectrum_sym(angle) for angle in linspace(-135*pi/180, 45*pi/180,100)])
        
        #db = db.clip(-100, 25)
        extent = [-135, 45, 1.45, 1.65]
        for to_plot, label in ((sq, 'trans'), (sqr, 'refl')):
            figure('PRX Regal fig2 ' + label, figsize=(12,8))
            subplot(3,1,index+1)
            db = 10*log10(abs(to_plot).T)
            imshow(db, aspect='auto', origin='lower', extent=extent)
            cb = colorbar()
            contour(db, levels=[0.], linewidths=3, colors='w', extent=extent)
            ylabel('$\Delta/2 \pi$')
            cb.set_label('Noise (dB)')
    for to_plot, label in ((sq, 'trans'), (sqr, 'refl')):
        xlabel('Homodyne phase $^\circ$')
    show()

    savefig(current_dir +  '/figures/regal_fig2.png')
    savefig(current_dir +  '/figures/regal_fig2.pdf')
