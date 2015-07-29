from pylab import *
import os.path as osp
from omr import OmR
from matplotlib import rcParams
current_dir = osp.join(osp.dirname(__file__))
if current_dir=='':
    current_dir = '.'


I_INTRA_FIG_DROITE = 480
TRANSMISSION_FIG_DROITE = 120

pilier = OmR()#0.304)
pilier.set_params(m=50.e-9,
        f_mech=4.e6,
        q=1e6,
        lambda_nm = 1064.,
        losses=30.e-6,
        transmission_input=30e-6,
        transmission_output=1.e-6,
        length=150e-6,
        i_incident=1000e-6,
        temp=.5,
        delta_hz=10,
        detection_losses=0.0)

pilier.freq = linspace(3.99e6,4.01e6, 1000)


figure('figure PF', figsize=(14, 6))

input_transmissions = linspace(1e-6,150e-6, 100)
intra_powers = linspace(0.05, 500, 100)

pilier.freq = array([4e6])
pilier.delta_hz = 0

input_transmissions = linspace(1e-6,150e-6, 100)
intra_powers = linspace(0.05, 500, 100)

rcParams['image.cmap'] = 'jet_r'

img = []
for p in intra_powers:
    pilier.i_intra = p
    line = []
    for t in input_transmissions:
        pilier.transmission_input = t
        line.append(pilier.left.output.squeezing.min())
        #print "argmin", pilier.left.output.squeezing.argmin()
    img.append(line)
img = -10*log10(abs(array(img)))
extent = [input_transmissions.min()*1e6, input_transmissions.max()*1e6,intra_powers.min(), intra_powers.max()]



subplot(121)
ax = gca()




imshow(img,
       origin='lower',
       extent=extent,
       aspect='auto')
cb = colorbar()
cb.set_label('Optimal squeezing (dB)')
cb.set_ticks([0,1,2,3,4,5])


label = ax.contour(img,
                   origin='lower',
                   levels=[5,4,3,2,1],
                   extent=extent,
                   aspect='auto')
xlabel("Coupler transmission (ppm)")
ylabel("Intracavity power (W)")

ax.clabel(label)

plot([TRANSMISSION_FIG_DROITE],[I_INTRA_FIG_DROITE], '+w', markersize=18, mew=2)



rcParams['image.cmap'] = 'jet'
subplot(122)
pilier.i_intra = I_INTRA_FIG_DROITE
pilier.transmission_input = TRANSMISSION_FIG_DROITE*1e-6

pilier.freq = linspace(3.9e6, 4.1e6, 1e3)
angle_min, angle_max = -90., 90.
sq = array([pilier.left.output.spectrum_sym(angle) for angle in linspace(angle_min*pi/180,
                                                                         angle_max*pi/180,100)])
extent = [angle_min, angle_max, pilier.freq.min()/1e6, pilier.freq.max()/1e6]
pilier.delta_hz = 0
db = 10*log10(abs(sq).T)
imshow(db, aspect='auto', origin='lower', extent=extent)
xlabel('Quadrature angle ($^\circ$)')
ylabel('Frequency (MHz)')
cb = colorbar()
contour(db,
        levels=[0.],
        linewidths=3, colors='w', extent=extent, origin='lower')
cb.set_label('Noise level (dB)')
cb.set_ticks([-5,0,5,10,15,20,25,30,35,40,45,50,55,60])
subplots_adjust(wspace=0.3)
savefig(current_dir +'/figures/figure_pf_wide_band.pdf')
