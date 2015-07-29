from pylab import *
import os.path as osp
from matplotlib.widgets import Slider, Button, RadioButtons

from omr import OmR
current_dir = osp.join(osp.dirname(__file__))
if current_dir=='':
    current_dir = '.'

#close("pilier")
#figure("pilier")
ax=gca()
N_DET = 3
ax.set_color_cycle([cm.jet(k) for k in linspace(0,1,N_DET)])

pilier = OmR()#0.304)
pilier.set_params(m=25.e-9,
        f_mech=4.e6,
        q=2e6,
        lambda_nm = 1064.,
        losses=15.e-6,
        transmission_input=30e-6,
        transmission_output=1.e-6,
        length=150e-6,
        i_incident=1000e-6,
        temp=1.,
        delta_hz=10,
        detection_losses=0.0)
pilier.freq = linspace(3.99e6,4.01e6, 1000)

"""
DETS = linspace(0,200e3, N_DET)
extent = [-90, 90, min(pilier.freq), max(pilier.freq)]
for index, delta_hz in enumerate(DETS):
    subplot(N_DET,1,index+1)  
    pilier.delta_hz = delta_hz
    img = array([10*log10(abs(pilier.left.output.spectrum_sym(phi))) for phi in linspace(-pi/2,pi/2)])
    img = img.T
    imshow(img, aspect='auto',origin='lower', extent=extent)
    cb = colorbar()
    contour(img, levels=[0.], linewidths=3, colors='w', extent=extent)
    ylabel('$\Delta/2 \pi$')
    cb.set_label('Noise (dB)')
#plot(pilier.freq, len(pilier.freq)*[1.], ':k')
show()
"""


figure('Optimize input coupler', figsize=(14, 10))
ax = gca()
input_transmissions = linspace(1e-6,150e-6, 100)
intra_powers = linspace(0.05, 500, 100)


#pilier.freq = linspace(3.99e6, 4.01e6, 100)
pilier.freq = array([4e6])

pilier.delta_hz = 0


plt.subplots_adjust(left=0.05, bottom=0.4)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axamp  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
axtemp  = plt.axes([0.25, 0.2, 0.65, 0.03], axisbg=axcolor)
axlength  = plt.axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)

temp = Slider(axtemp,
              'Temp.',
              0.,
              10.,
              valinit=1.,
              valfmt='%1.1f K')

slosses = Slider(axfreq,
                 'losses (ppm)',
                 1.,
                 100.0,
                 valinit=30.,
                 valfmt='%1.0f ppm')
sdl = Slider(axamp,
             'det. losses %',
             0.0,
             100.0,
             valinit=0.,
             valfmt='%1.0f')
length = Slider(axlength,
             'length',
             0.0,
             1000.0,
             valinit=150.,
             valfmt='%1.f $\mum')

for i in (temp, slosses, sdl, length):
    i.on_changed(lambda event:button.label.set_text('calculate'))


calcax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(calcax, 'save', color=axcolor, hovercolor='0.975')
#save_ax = plt.axes([0.4, 0.025, 0.1, 0.04])
#button_save = Button(save_ax, 'Save', color=axcolor, hovercolor='0.975')


plt.axes(ax)



save_id = 0
calc_id = 0
def save():
    pname = "length=" + str(pilier.length)+ "opt_losses="+str(pilier.losses)+", det_loss=" + str(pilier.detection_losses) + ", temp=" + str(pilier.temp)
    savefig(current_dir +  '/figures/pilier/' + pname + '.png')
    savefig(current_dir +  '/figures/pilier/' + pname + '.pdf')

def calc():
    ax.clear()
    pilier.detection_losses = float(sdl.valtext.get_text().rstrip(''))*1e-2
    pilier.losses = float(slosses.valtext.get_text().rstrip(' ppm'))*1e-6
    pilier.temp = float(temp.valtext.get_text().rstrip(' K'))
    pilier.length = float(length.valtext.get_text().rstrip(' $\mu$m'))*1e-6
    img = []
    extent = [input_transmissions.min()*1e6, input_transmissions.max()*1e6,intra_powers.min(), intra_powers.max()]
    for p in intra_powers:
        pilier.i_intra = p
        line = []
        for t in input_transmissions:
            pilier.transmission_input = t
            line.append(pilier.left.output.squeezing.min())
            #print "argmin", pilier.left.output.squeezing.argmin()
        img.append(line)
    img = abs(array(img))
    label = ax.contour(img, origin='lower', extent=extent)
    xlabel("Coupler transmission (ppm)")
    ylabel("Intracavity Power (W)")
    #title("Optimum squeezing, " + pname)
    ax.clabel(label)
    show()
#imshow()
## http://matplotlib.org/examples/widgets/slider_demo.html

def clicked(event):
    if button.label.get_text()=='calculate':
        calc()
        button.label.set_text('save')
    else:
        if button.label.get_text()=='save':
            save()
            button.label.set_text('calculate')

calc_id = button.on_clicked(clicked)
calc()
