from pylab import *
from scipy.constants import k,c, hbar

from input_output import InOutSpec, Field
from cavity import Cavity

class OmR(Cavity):
    """
    Optomechanical Resonator: base class to calculate covariance spectra
    """
    
    def init_specific(self,
                 #CAVITY PARAMETERS
                 lambda_nm = 1064., 
                 losses=30e-6,
                 transmission_input=30e-6,
                 transmission_output=30e-6,
                 length=100e-6,
                 delta_hz=1e6,
                 i_incident=1e-3,
                 detection_losses=0.,
                 #MECHANICAL PARAMETERS
                 temp=4., 
                 m=25e-9,
                 f_mech=4e6,
                 q=1e6,
                 #FERQ PARAMETERS
                 freq_start=3.8e6,
                 freq_stop=4.2e6,
                 n_freqs=1000,
                 ):
        ##==============================================
        ## BASE PROPERTIES MANUAL INITIALIZATION
        ##==============================================
        self._m = m
        self._f_mech = f_mech
        self._q = q
        self._temp = temp
        
        ##==============================================
        ## PARENT CLASS INITIALIZATION
        ##==============================================
        super(OmR, self).init_specific(lambda_nm=lambda_nm,
                                  length=length,
                                  i_incident=i_incident,
                                  transmission_input=transmission_input,
                                  transmission_output=transmission_output,
                                  losses=losses,
                                  delta_hz=delta_hz,
                                  detection_losses=detection_losses)
        

        ##================================================
        ##   INTRA FIELDS DEFINITIONS
        ##================================================
        self.add_intra_field('mech', index_x=2, index_y=3)

        ##================================================
        ##   PORTS DEFINITIONS
        ##================================================
        self.add_in_out_port('mech_bath',
                            index_x=2,
                            index_y=3,
                            coupling_hz=self.gamma_m_hz)#1./(self.gamma_m_hz*2*pi))

        


##================================================
##   EVOLUTION MATRIX
##================================================
    def get_evolution_mat(self):
        im2g = 2*self.g*sin(self.intra.mean_field_angle)
        re2g = 2*self.g*cos(self.intra.mean_field_angle)                     
        return array([[-self.kappa/2.,self.delta,im2g,0],
                      [-self.delta,-self.kappa/2.,-re2g,0],
                      [0,0,-self.gamma_m/2.,self.omega_m],
                      [-re2g,-im2g,-self.omega_m,-self.gamma_m/2.]])
               #array([[-self.kappa/2.,self.delta,0,0],
               #       [-self.delta,-self.kappa/2.,-2*self.g,0],
               #       [0,0,-self.gamma_m/2.,self.omega_m],
               #       [-2*self.g,0,-self.omega_m,-self.gamma_m/2.]])
               #array([[-self.kappa/2.,0,self.delta,0],
               #    [0,-self.gamma_m/2.,0,self.omega_m],
               #    [-self.delta,-2*self.g, -self.kappa/2.,0],
               #    [-2*self.g, -self.omega_m,0,-self.gamma_m/2.]])
    def get_evolution_mat_dc(self):
        """
        No DC radiation pressure... detuning becomes "effective detuning".
        """
        return array([[-self.kappa/2.,self.delta,0,0],
                      [-self.delta,-self.kappa/2.,0,0],
                      [0,0,-self.gamma_m/2.,self.omega_m],
                      [0,0,-self.omega_m,-self.gamma_m/2.]])

        
##================================================
##UPDATE FUNCTION
##================================================
    def update(self):
        self.mech_bath.coupling_hz = self.gamma_m_hz
        self.mech_bath.input.set_thermal_state(self.n_th)
        super(OmR, self).update()


##================================================
##  BASE PROPERTIES REQUIRING AN UPDATE :
##    lambda_nm, length, detection_losses,
##    delta, i_incident,transmission_input,
##    transmission_output, losses
##================================================
    @property
    def temp(self):
        """
        Thermal occupancy of the environnement.
        Setting a value affects the temperature
        """
        return self._temp

    @temp.setter
    def temp(self, val):
        """
        Takes care of also changing the bath property
        """
        
        self._temp = val
        self.update()
        return val

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, val):
        self._m = val
        self.update()
        return val

    @property
    def q(self):
        return self._q

    @q.setter
    def q(self, val):
        self._q = val
        self.update()

    @property
    def f_mech(self):
        return self._f_mech

    @f_mech.setter
    def f_mech(self, val):
        self._f_mech = val
        self.update()
        return val
##================================================
## SECONDARY PROPERTIES USING BASE PROPERTIES
##    (update() is called implicitly)
##================================================
    ## MECHANICAL
    @property
    def gamma_m_hz(self):
        """
        Mechanical damping. Setting a value affects q
        """
        
        return self.f_mech/self.q

    @gamma_m_hz.setter
    def gamma_m_hz(self, val):
        self.q = self.f_mech/val    
        return val

    @property
    def n_th(self):
        return self.temp*k/(hbar*self.omega_m)

    @n_th.setter
    def n_th(self, val):
        self.temp = val*(hbar*self.omega_m)/k
        return val

    ## OPTOMECHANICAL
    @property
    def rsf(self):
        """
        Resolved Sideband Factor, setting a value affects the cavity length 
        """
        return self.f_mech/self.bandwidth_hz

    @rsf.setter
    def rsf(self, val):
        self.bandwidth_hz = self.f_mech/val
        return val

    @property
    def cooperativity(self):
        """
        Multiphoton cooperativity. Setting a value affects the incident power.
        """
        return 8.*self.finesse*self.i_intra/(self.lambda_nm*1e-9*c*self.gamma_m*self.omega_m*self.m)

    @cooperativity.setter
    def cooperativity(self, val):
        self.i_intra = val*self.gamma_m*self.omega_m*self.lambda_nm*1e-9*c*self.m/(8*self.finesse)

##================================================
##  READ-ONLY PROPERTIES
##================================================
    ## MECHANICAL
    @property
    def omega_m(self):
        return 2*pi*self.f_mech

    @property
    def gamma_m(self):
        """
        Mechanical damping angular freq. Setting a value affects q
        """
        return self.omega_m/self.q


    ##OPTOMECHANICAL
    @property
    def xzpf(self):
        """
        zero-point fluctuations of mechanical mode
        """
        return sqrt(hbar/(2*self.m*self.omega_m))

    @property
    def g0(self):
        """
        Vacuum optomechanical coupling rate (angular frequency)
        """
        return self.xzpf*self.omega_c/self.length

    @property
    def g0_hz(self):
        """
        Vacuum optomechanical coupling rate (hertz)
        """
        return 1./(2*pi)*self.xzpf*self.omega_c/self.length

    @property
    def g(self):
        """
        multiphoton coupling rate.
        """
        return self.g0*sqrt(self.n_photons)
    
##================================================
##         USEFUL FUNCTIONS                          
##================================================
    def show_optical_spectrum(self):
        figure('Optical spectrum')
        upper_sideband = self.delta_hz + self.f_mech
        lower_sideband = self.delta_hz - self.f_mech
        bound = max(3*self.kappa_hz, self.f_mech*1.5)
        x = linspace(min(-3*self.kappa_hz, lower_sideband*1.2), max(3*self.kappa_hz, upper_sideband*1.2), 1000)
        plot(x, 1./(1 + (2*x/self.kappa_hz)**2), '--g')
        plot([self.delta_hz, self.delta_hz],[-0.1,1.1], '-k', label='laser frequency')
        xlabel('$\omega_{opt}/2 \pi$ (Hz)')
        ylabel('Optical mode density (arb.)')
        plot([upper_sideband, upper_sideband], [-0.1, 1.1], '-b', label='upper mech. sideband')
        
        plot([lower_sideband, lower_sideband], [-0.1, 1.1], '-r', label='lower mech. sideband')
        plot([x.min(), x.max()],[0,0], ':g')
        plot([x.min(), x.max()],[1,1], ':g')
        ylim((-0.1, 1.1))
        xlim([x.min(), x.max()])
        legend()
        show()

