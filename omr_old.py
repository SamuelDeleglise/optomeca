from pylab import *
from scipy.constants import k,c, hbar

from input_output_old import InOutSpec, Field

class OmR(InOutSpec):
    """
    Optomechanical Resonator: base class to calculate covariance spectra
    """
    
    def __init__(self,
                 m=25e-9,
                 f_mech=4e6,
                 q=1e6,
                 lambda_nm = 1064.,
                 losses=30e-6,
                 transmission_input=30e-6,
                 transmission_output=30e-6,
                 length=100e-6,
                 i_incident=1e-3,
                 temp=4.,
                 delta_hz=1e6,
                 freq_start=3.8e6,
                 freq_stop=4.2e6,
                 n_freqs=1000,
                 detection_losses=0.):
        self.m = m
        self.f_mech = f_mech
        self.q = q
        self.lambda_nm = lambda_nm
        self.losses = losses
        self.transmission_input= transmission_input
        self.transmission_output = transmission_output
        self.length = length
        self.i_incident = i_incident
        self.temp = temp
        self.delta_hz = delta_hz
        self.detection_losses = detection_losses
        super(OmR, self).__init__(freq=linspace(freq_start, freq_stop, n_freqs))
        

                
        self.add_in_out_port('right',
                                index_x=0,
                                index_y=2,
                                get_coupling_sqrthz=self._sqrt_kappa_right,
                                get_output_losses=lambda:self.detection_losses)
        
        
        
        self.add_in_out_port('left',
                                  index_x=0,
                                  index_y=2,
                                  get_coupling_sqrthz=self._sqrt_kappa_left,
                                  get_output_angle=self._reflected_angle,
                                  get_output_losses=lambda:self.detection_losses)


        self.add_in_out_port('bath',
                                index_x=0,
                                index_y=2,
                                get_coupling_sqrthz=self._sqrt_kappa_loss,
                                get_output_losses=lambda:self.detection_losses)
        


        self.add_in_out_port('mech_bath',
                            index_x=1,
                            index_y=3,
                            get_coupling_sqrthz=self._get_sqrt_gamma_m,
                            get_input_cov=self._get_mechanical_cov,
                            get_output_losses=lambda:self.detection_losses)

        self.add_intra_field('mech', index_x=1, index_y=3)



    def _sqrt_kappa_right(self):
            return 1./sqrt(self.tau_output)
    def _sqrt_kappa_left(self):
            return 1./sqrt(self.tau_input)
    def _reflected_angle(self):
        return arctan(2*self.delta/self.kappa)
    def _sqrt_kappa_loss(self):
        return 1./sqrt(self.tau_loss)
    def _get_mechanical_cov(self):
        return array([[2*self.n_th + 1, 0],[0, 2*self.n_th + 1]])
    def _get_sqrt_gamma_m(self):
        return sqrt(self.gamma_m)
########################################################
###         USEFUL FUNCIONS                          
########################################################
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


########################################################
###         OPTICAL PROPERTIES
########################################################
    @property
    def delta(self):
        return 2*pi*self.delta_hz

    @delta.setter
    def delta(self, val):
        self.delta_hz = val/(2*pi)
        return val



    @property
    def bandwidth_hz(self):
        """
        Half of optical linewidth, setting a value affects the cavity length 
        """
        return c/(4*self.length*self.finesse)

    @bandwidth_hz.setter
    def bandwidth_hz(self, val):
        self.length = c/(4*val*self.finesse)

    @property
    def finesse(self):
        return 2*pi/(self.losses + self.transmission_input + self.transmission_output)

    @property
    def i_intra(self):
        """
        Intracavity power (W), setting a value affects the incident power
        """
        return 2./pi*self.i_incident*self.finesse*(self.kappa**2)/(self.kappa**2 + 4*self.delta**2) ## to check

    @i_intra.setter
    def i_intra(self, val):
        self.i_incident = (self.kappa**2 + 4*self.delta**2)/(self.kappa**2)*pi*val/(2*self.finesse) #to check
        return val

    @property
    def tau_input(self):
        """
        Coupling limited lifetime, setting  a value affects the transmission 
        """
        return self.round_trip_time/self.transmission_input

    @tau_input.setter
    def tau_input(self, val):
        self.transmission_input = self.round_trip_time/val
        return val

    @property
    def tau_loss(self):
        """
        Internal loss limited lifetime, setting a value affects the losses
        """
        return self.round_trip_time/self.losses

    @tau_loss.setter
    def tau_loss(self, val):
        self.losses = self.round_trip_time/val
        return val

    @property
    def tau_output(self):
        """
        Internal loss limited lifetime, setting a value affects the losses
        """
        return self.round_trip_time/self.transmission_output

    @tau_output.setter
    def tau_output(self, val):
        self.transmission_output = self.round_trip_time/val

    @property
    def kappa(self):
        """
        Cavity linewidth
        """

        return 1./self.tau_input + 1./self.tau_loss + 1./self.tau_output

    @property
    def kappa_hz(self):
        return self.kappa/(2*pi)

    
    @property
    def omega_c(self):
        """
        Resonant optical angular frequency
        """
        return 2*pi*c/(self.lambda_nm*1.e-9)
    
    @property
    def round_trip_time(self):
        return 2*self.length/c
        
    @property
    def n_photons(self):
        """
        Number of intracavity photons, setting a value affects the incident power
        """
        return self.i_intra*self.round_trip_time/(hbar*self.omega_c)

    @n_photons.setter
    def n_photons(self, val):
        self.i_intra = val*hbar*self.omega_c/self.round_trip_time
        return val
########################################################
###         MECHANICAL PROPERTIES
########################################################

    @property
    def omega_m(self):
        return 2*pi*self.f_mech

    @property
    def gamma_m(self):
        """
        Mechanical damping angular freq. Setting a value affects q
        """
        return self.omega_m/self.q

    @gamma_m.setter
    def gamma_m(self, val):
        self.q = self.omega_m/val    
        return val
    
    @property
    def gamma_m_hz(self):
        """
        Mechanical damping. Setting a value affects q
        """
        return self.gamma_m/(2*pi)

    @gamma_m_hz.setter
    def gamma_m_hz(self, val):
        self.gamma_m = 2*pi*val    
        return val

    @property
    def n_th(self):
        """
        Thermal occupancy of the environnement.
        Setting a value affects the temperature
        """
        return k*self.temp/(hbar*self.omega_m)

    @n_th.setter
    def n_th(self, val):
        self.temp = val*(hbar*self.omega_m)/k
        return val

########################################################
###         OPTO-MECHANICAL PROPERTIES
########################################################
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


    @property
    def cooperativity(self):
        """
        Multiphoton cooperativity. Setting a value affects the incident power.
        """
        return 8.*self.finesse*self.i_intra/(self.lambda_nm*1e-9*c*self.gamma_m*self.omega_m*self.m)


    
    @cooperativity.setter
    def cooperativity(self, val):
        self.i_intra = val*self.gamma_m*self.omega_m*self.lambda_nm*1e-9*c*self.m/(8*self.finesse)



########################################################
###         FUNCTION OF PARENT CLASS
########################################################
    """
    def cov_input(self):
        return array([  [1.,0,0,0,0,0,0,0],
                        [0,1.,0,0,0,0,0,0],
                        [0,0,1.,0,0,0,0,0],
                        [0,0,0,2*self.n_th + 1.,0,0,0,0],
                        [0,0,0,0,1.,0,0,0],
                        [0,0,0,0,0,1.,0,0],
                        [0,0,0,0,0,0,1.,0],
                        [0,0,0,0,0,0,0,2*self.n_th + 1.]])
"""
    
    def evolution_mat(self):
        return array([[-self.kappa/2.,0,self.delta,0],
           [0,-self.gamma_m/2.,0,self.omega_m],
           [-self.delta,-2*self.g, -self.kappa/2.,0],
           [-2*self.g, -self.omega_m,0,-self.gamma_m/2.]])
    
    """ 
    def coupling_mat(self):
        return array([[sqrt(1./self.tau_loss),sqrt(1./self.tau_input),sqrt(1./self.tau_output),0,0,0,0,0],
                    [0,0,0,sqrt(self.gamma_m),0,0,0,0],
                    [0,0,0,0,sqrt(1./self.tau_loss),sqrt(1./self.tau_input),sqrt(1./self.tau_output),0],
                    [0,0,0,0,0,0,0,sqrt(self.gamma_m)]])
    """
