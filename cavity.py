from pylab import *
from scipy.constants import c, hbar

from input_output import InOutSpec


class Cavity(InOutSpec):
    """
    Class for a cavity
    """
    def init_specific(self,
                 lambda_nm=1064.,
                 losses=1e-6,
                 transmission_input=30e-6,
                 transmission_output=30e-6,
                 length=10e-3,
                 delta_hz=1e6,
                 detection_losses=0.,
                 i_incident=1e-3):
        ##==============================================
        ## BASE PROPERTIES MANUAL INITIALIZATION
        ##==============================================
        self._i_incident = i_incident
        self._lambda_nm = lambda_nm
        self._delta = delta_hz*(2*pi)
        self._length = length
        self._detection_losses = detection_losses
        self._losses = losses
        self._transmission_input = transmission_input
        self._transmission_output = transmission_output

        ##================================================
        ## PARENT CLASS INITIALIZATION
        ##================================================
        #super(Cavity, self).__init__(freq=freq)

        ##================================================
        ##   INTRA FIELDS DEFINITIONS
        ##================================================      
        self.add_intra_field('intra', index_x=0, index_y=1)

        ##================================================
        ##   PORTS DEFINITIONS
        ##================================================

        self.add_in_out_port('right',
                             index_x=0,
                             index_y=1,
                             coupling_hz=self.transmission_output/(2*pi*self.round_trip_time))
        self.add_in_out_port('left',
                             index_x=0,
                             index_y=1,
                             coupling_hz=self.transmission_input/(2*pi*self.round_trip_time),
                             mean_field_input=1.)
        self.add_in_out_port('lost',
                             index_x=0,
                             index_y=1,
                             coupling_hz=self.losses/(2*pi*self.round_trip_time))

        #self.update()
##================================================
##   EVOLUTION MATRIX
##================================================
    def get_evolution_mat(self):
        return array([[-self.kappa/2., self.delta],
                      [-self.delta, -self.kappa/2.]])

##================================================
##UPDATE FUNCTION
##================================================
    def update(self):
        self.left.coupling_hz = self.transmission_input/(2*pi*self.round_trip_time)
        self.right.coupling_hz = self.transmission_output/(2*pi*self.round_trip_time)
        self.lost.coupling_hz = self.losses/(2*pi*self.round_trip_time)

        for port in self.left, self.right, self.lost:
            port.output.losses = self.detection_losses
        super(Cavity, self).update()
##================================================
##  BASE PROPERTIES REQUIRING AN UPDATE :
##    lambda_nm, length, detection_losses,
##    delta, i_incident,transmission_input,
##    transmission_output, losses
##================================================
    @property
    def lambda_nm(self):
        return self._lambda_nm

    @lambda_nm.setter
    def lambda_nm(self, val):
        self._lambda_nm = val
        self.update()
        return val
    
    @property
    def length(self):
        return self._length
    
    @length.setter
    def length(self, val):
        self._length = val
        self.update()
        return val
    
    @property
    def detection_losses(self):
        return self._detection_losses

    @detection_losses.setter
    def detection_losses(self, val):
        self._detection_losses = val
        self.update()
        return val
    
    @property
    def delta(self):
        return self._delta

    @delta.setter
    def delta(self, val):
        self._delta = val
        self.update()
        return val
    
    @property
    def i_incident(self):
        return self._i_incident

    @i_incident.setter
    def i_incident(self, val):
        self._i_incident = val
        self.update()
        return val

    
    @property
    def transmission_input(self):
        return self._transmission_input

    @transmission_input.setter
    def transmission_input(self, val):
        self._transmission_input = val
        self.update()
        return val

    @property
    def transmission_output(self):
        return self._transmission_output

    @transmission_output.setter
    def transmission_output(self, val):
        self._transmission_output = val
        self.update()
        return val

    @property
    def losses(self):
        return self._losses

    @losses.setter
    def losses(self, val):
        self._losses = val
        self.update()
        return val

##================================================
## SECONDARY PROPERTIES USING BASE PROPERTIES
##    (update() is called implicitly)
##================================================
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
    def delta_hz(self):
        return self.delta/(2*pi)

    @delta_hz.setter
    def delta_hz(self, val):
        self.delta = 2*pi*val
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
        return val

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
    def n_photons(self):
        """
        Number of intracavity photons, setting a value affects the incident power
        """
        return self.i_intra*self.round_trip_time/(hbar*self.omega_c)

    @n_photons.setter
    def n_photons(self, val):
        self.i_intra = val*hbar*self.omega_c/self.round_trip_time
        return val
    
##================================================
##  READ-ONLY PROPERTIES
##================================================

    @property
    def finesse(self):
        return 2*pi/(self.losses + self.transmission_input + self.transmission_output)


    @property
    def kappa(self):
        """
        Cavity linewidth
        """

        return (self.losses + self.transmission_input + self.transmission_output)/self.round_trip_time
        #1./self.tau_input + 1./self.tau_loss + 1./self.tau_output

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
