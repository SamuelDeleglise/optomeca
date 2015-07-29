from pylab import *

class InOutSpec(object):
    def __init__(self):
        """
        This class solves in the frequency domain the Langevin problem:
                dX/dt = M.X + N.\\xi
        

        *cov_input is the covariance matrix of input noises:
            input_cov = <\\xi(-\\omega), \\xi\^T(\\omega)>.
        *evolution_mat = M is the Evolution matrix for modes X.
        *coupling_mat = N is the coupling matrix between langevin noises and modes X
          it has coefficients sqrt(kappa) in position i,j
          when the mode i is coupled to (loss) channel j.
          It can be a rectangular matrix if there are additional loss channels

        If an argument is not provided, the corresponding property should be defined
        in the derived class
        """

        self.freq = linspace(0.99e6, 1.01e6, 100)
        self.fields = []
        self.ports = []

        self.evolution_mat_dc = None
        self.current_omega = 0.
        self.init_specific()
        self.update()
        
        

    def set_params(self, **kwds):
        for k, v in kwds.iteritems():
            getattr(self, k) #raises an exception if property doesn't exist
            setattr(self, k, v)
##======================================================
##  TO BE IMPLEMENTED IN BASE CLASS
##======================================================
    def init_specific(self):
        raise NotImplementedError("This function should be implemented in base class")
    def get_evolution_mat(self):
        raise NotImplementedError("This function should be implemented in base class")
    def get_evolution_mat_dc(self):
        """
        Returns the evolution matrice at DC frequency. By default, returns self.evolution_mat().
        However, it can be beneficial to overwrite this default behaviour, for instance to remove
        DC radiation pressure effects and get rid of the autocoherent problem, bistability, and so on.
        """
        return self.get_evolution_mat()
    
##======================================================
##  UPDATE FUNCTION:
##      - TO BE CALLED WHEN PARAMETERS HAVE CHANGED.
##      - REQUIRES AN IMPLEMENTATION OF get_evolution_mat.
##      - get_coupling_mat AND get_cov_input ARE NOT NEEDED
##      SINCE THEY HAVE A DEFAULT IMPLEMENTATION
##      BASED ON PORT DEFINITIONS
##======================================================
    def update(self):
        """
        Function to be called when a parameter is changed
        (excluding current_frequency)
        """

        self.coupling_mat = self.get_coupling_mat()
        self.evolution_mat_dc = self.get_evolution_mat_dc()
        self.create_intermediates()
        self.calculate_mean_fields()
        self.evolution_mat = self.get_evolution_mat()
        self.cov_input = self.get_cov_input()
        
    def get_cov_input(self):
        mat = zeros((2*self.n_ports, 2*self.n_ports))
        for index, port in enumerate(self.ports):
            mat[2*index, 2*index] = port.input.cov[0,0]
            mat[2*index, 2*index + 1] = port.input.cov[0,1]
            mat[2*index + 1, 2*index] = port.input.cov[1,0]
            mat[2*index + 1, 2*index + 1] = port.input.cov[1,1]
        return mat

    def get_coupling_mat(self):
        mat = zeros((self.n_intra, 2*self.n_ports))
        for index, port in enumerate(self.ports):
            mat[port.index_x, 2*index] = sqrt(2*pi*port.coupling_hz)
            mat[port.index_y, 2*index + 1] = sqrt(2*pi*port.coupling_hz)
        return mat
##======================================================
##  ITERATOR OVER FREQUENCIES
##======================================================

    @property
    def current_freq(self):
        """
        Fourier angular frequencies to sample
        """
        return self.current_omega/(2*pi)

    
    def __iter__(self):
        for omega in self.omega:
            self.current_omega = omega
            yield self
##======================================================
##    PORTS/FIELDS DEFINITIONS
##======================================================   
    def add_in_out_port(self, name, **kwds):
        """
        Specify coupling to a channel.
        arguments:
          - name: name of the port (ex. left, right, bath...)
          - index_x: index of the x quadratue in the intracavity evolution matrix
          - index_y: index of the y quadratue in the intracavity evolution matrix
          - get_coupling_sqrthz: a function to calculate the coupling coefficient
          - get_cov: a function to calculate the covariance matrix
        """
        port = InOutPort(self, **kwds)
        self.__setattr__(name, port)
        self.ports.append(port)

    def add_intra_field(self, name, **kwds):
        field = IntraField(self, **kwds)
        self.__setattr__(name, field)
        self.fields.append(field)

    @property
    def n_intra(self):
        """
        number of 'intracavity' modes.
        """

        if self.evolution_mat_dc==None:
            self.evolution_mat_dc = self.get_evolution_mat_dc()
        return len(self.evolution_mat_dc) ## Exists before evolution_mat

    @property
    def n_ports(self):
        """
        Number of propagating modes
        """
        return len(self.ports)
       
##======================================================
##     FREQUENCY PROPERTIES
##======================================================
    @property
    def n_freqs(self):
        """
        Number of frequency steps. can be set.
        """
        
        return len(self.freq)

    @n_freqs.setter
    def n_freqs(self, val):
        min_freq, max_freq = self.freq_lim
        self.freq = linspace(min_freq, max_freq, val)
        return val
        
    @property
    def freq_lim(self):
        """
        A tuple of value with min and max fourier frequencies. These values can be set.
        """
        return min(self.freq), max(self.freq)

    @freq_lim.setter
    def freq_lim(self, val):
        min_freq, max_freq = val
        self.freq = linspace(min_freq, max_freq, self.n_freqs)
        return val

    @property
    def omega(self):
        """
        Fourier angular frequencies to sample
        """
        return 2*pi*self.freq
    
##======================================================
##     RESULTING COVARIANCE MATRICES
##======================================================

    @property
    def cov_intra(self):
        in_to_intra = self.in_to_intra_mat(self.current_omega)
        return in_to_intra.dot(self.cov_input).dot(conjugate(in_to_intra).T)
    
    def cov_out(self):
        in_to_out = self.in_to_out_mat(self.current_omega)
        return in_to_out.dot(self.cov_input).dot(conjugate(in_to_out).T)
    
##======================================================
##     CALCULATION INTERMEDIATES
##======================================================
    ##Once per spectrum
    def calculate_mean_fields(self):
        self.mean_fields_input = zeros(self.n_ports*2)
        for index, port in enumerate(self.ports):
            self.mean_fields_input[port.input.index_x] = real(port.input.mean_field)
            self.mean_fields_input[port.input.index_y] = imag(port.input.mean_field)
        self.left_DC = self.in_to_out_mat_dc()
        self.in_to_intra_DC = self.in_to_intra_mat_dc()

        self.mean_fields_output = self.left_DC.dot(self.mean_fields_input)
        self.mean_fields_intra = self.in_to_intra_DC.dot(self.mean_fields_input)
    
    def create_intermediates(self):
        self._trans_coupling = self.coupling_mat.T
        self._id_intra = identity(self.n_intra)
        self._id_output = identity(self.n_ports*2)
        
    def sub_cov(self, index_x, index_y):
        return self.cov_intra[np.ix_([index_x, index_y], [index_x, index_y])]#arr #ix should be cached

    def sub_cov_out(self, index_x, index_y):
        return self.cov_out()[np.ix_([index_x, index_y], [index_x, index_y])]


    ## For every omega_current
    def in_to_intra_mat(self, omega):
        """
        Matrix that transforms input noises into intra noises
        """
        
        return linalg.inv(1j*omega*self._id_intra - self.evolution_mat).dot(self.coupling_mat)

    def in_to_intra_mat_dc(self):
        """
        Matrix that transforms input noises into intra noises at DC frequency.
        """
        
        return linalg.inv(-self.evolution_mat_dc).dot(self.coupling_mat)   

    def in_to_out_mat(self, omega):
        """
        Matrix that transforms input noises into output noises
        """
        
        return (self._trans_coupling).dot(self.in_to_intra_mat(omega)) - self._id_output
        
    def in_to_out_mat_dc(self):
        """
        Matrix that transforms input noises into intra noises at DC frequency.
        """
        
        return (self._trans_coupling).dot(self.in_to_intra_mat_dc()) - self._id_output


    
class Field(object):
    """
    Object encapsulating 2-quadratures and a mean-field. 
    """

    def __init__(self, parent_inout_spec, index_x, index_y):
        self.parent = parent_inout_spec
        self.index_x = index_x
        self.index_y = index_y

    def __iter__(self):
        for parent in self.parent:
            yield self # to check

    @property
    def mean_field_angle(self):
        return angle(self.mean_field)
    
    def rotation_mat(self, angle):
        return array([[np.cos(angle), np.sin(angle)],[-np.sin(angle), np.cos(angle)]])


    @property
    def cov_rotated(self):
        """
        Covariance matrix in the frame of the mean field.
        """

        if self.mean_field_angle!=0:
            return self.do_rotate(self.cov_unrotated, self.mean_field_angle)
        return self.cov_unrotated
    
    def do_rotate(self, matrice, angle):
        rot = self.rotation_mat(angle)
        rott = rot.T
        return array(rot.dot(matrice).dot(rott))


    def spectrum_sym(self, quad_angle):
        """
        Returns the symmetrized spectrum along the quadrature quad_angle:
          - 0 is for intensity quadrature
          - pi/2 is for phase quadrature
        """
        
        rot = self.rotation_mat(-quad_angle)
        return array([rot.dot(obj.cov).dot(rot.T)[0,0] for obj in self])

    def spectrum_ss(self, side='upper'):
        """
        Returns the single-sided spectrum of the field
          - 'upper' is for omega>0 (blue sideband)
          - 'lower' is for omega<0 (red sideband)
        """

        if side=='upper':
            return array([0.5*(obj.cov_rotated[0,0] + obj.cov_rotated[1,1] + \
                               1j*(obj.cov_rotated[1,0] - obj.cov_rotated[0,1])) for obj in self])
        if side=='lower':
            return array([0.5*(obj.cov_rotated[0,0] + obj.cov_rotated[1,1] - \
                               1j*(obj.cov_rotated[1,0] - obj.cov_rotated[0,1])) for obj in self])
        raise(ValueError('side should be either "upper" (omega>0) or "lower" (omega<0)'))



    @property
    def squeezing(self):
        """
        Returns the lowest noise level for any quadrature (by diagonalizing the covariance matrix)
        """

        return array([min(linalg.eigvals(obj.cov)) for obj in self])


    #    @property
    #    def squeezing_angle(self):
    #        """
    #        Returns the lowest noise level for any quadrature (by diagonalizing the covariance matrix)
    #        """
    #
    #        vals, vects = linalg.eig(cav.left.output.cov)
    #        min_val = vals.min()
    #        vects[][1]
    #        return arctan2(real(vects[][1]),real(vects[vals.min()][0]))

    @property
    def anti_squeezing(self):
        """
        Returns the highest noise level for any quadrature (by diagonalizing the covariance matrix)
        """
        
        return array([max(linalg.eigvals(obj.cov)) for obj in self])


class IntraField(Field):
    def __init__(self, parent_inout_spec, index_x, index_y):
        super(IntraField, self).__init__(parent_inout_spec, index_x, index_y)

    @property
    def mean_field(self):
        return self.parent.mean_fields_intra[self.index_x] +\
               1j*self.parent.mean_fields_intra[self.index_y]
    @property
    def cov_unrotated(self):
        return self.parent.sub_cov(self.index_x, self.index_y)

    @property
    def cov(self):
        return self.cov_rotated


ID2 = identity(2)
class InOutField(Field):
    """
    A propagating field (either towards of from the system), adds the concept of losses.
    """
    def __init__(self, parent_inout_spec, index_x, index_y, losses=0.):
        super(InOutField, self).__init__(parent_inout_spec, index_x, index_y)
        self.losses = 0

    @property
    def cov(self):
        return (1 - self.losses)*self.cov_rotated + self.losses*ID2


class InField(InOutField):
    """
    A propagating field used as input. Usual states can be specified such as squeezing (freq. indep.)
    or thermal. In addition, one can provide a "mean_field" for DC Ccalculations
    """

    def __init__(self ,mean_field=0. + 0j, *args, **kwds):
        super(InField, self).__init__(*args, **kwds)
        self.cov_unrotated = identity(2)
        self.mean_field = mean_field
    
    def set_squeezed_state(self, level, angle=0.):
        """
        Sets this mode's covariance matrix as that of a squeezed state.
        """

        self.cov_unrotated = array([[level, 0],[0, 1./level]]) #to check (angle)
        #self.parent.update()

    def set_thermal_state(self, n_th=0.):
        self.cov_unrotated = array([[2*n_th + 1, 0],[0, 2*n_th + 1]])
        #self.parent.update()

    def set_classical_noise(self, n_quanta):
        raise NotImplementedError

    def set_shot_noise(self):
        """Default behaviour"""

        self.cov_unrotated = identity(2)
        #self.parent.update()

class OutField(InOutField):
    """
    The spectra are calculated according to the parent out_cov matrix
    """
        
    @property
    def cov_unrotated(self):
        return self.parent.sub_cov_out(self.index_x, self.index_y)

    @property
    def mean_field(self):
        return self.parent.mean_fields_output[self.index_x] +\
               1j*self.parent.mean_fields_output[self.index_y]

class InOutPort(object):
    """
    Object encapsulating two fields: an input field and an output field.
    """

    def __init__(self,
                 parent,
                 index_x,
                 index_y,
                 coupling_hz,
                 input_loss=0,
                 output_loss=0,
                 mean_field_input=0. + 0j):
        """
        Creates all the fields for you, just make sure index_x, index_y are the
        quadratures index corresponding to the field in parent.evolution_matrix
        """

        self.parent = parent
        
        self.index_x = index_x
        self.index_y = index_y
        self.index_x_out = 2*len(parent.ports)
        self.index_y_out = 2*len(parent.ports) + 1
        
        self.input = InField(parent_inout_spec=self.parent,
                             mean_field=mean_field_input,
                             index_x=self.index_x_out,
                             index_y=self.index_y_out,
                             losses=input_loss)
        self.output = OutField(self.parent, index_x=self.index_x_out,
                                  index_y=self.index_y_out,
                                  losses=output_loss)
        self.coupling_hz = coupling_hz
