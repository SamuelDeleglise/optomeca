from pylab import *

class InOutSpec(object):
    def __init__(self, evolution_mat=None, freq=linspace(0.99e6,1.01e6)):
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

#        if cov_input is not None:
#            self._cov_input = cov_input
        if evolution_mat is not None:
            self._evolution_mat        = evolution_mat
#        if coupling_mat is not None:
#            self._coupling_mat         = coupling_mat
        if freq is not None:
            self._freq      = freq

        self.fields = []
        self.ports = []
        self.current_freq = 0.

    
    def evolution_mat(self):
        return self._evolution_mat
########################################################
###     FUNCTIONS TO IMPLEMENT IN CHILD CLASS
########################################################
##    def cov_input(self):
##        return self._cov_input
##
##    def evolution_mat(self):
##        return self._evolution_mat

########################################################
###     PORTS/FIELDS DEFINITIONS
########################################################

    
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
        
        return len(self.evolution_mat())

    @property
    def n_ports(self):
        """
        Number of propagating modes
        """
        return len(self.ports)
        #return len(self.coupling_mat()[0,:])
########################################################
###     FREQUENCY PROPERTIES
########################################################
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


    @property
    def freq(self):
        return self._freq

    @freq.setter
    def freq(self, val):
        self._freq = val
        return val
########################################################
###     RESULTING COVARIANCE MATRICES
########################################################

    @property
    def cov_intra(self):
        ci = self.cov_input()
        return array([(left.dot(ci)).dot(right.T) \
                          for left, right in zip(self._left_intra(), self._right_intra())])

    def cov_input(self):
        mat = zeros((2*self.n_ports, 2*self.n_ports))
        for index, port in enumerate(self.ports):
            mat[2*index, 2*index] = port.input.cov[0,0]
            mat[2*index, 2*index + 1] = port.input.cov[0,1]
            mat[2*index + 1, 2*index] = port.input.cov[1,0]
            mat[2*index + 1, 2*index + 1] = port.input.cov[1,1]
        return mat

    def coupling_mat(self):
        mat = zeros((self.n_intra, 2*self.n_ports))
        for index, port in enumerate(self.ports):
            mat[port.index_x, 2*index] = port.coupling_sqrthz
            mat[port.index_y, 2*index + 1] = port.coupling_sqrthz
        return mat

    
    def cov_out(self):
        trans_coupling = self.coupling_mat().T
        idn = identity(self.n_ports*2)
        ci = self.cov_input()
        def _out_from_in(intra):
            return (trans_coupling).dot(intra) - idn
        
        return array([_out_from_in(left).dot(ci).dot((_out_from_in(right)).T)\
                      for left, right in zip(self._left_intra(), self._right_intra())])
    
########################################################
###     CALCULATION INTERMEDIATES
########################################################
    def sub_cov(self, index_x, index_y):
        parent_cov = self.cov_intra
        arr = array([[parent_cov[:,index_x, index_x], parent_cov[:,index_x, index_y]],
                    [parent_cov[:,index_y, index_x], parent_cov[:,index_y, index_y]]])
        return arr.transpose(2,0,1)

    def sub_cov_out(self, index_x, index_y):
        parent_cov = self.cov_out()
        arr = array([[parent_cov[:,index_x, index_x], parent_cov[:,index_x, index_y]],
                    [parent_cov[:,index_y, index_x], parent_cov[:,index_y, index_y]]])
        return arr.transpose(2,0,1)

    def _left_intra(self):
        left  =  1j*tensordot(self.omega, identity(self.n_intra),axes=0) - self.evolution_mat()
        return linalg.inv(left).dot(self.coupling_mat())

    def _right_intra(self):
        right = -1j*tensordot(self.omega, identity(self.n_intra),axes=0) - self.evolution_mat()
        return linalg.inv(right).dot(self.coupling_mat())
        


class Field(object):
    """
    Object encapsulating 2-quadratures and a mean-field. 
    """

    def __init__(self, get_cov=None, get_mean_field_angle=lambda:0., get_losses=lambda:0.):
        self._get_mean_field_angle = get_mean_field_angle
        if get_cov is not None:
            self._get_cov = get_cov
            self._get_squeezing = self._calculate_squeezing
            self._get_anti_squeezing = self._calculate_anti_squeezing
        else:
            self._squeezing = 1.
            self._squeezing_angle = 0.
            self._get_cov = self._get_arb_cov
        self._get_losses = get_losses

    def _get_arb_cov(self):
        return array([[self.squeezing, 0],[0, 1./self.squeezing]])


    @property
    def losses(self):
        return self._get_losses()

    @property
    def mean_field_angle(self):
        return self._get_mean_field_angle()
    
    def spectrum_sym(self, quad_angle):
        """
        Returns the symmetrized spectrum along the quadrature quad_angle:
          - 0 is for intensity quadrature
          - pi/2 is for phase quadrature
        """
        
        rot = self.rotation_mat(quad_angle - self.mean_field_angle)
        return array([rot.dot(cov).dot(rot.T)[0,0] for cov in self.cov])

    def rotation_mat(self, angle):
        return array([[np.cos(angle), -np.sin(angle)],[np.sin(angle), np.cos(angle)]])

    def spectrum_ss(self, side='upper'):
        """
        Returns the single-sided spectrum of the field
          - 'upper' is for omega>0 (blue sideband)
          - 'lower' is for omega<0 (red sideband)
        """

 #       rot = self.rotation_mat(-self.mean_field_angle)
#        cov_r = rot.dot(self.cov).dot(rot.T)

        if side=='upper':
            return array([0.5*(cov_r[0,0] + cov_r[1,1] + 1j*(cov_r[1,0] - cov_r[0,1])) for cov_r in self.cov_rotated])
        if side=='lower':
            return array([0.5*(cov_r[0,0] + cov_r[1,1] - 1j*(cov_r[1,0] - cov_r[0,1])) for cov_r in self.cov_rotated])
        raise(ValueError('side should be either "upper" (omega>0) or "lower" (omega<0)'))

    
    @property
    def cov(self):
        """
        Reimplement this property 
        """
        return (1 - self.losses)*self._get_cov() + identity(2)*self.losses

    @property
    def cov_rotated(self):
        """
        Covariance matrix in the frame of the mean field.
        """

        return self.do_rotate(self.cov, -self.mean_field_angle)
    
    def do_rotate(self, matrices, angle):
        rot = self.rotation_mat(angle)
        rott = rot.T
        return array([rot.dot(a).dot(rott) for a in matrices])

    
    def _calculate_squeezing(self):
        """
        Returns the lowest noise level for any quadrature (by diagonalizing the covariance matrix)
        """

        return array([min(linalg.eigvals(cov)) for cov in self.cov])
    
    def _calculate_anti_squeezing(self):
        """
        Returns the highest noise level for any quadrature (by diagonalizing the covariance matrix)
        """
        return array([max(linalg.eigvals(cov)) for cov in self.cov])


    def _get_squeezing(self):
        return self._squeezing

    def _get_anti_squeezing(self):
        return 1./self._squeezing
    
    @property
    def squeezing(self):
        return self._get_squeezing()

    @property
    def anti_squeezing(self):
        return self._get_anti_squeezing()

    @squeezing.setter
    def squeezing(self, val):
        if not hasattr(self, '_squeezing'):
            raise(ValueError('This field is not an input field: squeezing cannot be set'))
        self._squeezing = val
        return val
        
class IntraField(Field):
    def __init__(self, parent, index_x, index_y, **kwds):
        self.parent = parent
#        self.index_in_out = len(parent.ports)

        def get_cov():
            return parent.sub_cov(index_x, index_y)
        super(IntraField, self).__init__(get_cov=get_cov, **kwds)

    


        
class InOutPort(object):
    """
    Object encapsulating two fields.
    """

    def __init__(self,
                 parent,
                 index_x,
                 index_y,
                 get_input_angle=lambda:0.,
                 get_output_angle=lambda:0.,
                 get_coupling_sqrthz=lambda:0.,
                 get_input_cov=None,
                 get_output_losses=lambda:0.):
        """
        If getter functions are set, uses the function each time the property is needed.
        if get_input_cov is None, then uses the native built-in function of self.input
        """
        self.input = Field(get_cov=get_input_cov,
                           get_mean_field_angle=get_input_angle)
        self.output = Field(get_mean_field_angle=get_output_angle,
                            get_cov=self.get_output_cov,
                            get_losses=get_output_losses)
        self.index_x = index_x
        self.index_y = index_y
        self.parent = parent
        self.index_x_out = 2*len(parent.ports)
        self.index_y_out = 2*len(parent.ports) + 1
        self._get_coupling_sqrthz = get_coupling_sqrthz

    def get_output_cov(self):
        return self.parent.sub_cov_out(self.index_x_out, self.index_y_out)
        """arr = self.parent.sub_cov_out(self.index_x, self.index_y)
        return self.coupling_sqrthz*arr - self.input.cov
        """
    @property
    def coupling_sqrthz(self):
        """
        Typically sqrt(\kappa) where \kappa is the energy decay time via this port.
        """
        
        return self._get_coupling_sqrthz()

"""

    @property
    def cov(self):
        parent_cov = self.parent_cov
        index_x = self.index_x
        index_y = self.index_y
        arr = array([[parent_cov[:,index_x, index_x], parent_cov[:,index_x, index_y]],
                                    [parent_cov[:,index_y, index_x], parent_cov[:,index_y, index_y]]])
        
        arr = arr.transpose(2,0,1)
        if self.mean_field_angle!=0:
            arr = do_rotate(arr, self.mean_field_angle)
        return arr
    
"""

##class InOutPort(Field):
##    def __init__(self, *args):
##        super(MeasPort, self).__init__(*args)
##
##    def input_angle(self):
##        """
##        Phase of the input beam with respect to intra field.
##        To be reimplemented or monkey patched.
##        """
##        
##        return 0
##
##    def output_angle(self):
##        """
##        Phase of the input beam with respect to intra field.
##        To be reimplemented or monkey patched.
##        """
##
##        return 0
##
##    def coupling_sqrthz(self):
##        """
##        Coupling between the port and the intra field: 'sqrt(kappa)'
##        """
##        
##        return 0
##
##    @property
##    def cov_out(self):
##        parent_cov = self.parent.cov
##        index_x = self.index_x
##        index_y = self.index_y
##
##        arr = array([[parent_cov[:,index_x, index_x], parent_cov[:,index_x, index_y]],
##                    [parent_cov[:,index_y, index_x], parent_cov[:,index_y, index_y]]])
##        arr = arr.transpose(2,0,1)
##        
##        #return (1. - self.losses)*arr.transpose(2,0,1) + self.losses*identity(2)
##        if self.output_angle!=0:
##            arr = self.do_rotate(arr, self.mean_field_angle)
##        return (1. - self.losses)*arr + self.losses*identity(2)
###        return array([rot.dot(((1. - self.losses)*a + self.losses*identity(2))).dot(rot.T)]) 
##
##    @property
##    def parent_cov(self):
##        return self.parent.cov_out
