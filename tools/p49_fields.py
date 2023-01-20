import yt
import numpy as np

NH_Av = 1.8e21                  # [atoms/cm^2/mag] 
N_avg = 1000 * 4.6 * 3.09e18    # [atoms/cm^2 * code_length
verbose = False
def get_n0(factor):
    """    
    C = AvCutoff * N_H/Av / N_avg
      =  3 * 1.8e21 atoms/cm^2/mag * (1000 atoms cm^-3/code_density * 4.6pc/code_length)^{-1}
    AvCutoff = 3 mag - the visual extinction above which there is complete depolarization
    N_H/Av = 1.8e21 atoms/cm^2/mag - the ratio of column density of hydrogen nuclei to visual extinction
    N_avg = (1000 cm^-3/code_density * 4.6pc/code_length) - the average collumn density of the 
                                                            simulation in units of code_length^{-2} 
    1 pc = 3.09e18 cm 
    """
    AvCutoff = 3                    # [mag]
    NH_Av = 1.8e21                  # [atoms/cm^2/mag] 
    N_avg = 1000 * 4.6 * 3.09e18    # [atoms/cm^2 * code_length^2)]

    C = AvCutoff * NH_Av / N_avg    # C is approximately 3.8/code_length^2
     
    L_dense_pix = 1.                # Approximate size of a dense region in pixels (512x512 res)
    L_dense = L_dense_pix/512.      # Approximate size of a dense region in code_length = pixels/resolution
                                    # Dense regions ~1 pixel for 512 resolution.
    n0 = C/L_dense * factor           # Critical density [code_density], above which there is depolarization 
                                    # n0 is approximately 195 for L_dense_pix = 1 with factor = 1
    n0_str = str(int(round(n0)))   
    return n0, n0_str

def add_epsilon_old(axis, factor, this_thing=yt):
    """ 
    makes an epsilon field for yt.
    epsilon is a scaling factor that scales the stokes parameters 
    for high density regions where there is depolarization. 
    Above a visual extinction of 3 mag.
    """
    n0, n0_str = get_n0(factor)
    NH_Av = 1.8e21                  # [atoms/cm^2/mag] 
    N_avg = 1000 * 4.6 * 3.09e18    # [atoms/cm^2 * code_length^2)]

    def _epsilon_local(field,data):
        """ 
        Calculate epsilon.
        epsilon = code_density/(1.8e21 atoms/cm^2/mag) * 1000 atoms/cm^3/code_density * 4.6pc/code_length
                = code_density * N_avg/NH_Av
        [epsilon] = mag/code_length 
        -> Extinction, Av = Integral of epsilon along line of sight.
        """
        p = data.get_field_parameter('p')
        if p is None: p = 0
        n = data['density']

        epsilon = np.ones(data['density'].shape)
        epsilon[ n <= n0 ] = (n.v)[ n <= n0 ] * N_avg/NH_Av 
        epsilon[ n > n0 ]  = (n0**(1-p) * n.v**p)[ n > n0 ] * N_avg/NH_Av

        return epsilon

    this_thing.add_field('epsilon_n0-%s_%s'%(n0_str, axis), function=_epsilon_local, sampling_type='cell')
    
    return n0_str

def add_epsilon(n0,p, this_thing=yt): 
    def _eps_local(field,data):
        """This function calculates the Stokes Parameter "Q" along an axis x, y, or z.
        Makes use of the depolarization factor "epsilon" using a power exponent.
        """
        
        epsilon = np.ones(data['density'].shape)
        n = data['density']
        epsilon[ n <= n0 ] = (n.v)[ n <= n0 ]  
        epsilon[ n > n0 ]  = (n0**(1-p) * n.v**p)[ n > n0 ]   
  
        return  epsilon 
    
    if verbose: print( 'adding yt field epsilon_n0-%04d_p-%d'%(n0,p))
    this_thing.add_field('epsilon_n0-%04d_p-%d'%(n0,p), function=_eps_local, force_override=True,sampling_type='cell')
def add_stokes(axis, n0, p, this_thing=yt):
    """makes a stokes field for yt.
    axis should be x,y,z."""
    field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
    field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]
    add_epsilon(n0,p)
    

    def _Q_local(field,data):
        """This function calculates the Stokes Parameter "Q" along an axis x, y, or z.
        Makes use of the depolarization factor "epsilon" using a power exponent.
        """
        
        n = data['density']
        B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0

        #epsilon = np.ones(data['density'].shape)
        #epsilon[ n <= n0 ] = (n.v)[ n <= n0 ]  
        #epsilon[ n > n0 ]  = (n0**(1-p) * n.v**p)[ n > n0 ]   
        epsilon = n
        out = ( epsilon * (((data[field_horizontal])**2.0) - ((data[field_vertical])**2.0))/B_sq )
        out[B_sq == 0] = 0
  
        return out
    
    #Q_fname= 'Q%s_n0-%04d_p-%d'%(axis,n0,p)
    Q_fname= 'gas','Q%s'%(axis) #_n0-%04d_p-%d'%(axis,n0,p)
    if verbose: print( 'adding yt field %s'%Q_fname)
    this_thing.add_field(Q_fname, units='code_density', function=_Q_local, force_override=True,sampling_type='cell')

    def _U_local(field,data):
        """Makes stokes U."""
        
        n = data['density']
        B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0    

        #epsilon = np.ones(data['density'].shape)
        #epsilon[ n <= n0 ] = (n.v)[ n <= n0 ]  
        #epsilon[ n > n0 ]  = n #(n0**(1-p) * n.v**p)[ n > n0 ] 
        epsilon = n
        
        out = (2.0 * epsilon * ((data[field_horizontal]) * (data[field_vertical]))/B_sq)
        out[B_sq == 0 ] = 0
        return out 

    #U_fname = 'U%s_n0-%04d_p-%d'%(axis,n0,p)
    U_fname = 'enzo','U%s'%(axis)
    if verbose: print( 'adding yt field %s'%U_fname)
    this_thing.add_field(U_fname, units='g/cm**3', function=_U_local, force_override=True,sampling_type='cell')
    return Q_fname, U_fname


def add_unweighted_stokes(axis, this_thing=yt):
    """makes a stokes field for yt.
    axis should be x,y,z.
    These should test the projection and the rest of the pipeline.
    Second it's an interesting demonstration of what Stokes traces."""
    field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
    field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]

    def _unweighted_Q_local(field,data):
        """This function calculates the Stokes Parameter "Q"."""
        B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0

        return (data[field_horizontal]**2.0 - data[field_vertical]**2.0)/B_sq

    fieldname = 'unweighted_Q%s'%axis
    if verbose: print ("ADDING", fieldname)
    this_thing.add_field(fieldname, units='dimensionless', function=_unweighted_Q_local, sampling_type='cell')

    def _unweighted_U_local(field,data):
        """Makes stokes U."""
        B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0

        return 2*(data[field_horizontal]) * (data[field_vertical])/B_sq 

    fieldname = 'unweighted_U%s'%axis
    if verbose: print( "Adding", fieldname)
    this_thing.add_field(fieldname, units='dimensionless', function=_unweighted_U_local, sampling_type='cell')

def add_N2(axis, n0, p, this_thing=yt):
    """ Makes a field that when projected is a correction to the column density used
    in calculating the polarization fraction. """
    field_horizontal = {'x':'By','y':'Bz','z':'Bx'}[axis]
    field_vertical   = {'x':'Bz','y':'Bx','z':'By'}[axis]

    def _N2_local(field,data):
        """ Calculate n2 = n * (cos^2(gamma)/2 - 1/3) where gamma is the 
        inclination angle between the B field and the plane of the sky. 
        cos(gamma) = B_sky dot B / |B_sky|*|B| 
        e.g. Line of sight along z axis. B = (Bx,By,Bz); B_sky = (Bx,By,0)
        cos(gamma) = Bx^2 + By^2 / (sqrt(Bx^2 + By^2) * sqrt(Bx^2 + By^2 + Bz^2))
        cos(gamma)^2 = Bx^2 + By^2 / (Bx^2 + By^2 + Bz^2)

        This function returns n2 as a dimensionless value since epsilon was calculated as dimensionless"""

        B_sq = data['Bx']**2.0 + data['By']**2.0 + data['Bz']**2.0
        cos_gamma_sq = (data[field_horizontal]**2.0 + data[field_vertical]**2.0)/B_sq

        n = data['density'].in_units('code_density')

        if 0:
            n0 = 1.0
            p=1.0
            if hasattr(n0,'units'):
                n0=data.ds.quan(n0,'gram/cm**3')
            epsilon = np.ones_like(data['density'])
            epsilon[ n <= n0 ] = (n)[ n <= n0 ]  
            epsilon[ n > n0 ]  = (n0**(p) * n**p)[ n > n0 ]  
        else:
            epsilon = n #data.ds.quan(1.0,'gram/cm**3')
        return epsilon * (0.5*cos_gamma_sq - 1/3) 
       
    fieldname = 'N2%s_n0-%04d_p-%d'%(axis,n0,p)
    this_thing.add_field(fieldname, units='code_density', function=_N2_local, sampling_type='cell')    
    if verbose: print( "Added", fieldname)

# Add yt fields for stokes and n2 along each axis with 
# different cutoff density n0 and powerlaw index p
def add_QU(this_ds):
    for axis in ['x', 'y', 'z']:
#   for n0 in [19,39,1945]:
#       add_stokes(axis, n0, p=0)
#       add_N2(axis, n0, p=0)
        add_stokes(axis, n0=1, p=1, this_thing=this_ds)
        #add_N2(axis, n0=1, p=1, this_thing=this_ds)
        #add_unweighted_stokes(axis)




