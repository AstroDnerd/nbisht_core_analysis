from starter2 import *
import xtra_operators as xo
reload(xo)
kinetic_validators=[yt.ValidateSpatial(1,['x-velocity','y-velocity','z-velocity','kinetic_energy_density'])]
#std_validators = [yt.ValidateSpatial(1,['Bx','By','Bz', 'x-velocity','y-velocity','z-velocity'])]
std_validators_2 = [yt.ValidateSpatial(1,['magnetic_field_x','magnetic_field_y','magnetic_field_z', 'x-velocity','y-velocity','z-velocity'])]
mag_validators = [yt.ValidateSpatial(1,['magnetic_field_x','magnetic_field_y','magnetic_field_z'])]
kinetic_validators=[yt.ValidateSpatial(1,['x-velocity','y-velocity','z-velocity','kinetic_energy_density'])]
kinetic_validators=[yt.ValidateSpatial(1,['density'])]
pressure_validators=[yt.ValidateSpatial(1,['x-velocity','y-velocity','z-velocity','pressure'])]
work_units='erg/(cm**3*s)'
def grad(data,fieldname,direction):
    iM1 = slice(None,-2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=tuple([all]*3)
    dxi=1./(2*data.dds )
    out = np.zeros_like(data[fieldname]*dxi[0])
    Left = [all]*3
    Right = [all]*3
    Right[direction] = iP1
    Left[direction] = iM1
    Left=tuple(Left); Right=tuple(Right)
    out[all_all] = (data[fieldname][ Right ]- data[fieldname][ Left]) *dxi[direction]
    return out

def add_v_grad(obj):
    def dxvx(field,data):
        return grad( data, 'velocity_x',0)
    def dxvy(field,data):
        return grad( data, 'velocity_y',0)
    def dxvz(field,data):
        return grad( data, 'velocity_z',0)
    def dyvx(field,data):
        return grad( data, 'velocity_x',1)
    def dyvy(field,data):
        return grad( data, 'velocity_y',1)
    def dyvz(field,data):
        return grad( data, 'velocity_z',1)
    def dzvx(field,data):
        return grad( data, 'velocity_x',2)
    def dzvy(field,data):
        return grad( data, 'velocity_y',2)
    def dzvz(field,data):
        return grad( data, 'velocity_z',2)
    obj.add_field(YT_dxvx, dxvx, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dxvy, dxvy, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dxvz, dxvz, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dyvx, dyvx, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dyvy, dyvy, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dyvz, dyvz, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dzvx, dzvx, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dzvy, dzvy, units='1/s', sampling_type='cell', validators=std_validators_2)
    obj.add_field(YT_dzvz, dzvz, units='1/s', sampling_type='cell', validators=std_validators_2)

def add_b_over_rho(obj):
    def bx_hat(field,data):
        return data['magnetic_field_x']/data['magnetic_field_strength']
    obj.add_field(YT_bx_hat,bx_hat, units='dimensionless', sampling_type='cell')
    def by_hat(field,data):
        return data['magnetic_field_y']/data['magnetic_field_strength']
    obj.add_field(YT_by_hat,by_hat, units='dimensionless', sampling_type='cell')
    def bz_hat(field,data):
        return data['magnetic_field_z']/data['magnetic_field_strength']
    obj.add_field(YT_bz_hat,bz_hat, units='dimensionless', sampling_type='cell')
    def bx_over_rho(field,data):
        output=xo.AdotDel(data,[YT_bx_hat,YT_by_hat,YT_bz_hat],'velocity_x')
        return output
    obj.add_field(YT_bx_over_rho,bx_over_rho, units='1/s',sampling_type='cell',validators=std_validators_2)
    def by_over_rho(field,data):
        output=xo.AdotDel(data,[YT_bx_hat,YT_by_hat,YT_bz_hat],'velocity_y')
        return output
    obj.add_field(YT_by_over_rho,by_over_rho, units='1/s',sampling_type='cell',validators=std_validators_2)
    def bz_over_rho(field,data):
        output=xo.AdotDel(data,['bx_hat','by_hat','bz_hat'],'velocity_z')
        return output
    obj.add_field(YT_bz_over_rho,bz_over_rho, units='1/s',sampling_type='cell',validators=std_validators_2)

def add_gdotgradrho(obj):
    def gdotgradrho(field,data):
        drho_dx = xo.grad(data,YT_density,0)
        drho_dy = xo.grad(data,YT_density,1)
        drho_dz = xo.grad(data,YT_density,2)
        output = data[YT_acceleration_x]*drho_dx+ data[YT_acceleration_y]*drho_dy+ data[YT_acceleration_z]*drho_dz
        return data.ds.arr(output,'g/(cm**3*s**2)')
    obj.add_field(YT_gdotgradrho,gdotgradrho,units='g/(cm**3*s**2)', validators=[yt.ValidateSpatial(1,YT_density)], sampling_type='cell')

def add_grav_test(obj):
    def also_rho(field,data):
        #gx =data.ds.arr(data[YT_acceleration_x],'code_length/code_time**2')
        #gy =data.ds.arr(data[YT_acceleration_y],'code_length/code_time**2')
        #gz =data.ds.arr(data[YT_acceleration_z],'code_length/code_time**2')
        dgxdx = xo.grad(data,YT_acceleration_x,0)
        dgydy = xo.grad(data,YT_acceleration_y,1)
        dgzdz = xo.grad(data,YT_acceleration_z,2)
        four_pi_G = data.ds['GravitationalConstant']
        rho = -(dgxdx+dgydy+dgzdz)/(four_pi_G)
        rho = data.ds.arr(rho,'code_density')
        return rho
    obj.add_field(YT_also_rho,also_rho, units='code_density', sampling_type='cell', 
                  validators=[yt.ValidateSpatial(1,[YT_acceleration_x,YT_acceleration_y,YT_acceleration_z])])

def add_energies(obj):
    if obj.parameters['SelfGravity']:
        def grav_energy2(field,data):
            gx =data.ds.arr(data[YT_acceleration_x],'code_length/code_time**2')
            gy =data.ds.arr(data[YT_acceleration_y],'code_length/code_time**2')
            gz =data.ds.arr(data[YT_acceleration_z],'code_length/code_time**2')
            alpha = 1./(2*data.ds['GravitationalConstant']) #=1/8 pi G)
            alpha = data.ds.quan(alpha, '1/(code_length**3/code_time**2/code_mass)')
            return ( -(gx**2+gy**2+gz**2)*alpha )
        obj.add_field(YT_grav_energy_2,grav_energy2, units='code_mass*code_length**2/(code_time**2*code_length**3)', sampling_type='cell')
        def grav_energy(field,data):
            gx =data.ds.arr(grad(data,YT_potential_field,0),'code_length/code_time**2')
            gy =data.ds.arr(grad(data,YT_potential_field,1),'code_length/code_time**2')
            gz =data.ds.arr(grad(data,YT_potential_field,2),'code_length/code_time**2')
            alpha = 1./(2*data.ds['GravitationalConstant']) #=1/8 pi G)
            alpha = data.ds.quan(alpha, '1/(code_length**3/code_time**2/code_mass)')
            return ( -(gx**2+gy**2+gz**2)*alpha )
        obj.add_field(YT_grav_energy,grav_energy,validators=[ yt.ValidateSpatial(1, 'PotentialField')],
                 units='code_mass*code_length**2/(code_time**2*code_length**3)', sampling_type='cell')
#       def abs_grav_energy(field,data):
#           return np.abs( data[YT_grav_energy] )
#       obj.add_field(YT_abs_grav_energy,abs_grav_energy,validators=[yt.ValidateSpatial(1,'PotentialField')],
#                units='code_mass*code_length**2/(code_time**2*code_length**3)', sampling_type='cell')

    def therm_energy(field,data):
        sound_speed = data.ds.quan(1.,'code_velocity')
        rho_0 = data.ds.quan(1.,'code_mass/code_length**2')
        e = sound_speed**2*np.log( data[YT_density]/rho_0)
        therme = data[YT_density]*e
        return therme
    obj.add_field(YT_therm_energy,therm_energy,
                 units='code_mass*code_length**2/(code_time**2*code_length**3)', sampling_type='cell')

def add_gravity(obj):
    add_energies(obj)
    def density_log(field,data):
        return np.log(data['density'].v)
    obj.add_field(YT_density_log,density_log,sampling_type='cell')
    if obj.parameters['SelfGravity']:
        #def rho_x(field,data):
        #    FFF = data[YT_potential_field].v
        #    gi = xo.gradf( FFF, 0,data.dds)
        #    pdb.set_trace()
        #    #gi =-grad(data,YT_density_log,0)
        #    #gi+=-grad(data,YT_density_log,1)
        #    #gi+=-grad(data,YT_density_log,1)
        #    return gi
        #obj.add_field(YT_density_grad_x,rho_x,validators=[yt.ValidateSpatial(1,[YT_density_log])],
        #             units='1/code_length',sampling_type='cell')
        def grav_x(field,data):

            gi =-data.ds.arr(grad(data,'PotentialField',0),'code_length/code_time**2')
            return gi
        obj.add_field(YT_grav_x,grav_x,validators=[yt.ValidateSpatial(1,[YT_potential_field])],
                     units='code_length/(code_time**2)',sampling_type='cell')
        def grav_y(field,data):
            gi =-data.ds.arr(grad(data,'PotentialField',1),'code_length/code_time**2')
            return gi
        obj.add_field(YT_grav_y,grav_y,validators=[yt.ValidateSpatial(1,[YT_potential_field])],
                     units='code_length/(code_time**2)',sampling_type='cell')
        def grav_z(field,data):
            gi =-data.ds.arr(grad(data,'PotentialField',2),'code_length/code_time**2')
            return gi
        obj.add_field(YT_grav_z,grav_z,validators=[yt.ValidateSpatial(1,[YT_potential_field])],
                     units='code_length/(code_time**2)',sampling_type='cell')
        def gravity_work(field,data):
            try:
                gi =(data[YT_velocity_x]*data[YT_grav_x]+
                     data[YT_velocity_y]*data[YT_grav_y]+
                     data[YT_velocity_z]*data[YT_grav_z])
                gi = -gi*data[YT_density]
            except:
                return np.zeros_like(data[YT_density])
            
            return gi
        obj.add_field(YT_gravity_work,gravity_work,validators=pressure_validators,
                     units=work_units,sampling_type='cell', take_log=True)
        def grav_norm(field,data):
            return np.sqrt( data[YT_grav_x]**2 + data[YT_grav_y]**2 + data[YT_grav_z]**2)
        obj.add_field(YT_grav_norm,grav_norm, validators=pressure_validators,
                      units='code_length/code_time**2', sampling_type='cell')
        def geke(field, data):
            out = data.ds.arr(np.zeros_like( data[YT_grav_energy].v), 'dimensionless')
            ok = np.abs(data[YT_grav_energy])>0
            out[ok]=np.abs(data[YT_grav_energy][ok].v)/data[YT_kinetic_energy][ok].v
            return out.v
        obj.add_field(YT_ge_ke, function=geke, sampling_type='cell', validators=[yt.ValidateSpatial(1)])
        def geke2(field, data):
            out = data.ds.arr(np.zeros_like( data[YT_grav_energy_2].v), 'dimensionless')
            ok = np.abs(data[YT_grav_energy_2])>0
            out[ok]=np.abs(data[YT_grav_energy_2][ok].v)/data[YT_kinetic_energy][ok].v
            return out.v
        obj.add_field(YT_ge_ke_2, function=geke2, sampling_type='cell', validators=[yt.ValidateSpatial(1)])
def add_force_terms(obj):
    def gas_pressure(field,data):
        return data[YT_density]*data.ds.quan(1,'code_velocity')**2
    obj.add_field(YT_gas_pressure,gas_pressure,units='dyne/cm**2', sampling_type='cell')
    def pressure_work(field,data):
        try:
            gi =-1*xo.AdotDel(data, [YT_velocity_x, YT_velocity_y, YT_velocity_z], YT_gas_pressure) #-data.ds.arr(grad(data,'PotentialField',2),'code_length/code_time**2')
        except:
            return np.zeros_like(data[YT_density])
        
        return gi
    obj.add_field(YT_pressure_work,pressure_work,validators=pressure_validators,
                 units=work_units,sampling_type='cell', take_log=True)
    def dpdx(field,data):
        output = xo.grad(data,YT_gas_pressure,0)
        return output
    obj.add_field(YT_dpdx,dpdx,units='dyne/cm**3', validators=[yt.ValidateSpatial(1,YT_pressure)], sampling_type='cell')
    def dpdy(field,data):
        output = xo.grad(data,YT_gas_pressure,1)
        return output
    obj.add_field(YT_dpdy,dpdy,units='dyne/cm**3', validators=[yt.ValidateSpatial(1,YT_pressure)], sampling_type='cell')
    def dpdz(field,data):
        output = xo.grad(data,YT_gas_pressure,2)
        return output
    obj.add_field(YT_dpdz,dpdz,units='dyne/cm**3', validators=[yt.ValidateSpatial(1,YT_pressure)], sampling_type='cell')
    def pressure_accel_mag(field,data):
        return np.sqrt( data[YT_dpdx]**2+data[YT_dpdy]**2+data[YT_dpdz]**2)/data[YT_density]
    obj.add_field(YT_pressure_accel_mag,pressure_accel_mag,units='code_length/code_time**2', validators=[yt.ValidateSpatial(1,[YT_pressure])], sampling_type='cell')
    def momentum_flux(field,data):
        f1 = xo.gradf(0.5*data[YT_velocity_x]*data[YT_kinetic_energy],0,data.dds)
        f2 = xo.gradf(0.5*data[YT_velocity_y]*data[YT_kinetic_energy],1,data.dds)
        f3 = xo.gradf(0.5*data[YT_velocity_z]*data[YT_kinetic_energy],2,data.dds)
        return (f1+f2+f3)
    obj.add_field(YT_momentum_flux,momentum_flux, validators=kinetic_validators, sampling_type='cell',
                  units=work_units)
    current_units='gauss/cm'
if 0:
    def current_x(field,data):
        Bvector = [data['magnetic_field_%s'%s] for s in 'xyz']
        return xo.Curl(Bvector   , data.dds,component=0)
    obj.add_field('current_x',current_x, validators=mag_validators, units=current_units, sampling_type='cell')
    def current_y(field,data):
        Bvector = [data['magnetic_field_%s'%s] for s in 'xyz']
        return xo.Curl(Bvector   , data.dds,component=1)
    obj.add_field('current_y',current_y, validators=mag_validators, units=current_units, sampling_type='cell')
    def current_z(field,data):
        Bvector = [data['magnetic_field_%s'%s] for s in 'xyz']
        return xo.Curl(Bvector   , data.dds,component=2)
    obj.add_field('current_z',current_z, validators=mag_validators, units=current_units, sampling_type='cell')
    def mag_work(field,data):
        """v dot j cross B"""
        Bvector = [data['magnetic_field_%s'%s] for s in 'xyz']
        current = xo.Curl(Bvector, data.dds)
        output  = data['velocity_x']*(current[1]*Bvector[2]-current[2]*Bvector[1])
        output += data['velocity_y']*(current[2]*Bvector[0]-current[0]*Bvector[2])
        output += data['velocity_z']*(current[0]*Bvector[1]-current[1]*Bvector[2])
        output *= 1./(np.pi*4)
        return output
    obj.add_field('mag_work',mag_work, validators=std_validators_2, units=work_units, sampling_type='cell',take_log=True)

    def mag_force_norm(field,data):
        """||j cross B||"""
        Bvector = [data['magnetic_field_%s'%s] for s in 'xyz']
        current = xo.Curl(Bvector, data.dds)
        output  = (current[1]*Bvector[2]-current[2]*Bvector[1])**2
        output += (current[2]*Bvector[0]-current[0]*Bvector[2])**2
        output += (current[0]*Bvector[1]-current[1]*Bvector[2])**2
        output = np.sqrt(output)
        output *= 1./(np.pi*4)
        return output
    obj.add_field('mag_force_norm',mag_force_norm, validators=std_validators_2, units=work_units, sampling_type='cell',take_log=True)
    def mag_vel_angle(field,data):
        return data['mag_work']/(data['mag_force_norm']*data['velocity_magnitude'])
    obj.add_field('mag_vel_angle',mag_vel_angle,validators=std_validators_2, units=work_units, sampling_type='cell')


    #def pressure_force(field,data):
    #    output  =-data['x-velocity']*xo.gradf(data['gas_pressure'],0,data.dds)
    #    output +=-data['y-velocity']*xo.gradf(data['gas_pressure'],1,data.dds)
    #    output +=-data['z-velocity']*xo.gradf(data['gas_pressure'],2,data.dds)
    #    output = data.ds.arr(np.zeros(data['density'].shape), 'dyne/(cm**2*s)')
    #    return output
    #obj.add_field('pressure_force',pressure_force,validators=pressure_validators,
    #              units='dyne/(cm**2*s)')

