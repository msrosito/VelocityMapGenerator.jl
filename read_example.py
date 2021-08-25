import numpy as np
import h5py
import os
from argparse import ArgumentParser

#def read_header():
#    f = h5py.File('./data/snap 028 z000p000.0.hdf5', 'r')
#    a = f['Header'].attrs.get('Time')
#    h = f['Header'].attrs.get('HubbleParam')
#    boxsize = f['Header'].attrs.get('BoxSize')
#    f.close()
    
#    return a, h, boxsize 

def read_dataset(itype, att, nfiles=16):
    """ Read a selected dataset, itype is the PartType and att is the attribute name. """

    # Output array.
    data = []

    # Loop over each file and extract the data.
    for i in range(nfiles):
        
        f = h5py.File('./data/snap_028_z000p000.%i.hdf5'%i, 'r')
        tmp = f['PartType%i/%s'%(itype, att)][...]
        data.append(tmp)

        # Get conversion factors.
        cgs     = f['PartType%i/%s'%(itype, att)].attrs.get('CGSConversionFactor')
        aexp    = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
        hexp    = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')

        # Get expansion factor and Hubble parameter from the header.
        a       = f['Header'].attrs.get('Time')
        h       = f['Header'].attrs.get('HubbleParam')

        f.close()

    # Combine to a single array.
    if len(tmp.shape) > 1:
        data = np.vstack(data)
    else:
        data = np.concatenate(data)

    # Convert to physical.
    if data.dtype != np.int32 and data.dtype != np.int64:
        data = np.multiply(data, cgs * a**aexp * h**hexp, dtype='f8')

    return data
    
def read_dataset_dm_mass():
    """ Special case for the mass of dark matter particles. """

    f           = h5py.File('./data/snap_028_z000p000.0.hdf5', 'r')
    h           = f['Header'].attrs.get('HubbleParam')
    a           = f['Header'].attrs.get('Time')
    dm_mass     = f['Header'].attrs.get('MassTable')[1]
    n_particles = f['Header'].attrs.get('NumPart_Total')[1]

    # Create an array of length n_particles each set to dm_mass.
    m = np.ones(n_particles, dtype='f8') * dm_mass 

    # Use the conversion factors from the mass entry in the gas particles.
    cgs  = f['PartType0/Mass'].attrs.get('CGSConversionFactor')
    aexp = f['PartType0/Mass'].attrs.get('aexp-scale-exponent')
    hexp = f['PartType0/Mass'].attrs.get('h-scale-exponent')
    f.close()

    # Convert to physical.
    m = np.multiply(m, cgs * a**aexp * h**hexp, dtype='f8')

    return m


## cm2kpc(A):
#
## converts from cm to kpc
def cm2kpc(A):
    factor = 3.24e-22
    return factor * A

## gr210MSun(A):
#
## converts from gr to 1e10 MSun
def gr210MSun(A):
    factor = 5.03e-44
    return factor * A
    
## cms2kms(A):
#
## converts from cm/s to km/s
def cms2kms(A):
    factor = 1e-5
    return factor * A

      
## Nsub(nfiles):
#
## Returns the number of subgroups in the simulation
# debe devolver el nro de elementos diferentes. Dejar en suspenso.
def Nsub():
    
    sgns_g = read_dataset(0, 'SubGroupNumber')
    sgns_d = read_dataset(1, 'SubGroupNumber') 
    sgns_s = read_dataset(4, 'SubGroupNumber')
    
    return max(np.concatenate(sgns_g, sgns_d, sgns_s))
        
    
## gas_att(hf, sub_gas, nfiles= 16)):
#    
## Creates hdf5 file considering gas attributes and filtrers by sub_gas
def gas_att(hf, sub_gas, nfiles= 16):

    part0 = hf.create_group("PartType0")
    
    # positions
    coord = read_dataset(0, "Coordinates", nfiles=16)
    coord = np.transpose(coord)
    c0 = coord[0][sub_gas]
    c1 = coord[1][sub_gas]
    c2 = coord[2][sub_gas]
    Npart = len(c0)
#    coord = np.transpose(np.array([[c0],[c1],[c2]]))
    x0 = part0.create_dataset("x", data=cm2kpc(c0))
    x0.attrs["Units"] = "kpc"
    x0.attrs["PartNum"] = Npart
    y0 = part0.create_dataset("y", data=cm2kpc(c1))
    y0.attrs["Units"] = "kpc"
    y0.attrs["PartNum"] = Npart
    z0 = part0.create_dataset("z", data=cm2kpc(c2))
    z0.attrs["Units"] = "kpc"
    z0.attrs["PartNum"] = Npart 
    
    # velocities
    vel = read_dataset(0, "Velocity", nfiles=16)
    vel = np.transpose(vel)
    v0 = vel[0][sub_gas]
    v1 = vel[1][sub_gas]
    v2 = vel[2][sub_gas]
    vx0 = part0.create_dataset("vx", data=cms2kms(v0))
    vx0.attrs["Units"] = "km/s"
    vx0.attrs["PartNum"] = Npart
    vy0 = part0.create_dataset("vy", data=cms2kms(v1))
    vy0.attrs["Units"] = "km/s"
    vy0.attrs["PartNum"] = Npart
    vz0 = part0.create_dataset("vz", data=cms2kms(v2))
    vz0.attrs["Units"] = "km/s"
    vz0.attrs["PartNum"] = Npart        
    
    # masses
    mass = read_dataset(0, "Mass", nfiles=16)
    mass = mass[sub_gas]
    mass0 = part0.create_dataset("Mass", data=gr210MSun(mass))
    mass0.attrs["Units"] = "1e10 Msol"
    mass0.attrs["PartNum"] = Npart

## DM(hf, sub_gas, nfiles= 16):
#
# ## Creates hdf5 file considering dark matter attributes and filtrers by sub_dark  
def DM_att(hf, sub_dark, nfiles= 16):
    
    part1 = hf.create_group("PartType1")

    # positions
    coord = read_dataset(1, "Coordinates", nfiles=16)
    coord = np.transpose(coord)
    c0 = coord[0][sub_dark]
    c1 = coord[1][sub_dark]
    c2 = coord[2][sub_dark]
    Npart = len(c0)
#    r = np.multiply(c0,c0)+np.multiply(c1,c1)+np.multiply(c2,c2)
#    r = cm2kpc(np.sqrt(r))
#    print('rs = ', min(r), max(r), np.median(r))
    x1 = part1.create_dataset("x", data=cm2kpc(c0))
    x1.attrs["Units"] = "kpc"
    x1.attrs["PartNum"] = Npart
    y1 = part1.create_dataset("y", data=cm2kpc(c1))
    y1.attrs["Units"] = "kpc"
    y1.attrs["PartNum"] = Npart
    z1 = part1.create_dataset("z", data=cm2kpc(c2))
    z1.attrs["Units"] = "kpc"
    z1.attrs["PartNum"] = Npart  

    # velocities
    vel = read_dataset(1, "Velocity", nfiles=16)
    vel = np.transpose(vel)
    v0 = vel[0][sub_dark]
    v1 = vel[1][sub_dark]
    v2 = vel[2][sub_dark]
    vx1 = part1.create_dataset("vx", data=cms2kms(v0))
    vx1.attrs["Units"] = "km/s"
    vx1.attrs["PartNum"] = Npart
    vy1 = part1.create_dataset("vy", data=cms2kms(v1))
    vy1.attrs["Units"] = "km/s"
    vy1.attrs["PartNum"] = Npart
    vz1 = part1.create_dataset("vz", data=cms2kms(v2))
    vz1.attrs["Units"] = "km/s"
    vz1.attrs["PartNum"] = Npart      

    # masses
    
    mass = gr210MSun(read_dataset_dm_mass())
    mass = mass[sub_dark]
    mass1 = part1.create_dataset("Mass", data=mass)
    mass1.attrs["Units"] = "1e10 Msol"
    mass1.attrs["PartNum"] = Npart

## star(hf, sub_gas, nfiles= 16):
#
# ## Creates hdf5 file considering star attributes and filtrers by sub_star
def star_att(hf, sub_star, nfiles=16):

    part4 = hf.create_group("PartType4")
    
    # positions
    coord = read_dataset(4, "Coordinates", nfiles=16)
    coord = np.transpose(coord)
    c0 = coord[0][sub_star]
    c1 = coord[1][sub_star]
    c2 = coord[2][sub_star]
    Npart = len(c0)
    r = np.multiply(c0,c0)+np.multiply(c1,c1)+np.multiply(c2,c2)
    r = cm2kpc(np.sqrt(r))
    print('rs = ', min(r), max(r), np.median(r))
    print(Npart)
    x4 = part4.create_dataset("x", data=cm2kpc(c0))
    x4.attrs["Units"] = "kpc"
    x4.attrs["PartNum"] = Npart
    y4 = part4.create_dataset("y", data=cm2kpc(c1))
    y4.attrs["Units"] = "kpc"
    y4.attrs["PartNum"] = Npart
    z4 = part4.create_dataset("z", data=cm2kpc(c2))
    z4.attrs["Units"] = "kpc"
    z4.attrs["PartNum"] = Npart  
    
    # velocities
    vel = read_dataset(4, "Velocity", nfiles=16)
    vel = np.transpose(vel)
    v0 = vel[0][sub_star]
    v1 = vel[1][sub_star]
    v2 = vel[2][sub_star]
    print('vel ', v2.shape)
    vx4 = part4.create_dataset("vx", data=cms2kms(v0))
    vx4.attrs["Units"] = "km/s"
    vx4.attrs["PartNum"] = Npart
    vy4 = part4.create_dataset("vy", data=cms2kms(v1))
    vy4.attrs["Units"] = "km/s"
    vy4.attrs["PartNum"] = Npart
    vz4 = part4.create_dataset("vz", data=cms2kms(v2))
    vz4.attrs["Units"] = "km/s"
    vz4.attrs["PartNum"] = Npart      
    
    # masses
    mass = read_dataset(4, "Mass", nfiles=16)
    mass = mass[sub_star]
    print('mass ', mass.shape)
    mass4 = part4.create_dataset("Mass", data=gr210MSun(mass))
    mass4.attrs["Units"] = "1e10 Msol"
    mass4.attrs["PartNum"] = Npart   
    
    # initial masses
    
    print('hola')
    im = read_dataset(4, "InitialMass", nfiles=16)
    im = im[sub_star]
    imass4 = part4.create_dataset("InitialMass", data=gr210MSun(im))
    mass4.attrs["Units"] = "1e10 Msol"
    mass4.attrs["PartNum"] = Npart 

    # ages
    a = read_dataset(4, "StellarFormationTime", nfiles=16)
    a = a[sub_star]
    age4 = part4.create_dataset("StellarFormationTime", data=a)
    print(len(age4))
    print('hola2')
    age4.attrs["Units"] = "Expansion factor, a, where star particle forms"
    age4.attrs["PartNum"] = Npart
    
    # metallicities
    
    met = read_dataset(4, "Metallicity", nfiles=16)
    met = met[sub_star]    
    Z4 = part4.create_dataset("Metallicity", data=met)
    Z4.attrs["Units"] = "Smoothed mass fraction of elements heavier than Helium"
    Z4.attrs["PartNum"] = Npart
    print('star final')


### EAGLE does not use particle types Disc and Bulge

### MAIN
 
output_path = './subhalos/'
nfiles = 16
Subs_gas = read_dataset(0, 'SubGroupNumber', nfiles)
Subs_dark = read_dataset(1, 'SubGroupNumber', nfiles)
Subs_star = read_dataset(4, 'SubGroupNumber', nfiles)
Grs_gas = read_dataset(0, 'GroupNumber', nfiles)
Grs_dark = read_dataset(1, 'GroupNumber', nfiles)
Grs_star = read_dataset(4, 'GroupNumber', nfiles)

# test subhalo 0

# file
hf = h5py.File(output_path + 'gr_1_sub_0.hdf5', "w")

# gas particles
sub_gas = np.logical_and(Subs_gas == 0, Grs_gas == 1)
f = np.ones(len(Subs_gas))
N_gas = len(f[sub_gas])
print('hola gas ', N_gas)
if N_gas > 0:
    gas_att(hf, sub_gas, nfiles)
    
# dark matter particles
sub_dark = np.logical_and(Subs_dark == 0, Grs_dark == 1)
f = np.ones(len(Subs_dark))
N_dark = len(f[sub_dark])
print('hola dark ', N_dark)
if N_dark > 0:
    DM_att(hf, sub_dark, nfiles)
    
# star particles
sub_star = np.logical_and(Subs_star == 0, Grs_star == 1)
f = np.ones(len(Subs_star))
N_star = len(f[sub_star])
print('hola star ', N_star)
if N_star > 0:
    star_att(hf, sub_star, nfiles)
#    print('star final')

hf.close()
