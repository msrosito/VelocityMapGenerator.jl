using SimSpin, FITSIO, Plots

# Define telescope
telescope = SAMI(filter="r") 

# Define environment
redshift = 0.01
inclination = 70. # degrees
inc = trunc(Int, inclination)
virial_radius = 200 #kpc
environment = Environment(redshift, inclination, virial_radius) 
# mass-to-light ratio is an optional argument and should be used when ssp=false


input_dir="./galaxies/"
output_dir="./vLOS_RefL0012N0188/"

ngal=16 #number of files

# Load data
for k in 0:(ngal-1)
    if (k!=1) || (k!=7) || (k!=15) #replace for a more general condition
        data = sim_data(string(input_dir,"test_",k,".hdf5"), ssp=true)

        # Calculate datacube of fluxes regarding space and velocity bins
        datacube, observe = build_datacube(data, telescope, environment)
        #sim_FITS(datacube, observe, "test.fits")

        i = inclination * pi /180
        vz = [ d.vz for d in data ]
        vy = [ d.vy for d in data ]
        v_los = vz * cos(i) .- vy * sin(i) #computation of v_LOS
        vmax = maximum(v_los)
        vmin = minimum(v_los)
        nbins = size(datacube)[3]
        deltav = (vmax - vmin) / (nbins - 1)
        vels = [ vmin + (j-1) * deltav for j in 1:nbins]
        @show vels

        m = size(datacube)[1]
        n = size(datacube)[2]
        velmap = zeros(m, n)
        #for each pixel computes the flux weighted average v_LOS in the datacube
        for i in 1:m, j in 1:n
            velmap[i, j] = sum(datacube[i, j, :] .* vels) / sum(datacube[i, j, :])
        end
        heatmap(velmap) #, aspect_ratio=:equal)
        savefig(string(output_dir, "i",inc,"/gal_",k))
    end
end
