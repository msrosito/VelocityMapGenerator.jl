using SimSpin, FITSIO, Plots

telescope = SAMI()
environment = Environment(0.05, 70, 200, 1.)
data = sim_data()
datacube, observe = build_datacube(data, telescope, environment)
sim_FITS(datacube, observe, "SimSpin_Example_Observation.fits")
f = FITS("SimSpin_Example_Observation.fits")
data = read(f[1])
heatmap(data[:,:,2])
savefig("SimSpin_Example_Observation_2.jpg")
