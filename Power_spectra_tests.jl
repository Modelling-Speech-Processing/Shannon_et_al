# to see if power in stimulus band increases (evoked respone) or not (phase resetting).
using Pkg; Pkg.activate(".") 

#load example stimulus into global interpolator, with NSR=0.0

#zero out the first 5 seconds of noise

#run 20 simulations for each model. and plot mean spectrogram.
#should see power increase (or not) for the evoked and phase-resetting models in the 4Hz band.
#intereted in what the NGNMM response is.