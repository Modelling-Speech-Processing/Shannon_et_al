
println("script started pre loading project and packages")

using Pkg; Pkg.activate("$(pwd())")
println("project activated")
flush(stdout)
Pkg.instantiate()
println("project instantiated")
flush(stdout)
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0  # Disable automatic precompilation
println("precompilation disabled")
flush(stdout)
using NGNMM_NSP_paper_code
println("NGNMM_NSP_paper_code loaded")
flush(stdout)
using ComponentArrays, Plots, JSON, WAV, CSV, Statistics, DelimitedFiles, Random
println("packages loaded")
flush(stdout)
println("Threads: ",Base.Threads.nthreads())
flush(stdout)


drive_amplitude_ratio=0.1# will be reset to give max stimulus = 1.0. no longer sweeping this? at least for now. #,0.496,2.46,12.2,60.5,100.0]


noise_stimulus_ratio=0.1





phase_modulation_idx=parse(Int32,ARGS[1]) 
phase_modulations=collect(range(start=0.0,stop=1.0,length=11))
# phase_modulations=collect(range(start=0.0,stop=1.0,length=51))
phase_modulation=phase_modulations[phase_modulation_idx]

# save_path=ARGS[3]
save_path="./Results/Phase_Reset_Model/"

#running all noise cases as a single 60 trial batch
noisesets=["NoiseSequences60.csv","NoiseSequences60_new.csv","NoiseSequences60_3.csv"]
noise_ref="allnoises"

#path to time to peak derivative data
time2PD_path="./"
peak_deriviative_data=readdlm(time2PD_path*"mean_times_to_peak_deriv.csv",Float64)[:,1]
sorting_indices=get_sorting_indices(time2PD_path)

#load phoneme stimuli sequence envelopes
# phoneme_sequence_envelopes=JSON.parsefile("PhonemeEnvelopes_allconditions.json")
phoneme_sequence_envelopes=Vector{Float32}(undef,5) #dont need the envelopes anymore, but will keep in the code. (have saved the stimuli to disk instead of creating on the fly.)

# set parameters
u0=ComponentArray(θ=0.0, r=1.0)
c=1.0 #will be updated to 0.7*pi*(1/(maximum(interpolators_global[1])*drive_amplitude)) given a particular test condition,.
time_range=(0.0,10.0)
F=4.0
phoneme_sampling_rate=44100
drive_amplitude=4.0 #will be updated in the test function to make the max value of the stimulus envelope (across noise cases) equal to 1.0
#for rectified version
# q_normalisation=4/(2*pi-2*pi*phase_modulation+2*phase_modulation) #scaling factor to make absolute area under stimulus modulation curve = 4. 
#for non rectified
if phase_modulation<=0.5
    q_normalisation=2/((1-phase_modulation)*π) #scaling factor to make absolute area under stimulus modulation curve = 4. 
else
    m=phase_modulation
    q_normalisation=4/(4*m*sqrt(1-((m-1)/m)^2)+2*(m-1)asin((m-1)/m))
end

p=ComponentArray(F=F, c=c, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,modulation=phase_modulation,q=q_normalisation,noise_case_reference=1)
ITPCrange=(5.5,10.0)

t=time()
run_coupled_oscillator_modulated_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond("SNR05_high_c_60randinitcond_seeded_noise_normalised_nonrectified_modulated_phasemod_$(phase_modulation)_integralnormalised_coupled_oscillator_control_test_$(noise_ref)",[drive_amplitude_ratio],sorting_indices,peak_deriviative_data,p,u0,time_range,ITPCrange,noise_stimulus_ratio,phoneme_sequence_envelopes,save_path,noisesets)
println("Time taken: ",time()-t)