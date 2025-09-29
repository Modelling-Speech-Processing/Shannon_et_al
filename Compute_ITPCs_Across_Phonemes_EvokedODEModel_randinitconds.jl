
println("script started pre loading project and packages")
using Pkg; Pkg.activate("$(pwd())")
println("project activated")
flush(stdout)
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

drive_amplitude_ratio_index=parse(Int64,ARGS[1])
drive_amplitude_ratios=[0.1,0.496,1.0,2.46,12.2,60.5,100.0]
drive_amplitude_ratio=drive_amplitude_ratios[drive_amplitude_ratio_index]


noise_stimulus_ratio=parse(Float64,ARGS[2])

# save_path=ARGS[3]
save_path="./Results/Evoked_Model/"

#running all noise cases as a single 60 trial batch
noisesets=["NoiseSequences60.csv","NoiseSequences60_new.csv","NoiseSequences60_3.csv"]
noise_ref="allnoises"

#path to time to peak derivative data
time2PD_path="./"
peak_derivative_data=readdlm(time2PD_path*"mean_times_to_peak_deriv.csv",Float64)[:,1]
sorting_indices=get_sorting_indices(time2PD_path)

#load phoneme stimuli sequence envelopes
# phoneme_sequence_envelopes=JSON.parsefile("PhonemeEnvelopes_allconditions.json")
phoneme_sequence_envelopes=Vector{Float32}(undef,5) #dont need the envelopes anymore, but will keep in the code. (have saved the stimuli to disk instead of creating on the fly.)

# set parameters
u0=ComponentArray(x1=0.0, x2=0.0)
time_range=(0.0,10.0)
phoneme_sampling_rate=44100
drive_amplitude=4.0 #will be updated in the test function to make drive amplitude equal to DAR setting.
α=1/0.03
p=ComponentArray(α=α, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,noise_case_reference=1)
ITPCrange=(5.5,10.0)


#choose stimulus type
# stimulus_type="envelope"
# stimulus_type="derivative"
# stimulus_type="peakrate"
# stimulus_type="peakenvelope"

for stimulus_type in ["envelope","peakenvelope","derivative","peakrate"]
println("Running tests for stimulus type: ",stimulus_type)
    flush(stdout)
    t=time()
    run_evoked_model_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond("SNR05_high_c_60randinitcond_seeded_noise_evokedODE_stimulustype_$(stimulus_type)_test_$(noise_ref)",[drive_amplitude_ratio],sorting_indices,peak_derivative_data,p,u0,time_range,ITPCrange,noise_stimulus_ratio,phoneme_sequence_envelopes,save_path,noisesets,false,stimulus_type)
    println("Time taken: ",time()-t)
end