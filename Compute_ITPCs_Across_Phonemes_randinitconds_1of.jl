println("script started pre loading project and packages")

using Pkg; Pkg.activate("$(pwd())")
println("project activated")
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
drive_amplitude_ratios=[0.1,0.496,2.46,12.2,60.5,100.0]
drive_amplitude_ratio=drive_amplitude_ratios[drive_amplitude_ratio_index]

noise_stimulus_ratio=parse(Float64,ARGS[2])
# noise_stimulus_ratio=0.7

# save_path=ARGS[3]
save_path="./Results/NMM/"

# noise_filename="NoiseSequences60"
#individual noise cases
noise_filename=ARGS[4]
noise_filename=noise_filename*".csv"
#running all noise cases as a single 60 trial batch
noisesets=["NoiseSequences60.csv","NoiseSequences60_new.csv","NoiseSequences60_3.csv"]
noise_ref="allnoises"
#path to time to peak derivative data
time2PD_path="./"
peak_deriviative_data=readdlm(time2PD_path*"mean_times_to_peak_deriv.csv",Float64)[:,1]
sorting_indices=get_sorting_indices(time2PD_path)


#load phoneme stimuli sequence envelopes
# phoneme_sequence_envelopes=JSON.parsefile("PhonemeEnvelopes_allconditions.json")
phoneme_sequence_envelopes=Vector{Float32}(undef,5) #dont need the envelopes anymore, but will keep in the code. (have saved the stimuli to disk instead of creating online.)

#Slow (C0066D025) run
# set parameters
α=1/0.035 #1/a is 30ms
k=0.105
# C = 0.033#capacitance in Farads
C = 0.066#capacitance in Farads
vsyn=-10.0
# η_0=21.5
η_0=1.5
# Δ=0.5
Δ=0.25
α_D=1/0.0056
Π=15.0
#Drive 
phoneme_sampling_rate=44100
drive_amplitude=Π
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1) 
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,20.0)
ITPCrange=(5.5,10.0)


for stimulus_type in ["envelope"]
    println("Running tests for stimulus type: ",stimulus_type)
    flush(stdout)
    run_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond_1overf_noise("1overf_Randinitconds_c066D025LowFreqTest$(noise_ref)_stimulustype_$(stimulus_type)",[drive_amplitude_ratio],sorting_indices,peak_deriviative_data,p,u0,time_range,ITPCrange,noise_stimulus_ratio,phoneme_sequence_envelopes,save_path,noisesets,stimulus_type)
    t=time()
    println("Time taken: ",time()-t)
end

#fast run (Aine parameters)
# set parameters
α=1/0.035 #1/a is 30ms
k=0.105
C = 0.033#capacitance in Farads
# C = 0.066#capacitance in Farads
vsyn=-10.0
η_0=21.5
# η_0=1.5
Δ=0.5
# Δ=0.25
α_D=1/0.0056
Π=15.0
#Drive 
phoneme_sampling_rate=44100
drive_amplitude=Π
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1) 
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,20.0)
ITPCrange=(5.5,10.0)

for stimulus_type in ["envelope"]
    println("Running tests for stimulus type: ",stimulus_type)
    flush(stdout)
    run_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond_1overf_noise("1overf_Randinitconds_AineRerunTest$(noise_ref)_stimulustype_$(stimulus_type)",[drive_amplitude_ratio],sorting_indices,peak_deriviative_data,p,u0,time_range,ITPCrange,noise_stimulus_ratio,phoneme_sequence_envelopes,save_path,noisesets,stimulus_type)
    t=time()
    println("Time taken: ",time()-t)
end
