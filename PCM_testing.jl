using Pkg; Pkg.activate(".")
using NGNMM_NSP_paper_code
using ComponentArrays, Plots



save_path="./"
save_traj=false

# set parameters
u0=ComponentArray(x1=0.0, x2=0.0)
time_range=(0.0,10.0)
phoneme_sampling_rate=44100
drive_amplitude=20.0 #will be updated in the test function to make drive amplitude equal to DAR setting.
PCM_range=(5.5,10.0)
freq_range=(2.0,15.0) #Hz
stimulus_type="peakenvelope"

#looping over alpha to see how this affects the PCM vectors. It should, as the time scale of
#the evoked response will determine what happens at faster frequencies? /how much of the stimulus
#contains a response etc... overlap may occur, and phase difference may come out different.
#this is expected overall, the paper result shows this. and it is the key difference to the
#oscillator model.
alphas=[1/0.01,1/0.03,1/0.1,1/0.3,1/1.0]
all_pcm_vectors=Vector{Vector{ComplexF64}}()
for alpha in alphas
    α=alpha
    p=ComponentArray(α=α, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,noise_case_reference=1)
    PCM_vectors=get_evoked_model_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
    push!(all_pcm_vectors,PCM_vectors)
end

all_pcm_vectors_matrix=hcat(all_pcm_vectors...)
PCM_plots=[]
for i in 1:length(alphas)
    p=scatter(all_pcm_vectors_matrix[:,i],label="α=$(round(alphas[i],digits=2))s",legend=:outerright,xlims=(-1,1),ylims=(-1,1),aspect_ratio=1)
    push!(PCM_plots, p)
end
plot(PCM_plots...,layout=(2,3),size=(1800,800),title="Evoked Model PCM Vectors across 50 Frequencies, varying α")

#add some other alphas
lower_alphas=[1/0.005,1/0.002,1/0.05]
for alpha in lower_alphas
    α=alpha
    p=ComponentArray(α=α, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,noise_case_reference=1)
    PCM_vectors=get_evoked_model_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
    push!(all_pcm_vectors,PCM_vectors)
end
all_pcm_vectors_matrix=hcat(all_pcm_vectors...)
PCM_plots=[]
alphas=vcat(alphas,lower_alphas)
for i in 1:length(all_pcm_vectors)
    p=scatter(all_pcm_vectors_matrix[:,i],label="α=$(round(alphas[i],digits=3))s",legend=:outerright,xlims=(-1,1),ylims=(-1,1),aspect_ratio=1)
    push!(PCM_plots, p)
end
plot(PCM_plots...,layout=(3,3),size=(1800,1200))