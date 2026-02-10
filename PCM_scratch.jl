
using NGNMM_NSP_paper_code

frequencies=range(2,15,length=60)
data=reduce(hcat,[sin.(2π*freq.*(0:0.001:10)) for freq in frequencies])

#stimulus is the clean data plus offset that is constant fraction of a period to test PCM with constant phase lag.
#this is the case expected for a purely phase resetting oscillator.
stimuli=Matrix{Float64}(undef,10001,60)
for (i,freq) in enumerate(frequencies)
    period=1.0/freq
    offset=period*0.5
    stimuli[:,i]=sin.(2π*freq.*((0:0.001:10) .+ offset)) 
end

#add noise to data:
data=data.+0.25*randn(size(data))

num_trials=size(data,2)
sr=1/0.001
PCMrange=(5.0,10.0) #seconds
mean_phase_diff_vectors=Vector{ComplexF64}(undef,num_trials)
for trial_idx in 1:num_trials
        #get the stimulus for this trial:
        stimulus_interpolator=stimuli[:,trial_idx]
        #get the peak frequency of this stimulus:
        stimulus_times=range(0.0,10.0,length=Int(sr*10))
        stimulus_values=stimulus_interpolator
        peak_frequency=frequencies[trial_idx]
        flush(stdout)
        #filter and get instantaneous phase:
        filtered_response=gaussian_filter_and_hilbert_for_PCM(data[:,trial_idx][Int64(PCMrange[1]*sr):Int64(PCMrange[2]*sr)],sr,peak_frequency)
        filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(stimulus_values[Int64(PCMrange[1]*sr):Int64(PCMrange[2]*sr)],sr,peak_frequency)
        #get phases
        instantaneous_response_phase=angle.(filtered_response)
        instantaneous_stimulus_phase=angle.(filtered_stimulus)
        #get phase diff:
        phase_diff=instantaneous_response_phase.-instantaneous_stimulus_phase
        #convert to complex unit vector form:
        phase_diff_vectors=exp.(im.*phase_diff)
        #get average phase difference over time
        mean_phase_diff_vector=mean(phase_diff_vectors)
        mean_phase_diff_vectors[trial_idx] = mean_phase_diff_vector 
   end
#check resulting vectors:
scatter(mean_phase_diff_vectors,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal)
#works. e.g. for a constant 0.5period offset, all vectors cluster round -1+0im as expected.
#and for 0.25 period offset, they cluster round 0-1im as expected.

#now with a constant offset, un tuned for the period. This is the case expected for an evoked response.
constant_offset_stimuli=Matrix{Float64}(undef,10001,60)
for (i,freq) in enumerate(frequencies)
    offset=0.1 #seconds
    constant_offset_stimuli[:,i]=sin.(2π*freq.*((0:0.001:10) .+ offset)) 
end
mean_phase_diff_vectors_constant_offset=Vector{ComplexF64}(undef,num_trials)
for trial_idx in 1:num_trials
        #get the stimulus for this trial:
        stimulus_interpolator=constant_offset_stimuli[:,trial_idx]
        #get the peak frequency of this stimulus:
        stimulus_times=range(0.0,10.0,length=Int(sr*10))
        stimulus_values=stimulus_interpolator
        peak_frequency=frequencies[trial_idx]
        flush(stdout)
        #filter and get instantaneous phase:
        filtered_response=gaussian_filter_and_hilbert_for_PCM(data[:,trial_idx][Int64(PCMrange[1]*sr):Int64(PCMrange[2]*sr)],sr ,peak_frequency)
        filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(stimulus_values[Int64(PCMrange[1]*sr):Int64(PCMrange[2]*sr)],sr,peak_frequency)
        #get phases
        instantaneous_response_phase=angle.(filtered_response)
        instantaneous_stimulus_phase=angle.(filtered_stimulus)
        #get phase diff:
        phase_diff=instantaneous_response_phase.-instantaneous_stimulus_phase
        #convert to complex unit vector form:
        phase_diff_vectors=exp.(im.*phase_diff)
        #get average phase difference over time
        mean_phase_diff_vector=mean(phase_diff_vectors)
        mean_phase_diff_vectors_constant_offset[trial_idx] = mean_phase_diff_vector 
   end
#check resulting vectors:
scatter(mean_phase_diff_vectors_constant_offset,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal)
#looks good - they are evenly spaced around the unit circle.
#which is expected as the constant offset corresponds to a different phase lag for each frequency,
#so for uniformly spaced frequencies we get uniform spacing of phase differences, overlap aside.


save_path="./"
save_traj=true

# set parameters
u0=ComponentArray(x1=0.0, x2=0.0)
time_range=(0.0,10.0)
phoneme_sampling_rate=44100
drive_amplitude=20.0 #will be updated in the test function to make drive amplitude equal to DAR setting.
α=1/0.03
p=ComponentArray(α=α, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,noise_case_reference=1)
PCM_range=(5.5,10.0)
freq_range=(2.0,15.0) #Hz
stimulus_type="peakenvelope"

PPP=get_evoked_model_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)


scatter(PPP,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="Evkd - peak envelope stim",marker_z=frequencies)

traj_path="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/evoked_model_response_to_sine_stimulus_freq_range_2.0to15.0Hz_envelope_transformed.arrow.lz4"
impulse_traj_path="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/evoked_model_response_to_sine_stimulus_freq_range_2.0to15.0Hz_peakenvelope_transformed.arrow.lz4"
using Arrow, DataFrames
trajectories=Matrix(DataFrame(Arrow.Table(traj_path)))
impulse_driven_trajectories=Matrix(DataFrame(Arrow.Table(impulse_traj_path)))
#load sine drive 1 data from path:
sine_path="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_1.jld2"
using JLD2
sine_drive_data=jldopen(sine_path,"r")
sine_drive_interpolators=sine_drive_data["drives"]
sine_impulses=get_scaled_peakenv_impulses_scratch(sine_drive_interpolators,44100.0)

traj_times=range(5.5,10.0,length=size(trajectories,1))
plot(traj_times,trajectories[:,1], label="EVKD sine stim.")
plot!(traj_times,impulse_driven_trajectories[:,1]*3000, label="EVKD impulses stim. *3000") #scaled for visibility
traj_idx=range(5.5*44100,10.0*44100,step=1)
traj_times=traj_idx./44100.0
traj_times=range(5.5,10.0,length=size(trajectories,1))
stimulus_indexes=traj_times.*44100
stimulus_times=stimulus_indexes./44100
plot!(stimulus_times,sine_drive_interpolators[1](stimulus_indexes),label="sine stim drive")
plot!(xlims=(5.5,8.5))
#plot smoothed response phase:
plot!(traj_times,angle.(gaussian_filter_and_hilbert_for_PCM(trajectories[:,1],44100,frequencies[1])),label="EVKD sine stim phase")
plot!(traj_times,angle.(gaussian_filter_and_hilbert_for_PCM(impulse_driven_trajectories[:,1],44100,frequencies[1])),label="EVKD impulses stim phase")
#plot stimulus phase:
plot!(stimulus_times,angle.(gaussian_filter_and_hilbert_for_PCM(sine_drive_interpolators[1](stimulus_indexes),44100,frequencies[1])),label="sine stim drive phase")
plot!(legend=:outertopright,size=(1300,600),xlims=(5.5,8.5))

#same again but for the 15th frequency.
plot(traj_times,trajectories[:,15], label="EVKD sine stim.")
plot!(traj_times,impulse_driven_trajectories[:,15]*3000, label="EVKD impulses stim. *3000") #scaled for visibility
plot!(stimulus_times,sine_drive_interpolators[15](stimulus_indexes),label="sine stim drive")
plot!(stimulus_times,sine_impulses[15](stimulus_indexes),label="sine impulse stim drive",linewidth=3,color=:black)
plot!(xlims=(5.5,8.5))
#plot smoothed response phase:
plot!(traj_times,angle.(gaussian_filter_and_hilbert_for_PCM(trajectories[:,15],44100,frequencies[15])),label="EVKD sine stim phase")
plot!(traj_times,angle.(gaussian_filter_and_hilbert_for_PCM(impulse_driven_trajectories[:,15],44100,frequencies[15])),label="EVKD impulses stim phase")
#plot stimulus phase:
plot!(stimulus_times,angle.(gaussian_filter_and_hilbert_for_PCM(sine_drive_interpolators[15](stimulus_indexes),44100,frequencies[15])),label="sine stim drive phase")
plot!(legend=:outertopright,size=(1300,600),xlims=(5.5,8.5))
plot!(xlims=(5.5,7.75))
plot(sine_impulses[15](stimulus_indexes),label="sine impulse stim drive")

# phase differences between response and stimulus for sine driven case:
response_phase=angle.(gaussian_filter_and_hilbert_for_PCM(trajectories[:,1],44100,frequencies[1]))
stimulus_phase=angle.(gaussian_filter_and_hilbert_for_PCM(sine_drive_interpolators[1](stimulus_indexes),44100,frequencies[1]))
phase_diff=response_phase.-stimulus_phase
plot(traj_times,phase_diff,label="phase diff EVKD sine stim.")
# phase differences between response and stimulus for impulse driven case:
response_phase_impulse=angle.(gaussian_filter_and_hilbert_for_PCM(impulse_driven_trajectories[:,1],44100,frequencies[1]))
phase_diff_impulse=response_phase_impulse.-stimulus_phase
plot!(traj_times,phase_diff_impulse,label="phase diff EVKD impulses stim.")
#looks very similar.

#now between freq 1 and 15
response_phase15=angle.(gaussian_filter_and_hilbert_for_PCM(trajectories[:,15],44100,frequencies[15]))
stimulus_phase15=angle.(gaussian_filter_and_hilbert_for_PCM(sine_drive_interpolators[15](stimulus_indexes),44100,frequencies[15]))
phase_diff15=response_phase15.-stimulus_phase15
plot(traj_times,phase_diff15,label="phase diff EVKD sine stim. freq 15")
plot!(traj_times,phase_diff,label="phase diff EVKD sine stim. freq 1")

#computing unit vectors in complex plane for phase diffs:
phase_diff_vectors=exp.(im.*phase_diff)
phase_diff_15_vectors=exp.(im.*phase_diff15)
scatter(phase_diff_vectors,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal,label="freq 1")
scatter!(phase_diff_15_vectors,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal,label="freq 15")

#now for impulse case:
stimulus_phase_impulse15=angle.(gaussian_filter_and_hilbert_for_PCM(sine_impulses[15](stimulus_indexes),44100,frequencies[15]))
response_phase_impulse15=angle.(gaussian_filter_and_hilbert_for_PCM(impulse_driven_trajectories[:,15],44100,frequencies[15]))
phase_diff_impulse15=response_phase_impulse15.-stimulus_phase_impulse15
plot(traj_times,phase_diff_impulse15,label="phase diff EVKD impulses stim. freq 15")
plot!(traj_times,phase_diff_impulse,label="phase diff EVKD impulses stim. freq 1")

#computing unit vectors in complex plane for phase diffs:
phase_diff_impulse_vectors=exp.(im.*phase_diff_impulse)
phase_diff_impulse_15_vectors=exp.(im.*phase_diff_impulse15)
scatter(phase_diff_impulse_vectors,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal,label="freq 1")
scatter!(phase_diff_impulse_15_vectors,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal,label="freq 15")

#double checking PCM scatter plot vs manually reading the phase diffs.
mean_phase_diff_vectors_impulse=Vector{ComplexF64}(undef,20)
for i in 1:20
    response_phase=angle.(gaussian_filter_and_hilbert_for_PCM(impulse_driven_trajectories[:,i],44100,frequencies[i]))
    stimulus_phase=angle.(gaussian_filter_and_hilbert_for_PCM(sine_impulses[i](stimulus_indexes),44100,frequencies[i]))
    phase_diff= response_phase.-stimulus_phase
    phase_diff_vectors=exp.(im.*phase_diff)
    mean_phase_diff_vector=mean(phase_diff_vectors)
    mean_phase_diff_vectors_impulse[i]=mean_phase_diff_vector
    println("Freq $(frequencies[i]) Hz: PCM from manual calc: $(mean_phase_diff_vector), PCM from function: $(PPP[i])")
end
#compared to function output:
scatter(PPP,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal,label="from function",title= "impulse driven, evkd pcm",alpha=0.7)
scatter!(mean_phase_diff_vectors_impulse,xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal,label="manual calc",title= "impulse driven, evkd pcm",alpha=0.7)

#resample to match trajectory time points:
sine_drive_interpolators[1]
stimulus_times=range(5.5*44100,10.0*44100,length=Int(44100*(10.0-5.5)))
stimulus_values=sine_drive_interpolators[1](stimulus_times)
plot(stimulus_values)
stimulus_values[Int64(5.5):Int64(10.0*10000)]
plot(stimulus_values[Int64(5.5*10000):Int64(10.0*10000)])
scatter!(sine_drive_interpolators[1](1:1:10))
plot(stimulus_times,sine_drive_interpolators[1])
sine_drive_interpolators[1]

length(stimulus_times)
plot(traj_times,trajectories[:,1])
plot!(stimulus_times,sine_drive_interpolators[1])
using Interpolations
function get_scaled_peakenv_impulses_scratch(envelopes::Vector{Any}, sampling_rate::Float64; scale_range=(1.0, 1.0))
    sample_grid= envelopes[1].itp.itp.parentaxes[1][1:end]  # Get the time grid from the first envelope
    peakenv_impulses = Vector{Any}(undef, length(envelopes))
    for (i, envelope) in enumerate(envelopes)
        # Create a zero vector of the same length as the envelope
        peaks,locations = findpeaks(envelope[:], 0.1)  # Threshold set to 0.1

        # Scale the peak magnitudes to the specified range
        if !isempty(peaks)
            min_peak = minimum(peaks)
            max_peak = maximum(peaks)
            scaled_peaks = scale_range[1] .+ (peaks .- min_peak) .* (scale_range[2] - scale_range[1]) ./ (max_peak - min_peak)
        else
            scaled_peaks = Float64[]
        end
        impulse_vector = zeros(length(sample_grid))
        # Set the peak locations to 1.0
        for (peak, loc) in zip(scaled_peaks, locations)
            impulse_vector[loc] = peak  # Set peak rate impulse to the value of the scaled envelope amplitude at the peak location
        end
        # Interpolate the impulse vector
        peakenv_impulses[i] = linear_interpolation(StepRange(sample_grid), impulse_vector, extrapolation_bc=Line());
    end

    return peakenv_impulses
end

stim_impulses=get_scaled_peakenv_impulses_scratch(sine_drive_interpolators,44100.0)
plot!(stimulus_times,stim_impulses[1])


#plot the response and stimulus impulses for each of the 20 frequencies.
ps=[]
for freq_idx in 1:20
    p=plot(traj_times,trajectories[:,freq_idx],xlabel="Time (s)")
    plot!(p,stimulus_times,sine_drive_interpolators[freq_idx],label="Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))
    push!(ps,p)
end

plot(ps[1:2:end]...,layout=(5,2),size=(1200,800))
sine_drive_interpolators
#plot the same but filtered rates and signal, to see what the PCM was calculated from
sine_drive_interpolators[1](stimulus_times)
trajectories[:,1]
ps_filtered=[]
for freq_idx in 1:20
    filtered_rate=gaussian_filter_and_hilbert_for_PCM(trajectories[:,freq_idx],1/0.0001,frequencies[freq_idx])
    filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(sine_drive_interpolators[freq_idx](stimulus_times),44100,frequencies[freq_idx])
    p=plot(traj_times,abs.(filtered_rate),xlabel="Time (s)")
    plot!(p,stimulus_times./44100,abs.(filtered_stimulus),label="Filtered Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))   
    push!(ps_filtered,p)
end

plot(ps_filtered[1:2:end]...,layout=(5,2),size=(1200,800))


#and the phase angle of each signal. amplitude is as expected. pure sine wave gives abs value of 1
#and the pure sine wave plus noise gets a higher abs value, as noise adds power around the peak frequency
#this is still picked up by the gauss_kernel.
ps_phase=[]
phase_diffs=[]

for freq_idx in 1:20
    filtered_rate=gaussian_filter_and_hilbert_for_PCM(trajectories[:,freq_idx],1/0.0001,frequencies[freq_idx])
    filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(sine_drive_interpolators[freq_idx](stimulus_times),44100,frequencies[freq_idx])
    phase_diff=angle.(filtered_rate)-angle.(filtered_stimulus)
    push!(phase_diffs,phase_diff)
    p=plot(traj_times,angle.(filtered_rate),xlabel="Time (s)")
    plot!(p,stimulus_times./44100,angle.(filtered_stimulus),label="Filtered Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))   
    push!(ps_phase,p)
end

plot(ps_phase[1:2:end]...,layout=(5,2),size=(1200,800))
#looks like constant phase lag?

#cant get the phase diffs here as different sample count.
#so re sampled the stimulus before filtering to match the trajectory time points.
#this means taking its long 44100 hz signal, and using the interpolator to pull out values at the trajectory time points.
traj_times=range(5.5,10.0,length=size(trajectories,1))
stim_idxs=traj_times.*44100
ps_phase_resampled=[]
phase_diffs_resampled=[]
for freq_idx in 1:50
    filtered_rate=gaussian_filter_and_hilbert_for_PCM(trajectories[:,freq_idx],1/0.0001,frequencies[freq_idx])
    stim_values_resampled=sine_drive_interpolators[freq_idx](stim_idxs)
    filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(stim_values_resampled,1/0.0001,frequencies[freq_idx])
    phase_diff=angle.(filtered_rate)-angle.(filtered_stimulus)
    push!(phase_diffs_resampled,phase_diff)
    p=plot(traj_times,angle.(filtered_rate),xlabel="Time (s)")
    plot!(p,traj_times,angle.(filtered_stimulus),label="Filtered Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))   
    push!(ps_phase_resampled,p)
end

plot(ps_phase_resampled[1:2:end]...,layout=(5,2),size=(1200,800))
plot(phase_diffs_resampled[1:2:end],xlabel="Time (s)",layout=(5,2),xlims=(5.5,10.0),size=(1200,800))
phase_diffs_array=hcat(phase_diffs_resampled...)
mean_phase_diffs=mean(phase_diffs_array,dims=1)
bar(mean_phase_diffs')


#to compare to scatter produced by my main PCM function, will get phase diff complex vectors,
#then get circular average, as the mean above is incorrect, (not circular)
phase_diff_exponential=exp.(im.*phase_diffs_array)
mean_phase_diff_vectors_resampled=mean(phase_diff_exponential,dims=1)
scatter(mean_phase_diff_vectors_resampled[1,:],xlims=(-1.0,1.0),ylims=(-1.0,1.0),aspect_ratio=:equal)  
PCM=abs.(mean(mean_phase_diff_vectors_resampled[1,:]))



frequencies[1]
filtered_rate=gaussian_filter_and_hilbert_for_PCM(trajectories[:,1],1/0.0001,frequencies[1])
plot(traj_times,abs.(filtered_rate),xlabel="Time (s)")
plot!(traj_times,trajectories[:,1],xlabel="Time (s)")
filtered_stimulus=abs.(gaussian_filter_and_hilbert_for_PCM(sine_drive_interpolators[1](stimulus_idxs),1/0.0001,frequencies[1]))
plot(stimulus_times,filtered_stimulus)
plot!(stimulus_times,sine_drive_interpolators[1](stimulus_idxs))

plot(sine_drive_interpolators[1])

stimulus_times

stimulus_idxs=stimulus_times.*44100



#### checking gaussian kernel function:
sinewave=sin.(2π*5.0.*(0:0.0001:10.0)) #at 10000Hz sampling rate
filtered_sinewave=gaussian_filter_and_hilbert_for_PCM(sinewave,10000.0,1.0)
plot(sinewave)
plot!(abs.(filtered_sinewave))
plot!(angle.(filtered_sinewave))


#double checking impulses from the sine wave:
sine_impulses=get_scaled_peakenv_impulses_scratch(sine_drive_interpolators,44100.0)
plot(sine_impulses[20])



##checking evoked model response to this impulse train.

global interpolators_global=sine_impulses
prob_func=vary_noise_and_initial_conditions_evokedmodel
p=ComponentArray(α=1/0.03, drive_amplitude=40.0, noise_selector=1, sampling_rate=44100.0,noise_case_reference=1)
u_init=ComponentArray(x1=0.0, x2=0.0)
time_range=(0.0,10.0)
saveat=0.0001
trialdata=Ensemble_EvokedModel(prob_func,time_range,p,u_init,20,saveat)
for i in 1:20
    display(plot(trialdata[i][1,:],ylims=(-0.01,0.01)))
end





###testing the coupled oscillator
u0=ComponentArray(θ=0.0, r=1.0)
c=1.0 #will be updated to 0.7*pi*(1/(maximum(interpolators_global[1])*drive_amplitude)) given a particular test condition,.
time_range=(0.0,10.0)
F=4.0
phoneme_sampling_rate=44100
drive_amplitude=4.0
phase_modulation=1.0
if phase_modulation<=0.5
    q_normalisation=2/((1-phase_modulation)*π) #scaling factor to make absolute area under stimulus modulation curve = 4. 
else
    m=phase_modulation
    q_normalisation=4/(4*m*sqrt(1-((m-1)/m)^2)+2*(m-1)asin((m-1)/m))
end

#pass in frequencies used by the stimulus ,so the oscillator can be set to match them
frequencies=range(2,15,length=60)
p=ComponentArray(F=F, c=c, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,modulation=phase_modulation,q=q_normalisation,noise_case_reference=1,frequencies=frequencies)
PCM_range=(5.5,10.0)


save_path="./"
save_traj=true

freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

PCMs_coupled_oscillator=get_coupled_oscillator_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(PCMs_coupled_oscillator,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="sine wave stim - coupled oscillator",marker_z=frequencies)

frequencies=range(2,15,length=60)
frequencies[10]
abs(PCMs_coupled_oscillator[9])
abs(PCMs_coupled_oscillator[10])
#the 10th freq is about 4, it is the only one that syncrhonises well enough to give a high PCM.
abs(PCMs_coupled_oscillator[11])

#does the paper adjust the models frequency with the stimulus frequency?
#yes. 

#so the synchronisation is fixed with the prob func change. 
#but now the phase lag is all around the circle as the frequency changeds.

#so inspect the trajectories:
traj_path_coupled_oscillator="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/oscillator_response_to_sine_stimulus_freq_range_2.0to15.0Hz_peakenvelope_transformed.arrow.lz4"
high_c_traj_path_coupled_oscillator="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/highccheck_2_oscillator_response_to_sine_stimulus_freq_range_2.0to15.0Hz_peakenvelope_transformed.arrow.lz4"
using Arrow, DataFrames
trajectories_coupled_oscillator=Matrix(DataFrame(Arrow.Table(traj_path_coupled_oscillator)))
high_c_trajectories_coupled_oscillator=Matrix(DataFrame(Arrow.Table(high_c_traj_path_coupled_oscillator)))
#load sine drive 1 data from path:
sine_path1="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_1.jld2"
sine_path2="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_2.jld2"
sine_path3="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_3.jld2"
using JLD2
sine_drive_data1=jldopen(sine_path1,"r")
sine_drive_interpolators1=sine_drive_data1["drives"]
sine_drive_data2=jldopen(sine_path2,"r")
sine_drive_interpolators2=sine_drive_data2["drives"]
sine_drive_data3=jldopen(sine_path3,"r")
sine_drive_interpolators3=sine_drive_data3["drives"]
sine_drive_interpolators=vcat(sine_drive_interpolators1,sine_drive_interpolators2,sine_drive_interpolators3)
sine_impulses=get_unitary_peakenv_impulses(sine_drive_interpolators,44100.0) 
traj_times=range(5.5,10.0,length=size(trajectories_coupled_oscillator,1))
stimulus_indexes=traj_times.*44100
stimulus_times=stimulus_indexes./44100
#plot the response and stimulus for each of the 20 frequencies.
ps_coupled=[]
for freq_idx in 1:5:50
    p=plot(traj_times,trajectories_coupled_oscillator[:,freq_idx],xlabel="Time (s)")
    plot!(p,stimulus_times,sine_impulses[freq_idx](stimulus_indexes),label="Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))
    push!(ps_coupled,p)
end
plot(ps_coupled...,layout=(5,2),size=(1200,800))


ps_highc_coupled=[]
for freq_idx in 1:5:50
    p=plot(traj_times,high_c_trajectories_coupled_oscillator[:,freq_idx],xlabel="Time (s)")
    plot!(p,stimulus_times,sine_impulses[freq_idx](stimulus_indexes),label="Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))
    push!(ps_highc_coupled,p)
end
plot(ps_highc_coupled...,layout=(5,2),size=(1200,800),plot_title="high c coupled oscillator")


#and inspect the trajectories for the sine wave driven oscillator:
traj_path_coupled_oscillator_sine="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/oscillator_response_to_sine_stimulus_freq_range_2.0to15.0Hz_envelope_transformed.arrow.lz4"
trajectories_coupled_oscillator_sine=Matrix(DataFrame(Arrow.Table(traj_path_coupled_oscillator_sine)))
#plot the response and stimulus for each of the 20 frequencies.
ps_coupled_sine=[]
for freq_idx in 1:5:50
    p=plot(traj_times,trajectories_coupled_oscillator_sine[:,freq_idx],xlabel="Time (s)")
    plot!(p,stimulus_times,sine_drive_interpolators[freq_idx](stimulus_indexes),label="Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))
    push!(ps_coupled_sine,p)
end
plot(ps_coupled_sine...,layout=(5,2),size=(1200,800),plot_title="sine driven coupled oscillator")



### testing NGNMM with sine wave drives: of interest as we are not changing the frequency of the model.
## slow parameters (4Hz)
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
time_range=(0.0,10.0)
PCM_range=(5.5,10.0)

save_path="./"
save_traj=true

freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

PCMs_NGNMM=get_NGNMM_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(PCMs_NGNMM,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="sine wave stim - 4Hz NGNMM",marker_z=frequencies)

#and for the fast 19Hz NGNMM:
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
PCM_range=(5.5,10.0)

save_path="./"
save_traj=true

freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

PCMs_NGNMM_fast=get_NGNMM_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(PCMs_NGNMM_fast,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="sine wave stim - 19Hz NGNMM",marker_z=frequencies)



#plot all four cases on one figure: PPP (evoked impulse stim), coupled oscillator sine stim, NGNMM slow sine stim, NGNMM fast sine stim
p1=scatter(PPP,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="Evoked Model - Impulse Stim",marker_z=frequencies)
p2=scatter(PCMs_coupled_oscillator,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="Coupled Oscillator - Sine Stim",marker_z=frequencies)
p3=scatter(PCMs_NGNMM,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="NGNMM Slow - Sine Stim",marker_z=frequencies)
p4=scatter(PCMs_NGNMM_fast,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="NGNMM Fast - Sine Stim",marker_z=frequencies)
plot(p1,p2,p3,p4,layout=(2,2),size=(1000,800))


### run coupled oscillator case without frequency adjustment:
u0=ComponentArray(θ=0.0, r=1.0)
c=0.7*pi*(1/(maximum(sine_drive_interpolators[1](
1:441000))*drive_amplitude)) #will be updated to make sure coupling is strong enough to get synchronisation at least at one frequency.
time_range=(0.0,10.0)
phoneme_sampling_rate=44100
drive_amplitude=4.0
phase_modulation=1.0
if phase_modulation<=0.5
    q_normalisation=2/((1-phase_modulation)*π) #scaling factor to make absolute area under stimulus modulation curve = 4. 
else
    m=phase_modulation
    q_normalisation=4/(4*m*sqrt(1-((m-1)/m)^2)+2*(m-1)asin((m-1)/m))
end
p=ComponentArray(F=4.0, c=c, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,modulation=phase_modulation,q=q_normalisation,noise_case_reference=1,frequencies=[F for i in 1:60])
PCM_range=(5.5,10.0)    
save_path="./"
save_traj=true
freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

PCMs_coupled_oscillator_nofreqadj=get_coupled_oscillator_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(PCMs_coupled_oscillator_nofreqadj,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="sine wave stim - coupled oscillator no freq adj",marker_z=frequencies)

#double check trajectories:
traj_path_coupled_oscillator_nofreqadj="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/nofreq_adj_oscillator_response_to_sine_stimulus_freq_range_2.0to15.0Hz_envelope_transformed.arrow.lz4"
using Arrow, DataFrames
trajectories_coupled_oscillator_nofreqadj=Matrix(DataFrame(Arrow.Table(traj_path_coupled_oscillator_nofreqadj)))
#plot the response and stimulus for each of the 20 frequencies.
ps_coupled_nofreqadj=[]
for freq_idx in 1:5:50
    p=plot(traj_times,trajectories_coupled_oscillator_nofreqadj[:,freq_idx],xlabel="Time (s)")
    plot!(p,stimulus_times,sine_drive_interpolators[freq_idx](stimulus_indexes),label="Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))
    push!(ps_coupled_nofreqadj,p)
end
plot(ps_coupled_nofreqadj...,layout=(5,2),size=(1200,800),plot_title="coupled oscillator no freq adj")