using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using Statistics
using ComponentArrays, OrdinaryDiffEq, Plots, Parameters, JLD2
using LaTeXStrings
using Interpolations, DSP
Plots.default(titlefontsize=16,legendfontsize=12,tickfontsize=9,guidefontsize=11)


"""
    get_scaled_peakrate_impulses(envelopes::Vector{T}, sampling_rate::Float64; scale_range=(0.5,1.0)) where T <: Interpolations.Extrapolation
    Takes in  envelopes, turns them into smooth derivatives, and then finds the peaks in the derivatives.
    Returns a vector of Interpolations.Extrapolation of impulse peak-rate events. originally of magnitude = rate.
    but then scaled to [0.5,1.0] as in "Oganian Y et al. Phase Alignment of Low-Frequency Neural Activity to the Amplitude Envelope of Speech Reflects Evoked Responses to Acoustic Edges, Not Oscillatory Entrainment. J Neurosci. 2023. doi: 10.1523/JNEUROSCI.1663-22.2023"
"""
function get_scaled_peakrate_impulses_and_locations(envelopes::Vector{T}, sampling_rate::Float64;scale_range=(0.0,1.0)) where T <: Interpolations.Extrapolation
    smooth_derivatives = get_smooth_derivative(envelopes, sampling_rate)
    sample_grid= envelopes[1].itp.itp.parentaxes[1][1:end]  # Get the time grid from the first envelope
    peakrate_impulses = Vector{T}(undef, length(envelopes))
    lcoations_list=Vector{Vector{Int64}}(undef,length(envelopes))
    for (i, derivative) in enumerate(smooth_derivatives)
        # Create a zero vector of the same length as the envelope
        derivative_samples=derivative(sample_grid)
        peaks,locations = findpeaks(derivative_samples, 0.1)  # Threshold set to 0.1

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
            impulse_vector[loc] = peak  # Set peak rate impulse to the value of the derivative at the peak location
        end
        # Interpolate the impulse vector
        peakrate_impulses[i] = linear_interpolation(StepRange(sample_grid), impulse_vector, extrapolation_bc=Line())
        lcoations_list[i]=locations
    end

    return peakrate_impulses, lcoations_list
end

function segment_about_phoneme_onsets(phase_data,stim_times,response_window,response_sr)
responses_2D=Array{Float64,2}(undef,length(stim_times),Int64(response_sr*response_window*2)+1)
    for i in 1:length(stim_times)
        responses_2D[i,:]=phase_data[Int64(round((stim_times[i]-response_window)*response_sr)):Int64(round((stim_times[i]+response_window)*response_sr))]
    end
    return responses_2D
end

function logspace(start, stop, num)
    return exp.(range(log(start), log(stop), length=num))
end

condition="b" #which phoneme condition to use for the drive
idx=1 #which test index to use
noisestimratio=0.0 #which noise to stimulus ratio to use
noise_filename="NoiseSequences60.csv" #which noise file to use

global NGNMM_NSP_paper_code.interpolators_global = jldopen("./Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]

#get exact peak rate event times from the drive interpolators.
peaks,locs=get_scaled_peakrate_impulses_and_locations(NGNMM_NSP_paper_code.interpolators_global,44100.0)

stim_times=collect(range(5.0,9.0,step=0.25))


#fast NGNMM model's response to each phoneme:
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
DAR=12.2
phoneme_sampling_rate=44100
drive_amplitude=η_0*DAR
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=drive_amplitude,noise_selector=1,noise_case_reference=1) 
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,10.0)
saveat=0.0001   

prob_func=vary_noise_and_initial_conditions_NGNMM

results=Ensemble_NoisyPhoneme(prob_func,time_range,p,u0,20,saveat)


## with 1/f noise on the firing rates rather than abs(Z).
noisy_rates = Vector{Vector{Float64}}(undef, 20)
for i in 1:20
    rates=get_firing_rate_NMM(results[i],p.C,p.vsyn)[1]
    
    trialwise_seed=1000+i #different seed for each trial.
    noise=seeded_noise(trialwise_seed+i, 1.0, 0.0, length(rates)) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
    
    noise_power=var(noise)
    rate_power=var(rates)
    
    desired_signal_to_noise_ratio=0.1
    noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
    
    scaled_noise=noise.*noise_scaling_factor
    noisy_rates[i]=rates+scaled_noise
end

# plot(noisy_rates[1])
# plot!(get_firing_rate_NMM(results[1],p.C,p.vsyn)[1])


## now filter the noisy trajectories into frequency bands.
freq_range=range(0.67,35,length=100)
trajectories=noisy_rates

filtered_trajectories=Array{Float64,3}(undef,length(trajectories),length(freq_range),length(results[1][3,:]))
filters=Vector{Any}(undef,length(freq_range))
for (i,freq) in enumerate(freq_range)
   if freq<1
    filters[i]=digitalfilter(Bandpass(0.1,freq+1),Butterworth(2);fs=10000)
   else
    filters[i]=digitalfilter(Bandpass(freq-1,freq+1),Butterworth(2);fs=10000)
   end
end

for i in 1:length(freq_range)
    for j in 1:length(trajectories)
         filtered_trajectories[j,i,:]=filtfilt(filters[i],trajectories[j])
    end
end

hilbert_transforms=Array{ComplexF64,3}(undef,size(filtered_trajectories))
for i in 1:size(filtered_trajectories,2)
    for j in 1:size(filtered_trajectories,1)
        hilbert_transforms[j,i,:]=hilbert(filtered_trajectories[j,i,:])
    end
end

phases=angle.(hilbert_transforms)

response_window=0.25
response_sr=Int64(1/saveat)

segmented_phase_data=Array{Float64,4}(undef,size(phases,1),size(phases,2),length(stim_times),Int64(response_sr*response_window*2)+1);
for i in 1:size(phases,1) #for each trajectory
    for j in 1:size(phases,2) #for each frequency
        locs_based_stim_times=locs[i]./44100
        locs_based_stim_times_in_5_to_9=filter(x->x>=5 && x<=9, locs_based_stim_times) #interested in onsets occuring between 5 and 9 seconds (at which the windows are centred)
        segmented_phase_data[i,j,:,:]=segment_about_phoneme_onsets(phases[i,j,:],locs_based_stim_times_in_5_to_9,response_window,response_sr)
    end
end
#check segmented phase data looks reasonable:
# plot(segmented_phase_data[2,14,:,:]')

#stack the segments of all trajectories into one dimension. that is 20x17=340 segments per frequency, each segment is 5001 long.
segmented_phase_data_reshaped=Array{Float64,3}(undef,size(segmented_phase_data,2),size(segmented_phase_data,3)*size(segmented_phase_data,1),size(segmented_phase_data,4))
for i in 1:size(segmented_phase_data,2) #foser each frequency
    segmented_phase_data_reshaped[i,:,:]=reshape(segmented_phase_data[:,i,:,:],size(segmented_phase_data,1)*size(segmented_phase_data,3),size(segmented_phase_data,4))
end

##now calculate the ITPC across events (every segment is an event, every trajectory has multiple segments.) for each time point and frequency:
ITPCs=Array{Float64,2}(undef,size(segmented_phase_data_reshaped,1),size(segmented_phase_data_reshaped,3))
for i in 1:size(segmented_phase_data_reshaped,1) #for each frequency
    for j in 1:size(segmented_phase_data_reshaped,3) #for each time point
        ITPCs[i,j]=abs(mean(exp.(im*segmented_phase_data_reshaped[i,:,j])))
    end
end


####
#PLOT FIGURE 9
####
#plot the ITPCs across frequencies and time points:
# p_ITPCs=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs,color=:viridis, xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}", ylabel=L"\textrm{frequency~(Hz)}", title=L"\textrm{ITPC~across~phoneme~onsets}",dpi=300)
p_ITPCs=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs,color=:viridis, xlabel="", ylabel=L"\textrm{frequency~(Hz)}", title="",dpi=300,cbar_title=L"$\operatorname{IEPC}$")

#remove baseline ITPC by subtracting the mean ITPC in the -0.2 to -0.1s window from each time point, to see the change in ITPC relative to baseline:
baseline_ITPCs_all_phonemes=mean(ITPCs[:,Int64(0.05*response_sr):Int64(0.15*response_sr)],dims=2)
ITPCs_minus_baseline_all_phonemes=ITPCs .- baseline_ITPCs_all_phonemes
#threshold to 0, to remove negative values:
ITPCs_minus_baseline_all_phonemes[ITPCs_minus_baseline_all_phonemes .< 0.0] .= 0.0
# p_ITPCs_minus_baseline_all_phonemes=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs_minus_baseline_all_phonemes,color=:viridis, xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}", ylabel=L"\textrm{frequency~(Hz)}", title=L"\textrm{ITPC~minus~baseline}",dpi=300)
p_ITPCs_minus_baseline_all_phonemes=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs_minus_baseline_all_phonemes,color=:viridis, xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}", ylabel=L"\textrm{frequency~(Hz)}", title="",dpi=300,cbar_title=L"$\operatorname{IEPC}$")

#plot all phoneme case and baseline removed all phoneme case together.
plot(p_ITPCs,p_ITPCs_minus_baseline_all_phonemes,layout=grid(2,1),size=figure_size_tuple(1,aspect_ratio=0.7),dpi=300,margin=2Plots.mm)
savefig("IEPC_all_phonemes.pdf")

###
#OTHER PLOTS/Tests not in paper:
###

## caclulate ITPCs from just the very first phoneme to see transient response from noise to stimulus application.
first_phoneme_phase_data=segmented_phase_data_reshaped[:,1:17:end,:] #first 20 segments are the first phoneme across the 20 trajectories.

ITPCs_first_phoneme=Array{Float64,2}(undef,size(first_phoneme_phase_data,1),size(first_phoneme_phase_data,3))
for i in 1:size(first_phoneme_phase_data,1) #for each frequency
    for j in 1:size(first_phoneme_phase_data,3) #for each time point
        ITPCs_first_phoneme[i,j]=abs(mean(exp.(im*first_phoneme_phase_data[i,:,j])))
    end
end 

p_ITPCs_first_phoneme=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs_first_phoneme,color=:viridis, xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}", ylabel=L"\textrm{frequency~(Hz)}", title=L"\textrm{ITPC~first~phoneme}",dpi=300)
p_ITPCs_first_phoneme=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs_first_phoneme,color=:viridis, xlabel="", ylabel=L"\textrm{frequency~(Hz)}",dpi=300)

#plot ITPC minus baseline, where baseline is the mean ITPC in the -0.25 to -0.15s window:
baseline_ITPCs=mean(ITPCs_first_phoneme[:,1:Int64(0.10*response_sr)],dims=2)
ITPCs_minus_baseline=ITPCs_first_phoneme .- baseline_ITPCs
#threshold to 0, to remove negative values
ITPCs_minus_baseline[ITPCs_minus_baseline .< 0.0] .= 0.0
p_ITPCs_minus_baseline=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs_minus_baseline,color=:viridis, xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}", ylabel=L"\textrm{frequency~(Hz)}", title=L"\textrm{ITPC~first~phoneme~minus~baseline}",dpi=300)
p_ITPCs_minus_baseline=heatmap(range(-0.25,0.25,step=1/response_sr),freq_range,ITPCs_minus_baseline,color=:viridis, xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}", ylabel=L"\textrm{frequency~(Hz)}", dpi=300)



## checking logic to get the first phoneme trajectories.
#right now segment the drive again using the times from locs, and see how well they align.
drive_sr=44100
new_drive_segments=Array{Float64,2}(undef,340,Int64(drive_sr*response_window*2)+1)
for i in 1:20
    new_stim_times=locs[i]./44100
    new_stim_times_in_5_to_9=filter(x->x>=5 && x<=9, new_stim_times)
    new_drive_segments[i*17-16:i*17,:]=segment_about_phoneme_onsets(NGNMM_NSP_paper_code.interpolators_global[i],new_stim_times_in_5_to_9,response_window,drive_sr)
end
plot(new_drive_segments[1:2:340,:]',label=nothing,dpi=300)
#add line at t=0 to check alignment:
plot!([response_window*drive_sr,response_window*drive_sr],[0.0,0.3],color=:black,linewidth=1,linestyle =:dash)

first_phoneme_drive_data=new_drive_segments[1:17:end,:] #first 20 segments are the first phoneme across the 20 trajectories.
plot(first_phoneme_drive_data')

first_phoneme_trajectories=Array{Float64,2}(undef,340,5001)
for i in 1:20
    locs_based_stim_times=locs[i]./44100
    locs_based_stim_times_in_5_to_9=filter(x->x>=5 && x<=9, locs_based_stim_times)
    first_phoneme_trajectories[i*17-16:i*17,:]=segment_about_phoneme_onsets(abs.(results[i][3,:]),locs_based_stim_times_in_5_to_9,response_window,response_sr)
end
first_phoneme_trajectories=first_phoneme_trajectories[1:17:end,:] #first 20 segments are the first phoneme across the 20 trajectories.
plot(first_phoneme_trajectories')




