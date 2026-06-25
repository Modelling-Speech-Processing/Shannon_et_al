### code to run Phase concentration metric tests on each model.
### needs sine wave stimuli, from which impulse trains of a range of frequencies can be extracted from the peaks.

using FFTW, DSP, Statistics


### reformed drive interpolators generator to create sine waves of varying frequency for PCM tests.
"""
Will make an "interpolators" array that contains interpolations of 20 sine waves of different frequency and put it in global variable Interpolators.
This can then be used in the ensemble simulation via the interpolation selector variable in a prob_func
"""
function generate_sine_drive_interpolators_for_saving(sampling_rate::Float64,frequency_range::Tuple{Float64,Float64},stimulus_set_number::Int)
    sr=sampling_rate
    stimlength=sr*10 #code expects 10 seconds, 5 seconds of noise, 5 seconds of stimulus. 
    times=range(0.0,10.0,length=Int(stimlength))
    frequencies=range(frequency_range[1],frequency_range[2],length=20*3) #20 frequencies per stimulus set, 3 stimulus sets. 60 freqs.
    sine_drives=Vector{Any}(undef,20)

    for trial in 1:20   
        #will go through the noise samples in triplets: 1 is the before stimulus noise, 1 is during stimulus and 1 is after stimulus
        #will loop over this, and create a single interpolator for each set of noise samples, with the stimulus added to the middle and appropriate scale applied
        #drive is [noise,stimulus,noise,zeroes] (in case it runs over, making extrapolation to be 0 only)
        trial_drive_input=sin.(2.0*pi.*frequencies[trial+20*(stimulus_set_number-1)].*times) #only for the appropriate 20 frequencies of this stimulus set
        
        xs=1:1:length(trial_drive_input)
       
        sine_drives[trial]=linear_interpolation(xs,trial_drive_input,extrapolation_bc=Line())
        #interpolators[trial]=interpolate(Vector(trial_drive_input),BSpline(Linear()))#,extrapolation_bc=Line())
    end
    return sine_drives #look into typed global
end

# new prob func for the oscillator. this will update its frequency to match the stimulus frequency.
function vary_noise_and_initial_conditions_and_natfreq(prob,i,repeat)
    prob.p.noise_selector=i
    freqs=prob.p.frequencies
    prob.p.F=freqs[Int64(i+20*(prob.p.noise_case_reference-1))] #set natural frequency to match stimulus frequency.
    prob.u0[1]=rand(MersenneTwister(Int64(123+i+20*(prob.p.noise_case_reference-1))))*2*π-π #random phase, but keeping r=1.0. indexed by noise selector i = 1:1:60.
    prob
end

function get_coupled_oscillator_PCM_across_60freqs_randinitcond(p,u_init,time_range,PCMrange,freq_range,save_traj,savepath,stimulus_type)
    condition="PCM"
    println("testing: ",condition)
    flush(stdout)
    saveat=1/44100
    response_sr=Int64(1/saveat)
    frequencies=range(freq_range[1],freq_range[2],length=60)#[1:50] #cut to 50 to avoid failed simulations at high freq end.
    all_60_trials=Vector{Any}()
    all_stimuli=Vector{Any}()
    for stimulus_idx in 1:3 #3 sets of 20 stimuli to cover 60 frequencies, structured this way to match prior phoneme drive generation.
        println("stimulus set: ",stimulus_idx)
        p.noise_case_reference=stimulus_idx
        p.noise_selector=1
        flush(stdout)
        #set stimulus based on type: either normal envelopes, or derivative or event based impulse trains (peak rate or peak envelope)
        if stimulus_type=="envelope"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
        elseif stimulus_type=="derivative"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_smooth_derivative(interpolators_global,Float64(response_sr)) 
        elseif stimulus_type=="peakrate"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_scaled_peakrate_impulses(interpolators_global,Float64(response_sr))
        elseif stimulus_type=="peakenvelope"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_unitary_peakenv_impulses(interpolators_global, Float64(response_sr))
        else
            error("Unknown stimulus type: $(stimulus_type), must be one of 'envelope', 'derivative', 'peakrate', or 'peakenvelope'.")
        end

        #drive scaling as per phoneme case workflow.
        max_stimulus_amplitude=maximum(maximum.(interpolators_global))
        p.drive_amplitude=20000/mean(sum.([i[5.0*response_sr:10.0*response_sr] for i in interpolators_global]))
        p.c=1.0*pi*(1/(max_stimulus_amplitude*p.drive_amplitude*p.q)) #including q. this is the final maximumum stimulus amplitude after all normalisation. 

        println("running sim")
        flush(stdout)
        if stimulus_idx==3 #to match the cut to 50 frequencies above. #no longer cut.
            trialdata=Ensemble_CoupledOscillators_modulated_PCM(vary_noise_and_initial_conditions_and_natfreq,time_range,p,u_init,20,saveat)        
        else
            trialdata=Ensemble_CoupledOscillators_modulated_PCM(vary_noise_and_initial_conditions_and_natfreq,time_range,p,u_init,20,saveat)
        end
        push!(all_60_trials,trialdata...)
        push!(all_stimuli,interpolators_global...)
    end
    
    PCM_vectors,rates=calculate_PCM_oscillator_noisyrates(all_60_trials,all_stimuli,frequencies,response_sr,time_range,PCMrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
    #save the rates for each trial, for comparison across natural frequencies, and computation of conglomerate activity if desired.
    if save_traj
        #check dir exists:
        if !isdir(savepath*"Trajectories/")
            mkdir(savepath*"Trajectories/")
        end
        local name="Trajectories/$(response_sr)sr_oscillator_response_to_sine_stimulus_freq_range_$(freq_range[1])to$(freq_range[2])Hz_$(stimulus_type)_transformed"
        CSV.write(savepath*name*".csv",DataFrame(rates,:auto),writeheader=true)
        generate_arrow(name,savepath)
        rm(savepath*name*".csv")
    end
    return(PCM_vectors)
end


function get_evoked_model_PCM_across_60freqs_randinitcond(p,u_init,time_range,PCMrange,freq_range,save_traj,savepath,stimulus_type)
    condition="PCM"
    println("testing: ",condition)
    flush(stdout)
    saveat=1/44100
    response_sr=Int64(1/saveat)
    frequencies=range(freq_range[1],freq_range[2],length=60) #cut to 50 to avoid failed simulations at high freq end.
    all_60_trials=Vector{Any}()
    all_stimuli=Vector{Any}()
    for stimulus_idx in 1:3 #3 sets of 20 stimuli to cover 60 frequencies, structured this way to match prior phoneme drive generation.
        println("stimulus set: ",stimulus_idx)
        p.noise_case_reference=stimulus_idx
        p.noise_selector=1
        flush(stdout)
        #set stimulus based on type: either normal envelopes, or derivative or event based impulse trains (peak rate or peak envelope)
        if stimulus_type=="envelope"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
        elseif stimulus_type=="derivative"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_smooth_derivative(interpolators_global,Float64(response_sr)) 
        elseif stimulus_type=="peakrate"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_scaled_peakrate_impulses(interpolators_global,Float64(response_sr))
        elseif stimulus_type=="peakenvelope"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_unitary_peakenv_impulses(interpolators_global,Float64(response_sr))
        else
            error("Unknown stimulus type: $(stimulus_type), must be one of 'envelope', 'derivative', 'peakrate', or 'peakenvelope'.")
        end

        println("running sim")
        flush(stdout)
        if stimulus_idx==3 #to match the cut to 50 frequencies above. #no longer cut.
            trialdata=Ensemble_EvokedModel(vary_noise_and_initial_conditions_evokedmodel,time_range,p,u_init,20,saveat)
        else
            trialdata=Ensemble_EvokedModel(vary_noise_and_initial_conditions_evokedmodel,time_range,p,u_init,20,saveat)
        end
        push!(all_60_trials,trialdata...)
        push!(all_stimuli,interpolators_global...)
    end

    PCM_vectors,rates=calculate_PCM_EvokedModel_noisyrates(all_60_trials,all_stimuli,frequencies,response_sr,time_range,PCMrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
    #save the rates for each trial, for comparison across natural frequencies, and computation of congolmerate activity if desired.
    if save_traj
        #check dir exists:
        if !isdir(savepath*"Trajectories/")
            mkdir(savepath*"Trajectories/")
        end
        local name="Trajectories/evoked_model_response_to_sine_stimulus_freq_range_$(freq_range[1])to$(freq_range[2])Hz_$(stimulus_type)_transformed"
        CSV.write(savepath*name*".csv",DataFrame(rates,:auto),writeheader=true)
        generate_arrow(name,savepath)
        rm(savepath*name*".csv")
    end
    return(PCM_vectors)
end

function get_NGNMM_PCM_across_60freqs_randinitcond(p,u_init,time_range,PCMrange,freq_range,save_traj,savepath,stimulus_type)
    condition="PCM"
    println("testing: ",condition)
    flush(stdout)
    saveat=1/44100
    response_sr=Int64(1/saveat)
    C=p.C
    vsyn=p.vsyn
    frequencies=range(freq_range[1],freq_range[2],length=60)#[1:50] #cut to 50 to avoid failed simulations at high freq end.
    @info "Frequencies being tested: $(frequencies)"
    
    all_60_trials=Vector{Any}()
    all_stimuli=Vector{Any}()
    for stimulus_idx in 1:3 #3 sets of 20 stimuli to cover 60 frequencies, structured this way to match prior phoneme drive generation.
        println("stimulus set: ",stimulus_idx)
        p.noise_case_reference=stimulus_idx
        p.noise_selector=1
        flush(stdout)
        #set stimulus based on type: either normal envelopes, or derivative or event based impulse trains (peak rate or peak envelope)
        if stimulus_type=="envelope"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]

        elseif stimulus_type=="derivative"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_smooth_derivative(interpolators_global,Float64(response_sr)) 
        elseif stimulus_type=="peakrate"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_scaled_peakrate_impulses(interpolators_global,Float64(response_sr))
        elseif stimulus_type=="peakenvelope"
            global interpolators_global = jldopen("./Sine_Drives/drive_interpolators_$(stimulus_idx).jld2","r")["drives"]
            interpolators_global=get_unitary_peakenv_impulses(interpolators_global,Float64(response_sr))
        else
            error("Unknown stimulus type: $(stimulus_type), must be one of 'envelope', 'derivative', 'peakrate', or 'peakenvelope'.")
        end
        println("running sim")
        flush(stdout)
        if stimulus_idx==3 #to match the cut to 50 frequencies above. #no longer cut.
                trialdata=Ensemble_NoisyPhoneme_PCM(vary_noise_and_initial_conditions_NGNMM,time_range,p,u_init,20,saveat)
        else
                trialdata=Ensemble_NoisyPhoneme_PCM(vary_noise_and_initial_conditions_NGNMM,time_range,p,u_init,20,saveat)
        end
        push!(all_60_trials,trialdata...)
        push!(all_stimuli,interpolators_global...)
    end
   
    PCM_vectors,rates=calculate_PCM_NGNMM_noisyrates(all_60_trials,all_stimuli,frequencies,response_sr,time_range,C,vsyn,PCMrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
    #save the rates for each trial, for comparison across natural frequencies, and computation of congolmerate activity if desired.
    if save_traj
        #check dir exists:
        if !isdir(savepath*"Trajectories/")
            mkdir(savepath*"Trajectories/")
        end
        local name="Trajectories/fast_NGNMM_response_to_abs_sine_stimulus_freq_range_$(freq_range[1])to$(freq_range[2])Hz_$(stimulus_type)_transformed"
        CSV.write(savepath*name*".csv",DataFrame(rates,:auto),writeheader=true)
        generate_arrow(name,savepath)
        rm(savepath*name*".csv")
    end
    return(PCM_vectors)
end


function gaussian_filter_and_hilbert_for_PCM(signal::Vector{Float64},sampling_rate,peak_frequency::Float64)
    num_samples=length(signal)
    #put in frequency domain
    freqs=fftfreq(num_samples,sampling_rate)
    signal_fft=fft(signal)

    #gaussian filter parameters
    center_frequency=peak_frequency
    filter_width=center_frequency/2.0

    #gaussian filter:
    gauss_kernel=exp.(-0.5.*((abs.(freqs).-center_frequency)./(filter_width)).^2)
    #apply in freq domain;
    filtered_signal_fft=signal_fft.*gauss_kernel
    #return to time domain:
    filtered_signal=ifft(filtered_signal_fft)
    #return envelope via hilbert
    analytic_signal=hilbert(real.(filtered_signal))
    return analytic_signal
end

function calculate_PCM_oscillator_noisyrates(vector_of_solutions,stimuli,frequencies,response_sampling_rate,timerange,PCMrange,noise_seed=123)
    num_trials=length(vector_of_solutions)

    rates_sr=response_sampling_rate
    #add noise to the model response ('firing rates'):
    rates=Array{Vector}(undef,num_trials)
    for (idx,data) in enumerate(vector_of_solutions)
        rates[idx]= coupled_oscillator_activity(data)[Int64(PCMrange[1]*rates_sr):Int64(PCMrange[2]*rates_sr)]#cos(theta).*r for coupled oscillator activity
        trialwise_seed=noise_seed+idx
        noises=seeded_noise(trialwise_seed, 1.0, 0.0, length(rates[idx])) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
        noise_power=mean(noises.^2)
        rate_power=mean(rates[idx].^2)
        # desired_signal_to_noise_ratio=0.5 #for high SNR test.
        desired_signal_to_noise_ratio=100.0 #for clean test.
        noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
        scaled_noise=noises.*noise_scaling_factor
        noisy_rates=(rates[idx]).+(scaled_noise)
        rates[idx]=noisy_rates
    end
    #compute PCM vectors.
    mean_phase_diff_vectors=Vector{ComplexF64}(undef,num_trials)
    for trial_idx in 1:num_trials
        #get the stimulus for this trial:
        stimulus_interpolator=stimuli[trial_idx]
        #get the peak frequency of this stimulus:
        peak_frequency=frequencies[trial_idx]
        #turn stimulus interpolator into values over time range matching the response sampling rate:
        # response_idx=range(PCMrange[1]*44100,PCMrange[2]*44100,step=1)
        response_idx=range(PCMrange[1]*rates_sr,PCMrange[2]*rates_sr,step=1)
        # response_times=response_idx./44100.0
        response_times=response_idx./rates_sr
        # stimulus_idxs=response_times.*44100.0 #stimulus at 44100Hz sampling rate. (hang over from phoneme drive interpolator.)
        stimulus_idxs=response_times.*rates_sr #stimulus at 44100Hz sampling rate. (hang over from phoneme drive interpolator.)
        stimulus_values=1.0 .+ stimulus_interpolator(stimulus_idxs)
        #filter and get instantaneous phase:
        if trial_idx==1
        @info size(rates[trial_idx])
        @info size(stimulus_values)
        end
        filtered_response=gaussian_filter_and_hilbert_for_PCM(rates[trial_idx],rates_sr,peak_frequency)
        filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(stimulus_values,rates_sr,peak_frequency)
        #get phases
        instantaneous_response_phase=angle.(filtered_response)
        instantaneous_stimulus_phase=angle.(filtered_stimulus)
        #get phase diff:
        phase_diff=instantaneous_response_phase.-instantaneous_stimulus_phase
        #convert to complex unit vector form:
        phase_diff_vectors=exp.(im.*phase_diff)
        #get average phase difference over time
        mean_phase_diff_vectors[trial_idx]=mean(phase_diff_vectors)
    end
    PCM_vectors=mean_phase_diff_vectors
    return PCM_vectors,rates
end

function calculate_PCM_EvokedModel_noisyrates(vector_of_solutions,stimuli,frequencies,response_sampling_rate,timerange,PCMrange,noise_seed=123)
    num_trials=length(vector_of_solutions)

    rates_sr=response_sampling_rate
    #add noise to the model response ('firing rates'):
    rates=Array{Vector}(undef,num_trials)
    for (idx,data) in enumerate(vector_of_solutions)
        rates[idx]= data[1,Int64(PCMrange[1]*rates_sr):Int64(PCMrange[2]*rates_sr)] #just the first component of the solution, the 'activity' of the evoked model, similar to global conductance in NGNMM.
        trialwise_seed=noise_seed+idx
        noises=seeded_noise(trialwise_seed, 1.0, 0.0, length(rates[idx])) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
        noise_power=mean(noises.^2)
        rate_power=mean(rates[idx].^2)
        # desired_signal_to_noise_ratio=0.5 #for high SNR test.
        desired_signal_to_noise_ratio=100.0 #for clean test.
        noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
        scaled_noise=noises.*noise_scaling_factor
        noisy_rates=(rates[idx]).+(scaled_noise)
        rates[idx]=noisy_rates
    end
    #compute PCM vectors.
    mean_phase_diff_vectors=Vector{ComplexF64}(undef,num_trials)
    for trial_idx in 1:num_trials
        #get the stimulus for this trial:
        stimulus_interpolator=stimuli[trial_idx]
        #get the peak frequency of this stimulus:
        peak_frequency=frequencies[trial_idx]
        #turn stimulus interpolator into values over time range matching the response sampling rate:
        response_idx=range(PCMrange[1]*rates_sr,PCMrange[2]*rates_sr,step=1)
        response_times=response_idx./rates_sr
        stimulus_idxs=response_times.*rates_sr #stimulus at 44100Hz sampling rate. (hang over from phoneme drive interpolator.)
        stimulus_values=stimulus_interpolator(stimulus_idxs)
        #filter and get instantaneous phase:
        if trial_idx==1
        @info size(rates[trial_idx])
        @info size(stimulus_values)
        end
        filtered_response=gaussian_filter_and_hilbert_for_PCM(rates[trial_idx],rates_sr,peak_frequency)
        filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(stimulus_values,rates_sr,peak_frequency)
        #get phases
        instantaneous_response_phase=angle.(filtered_response)
        instantaneous_stimulus_phase=angle.(filtered_stimulus)
        #get phase diff:
        phase_diff=instantaneous_response_phase.-instantaneous_stimulus_phase
        #convert to complex unit vector form:
        phase_diff_vectors=exp.(im.*phase_diff)
        #get average phase difference over time
        mean_phase_diff_vectors[trial_idx]=mean(phase_diff_vectors)
    end
    PCM_vectors=mean_phase_diff_vectors
    return PCM_vectors,rates
end

function calculate_PCM_NGNMM_noisyrates(vector_of_solutions,stimuli,frequencies,response_sampling_rate,timerange,C,vsyn,PCMrange,noise_seed=123)
    num_trials=length(vector_of_solutions)

    rates_sr=response_sampling_rate
    #add noise to the model response ('firing rates'):
    rates=Array{Vector}(undef,num_trials)
    for (idx,data) in enumerate(vector_of_solutions)
        rates[idx]=get_firing_rate_NMM(data,C,vsyn)[1][Int64(PCMrange[1]*rates_sr):Int64(PCMrange[2]*rates_sr)]#
        trialwise_seed=noise_seed+idx
        noises=seeded_noise(trialwise_seed, 1.0, 0.0, length(rates[idx])) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
        noise_power=mean(noises.^2)
        rate_power=mean(rates[idx].^2)
        # desired_signal_to_noise_ratio=0.5 #for high SNR test.
        desired_signal_to_noise_ratio=100.0 #for clean test.
        noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
        scaled_noise=noises.*noise_scaling_factor
        noisy_rates=(rates[idx]).+(scaled_noise)
        rates[idx]=noisy_rates
    end
    #compute PCM vectors.
    mean_phase_diff_vectors=Vector{ComplexF64}(undef,num_trials)
    for trial_idx in 1:num_trials
        #get the stimulus for this trial:
        stimulus_interpolator=stimuli[trial_idx]
        #get the peak frequency of this stimulus:
        peak_frequency=frequencies[trial_idx]
        if trial_idx==1
            @info "peak frequency for trial 1: $(peak_frequency)" 
        end
        #turn stimulus interpolator into values over time range matching the response sampling rate:
        response_idx=range(PCMrange[1]*rates_sr,PCMrange[2]*rates_sr,step=1)
        response_times=response_idx./rates_sr
        stimulus_idxs=response_times.*rates_sr #stimulus at 44100Hz sampling rate. (hang over from phoneme drive interpolator.)
        # stimulus_values=stimulus_interpolator(stimulus_idxs)
        stimulus_values=1.0 .+(stimulus_interpolator(stimulus_idxs)) #set to abs to match what the NGNMM recieved.
        #filter and get instantaneous phase:
        if trial_idx==1
        @info size(rates[trial_idx])
        @info size(stimulus_values)
        end
        filtered_response=gaussian_filter_and_hilbert_for_PCM(rates[trial_idx],rates_sr,peak_frequency)
        filtered_stimulus=gaussian_filter_and_hilbert_for_PCM(stimulus_values,rates_sr,peak_frequency)
        #get phases
        instantaneous_response_phase=angle.(filtered_response)
        instantaneous_stimulus_phase=angle.(filtered_stimulus)
        #get phase diff:
        phase_diff=instantaneous_response_phase.-instantaneous_stimulus_phase
        #convert to complex unit vector form:
        phase_diff_vectors=exp.(im.*phase_diff)
        #get average phase difference over time
        mean_phase_diff_vectors[trial_idx]=mean(phase_diff_vectors)
    end
    PCM_vectors=mean_phase_diff_vectors
    return PCM_vectors,rates
end
