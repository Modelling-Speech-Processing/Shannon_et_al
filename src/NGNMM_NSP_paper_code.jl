module NGNMM_NSP_paper_code

# __precompile__(false)z
using FFTW, DSP, CSV, JSON, WAV, LinearAlgebra, Random, Interpolations, Statistics, WAV, Parameters, ComponentArrays, OrdinaryDiffEq, DelimitedFiles, Distributions, DataFrames, Arrow, JLD2, PowerLawNoise, Plots

#Dynamical Systems Models
export NMM_CA,NMM_PhonemeDrive_Noisy,Kuramoto_CA, Kuramoto_PhonemeDrive_Noisy, kuramoto_phase_reset_condition, kuramoto_phase_reset_effect!, cartesian_kuramoto
#Stimuli Processing
export narrowband_envelopes, store_envelopes!,interpolators, generate_drive_interpolators_specify_noise2stimulus_ratio, generate_drive_interpolators_specify_noise2stimulus_ratio_forsaving, select_and_modify_20_random_phonemes
#Experiment Running Functions
export vary_noise, Ensemble_NoisyPhoneme, get_ITPC_t2pd_correlation_varynoisestimratio_150Hz!, run_noise_tests_serial_varydriveamplituderatio
#Data Analysis Functions
export calculate_ITPC, get_synaptic_current_NMM, get_firing_rate_NMM, get_frequencies, get_sorting_indices, calculate_ITPC_1overf_noise
export cos_activity

export generate_arrow
export coupled_oscillator!, oscillator_phase_reset_condition_up, oscillator_phase_reset_effect_up!, oscillator_phase_reset_condition_down, oscillator_phase_reset_effect_down!
export Ensemble_CoupledOscillators, calculate_ITPC_CoupledOscillators_noisyrates, run_coupled_oscillator_noise_tests_serial_varydriveamplituderatio, run_coupled_oscillator_noise_tests_serial_varynoisestimratio, coupled_oscillator_activity
#tests across all 60 trials at once:
export run_coupled_oscillator_noise_tests_serial_varydriveamplituderatio_60trials, run_noise_tests_serial_varydriveamplituderatio_60trials, get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials!, get_coupled_oscillator_activity_varynoisestimratio_60trials
export run_coupled_oscillator_modulated_noise_tests_serial_varydriveamplituderatio_60trials, get_coupled_oscillator_modulated_activity_varynoisestimratio_60trials!, Ensemble_CoupledOscillators_modulated, coupled_oscillator_modulated!
export seeded_noise
#random initial conditions per noise set:
export run_coupled_oscillator_modulated_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond, get_coupled_oscillator_modulated_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!
export run_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond, get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!
export run_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond_1overf_noise, get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond_1overf_noise!
#evoked models
export alpha_kernel_ODE!, findpeaks, get_smooth_derivative, get_scaled_peakenv_impulses, get_scaled_peakrate_impulses
export run_evoked_model_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond, get_evoked_model_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!, calculate_ITPC_EvokedModel_noisyrates, Ensemble_EvokedModel, vary_noise_and_initial_conditions_evokedmodel
## Dynamical Systems Models
function NMM_CA(D,u,p,t)
    @unpack α,k,C,Δ,η_0,vsyn,drive_switch,α_D,Π = p
    @unpack g_dot, g, Z, A_dot, A = u #A is drive 

    #u[1] is g', 
    D.g_dot = (α^2)*(k/(C*pi)) * ((1-(abs(Z)^2))/(1+Z+conj(Z)+(abs(Z)^2))) - (2*α)*g_dot - (α^2)*g

    #u[2] is g
    D.g = g_dot

    #u[3] is Z
    D.Z = (1/C) * (-im*(((Z-1)^2)/2) + (((Z+1)^2)/2)*(-Δ+im*(η_0+A)+im*vsyn*g) - (((Z^2)-1)/2)*g)

    D.A_dot=(α_D^2)*Π*drive_switch-2*α_D*A_dot-(α_D^2)*A  #I believe swapping Π for Π(t-drivestarttime) is the way to incorporate ∀ drive functions (of time)

    D.A=A_dot

end
function NMM_PhonemeDrive_Noisy(D,u,p,t)
    @unpack α,k,C,Δ,η_0,vsyn,α_D, sampling_rate, drive_amplitude, noise_selector, noise_case_reference = p
    @unpack g_dot, g, Z, A_dot, A = u #A is drive 
    
    #u[1] is g', 
    D.g_dot = (α^2)*(k/(C*pi)) * ((1-(abs(Z)^2))/(1+Z+conj(Z)+(abs(Z)^2))) - (2*α)*g_dot - (α^2)*g
    
    #u[2] is g
    D.g = g_dot
       
    #u[3] is Z
    D.Z = (1/C) * (-im*(((Z-1)^2)/2) + (((Z+1)^2)/2)*(-Δ+im*(η_0+A)+im*vsyn*g) - (((Z^2)-1)/2)*g)

    #interpolators[NoiseSelector] chooses a given interpolator out of the previously constructed global interpolator vector (containing all the noise conditions and the given stimulus)
    D.A_dot=(α_D^2).*interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude-2*α_D*A_dot-(α_D^2)*A 
    
    D.A=A_dot
    
end

#coupled oscillator control case:
function oscillator_phase_reset_condition_up(u,t,integrator)
    u[1] - pi
end

function oscillator_phase_reset_effect_up!(integrator)
    # N=length(integrator.u.θs)
    integrator.u[1]=-pi
end

function oscillator_phase_reset_condition_down(u,t,integrator)
    u[1] - -pi 
end

function oscillator_phase_reset_effect_down!(integrator)
    # N=length(integrator.u.θs)
    integrator.u[1]=pi
end

function coupled_oscillator!(du, u, p, t)
    @unpack F, c, drive_amplitude, noise_selector, sampling_rate,q = p
    @unpack θ, r = u
    du.θ = 2*pi*F - c * q * ((interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude)/r) * sin(θ) 
    # du.θ = 2*pi*F - c * ((interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude)/r) * sin(θ/2) #for sin(1/2theta) run
    # du.θ = 2*pi*F - c * q * ((interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude)/r)# * sin(θ) #for no sin run
    # du.θ = 2*pi*F + c * ((interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude)/r)# * sin(θ) #for positive correction run
    du.r = r*(1-r^2) + c * q * (interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude) * cos(θ)
    # du.r = r*(1-r^2) + c * q * (interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude) #* cos(θ) #for no sin no cos run
    # du.r = r*(1-r^2) - c * (interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude) * cos(θ) #for positivenegative correction run
    # du.r = r*(1-r^2) + c * (interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude) * cos(θ/2) #for sin(1/2theta) run 
    return nothing
end

# has a discontinuity when modulation makes the stimulus correction always positive.
function coupled_oscillator_modulated!(du, u, p, t)
    @unpack F, c, drive_amplitude, noise_selector, sampling_rate,modulation,q, noise_case_reference = p
    @unpack θ, r = u
    du.θ = 2*pi*F - c * q * ((interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude)/r) * (1*(1-modulation)+modulation*sin(θ)) #modulation of 1.0 means full phase modulation. 0.0 = no phase modulation.
    du.r = r*(1-r^2) + c * q * (interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude) * (1*(1-modulation)+modulation*cos(θ)) 
    return nothing
end

# to avoid the discontinuity, have positive only corrections. and scale those.
# function coupled_oscillator_modulated!(du, u, p, t)
#     @unpack F, c, drive_amplitude, noise_selector, sampling_rate,modulation,q, noise_case_reference = p
#     @unpack θ, r = u
#     du.θ = 2*pi*F - c * q * ((interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude)/r) * (1*(1-modulation)+modulation*maximum([sin(θ),0.0]))
#     du.r = r*(1-r^2) + c * q * (interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude) * (1*(1-modulation)+modulation*maximum([cos(θ),0.0])) 
#     return nothing
# end

"""
a bank of 20 modulated coupled oscillators with natural frequencies varying from 0.1 to 10Hz. t
They are not coupled to each other, but are driven by the same stimulus.
"""
function coupled_oscillator_bank_modulated!(du,u,p,t)
    @unpack F, c, drive_amplitude, noise_selector, sampling_rate,modulation,q = p
    @unpack θs, rs = u
    N=Int64(length(θs))
    for i in Int64(1):Int64(N)
        du.θs[i] = 2*pi*F - c * q * ((interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude)/rs[i]) * (1*(1-modulation)+modulation*maximum([sin(θs[i]),0.0]))
        du.rs[i] = rs[i]*(1-rs[i]^2) + c * q * (interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude) * (1*(1-modulation)+modulation*maximum([cos(θs[i]),0.0])) 
    end
    return nothing
end


"""
Equivalent to alpha-kernel convolution, ODE form, using the synaptic filter from the NGNMM. 2nd Order ODE.
"""
function alpha_kernel_ODE!(du,u,p,t)
    @unpack α,drive_amplitude, noise_selector,sampling_rate, noise_case_reference = p
    @unpack x1,x2 = u
    #x1 = synaptic current, x2 = synaptic current derivative
    du.x1 = x2
    du.x2 = (α^2).*interpolators_global[Int64(noise_selector)](t*sampling_rate+1.0)*drive_amplitude-2*α*x2-(α^2)*x1 
    return nothing
end


## Stimuli Processing
function narrowband_envelopes(signal,sampling_rate,freq_low=80,freq_high=8000,num_bands=32)
    logspace=range(log10(freq_low), stop=log10(freq_high), length=num_bands+1)
    frequency_bands=10 .^ logspace
    for i in 1:length(frequency_bands)
        if frequency_bands[i] < 100
            frequency_bands[i]=round(frequency_bands[i],digits=-1)
        elseif frequency_bands[i] <1000
            frequency_bands[i]=round(frequency_bands[i],digits=-1)
        else
            frequency_bands[i]=round(frequency_bands[i],digits=-2)
        end
    end

    Filtered_sounds=Array{Float64}(undef, length(signal[:,1]),length(frequency_bands)-1)
    
    for i in 1:(length(frequency_bands)-1)
        freq_band=[frequency_bands[i],frequency_bands[i+1]]
        freqrad=freq_band
        # println(freqrad)
        a=digitalfilter(Bandpass(freqrad[1],freqrad[2];fs=sampling_rate),Butterworth(2))
        Filtered_sounds[:,i]=filt(a,signal[:,1])
    end

    Filtered_Envelopes=Array{Complex}(undef, length(signal[:,1]),length(frequency_bands)-1)

    for i in 1:length(Filtered_sounds[1,:])
        # println(i)
        Filtered_Envelopes[:,i]=hilbert(Filtered_sounds[:,i])
    end
    return Filtered_Envelopes
end

function store_envelopes!(d,condition,stim_path_dict)
   println("storing: ",condition)
   get!(d,condition) do 
       for_storage=Vector{Vector{Float64}}(undef,3) #3 stimuli per condition
       paths=stim_path_dict[condition] 
       for (idx,stimuluspath) in enumerate(paths)
           println("generating envelopes")
           stimulus, fs, _, _=wavread(join(["./StimuliNorm/",stimuluspath]))
           stimulus_envelopes=abs.(narrowband_envelopes(stimulus,fs))
           stimulus_envelope_sum=sum(stimulus_envelopes,dims=2)[:]
           for_storage[idx]=stimulus_envelope_sum #store entire envelope into vector into dictionary under the condition label key.
       end
       for_storage
   end
end

const interpolators = Ref(Vector{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Linear{Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Linear{Throw{OnGrid}}}, Tuple{StepRange{Int64, Int64}}}, BSpline{Linear{Throw{OnGrid}}}, Line{Nothing}}}(undef,20))
global interpolators_global::Vector{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Linear{Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Linear{Throw{OnGrid}}}, Tuple{StepRange{Int64, Int64}}}, BSpline{Linear{Throw{OnGrid}}}, Line{Nothing}}}
"""
Will make an "interpolators" array that contains 20 interpolations of the given stimulus with different noise before/after and during the stimulus. and put it in global variable Interpolators.
This can then be used in the ensemble simulation via the interpolation selector variable in a prob_func

before the stimulus is applied, noise is applied of amplitude scale 1.0 (which will be multiplied by the drive_amplitude parameter of the ODE model later)
the noisestimratio parameter sets the ratio of noise to stimulus amplitude to use when the stimulus is applied. for instance, if a ratio of 0.4 noise to 0.6 stimulus is required
just set noisestimratio to 0.4. This results in a constant amplitude drive being applied to the model of format (1.0*noise, (0.4*noise + 0.6*stimulus), 1.0*noise).
"""
function generate_drive_interpolators_specify_noise2stimulus_ratio(stimulus_envelope_sum,noisestimratio,noise_filename)
    stimlength=length(stimulus_envelope_sum)
    noise=readdlm("./$(noise_filename)",',')
    if typeof(noise[1])==SubString{String}
        noise=readdlm("./$(noise_filename)")
    end
    all_noise_vector=[noise[i,:] for i in 1:size(noise,1)]
    noise_length=5*44100
    out_of_stim_noise_scale=1.0
    during_stim_noise_scale=noisestimratio #scales size of amplitude for the during-stimulus noise.

    for trial in 1:20   
        #will go through the noise samples in triplets: 1 is the before stimulus noise, 1 is during stimulus and 1 is after stimulus
        #will loop over this, and create a single interpolator for each set of noise samples, with the stimulus added to the middle and appropriate scale applied
        #drive is [noise,stimulus,noise,zeroes] (in case it runs over, making extrapolation to be 0 only)
        trial_drive_input=[all_noise_vector[1+(trial-1)*3]*out_of_stim_noise_scale;(stimulus_envelope_sum.*(1-during_stim_noise_scale)).+(all_noise_vector[2+(trial-1)*3][noise_length-stimlength+1:end]*during_stim_noise_scale);all_noise_vector[3+(trial-1)*3]*out_of_stim_noise_scale;zeros(length(all_noise_vector[1]))]
        
        xs=1:1:length(trial_drive_input)
       
        interpolators[][trial]=linear_interpolation(xs,Vector(trial_drive_input),extrapolation_bc=Line())
        #interpolators[trial]=interpolate(Vector(trial_drive_input),BSpline(Linear()))#,extrapolation_bc=Line())
    end
    return interpolators[] #look into typed global
end


function generate_drive_interpolators_specify_noise2stimulus_ratio_forsaving(stimulus_envelope_sum,noisestimratio,noise_filename)
    stimlength=length(stimulus_envelope_sum)
    noise=readdlm("./$(noise_filename)",',')
    if typeof(noise[1])==SubString{String}
        noise=readdlm("./$(noise_filename)")
    end
    all_noise_vector=[noise[i,:] for i in 1:size(noise,1)]
    noise_length=5*44100
    out_of_stim_noise_scale=1.0
    during_stim_noise_scale=noisestimratio #scales size of amplitude for the during-stimulus noise.
    drives=Vector{Any}(undef,20)
    for trial in 1:20   
        #will go through the noise samples in triplets: 1 is the before stimulus noise, 1 is during stimulus and 1 is after stimulus
        #will loop over this, and create a single interpolator for each set of noise samples, with the stimulus added to the middle and appropriate scale applied
        #drive is [noise,stimulus,noise,zeroes] (in case it runs over, making extrapolation to be 0 only)
        trial_drive_input=[all_noise_vector[1+(trial-1)*3]*out_of_stim_noise_scale;(stimulus_envelope_sum.*(1-during_stim_noise_scale)).+(all_noise_vector[2+(trial-1)*3][noise_length-stimlength+1:end]*during_stim_noise_scale);all_noise_vector[3+(trial-1)*3]*out_of_stim_noise_scale;zeros(length(all_noise_vector[1]))]
        println(typeof(trial_drive_input))
        println(typeof(trial_drive_input[1]))
        xs=1:1:length(trial_drive_input)
       
        drives[trial]=linear_interpolation(xs,Vector(trial_drive_input),extrapolation_bc=Line())
        #interpolators[trial]=interpolate(Vector(trial_drive_input),BSpline(Linear()))#,extrapolation_bc=Line())
    end
    return drives #look into typed global
end

"""
Selects 20 random phoneme envelopes from the individual phoneme envelopes csv file, stretches/squeezes them to some degree. \n 
Returns a 5 second (at 44100Hz sampling rate) sequence to be used as noise stimulus before \n
test stimulus is applied - save a set of these so later will only need to read a csv file rather than construct them every time.
"""
function select_and_modify_20_random_phonemes(phonemes,rng)
    selection_indices=rand(rng,eachindex(phonemes),20)
    chosen_phonemes=phonemes[selection_indices]
    adjusted_phonemes=similar(chosen_phonemes)
    #if remove half of the samples, the phoneme is half the size.
    #if duplicate half the samples, the phoneme is 1.5 times the size. so keep correct total size by matching the number of removed and duplicated samples
    #squeeze 10 of the phonemes
    for i in 1:2
        for (idx,squeeze) in enumerate([3,6,10,15,20]) #remove either every 3rd, 6th, 10th, 15th or 20th sample 
            adjusted_phonemes[idx+(i-1)*5]=chosen_phonemes[idx+(i-1)*5][[k for k in eachindex(chosen_phonemes[idx+(i-1)*5]) if k%squeeze!=0]]
        end
    end
    #stretch 10 phonemes
    for i in 3:4
        for (idx,stretch) in enumerate([3,5,8,10,12]) #duplicate either every 3th, 5th, 8th, 10th or 12th sample, more are duplicated than were removed, to ensure sequence is long enough, will be cut from the start to ensure the end is clean.
            all_entries_duplicated=repeat(chosen_phonemes[idx+(i-1)*5],inner=2)
            adjusted_phonemes[idx+(i-1)*5]=all_entries_duplicated[[k for k in eachindex(all_entries_duplicated) if k%2!=0 || k%(2*stretch)==0]] #k%!=0 keeps all original samples (the odd ones). k%(2*stretch)==0 keeps all new samples that were generated from samples at indexes that are mulitples of k.
        end
    end
    #5 seconds at 44100Hz sampling rate is 220500 samples
    n_samples=44100*5
    phoneme_seq=reduce(vcat,shuffle(rng,adjusted_phonemes)) #shuffle order of modified phonemes so stretched and squeezed in random order, then concatenate into single vector
    num_samples_to_cut=length(phoneme_seq)-n_samples
    return phoneme_seq[num_samples_to_cut+1:end]
end


## Experiment Running Functions

function vary_noise(prob,i,repeat)
    prob.p.noise_selector=i
    prob
end

function vary_noise_and_initial_conditions(prob,i,repeat)
    prob.p.noise_selector=i
    prob.u0[1]=rand(MersenneTwister(Int64(123+i+20*(prob.p.noise_case_reference-1))))*2*π-π #random phase, but keeping r=1.0. indexed by noise selector i = 1:1:60.
    prob
end


function vary_noise_and_initial_conditions_NGNMM(prob,i,repeat)
    prob.p.noise_selector=i
    #random Z in unit circle, with radius buffer so it definitely will not be on or too close to the boundary.
    random_argument=rand(MersenneTwister(Int64(123+i+20*(prob.p.noise_case_reference-1))))*2*π-π
    random_radius=rand(MersenneTwister(Int64(223+i+20*(prob.p.noise_case_reference-1))))*0.9
    prob.u0[3]=random_radius*exp(im*random_argument)
    #A, A_dot remains at 0 to start, it is the stimulus
    #g g_dot set at 0 to avoid discontinuous transient, or unrealistic stimulus, ensuring Z stays in unit circle.
    prob
end

function Ensemble_NoisyPhoneme(prob_func,timerange,p,u0,trajectories,saveat)
    model=NMM_PhonemeDrive_Noisy
    #create ensemble ODE problem using NMM_PhonemeDrive as the model & the given prob_func.
    prob=ODEProblem(model,u0,timerange,p)
    EnsembleProb=EnsembleProblem(prob,prob_func=prob_func)

    #solve for given number of trajectories
    results=solve(EnsembleProb,EnsembleThreads(),trajectories=trajectories,saveat=saveat)
    return(results)
end

function get_ITPC_t2pd_correlation_varynoisestimratio_150Hz!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noise_filename)
    println("testing: ",condition)
    C=p.C
    vsyn=p.vsyn
    saveat=0.0001
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        envelopes=Phoneme_Envelopes[condition]
        for (idx,stimulus_envelope) in enumerate(envelopes)
            println("NSR: ",noisestimratio," Test: ",idx)
            generate_drive_interpolators_specify_noise2stimulus_ratio(stimulus_envelope,noisestimratio,noise_filename)
            println("running sim")
            trialdata=Ensemble_NoisyPhoneme(vary_noise,time_range,p,u_init,20,saveat)
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,_=calculate_ITPC(trialdata,sr,ITPCrange,C,vsyn,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
        end
        for_storage
    end
end

"""
Will run a series of 20 trial tests for ITPC and correlation to time 2 peak derivative with differing Π/η_0 ratios, on both "Aines" parameters and "my" parameter set. THis is done by changing Π accordingly for fixed η_0. 
It will do so for all the noisestimratio conditions also. Each run through 6 noise ratios takes around an hour.
Π/Π for initial test with aines parameters was: 15/21.5 = 0.70
Π/Π for initial test with my parameters was: 4.0/0.03285 = 121.77
Here I will test 6 ratios between 0.1 to 300 logarithmically spaced as there is an order 100 scale change.
Also will save the ITPC values up to 150Hz rather than just 32.
"""
function run_noise_tests_serial_varydriveamplituderatio(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,save_path,noise_filename)
    
    freq_indexes=[19,37,55]
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        correlations=Dict{String,Any}()
        println("testing DAR = ",DAR)
        η_0=p.η_0
        p.drive_amplitude=η_0*DAR

        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                get_ITPC_t2pd_correlation_varynoisestimratio_150Hz!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noise_filename)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end

        #make the directory exist, 
        mkpath(save_path)
        tosave=json(correlations)
        open(string(save_path,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

function get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noisesets)
    println("testing: ",condition)
    C=p.C
    vsyn=p.vsyn
    saveat=0.0001
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        for idx in 1:3
            all_60_trials=Vector{Any}()
            for noise_filename in noisesets
                p.noise_selector=1
                println("NSR: ",noisestimratio,"Noiseset: ",noise_filename," Test: ",idx)
                # generate_drive_interpolators_specify_noise2stimulus_ratio(stimulus_envelope,noisestimratio,noise_filename)
                global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                println("running sim")
                trialdata=Ensemble_NoisyPhoneme(vary_noise,time_range,p,u_init,20,saveat)
                push!(all_60_trials,trialdata...)
            end
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,_=calculate_ITPC(all_60_trials,sr,ITPCrange,C,vsyn,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
        end
        for_storage
    end
end

"""
Will run a series of 20 trial tests for ITPC and correlation to time 2 peak derivative with differing Π/η_0 ratios, on both "Aines" parameters and "my" parameter set. THis is done by changing Π accordingly for fixed η_0. 
It will do so for all the noisestimratio conditions also. Each run through 6 noise ratios takes around an hour.
Π/Π for initial test with aines parameters was: 15/21.5 = 0.70
Π/Π for initial test with my parameters was: 4.0/0.03285 = 121.77
Here I will test 6 ratios between 0.1 to 300 logarithmically spaced as there is an order 100 scale change.
Also will save the ITPC values up to 150Hz rather than just 32.
"""
function run_noise_tests_serial_varydriveamplituderatio_60trials(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,save_path,noisesets)
    
    freq_indexes=[19,37,55] #4, 8 and 12 Hz indexes
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        correlations=Dict{String,Any}()
        println("testing DAR = ",DAR)
        η_0=p.η_0
        p.drive_amplitude=η_0*DAR

        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noisesets)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sqr_correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j].^2) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end

        #make the directory exist, 
        mkpath(save_path)
        tosave=json(correlations)
        open(string(save_path,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

function get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noisesets)
    println("testing: ",condition)
    C=p.C
    vsyn=p.vsyn
    saveat=0.0001
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        for idx in 1:3
            all_60_trials=Vector{Any}()
            for noise_filename in noisesets
                p.noise_selector=1
                println("NSR: ",noisestimratio,"Noiseset: ",noise_filename," Test: ",idx)
                # generate_drive_interpolators_specify_noise2stimulus_ratio(stimulus_envelope,noisestimratio,noise_filename)
                global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                println("running sim")
                trialdata=Ensemble_NoisyPhoneme(vary_noise_and_initial_conditions_NGNMM,time_range,p,u_init,20,saveat)
                push!(all_60_trials,trialdata...)
            end
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,_=calculate_ITPC(all_60_trials,sr,ITPCrange,C,vsyn,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
        end
        for_storage
    end
end

"""
Will run a series of 20 trial tests for ITPC and correlation to time 2 peak derivative with differing Π/η_0 ratios, on both "Aines" parameters and "my" parameter set. THis is done by changing Π accordingly for fixed η_0. 
It will do so for all the noisestimratio conditions also. Each run through 6 noise ratios takes around an hour.
Π/Π for initial test with aines parameters was: 15/21.5 = 0.70
Π/Π for initial test with my parameters was: 4.0/0.03285 = 121.77
Here I will test 6 ratios between 0.1 to 300 logarithmically spaced as there is an order 100 scale change.
Also will save the ITPC values up to 150Hz rather than just 32.
"""
function run_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,save_path,noisesets)
    
    freq_indexes=[19,37,55] #4, 8 and 12 Hz indexes
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        correlations=Dict{String,Any}()
        println("testing DAR = ",DAR)
        η_0=p.η_0
        p.drive_amplitude=η_0*DAR

        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noisesets)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sqr_correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j].^2) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end

        #make the directory exist, 
        mkpath(save_path)
        tosave=json(correlations)
        open(string(save_path,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

function get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond_1overf_noise!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noisesets,stimulus_type)
    println("testing: ",condition)
    C=p.C
    vsyn=p.vsyn
    saveat=0.0001
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        for idx in 1:3
            all_60_trials=Vector{Any}()
            for noise_filename in noisesets
                p.noise_selector=1
                println("NSR: ",noisestimratio,"Noiseset: ",noise_filename," Test: ",idx)
                flush(stdout)
                # generate_drive_interpolators_specify_noise2stimulus_ratio(stimulus_envelope,noisestimratio,noise_filename)
                # global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                if stimulus_type=="envelope"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                elseif stimulus_type=="derivative"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_smooth_derivative(interpolators_global,44100.0) 
                elseif stimulus_type=="peakrate"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_scaled_peakrate_impulses(interpolators_global,44100.0;scale_range=(0.5,1.0))
                elseif stimulus_type=="peakenvelope"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_scaled_peakenv_impulses(interpolators_global,44100.0;scale_range=(0.5,1.0))
                else
                    error("Unknown stimulus type: $(stimulus_type), must be one of 'envelope', 'derivative', 'peakrate', or 'peakenvelope'.")
                end
                println("running sim")
                flush(stdout)
                trialdata=Ensemble_NoisyPhoneme(vary_noise_and_initial_conditions_NGNMM,time_range,p,u_init,20,saveat)
                push!(all_60_trials,trialdata...)
            end
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,_=calculate_ITPC_1overf_noise(all_60_trials,sr,ITPCrange,C,vsyn,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
        end
        for_storage
    end
end

"""
Will run a series of 20 trial tests for ITPC and correlation to time 2 peak derivative with differing Π/η_0 ratios, on both "Aines" parameters and "my" parameter set. THis is done by changing Π accordingly for fixed η_0. 
It will do so for all the noisestimratio conditions also. Each run through 6 noise ratios takes around an hour.
Π/Π for initial test with aines parameters was: 15/21.5 = 0.70
Π/Π for initial test with my parameters was: 4.0/0.03285 = 121.77
Here I will test 6 ratios between 0.1 to 300 logarithmically spaced as there is an order 100 scale change.
Also will save the ITPC values up to 150Hz rather than just 32.
"""
function run_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond_1overf_noise(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,save_path,noisesets,stimulus_type="envelope")
    
    freq_indexes=[19,37,55] #4, 8 and 12 Hz indexes
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        correlations=Dict{String,Any}()
        println("testing DAR = ",DAR)
        flush(stdout)
        η_0=p.η_0
        p.drive_amplitude=η_0*DAR

        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                flush(stdout)
                get_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond_1overf_noise!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noisesets,stimulus_type)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            flush(stdout)
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sqr_correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j].^2) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end

        #make the directory exist, 
        mkpath(save_path)
        tosave=json(correlations)
        open(string(save_path,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end


### for coupled oscillator control case:
function Ensemble_CoupledOscillators(prob_func,timerange,p,u0,trajectories,saveat)
    model=coupled_oscillator!
    #create ensemble ODE problem using NMM_PhonemeDrive as the model & the given prob_func.
    callback_up=ContinuousCallback(oscillator_phase_reset_condition_up, oscillator_phase_reset_effect_up!)
    callback_down=ContinuousCallback(oscillator_phase_reset_condition_down, oscillator_phase_reset_effect_down!)
    cb=CallbackSet(callback_up,callback_down)
    prob=ODEProblem(model,u0,timerange,p)
    EnsembleProb=EnsembleProblem(prob,prob_func=prob_func)
    #solve for given number of trajectories
    results=solve(EnsembleProb,EnsembleThreads(),alg=Tsit5(),trajectories=trajectories,saveat=saveat,callback=cb)
    return(results)
end
function Ensemble_CoupledOscillators_modulated(prob_func,timerange,p,u0,trajectories,saveat)
    model=coupled_oscillator_modulated!
    #create ensemble ODE problem using NMM_PhonemeDrive as the model & the given prob_func.
    callback_up=ContinuousCallback(oscillator_phase_reset_condition_up, oscillator_phase_reset_effect_up!)
    callback_down=ContinuousCallback(oscillator_phase_reset_condition_down, oscillator_phase_reset_effect_down!)
    cb=CallbackSet(callback_up,callback_down)
    prob=ODEProblem(model,u0,timerange,p)
    EnsembleProb=EnsembleProblem(prob,prob_func=prob_func)
    #solve for given number of trajectories
    results=solve(EnsembleProb,EnsembleThreads(),alg=Tsit5(),trajectories=trajectories,saveat=saveat,callback=cb)
    return(results)
end

function get_coupled_oscillator_ITPC_t2pd_correlation_varynoisestimratio_150Hz!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noise_filename,log_save_name)
    println("testing: ",condition)
    saveat=0.0001 
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        # for (idx,stimulus_envelope) in enumerate(envelopes)
        for idx in 1:3
            println("NSR: ",noisestimratio," Test: ",idx)
            # interpolators[] = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
            global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
            # generate_drive_interpolators_specify_noise2stimulus_ratio(stimulus_envelope,noisestimratio,noise_filename)
            
            #update parameters for this noise condition: DAR so that max stimulus amplitude is 1.0. and set c so that max phase correction is 0.7π
            max_stimulus_amplitude=maximum(maximum.(interpolators_global))
            # p.drive_amplitude=1.0/max_stimulus_amplitude

            #DAR so that stimulus integral is = 20000 (unmodified first vowel stimulus was 20315, so approximating that.)
            # println("mean sum of phoneme envelope", mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global])))
            p.drive_amplitude=20000/mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global]))
            
            p.c=0.7*pi*(1/(max_stimulus_amplitude*p.drive_amplitude))

            println("running sim")
            flush(stdout)
            trialdata=Ensemble_CoupledOscillators(vary_noise,time_range,p,u_init,20,saveat)
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,_=calculate_ITPC_CoupledOscillators_noisyrates(trialdata,sr,ITPCrange,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
        end
        for_storage
    end
end

function run_coupled_oscillator_noise_tests_serial_varydriveamplituderatio(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,savepath,noise_filename)
    
    freq_indexes=[19,37,55]
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    # Condition_keys=["vowel"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        log_save_name=string(parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json")
        correlations=Dict{String,Any}()
        println("testing DAR = N/A, fixed to give max stimulus = 1.0",)
        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                get_coupled_oscillator_ITPC_t2pd_correlation_varynoisestimratio_150Hz!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noise_filename,log_save_name)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end
        tosave=json(correlations)
        open(string(savepath,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

function get_coupled_oscillator_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noisesets)
    println("testing: ",condition)
    saveat=0.0001 
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        # for (idx,stimulus_envelope) in enumerate(envelopes)
        for idx in 1:3
            all_60_trials=Vector{Any}()
            for noise_filename in noisesets
                p.noise_selector=1
                println("NSR: ",noisestimratio,"Noise set: ",noise_filename," Test: ",idx)
                global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                
                #update parameters for this noise condition: DAR so that max stimulus amplitude is 1.0. and set c so that max phase correction is 0.7π
                max_stimulus_amplitude=maximum(maximum.(interpolators_global))
                # p.drive_amplitude=1.0/max_stimulus_amplitude

                #DAR so that stimulus integral is = 20000 (unmodified first vowel stimulus was 20315, so approximating that.)
                # println("mean sum of phoneme envelope", mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global])))
                p.drive_amplitude=20000/mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global]))
                
                p.c=0.7*pi*(1/(max_stimulus_amplitude*p.drive_amplitude*p.q))

                println("running sim")
                flush(stdout)
                trialdata=Ensemble_CoupledOscillators(vary_noise,time_range,p,u_init,20,saveat)
                push!(all_60_trials,trialdata...)
            end
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,_=calculate_ITPC_CoupledOscillators_noisyrates(all_60_trials,sr,ITPCrange,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
        end
        for_storage
    end
end

function run_coupled_oscillator_noise_tests_serial_varydriveamplituderatio_60trials(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,savepath,noisesets)
    
    freq_indexes=[19,37,55]
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    # Condition_keys=["vowel"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        correlations=Dict{String,Any}()
        println("testing DAR = N/A, fixed to give max stimulus = 1.0",)
        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                get_coupled_oscillator_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noisesets)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end
        tosave=json(correlations)
        open(string(savepath,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

function get_coupled_oscillator_modulated_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noisesets)
    println("testing: ",condition)
    saveat=0.0001 
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        # for (idx,stimulus_envelope) in enumerate(envelopes)
        for idx in 1:3
            all_60_trials=Vector{Any}()
            for noise_filename in noisesets
                p.noise_selector=1
                println("NSR: ",noisestimratio,"Noise set: ",noise_filename," Test: ",idx)
                global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                
                #update parameters for this noise condition: DAR so that max stimulus amplitude is 1.0. and set c so that max phase correction is 0.7π
                max_stimulus_amplitude=maximum(maximum.(interpolators_global))
                # p.drive_amplitude=1.0/max_stimulus_amplitude

                #DAR so that stimulus integral is = 20000 (unmodified first vowel stimulus was 20315, so approximating that.)
                # println("mean sum of phoneme envelope", mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global])))
                p.drive_amplitude=20000/mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global]))
                
                p.c=0.7*pi*(1/(max_stimulus_amplitude*p.drive_amplitude*p.q)) #including q. this is the final maximumum stimulus amplitude after all normalisation. 
                                                                            #s(t) normalised to 20000, then s(t)*phasemod normalised to 4.0*20000, now set so max phase correction after all this is 0.7π.

                println("running sim")
                flush(stdout)
                trialdata=Ensemble_CoupledOscillators_modulated(vary_noise,time_range,p,u_init,20,saveat)
                push!(all_60_trials,trialdata...)
            end
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,_=calculate_ITPC_CoupledOscillators_noisyrates(all_60_trials,sr,ITPCrange,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
        end
        for_storage
    end
end

function run_coupled_oscillator_modulated_noise_tests_serial_varydriveamplituderatio_60trials(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,savepath,noisesets)
    
    freq_indexes=[19,37,55]
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    # Condition_keys=["vowel"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        correlations=Dict{String,Any}()
        println("testing DAR = N/A, fixed to give max stimulus = 1.0",)
        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                get_coupled_oscillator_modulated_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noisesets)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end
        tosave=json(correlations)
        open(string(savepath,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

function get_coupled_oscillator_modulated_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noisesets,save_traj,savepath,stimulus_type)
    println("testing: ",condition)
    flush(stdout)
    saveat=0.0001 
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        # for (idx,stimulus_envelope) in enumerate(envelopes)
        for idx in 1:3
            all_60_trials=Vector{Any}()
            for (noise_idx,noise_filename) in enumerate(noisesets)
                p.noise_case_reference=noise_idx
                p.noise_selector=1
                println("NSR: ",noisestimratio,"Noise set: ",noise_filename," Test: ",idx)
                flush(stdout)
                if stimulus_type=="envelope"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                elseif stimulus_type=="derivative"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_smooth_derivative(interpolators_global,44100.0) 
                elseif stimulus_type=="peakrate"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_scaled_peakrate_impulses(interpolators_global,44100.0;scale_range=(0.5,1.0))
                elseif stimulus_type=="peakenvelope"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_scaled_peakenv_impulses(interpolators_global,44100.0;scale_range=(0.5,1.0))
                else
                    error("Unknown stimulus type: $(stimulus_type), must be one of 'envelope', 'derivative', 'peakrate', or 'peakenvelope'.")
                end
                #update parameters for this noise condition: DAR so that max stimulus amplitude is 1.0. and set c so that max phase correction is 0.7π
                max_stimulus_amplitude=maximum(maximum.(interpolators_global))
                # p.drive_amplitude=1.0/max_stimulus_amplitude

                #DAR so that stimulus integral is = 20000 (unmodified first vowel stimulus was 20315, so approximating that.)
                # println("mean sum of phoneme envelope", mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global])))
                p.drive_amplitude=20000/mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global]))
                
                # p.c=0.7*pi*(1/(max_stimulus_amplitude*p.drive_amplitude*p.q)) #including q. this is the final maximumum stimulus amplitude after all normalisation. 
                #s(t) normalised to 20000, then s(t)*phasemod normalised to 4.0*20000, now set so max phase correction after all this is 0.7π.
                #for high c simulation, a max reset that is the full half circle
                p.c=1.0*pi*(1/(max_stimulus_amplitude*p.drive_amplitude*p.q)) #including q. this is the final maximumum stimulus amplitude after all normalisation. 
                
                println("running sim")
                flush(stdout)
                trialdata=Ensemble_CoupledOscillators_modulated(vary_noise_and_initial_conditions,time_range,p,u_init,20,saveat)
                push!(all_60_trials,trialdata...)
            end
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,rates=calculate_ITPC_CoupledOscillators_noisyrates(all_60_trials,sr,ITPCrange,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
            #save the rates for each trial, for comparison across natural frequencies, and computation of congolmerate activity if desired.
            if save_traj
                local name="Trajectories/condition_$(condition)_trial_$(idx)_F$(p.F)Hz_SNR05_c1_60randinitcond_seeded_noise_normalised_nonrectified_modulated_phasemod_$(p.modulation)_integralnormalised_coupled_oscillator_control_test_allnoises_NSR$(noisestimratio)_trajectories"
                CSV.write(savepath*name*".csv",DataFrame(rates,:auto),writeheader=true)
                generate_arrow(name,savepath)
                rm(savepath*name*".csv")
            end
        end
        for_storage
    end
end

function run_coupled_oscillator_modulated_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,savepath,noisesets,save_traj=false,stimulus_type="envelope")
    
    freq_indexes=[19,37,55]
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    # Condition_keys=["vowel"]
    p=p
    #loop over driveamplituderatio
    for DAR in driveamplituderatios
        correlations=Dict{String,Any}()
        println("testing DAR = N/A, fixed to give max stimulus = 1.0",)
        flush(stdout)
        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                flush(stdout)
                get_coupled_oscillator_modulated_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noisesets,save_traj,savepath,stimulus_type)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            flush(stdout)
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end
        tosave=json(correlations)
        open(string(savepath,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

function get_evoked_model_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!(d,condition,p,u_init,time_range,ITPCrange,noisestimratio,Phoneme_Envelopes,noisesets,save_traj,savepath,stimulus_type)
    println("testing: ",condition)
    flush(stdout)
    saveat=0.0001 
    get!(d,condition) do 
        for_storage=Vector{Vector{Float64}}(undef,3)
        # envelopes=Phoneme_Envelopes[condition]
        # for (idx,stimulus_envelope) in enumerate(envelopes)
        for idx in 1:3
            all_60_trials=Vector{Any}()
            for (noise_idx,noise_filename) in enumerate(noisesets)
                p.noise_case_reference=noise_idx
                p.noise_selector=1
                println("NSR: ",noisestimratio,"Noise set: ",noise_filename," Test: ",idx)
                flush(stdout)

                #set stimulus based on type: either normal envelopes, or derivative or event based impulse trains (peak rate or peak envelope)
                if stimulus_type=="envelope"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                elseif stimulus_type=="derivative"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_smooth_derivative(interpolators_global,44100.0) 
                elseif stimulus_type=="peakrate"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_scaled_peakrate_impulses(interpolators_global,44100.0;scale_range=(0.5,1.0))
                elseif stimulus_type=="peakenvelope"
                    global interpolators_global = jldopen("/user/work/as15635/input_data/Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
                    interpolators_global=get_scaled_peakenv_impulses(interpolators_global,44100.0;scale_range=(0.5,1.0))
                else
                    error("Unknown stimulus type: $(stimulus_type), must be one of 'envelope', 'derivative', 'peakrate', or 'peakenvelope'.")
                end

                println("running sim")
                flush(stdout)
                trialdata=Ensemble_EvokedModel(vary_noise_and_initial_conditions_evokedmodel,time_range,p,u_init,20,saveat)
                push!(all_60_trials,trialdata...)
            end
            sr=Int64(1/saveat)
            ITPC,_,_,frequencies,rates=calculate_ITPC_EvokedModel_noisyrates(all_60_trials,sr,ITPCrange,ITPCrange) #ITPCrange can now exclude first 500ms of stimulus as per O.Cucu paper.
            for_storage[idx]=ITPC[1:676] #frequencies up to 150Hz
            #save the rates for each trial, for comparison across natural frequencies, and computation of congolmerate activity if desired.
            if save_traj
                local name="Trajectories/condition_$(condition)_trial_$(idx)_F$(p.F)Hz_SNR05_c1_60randinitcond_seeded_noise_normalised_nonrectified_modulated_phasemod_$(p.modulation)_integralnormalised_coupled_oscillator_control_test_allnoises_NSR$(noisestimratio)_trajectories"
                CSV.write(savepath*name*".csv",DataFrame(rates,:auto),writeheader=true)
                generate_arrow(name,savepath)
                rm(savepath*name*".csv")
            end
        end
        for_storage
    end
end
function vary_noise_and_initial_conditions_evokedmodel(prob,i,repeat)
    prob.p.noise_selector=i
    prob.u0[1]=rand(MersenneTwister(Int64(123+i+20*(prob.p.noise_case_reference-1)))) #random activity x1 in [0,1], but keeping derivative x2 at 0.0 initially.
    prob
end

function Ensemble_EvokedModel(prob_func,timerange,p,u0,trajectories,saveat)
    model=alpha_kernel_ODE!
    prob=ODEProblem(model,u0,timerange,p)
    EnsembleProb=EnsembleProblem(prob,prob_func=prob_func)
    #solve for given number of trajectories
    sr=p.sampling_rate
    results=solve(EnsembleProb,EnsembleThreads(),alg=Tsit5(),trajectories=trajectories,saveat=saveat,tstops=timerange[1]:1/sr:timerange[2])
    return(results)
end

function run_evoked_model_noise_tests_serial_varydriveamplituderatio_60trials_randinitcond(parametersetname,driveamplituderatios,sortingIndexes,time2PDdata,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,savepath,noisesets,save_traj=false,stimulus_type="envelope")
    
    freq_indexes=[19,37,55]
    Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
    # Condition_keys=["vowel"]
    p=p
    for DAR in driveamplituderatios
        p.drive_amplitude=DAR 
        correlations=Dict{String,Any}()
        println("testing DAR = $(DAR), fixed to give max stimulus = 1.0",)
        flush(stdout)
        ITPCs=Dict{String,Vector{Vector{Float64}}}()
        ITPC4812s=Array{Float64}(undef,15,3)
        for_storage=Dict{String,Any}()
        get!(correlations,string("noisestimratio: ",NSR)) do
            for condition in Condition_keys
                println("getting ITPCs for: ",condition)
                flush(stdout)
                get_evoked_model_ITPC_t2pd_correlation_varynoisestimratio_150Hz_60trials_randinitcond!(ITPCs,condition,p,u_init,time_range,ITPCrange,NSR,Phoneme_Envelopes,noisesets,save_traj,savepath,stimulus_type)
            end
            for (idx,condition) in enumerate(Condition_keys)
                for j in 1:3
                    ITPC4812s[idx,j]=mean(ITPCs[condition])[freq_indexes[j]]
                end
            end
            ITPC4812s_and_peakderivtime=hcat(ITPC4812s,time2PDdata)
            sortedITPCs_and_pds=ITPC4812s_and_peakderivtime[sortingIndexes,:]
            println("storing results for NSR=",NSR," test.")
            flush(stdout)
            for_storage=Dict("correlations"=>[cor(sortedITPCs_and_pds[:,4],sortedITPCs_and_pds[:,j]) for j in 1:3],"sortedITPCS_t2PD"=>sortedITPCs_and_pds,"ITPCs"=>ITPCs)
            for_storage
        end
        tosave=json(correlations)
        open(string(savepath,parametersetname,"DAR",round(DAR,sigdigits=4),"NSR",NSR,"data.json"), "w") do f
            write(f,tosave)
        end

    end
end

## Data Analysis Functions

function calculate_ITPC(vector_of_solutions,sampling_rate,timerange,C,vsyn,ITPCrange)
    num_trials=length(vector_of_solutions)
    sr=sampling_rate
    wl=(timerange[2]-timerange[1])*sr
    rates=Array{Vector}(undef,num_trials)
    fourier_coeffs_alltrials=Array{Vector}(undef,num_trials)
    phase_unit_vectors=Array{Vector}(undef,num_trials)
    ITPC=Array{Float64}(undef,Int64((wl/sr)*(sr/2)))
    freqs=collect(0:1.0/(wl/sr):(sr/2))
    for (idx,data) in enumerate(vector_of_solutions)
        rates[idx]=get_firing_rate_NMM(data,C,vsyn)[1][Int64(ITPCrange[1]*sr):Int64(ITPCrange[2]*sr)]#
        fourier_coeffs_alltrials[idx]=rfft(abs.(rates[idx]))
        phase_unit_vectors[idx]=fourier_coeffs_alltrials[idx]./abs.(fourier_coeffs_alltrials[idx])
    end
    #ITPC=abs.(sum(phase_unit_vectors,dims=2))./num_trials
    ITPC=abs.(sum(phase_unit_vectors,dims=1)[1])./num_trials
    return ITPC,fourier_coeffs_alltrials,phase_unit_vectors,freqs,rates
end

function calculate_ITPC_1overf_noise(vector_of_solutions,sampling_rate,timerange,C,vsyn,ITPCrange,noise_seed=123)
    num_trials=length(vector_of_solutions)
    sr=sampling_rate
    wl=(timerange[2]-timerange[1])*sr
    rates=Array{Vector}(undef,num_trials)
    fourier_coeffs_alltrials=Array{Vector}(undef,num_trials)
    phase_unit_vectors=Array{Vector}(undef,num_trials)
    ITPC=Array{Float64}(undef,Int64((wl/sr)*(sr/2)))
    freqs=collect(0:1.0/(wl/sr):(sr/2))
    for (idx,data) in enumerate(vector_of_solutions)
        rates[idx]=get_firing_rate_NMM(data,C,vsyn)[1][Int64(ITPCrange[1]*sr):Int64(ITPCrange[2]*sr)]#

        trialwise_seed=noise_seed+idx
        # noises=PowerLawNoise.noise(1.0, 0.0, length(rates[idx])) 
        noises=seeded_noise(trialwise_seed, 1.0, 0.0, length(rates[idx])) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
        noise_power=mean(noises.^2)
        rate_power=mean(rates[idx].^2)
        # desired_signal_to_noise_ratio=0.1
        desired_signal_to_noise_ratio=0.5 #for high SNR test.
        noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
        scaled_noise=noises.*noise_scaling_factor
        noisy_rates=(rates[idx]).+(scaled_noise)
        rates[idx]=noisy_rates

        fourier_coeffs_alltrials[idx]=rfft(abs.(rates[idx]))
        phase_unit_vectors[idx]=fourier_coeffs_alltrials[idx]./abs.(fourier_coeffs_alltrials[idx])
    end
    #ITPC=abs.(sum(phase_unit_vectors,dims=2))./num_trials
    ITPC=abs.(sum(phase_unit_vectors,dims=1)[1])./num_trials
    return ITPC,fourier_coeffs_alltrials,phase_unit_vectors,freqs,rates
end


function calculate_ITPC_CoupledOscillators_noisyrates(vector_of_solutions,sampling_rate,timerange,ITPCrange,noise_seed=123)
    num_trials=length(vector_of_solutions)
    sr=sampling_rate
    wl=(timerange[2]-timerange[1])*sr
    rates=Array{Vector}(undef,num_trials)
    fourier_coeffs_alltrials=Array{Vector}(undef,num_trials)
    phase_unit_vectors=Array{Vector}(undef,num_trials)
    ITPC=Array{Float64}(undef,Int64((wl/sr)*(sr/2)))
    freqs=collect(0:1.0/(wl/sr):(sr/2))
    for (idx,data) in enumerate(vector_of_solutions)
        # rates[idx]=get_firing_rate_NMM(data,C,vsyn)[1][Int64(ITPCrange[1]*sr):Int64(ITPCrange[2]*sr)]# #from NGNMM implementation
        # rates[idx]=(mean(cos_activity(xytophase(data)),dims=1)[1,:].^2)[Int64(ITPCrange[1]*sr):Int64(ITPCrange[2]*sr)] #squared cos activity, mean across oscillators. the effective firing rate. a vector. #for kuramoto oscillators
        rates[idx]= coupled_oscillator_activity(data)[Int64(ITPCrange[1]*sr):Int64(ITPCrange[2]*sr)]#cos(theta).*r for coupled oscillator activity
        
        trialwise_seed=noise_seed+idx
        # noises=PowerLawNoise.noise(1.0, 0.0, length(rates[idx])) 
        noises=seeded_noise(trialwise_seed, 1.0, 0.0, length(rates[idx])) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
        noise_power=mean(noises.^2)
        rate_power=mean(rates[idx].^2)
        # desired_signal_to_noise_ratio=0.1
        desired_signal_to_noise_ratio=0.5 #for high SNR test.
        noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
        scaled_noise=noises.*noise_scaling_factor
        noisy_rates=(rates[idx]).+(scaled_noise)
        rates[idx]=noisy_rates
        fourier_coeffs_alltrials[idx]=rfft(rates[idx])
        phase_unit_vectors[idx]=fourier_coeffs_alltrials[idx]./abs.(fourier_coeffs_alltrials[idx])
    end
    #ITPC=abs.(sum(phase_unit_vectors,dims=2))./num_trials
    ITPC=abs.(sum(phase_unit_vectors,dims=1)[1])./num_trials
    return ITPC,fourier_coeffs_alltrials,phase_unit_vectors,freqs,rates
end

function calculate_ITPC_EvokedModel_noisyrates(vector_of_solutions,sampling_rate,timerange,ITPCrange,noise_seed=123)
    num_trials=length(vector_of_solutions)
    sr=sampling_rate
    wl=(timerange[2]-timerange[1])*sr
    rates=Array{Vector}(undef,num_trials)
    fourier_coeffs_alltrials=Array{Vector}(undef,num_trials)
    phase_unit_vectors=Array{Vector}(undef,num_trials)
    ITPC=Array{Float64}(undef,Int64((wl/sr)*(sr/2)))
    freqs=collect(0:1.0/(wl/sr):(sr/2))
    for (idx,data) in enumerate(vector_of_solutions)
        rates[idx]= data[1,Int64(ITPCrange[1]*sr):Int64(ITPCrange[2]*sr)] #just the first component of the solution, the 'activity' of the evoked model, similar to global conductance in NGNMM.
        trialwise_seed=noise_seed+idx
        noises=seeded_noise(trialwise_seed, 1.0, 0.0, length(rates[idx])) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
        noise_power=mean(noises.^2)
        rate_power=mean(rates[idx].^2)
        desired_signal_to_noise_ratio=0.5 #for high SNR test.
        noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
        scaled_noise=noises.*noise_scaling_factor
        noisy_rates=(rates[idx]).+(scaled_noise)
        rates[idx]=noisy_rates
        fourier_coeffs_alltrials[idx]=rfft(rates[idx])
        phase_unit_vectors[idx]=fourier_coeffs_alltrials[idx]./abs.(fourier_coeffs_alltrials[idx])
    end
    ITPC=abs.(sum(phase_unit_vectors,dims=1)[1])./num_trials
    return ITPC,fourier_coeffs_alltrials,phase_unit_vectors,freqs,rates
end

function get_synaptic_current_NMM(sol,C,vsyn)
    Z=sol[3,:]
    W=(1 .-conj.(Z))./(1 .+conj.(Z))
    time=sol.t[:]
    g=real.(sol[2,:])
    I_syn=g.*(vsyn.-imag.(W))
    return I_syn,time
end

function get_firing_rate_NMM(sol,C,vsyn)
    Z=sol[3,:]
    W=(1 .-conj.(Z))./(1 .+conj.(Z))
    time=sol.t[:]   
    rate=real.(W)/(pi*C)
    return rate,time
end

function get_frequencies(window_length,sampling_rate)
    return collect(0:1.0/(window_length/sampling_rate):(sampling_rate/2)) #window length in number of samples.
end

"""
Reads the mean_times_to_peak_deriv.csv file and returns a set of indices that order them in ascending order.
To be applied to the condition labels, and the ITPC's accordingly (in run_noise_tests_serial_varydriveamplituderatio) 
"""
function get_sorting_indices(path_to_mean_times_to_peak_deriv)
    time2PDdata=readdlm(path_to_mean_times_to_peak_deriv*"mean_times_to_peak_deriv.csv")
    Q=collect(zip(time2PDdata,1:1:(length(time2PDdata))))
    sortedQ=sort(Q, by=first)
    sortingIndexes=[x[2] for x in sortedQ]
    return sortingIndexes
end

## For evoked model stimuli variations:
"""
    findpeaks(y::Vector{T},threshold) where T <: Real
    peak finder with threshold. returns tuple of (Vector(peaks), Vector(locations)).
"""
function findpeaks(y::Vector{T},threshold) where T
    # Find peaks in the vector y
    peaks = Vector{T}()
    locs = Int[]
    n = length(y)
    for i in eachindex(y[2:n-1])  # Avoid first and last points
        i=i+1  # Adjust index to match original vector
        if y[i] > y[i-1] && y[i] > y[i+1] && y[i] > threshold  # Ensure peak is above threshold
            push!(peaks, y[i])
            push!(locs, i)  # +1 because skipping the first point
        end
    end
    return peaks, locs
end
"""
    get_smooth_derivative(envelopes::Vector{T},sampling_rate::Float64) where T <: Interpolations.Extrapolation
    Calculates the derivative of a set of envelopes after applying a 10Hz low-pass filter. Returns a vector of Interpolations.Extrapolation.
"""
function get_smooth_derivative(envelopes::Vector{T},sampling_rate::Float64) where T <: Interpolations.Extrapolation
    fc = 10 / (sampling_rate / 2)  # normalized cutoff frequency
    b = digitalfilter(Lowpass(fc), Butterworth(4))
    samples=envelopes[1].itp.itp.parentaxes[1][1:end]  
    derivatives = Vector{T}(undef, length(envelopes))
     for (i, env) in enumerate(envelopes)
        #Sample the envelope
        y = env(samples)
        #smooth the envelope using a low-pass filter
        y_smooth = filtfilt(b, y)
        #Compute the discrete derivative (multiply by sr for units/s)
        dy = diff(y_smooth) * sampling_rate
        # Interpolate the derivative (grid is now 1:n_samples-1)
        derivatives[i] = linear_interpolation(StepRange(samples[1:end-1]),Vector(dy),extrapolation_bc=Line())
    end
    return derivatives
end

"""
    get_scaled_peakrate_impulses(envelopes::Vector{T}, sampling_rate::Float64; scale_range=(0.5,1.0)) where T <: Interpolations.Extrapolation
    Takes in  envelopes, turns them into smooth derivatives, and then finds the peaks in the derivatives.
    Returns a vector of Interpolations.Extrapolation of impulse peak-rate events. originally of magnitude = rate.
    but then scaled to [0.5,1.0] as in "Oganian Y et al. Phase Alignment of Low-Frequency Neural Activity to the Amplitude Envelope of Speech Reflects Evoked Responses to Acoustic Edges, Not Oscillatory Entrainment. J Neurosci. 2023. doi: 10.1523/JNEUROSCI.1663-22.2023"
"""
function get_scaled_peakrate_impulses(envelopes::Vector{T}, sampling_rate::Float64;scale_range=(0.5,1.0)) where T <: Interpolations.Extrapolation
    smooth_derivatives = get_smooth_derivative(envelopes, sampling_rate)
    sample_grid= envelopes[1].itp.itp.parentaxes[1][1:end]  # Get the time grid from the first envelope
    peakrate_impulses = Vector{T}(undef, length(envelopes))
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
    end

    return peakrate_impulses
end

"""
    get_scaled_peakenv_impulses(envelopes::Vector{T}, sampling_rate::Float64; scale_range=(0.5, 1.0)) where T <: Interpolations.Extrapolation
    Takes in envelopes, smooths them with a low-pass filter, finds the peaks in the smoothed envelopes,
    and returns a vector of Interpolations.Extrapolation of impulse peak-envelope events.
    The peak magnitudes are scaled to the range [0.5, 1.0].
"""
function get_scaled_peakenv_impulses(envelopes::Vector{T}, sampling_rate::Float64; scale_range=(0.5, 1.0)) where T <: Interpolations.Extrapolation

    fc = 10 / (sampling_rate / 2)  # normalized cutoff frequency
    b = digitalfilter(Lowpass(fc), Butterworth(4))
    sample_grid= envelopes[1].itp.itp.parentaxes[1][1:end]  # Get the time grid from the first envelope
    smooth_envelopes=[filtfilt(b, env(sample_grid)) for env in envelopes]  # Smooth the envelopes using a low-pass filter
    peakenv_impulses = Vector{T}(undef, length(envelopes))
    for (i, envelope) in enumerate(smooth_envelopes)
        # Create a zero vector of the same length as the envelope
        peaks,locations = findpeaks(envelope, 0.1)  # Threshold set to 0.1

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
        peakenv_impulses[i] = linear_interpolation(StepRange(sample_grid), impulse_vector, extrapolation_bc=Line());
    end

    return peakenv_impulses
end


#for control cases
function cos_activity(phase_data)
    return cos.(phase_data)
end
function coupled_oscillator_activity(solution)
    return cos.(solution[1,:]).*solution[2,:]
end


"""
    seeded_noise(seed, β, ν₀, dims...)
    
A wrapper around PowerLawNoise.noise that uses a specific random seed
for reproducible noise generation.

Parameters:
- seed: Integer seed for reproducible noise generation
- β: Power law exponent (1.0 for pink noise, 0 for white noise, etc.)
- ν₀: Minimum frequency
- dims: Dimensions of the noise array

Returns:
- An array of noise values with the specified power law distribution
~-Claude 3.7 Sonnet.
"""
function seeded_noise(seed, β, ν₀, dims...)

    # Set the seed for reproducibility
    Random.seed!(seed)
    # Generate the noise with the seeded RNG
    noise_data = PowerLawNoise.noise(β, ν₀, dims...)
    #reset seed
    Random.seed!()
    return noise_data
end

#From David Barton:
"""
Convert a csv file to compressed arrow format. Courtesy of Professor David Barton, Engineering Maths
"""
function generate_arrow(name, data_path; force = false, compress = true)
    filename = joinpath(data_path, name)
    arrow_file = filename * ".arrow"
    if compress
        arrow_file = arrow_file * ".lz4"
    end
    csv_file = filename * ".csv"
    if !isfile(arrow_file) || force
        data = CSV.read(csv_file, DataFrame; normalizenames = true)
        if compress
            Arrow.write(arrow_file, data, compress = :lz4)
        else
            Arrow.write(arrow_file, data)
        end
    end
    return nothing
end

end