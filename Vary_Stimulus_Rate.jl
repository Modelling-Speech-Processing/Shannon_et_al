using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using ComponentArrays, OrdinaryDiffEq, Plots, Parameters, JLD2, ProgressLogging
#to vary the stimulus rate I need a new version of each model that scales time in its 
#call of the stimulus interpolation object.
condition="b" #which phoneme condition to use for the drive
idx=1 #which test index to use
noisestimratio=0.3 #which noise to stimulus ratio to use
DAR=12.2
noise_filename="NoiseSequences60.csv" #which noise file to use

global interpolators_global = jldopen("./Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
#testing that here first:

function NMM_PhonemeDrive_Noisy_VaryStimRate(D,u,p,t)
    @unpack α,k,C,Δ,η_0,vsyn,α_D, sampling_rate, drive_amplitude, noise_selector, noise_case_reference, τ = p
    @unpack g_dot, g, Z, A_dot, A = u #A is drive 
    
    #u[1] is g', 
    D.g_dot = (α^2)*(k/(C*pi)) * ((1-(abs(Z)^2))/(1+Z+conj(Z)+(abs(Z)^2))) - (2*α)*g_dot - (α^2)*g
    
    #u[2] is g
    D.g = g_dot
       
    #u[3] is Z
    D.Z = (1/C) * (-im*(((Z-1)^2)/2) + (((Z+1)^2)/2)*(-Δ+im*(η_0+A)+im*vsyn*g) - (((Z^2)-1)/2)*g)

    #interpolators[NoiseSelector] chooses a given interpolator out of the previously constructed global interpolator vector (containing all the noise conditions and the given stimulus)
    # τ is the stimulus rate scaling factor.
    D.A_dot=(α_D^2).*interpolators_global[Int64(noise_selector)](τ * t * sampling_rate+1.0)*drive_amplitude-2*α_D*A_dot-(α_D^2)*A 
    
    D.A=A_dot
end

function NMM_PhonemeDrive_Noisy_VaryStimRate(D,u,p,t)
    @unpack α,k,C,Δ,η_0,vsyn,α_D, sampling_rate, drive_amplitude, noise_selector, noise_case_reference, τ = p
    @unpack g_dot, g, Z, A_dot, A = u 
    
    # Define the duration of each phase in the original data (5 seconds)
    T_orig = 5.0

    if t < 5.0
        # --- PHASE 1: Noise (0.0 to 5.0 in interpolator) ---
        # Scale time by τ and wrap within the [0, 5] range
        t_lookup = mod(τ * t, T_orig)
    else
        # --- PHASE 2: Stimulus (5.0 to 10.0 in interpolator) ---
        # Calculate time elapsed since the start of Phase 2
        Δt = t - 5.0
        # Scale elapsed time by τ and wrap within the [0, 5] range
        # Then shift by 5.0 to point to the Stimulus section of the interpolator
        t_lookup = 5.0 + mod(τ * Δt, T_orig)
    end

    # Convert to index (preserving your +1.0 offset convention)
    idx = t_lookup * sampling_rate + 1.0
    
    # --- ODE Dynamics ---
    # (Equations for g'', g', and Z remain unchanged)
    D.g_dot = (α^2)*(k/(C*pi)) * ((1-(abs(Z)^2))/(1+Z+conj(Z)+(abs(Z)^2))) - (2*α)*g_dot - (α^2)*g
    D.g = g_dot
    D.Z = (1/C) * (-im*(((Z-1)^2)/2) + (((Z+1)^2)/2)*(-Δ+im*(η_0+A)+im*vsyn*g) - (((Z^2)-1)/2)*g)

    # --- Drive Dynamics ---
    drive_input = interpolators_global[Int64(noise_selector)](idx)
    D.A_dot = (α_D^2) * drive_input * drive_amplitude - 2*α_D*A_dot - (α_D^2)*A 
    D.A = A_dot
end


function Ensemble_NoisyPhoneme_VaryStimRate(prob_func,timerange,p,u0,trajectories,saveat)
    model=NMM_PhonemeDrive_Noisy_VaryStimRate
    #create ensemble ODE problem using NMM_PhonemeDrive as the model & the given prob_func.
    prob=ODEProblem(model,u0,timerange,p)
    EnsembleProb=EnsembleProblem(prob,prob_func=prob_func)

    #solve for given number of trajectories
    results=solve(EnsembleProb,EnsembleThreads(),trajectories=trajectories,saveat=saveat)
    return(results)
end

prob_func=vary_noise_and_initial_conditions_NGNMM


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
drive_amplitude=η_0*DAR
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,10.0)
saveat=0.0001   



taus=collect(range(0.25,5.0,length=20))
stim_rates=[4*tau for tau in taus] #in Hz
freq_idxs=Vector{Int64}(undef,length(stim_rates)) #to store after first ITPC calculation. freq per tau. 
stimrate_ITPCs=Vector{Float64}(undef,length(taus))
fourHz_ITPCs=Vector{Float64}(undef,length(taus))
fourHzidx=19
@progress for (tidx,τ) in enumerate(taus)
    p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1, τ=τ) 
    results=Ensemble_NoisyPhoneme_VaryStimRate(prob_func,time_range,p,u0,20,saveat)
    ITPC,_,_,freqs,_=calculate_ITPC_1overf_noise(results,1/saveat,(5.5,10.0),C,vsyn,(5.5,10.0))
    if tidx==1
        for j in eachindex(stim_rates)
            freq_idxs[j]=Int64(findmin(abs.(freqs .- stim_rates[j]))[2])
        end
    end
    stimrate_ITPCs[tidx]=ITPC[freq_idxs[tidx]].^2
    fourHz_ITPCs[tidx]=ITPC[fourHzidx].^2
end
slow_itpc_plot=bar(stim_rates,stimrate_ITPCs,xlabel="Stimulus Rate (Hz)",ylabel="ITPC",label="stimrate ITPC",title="4Hz NGNMM 'b'-Stim",size=(800,600),alpha=0.5);
bar!(slow_itpc_plot,stim_rates,fourHz_ITPCs,xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5);

### and for the fast NGNMM too!
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
drive_amplitude=η_0*DAR
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,10.0)
ITPCrange=(5.5,10.0)

stimrate_ITPCs_fastNGNMM=Vector{Float64}(undef,length(taus))
fourHz_ITPCs_fastNGNMM=Vector{Float64}(undef,length(taus))

@progress for (tidx,τ) in enumerate(taus)
    p_fast_NGNMM=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1, τ=τ) 
    results=Ensemble_NoisyPhoneme_VaryStimRate(prob_func,time_range,p_fast_NGNMM,u0,20,saveat)
    ITPC,_,_,freqs,_=calculate_ITPC_1overf_noise(results,1/saveat,(5.5,10.0),C,vsyn,(5.5,10.0))
    if tidx==1
        for j in eachindex(stim_rates)
            freq_idxs[j]=Int64(findmin(abs.(freqs .- stim_rates[j]))[2])
        end
    end
    stimrate_ITPCs_fastNGNMM[tidx]=ITPC[freq_idxs[tidx]].^2
    fourHz_ITPCs_fastNGNMM[tidx]=ITPC[fourHzidx].^2
end


fast_itpc_plot=bar(stim_rates,stimrate_ITPCs_fastNGNMM,xlabel="Stimulus Rate (Hz)",ylabel="ITPC",label="Stimulus Rate ITPC",title="15Hz NGNMM 'b'-Stim",size=(800,600),alpha=0.5);
bar!(fast_itpc_plot,stim_rates,fourHz_ITPCs_fastNGNMM,xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5);
plot(slow_itpc_plot,fast_itpc_plot,layout=(2,1),size=(900,1200),ylims=(0.0,1.0),ylabel="ITPC^2",margin=5Plots.mm)

#now the same for the evoked and phase-resetting models too:
# has a discontinuity when modulation makes the stimulus correction always positive.
function coupled_oscillator_modulated_varystimrate!(du, u, p, t)
    @unpack F, c, drive_amplitude, noise_selector, sampling_rate,modulation,q, noise_case_reference, τ = p
    @unpack θ, r = u

    # Define the duration of each phase in the original data (5 seconds)
    T_orig = 5.0

    if t < 5.0
        # --- PHASE 1: Noise (0.0 to 5.0 in interpolator) ---
        # Scale time by τ and wrap within the [0, 5] range
        t_lookup = mod(τ * t, T_orig)
    else
        # --- PHASE 2: Stimulus (5.0 to 10.0 in interpolator) ---
        # Calculate time elapsed since the start of Phase 2
        Δt = t - 5.0
        # Scale elapsed time by τ and wrap within the [0, 5] range
        # Then shift by 5.0 to point to the Stimulus section of the interpolator
        t_lookup = 5.0 + mod(τ * Δt, T_orig)
    end

    # Convert to index (preserving your +1.0 offset convention)
    idx = t_lookup * sampling_rate + 1.0

    drive_input = interpolators_global[Int64(noise_selector)](idx)

    du.θ = 2*pi*F - c * q * ((drive_input*drive_amplitude)/r) * (1*(1-modulation)+modulation*sin(θ)) #modulation of 1.0 means full phase modulation. 0.0 = no phase modulation.
    du.r = r*(1-r^2) + c * q * (drive_input*drive_amplitude) * (1*(1-modulation)+modulation*cos(θ)) 
    return nothing
end


"""
Equivalent to alpha-kernel convolution, ODE form, using the synaptic filter from the NGNMM. 2nd Order ODE.
"""
function alpha_kernel_ODE_varystimrate!(du,u,p,t)
    @unpack α,drive_amplitude, noise_selector,sampling_rate, noise_case_reference, τ = p
    @unpack x1,x2 = u

    # Define the duration of each phase in the original data (5 seconds)
    T_orig = 5.0

    if t < 5.0
        # --- PHASE 1: Noise (0.0 to 5.0 in interpolator) ---
        # Scale time by τ and wrap within the [0, 5] range
        t_lookup = mod(τ * t, T_orig)
    else
        # --- PHASE 2: Stimulus (5.0 to 10.0 in interpolator) ---
        # Calculate time elapsed since the start of Phase 2
        Δt = t - 5.0
        # Scale elapsed time by τ and wrap within the [0, 5] range
        # Then shift by 5.0 to point to the Stimulus section of the interpolator
        t_lookup = 5.0 + mod(τ * Δt, T_orig)
    end

    # Convert to index (preserving your +1.0 offset convention)
    idx = t_lookup * sampling_rate + 1.0
    
    drive_input = interpolators_global[Int64(noise_selector)](idx)


    #x1 = synaptic current, x2 = synaptic current derivative
    du.x1 = x2
    du.x2 = (α^2).*drive_input*drive_amplitude-2*α*x2-(α^2)*x1 
    return nothing
end

function Ensemble_CoupledOscillators_modulated_varystimrate(prob_func,timerange,p,u0,trajectories,saveat)
    model=coupled_oscillator_modulated_varystimrate!
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

function Ensemble_EvokedModel_varystimrate(prob_func,timerange,p,u0,trajectories,saveat)
    model=alpha_kernel_ODE_varystimrate!
    prob=ODEProblem(model,u0,timerange,p)
    EnsembleProb=EnsembleProblem(prob,prob_func=prob_func)
    #solve for given number of trajectories
    sr=p.sampling_rate
    results=solve(EnsembleProb,EnsembleThreads(),alg=Tsit5(),trajectories=trajectories,saveat=saveat,tstops=timerange[1]:1/sr:timerange[2])
    return(results)
end



# phase resetting oscillator case:
phase_modulation=1 
# set parameters
u0=ComponentArray(θ=0.0, r=1.0)
max_stim=maximum(maximum.(interpolators_global))
drive_amplitude=20000/mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global]))
c=0.7*pi*(1/(max_stim*drive_amplitude))

time_range=(0.0,10.0)
F=4.0
phoneme_sampling_rate=44100
#for rectified version
# q_normalisation=4/(2*pi-2*pi*phase_modulation+2*phase_modulation) #scaling factor to make absolute area under stimulus modulation curve = 4. 
#for non rectified
if phase_modulation<=0.5
    q_normalisation=2/((1-phase_modulation)*π) #scaling factor to make absolute area under stimulus modulation curve = 4. 
else
    m=phase_modulation
    q_normalisation=4/(4*m*sqrt(1-((m-1)/m)^2)+2*(m-1)asin((m-1)/m))
end
ITPCrange=(5.5,10.0)

stimrate_ITPCs_phasereset=Vector{Float64}(undef,length(taus))
fourHz_ITPCs_phasereset=Vector{Float64}(undef,length(taus))
prob_func=vary_noise_and_initial_conditions
@progress for (tidx,τ) in enumerate(taus)
    @info tidx
    p=ComponentArray(F=F, c=c, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,modulation=phase_modulation,q=q_normalisation,noise_case_reference=1,τ=τ)
    results=Ensemble_CoupledOscillators_modulated_varystimrate(prob_func,time_range,p,u0,20,saveat)
    ITPC,_,_,freqs,_=calculate_ITPC_CoupledOscillators_noisyrates(results,1/saveat,(5.5,10.0),(5.5,10.0))
    if tidx==1
        for j in eachindex(stim_rates)
            freq_idxs[j]=Int64(findmin(abs.(freqs .- stim_rates[j]))[2])
        end
    end
    stimrate_ITPCs_phasereset[tidx]=ITPC[freq_idxs[tidx]].^2
    fourHz_ITPCs_phasereset[tidx]=ITPC[fourHzidx].^2
end

phasereset_itpc_plot=bar(stim_rates,stimrate_ITPCs_phasereset,xlabel="Stimulus Rate (Hz)",ylabel="ITPC",label="Stimulus Rate ITPC",title="UnTuned Phase Resetting Model 'b'-Stim",size=(800,600),alpha=0.5);
bar!(phasereset_itpc_plot,stim_rates,fourHz_ITPCs_phasereset,xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,ylims=(0.0,1.0));

#and with frequency tuning:
stimrate_ITPCs_tunedphasereset=Vector{Float64}(undef,length(taus))
fourHz_ITPCs_tunedphasereset=Vector{Float64}(undef,length(taus))
prob_func=vary_noise_and_initial_conditions
@progress for (tidx,τ) in enumerate(taus)
    @info tidx
    p=ComponentArray(F=F, c=c, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,modulation=phase_modulation,q=q_normalisation,noise_case_reference=1,τ=τ)
    p.F=stim_rates[tidx] #tune each oscillator to the stimulus rate
    results=Ensemble_CoupledOscillators_modulated_varystimrate(prob_func,time_range,p,u0,20,saveat)
    ITPC,_,_,freqs,_=calculate_ITPC_CoupledOscillators_noisyrates(results,1/saveat,(5.5,10.0),(5.5,10.0))
    if tidx==1
        for j in eachindex(stim_rates)
            freq_idxs[j]=Int64(findmin(abs.(freqs .- stim_rates[j]))[2])
        end
    end
    stimrate_ITPCs_tunedphasereset[tidx]=ITPC[freq_idxs[tidx]].^2
    fourHz_ITPCs_tunedphasereset[tidx]=ITPC[fourHzidx].^2
end

tunedphasereset_itpc_plot=bar(stim_rates,stimrate_ITPCs_tunedphasereset,xlabel="Stimulus Rate (Hz)",ylabel="ITPC",label="Stimulus Rate ITPC",title="Tuned Phase Resetting Model 'b'-Stim",size=(800,600),alpha=0.5);
bar!(tunedphasereset_itpc_plot,stim_rates,fourHz_ITPCs_tunedphasereset,xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,ylims=(0.0,1.0));


#evoked model:
u0=ComponentArray(x1=0.0, x2=0.0)
time_range=(0.0,10.0)
phoneme_sampling_rate=44100
drive_amplitude=4.0 #will be updated in the test function to make drive amplitude equal to DAR setting.
α=1/0.03
ITPCrange=(5.5,10.0)

stimrate_ITPCs_evoked=Vector{Float64}(undef,length(taus))
fourHz_ITPCs_evoked=Vector{Float64}(undef,length(taus))
prob_func=vary_noise_and_initial_conditions_evokedmodel
@progress for (tidx,τ) in enumerate(taus)
    @info tidx
    p=ComponentArray(α=α, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,noise_case_reference=1, τ=τ)
    results=Ensemble_EvokedModel_varystimrate(prob_func,time_range,p,u0,20,saveat)
    ITPC,_,_,freqs,_=calculate_ITPC_EvokedModel_noisyrates(results,1/saveat,(5.5,10.0),(5.5,10.0))
    if tidx==1
        for j in eachindex(stim_rates)
            freq_idxs[j]=Int64(findmin(abs.(freqs .- stim_rates[j]))[2])
        end
    end
    stimrate_ITPCs_evoked[tidx]=ITPC[freq_idxs[tidx]].^2
    fourHz_ITPCs_evoked[tidx]=ITPC[fourHzidx].^2
end

evoked_itpc_plot=bar(stim_rates,stimrate_ITPCs_evoked,xlabel="Stimulus Rate (Hz)",ylabel="ITPC",label="Stimulus Rate ITPC",title="Evoked Model 'b'-Stim",size=(800,600),alpha=0.5);
bar!(evoked_itpc_plot,stim_rates,fourHz_ITPCs_evoked,xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,ylims=(0.0,1.0));
plot(phasereset_itpc_plot,tunedphasereset_itpc_plot,evoked_itpc_plot,layout=(3,1),size=(900,1800),ylabel="ITPC^2",margin=5Plots.mm)

## all plotted together:
plot(slow_itpc_plot,fast_itpc_plot,phasereset_itpc_plot,tunedphasereset_itpc_plot,evoked_itpc_plot,layout=(3,2),size=(1500,1500),ylabel="ITPC^2",margin=5Plots.mm,plot_title="0.3 noise to stim ratio")


#and scaled so that each plot has its maximum at 1.0 for the stimulus rate case, use the same factor for the 4Hz Case too.
slow_itpc_plot_norm=bar(stim_rates,stimrate_ITPCs./maximum(stimrate_ITPCs),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="stimrate ITPC",title="4Hz NGNMM 'b'-Stim",size=(800,600),alpha=0.5);
bar!(slow_itpc_plot_norm,stim_rates,fourHz_ITPCs./maximum(stimrate_ITPCs),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5);
fast_itpc_plot_norm=bar(stim_rates,stimrate_ITPCs_fastNGNMM./maximum(stimrate_ITPCs_fastNGNMM),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="15Hz NGNMM 'b'-Stim",size=(800,600),alpha=0.5);
bar!(fast_itpc_plot_norm,stim_rates,fourHz_ITPCs_fastNGNMM./maximum(stimrate_ITPCs_fastNGNMM),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,ylims=(0.0,1.0));
phasereset_itpc_plot_norm=bar(stim_rates,stimrate_ITPCs_phasereset./maximum(stimrate_ITPCs_phasereset),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="UnTuned Phase Resetting Model 'b'-Stim",size=(800,600),alpha=0.5);
bar!(phasereset_itpc_plot_norm,stim_rates,fourHz_ITPCs_phasereset./maximum(stimrate_ITPCs_phasereset),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,ylims=(0.0,1.0));
tunedphasereset_itpc_plot_norm=bar(stim_rates,stimrate_ITPCs_tunedphasereset./maximum(stimrate_ITPCs_tunedphasereset),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="Tuned Phase Resetting Model 'b'-Stim",size=(800,600),alpha=0.5);
bar!(tunedphasereset_itpc_plot_norm,stim_rates,fourHz_ITPCs_tunedphasereset./maximum(stimrate_ITPCs_tunedphasereset),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,ylims=(0.0,1.0));
evoked_itpc_plot_norm=bar(stim_rates,stimrate_ITPCs_evoked./maximum(stimrate_ITPCs_evoked),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="Evoked Model 'b'-Stim",size=(800,600),alpha=0.5);
bar!(evoked_itpc_plot_norm,stim_rates,fourHz_ITPCs_evoked./maximum(stimrate_ITPCs_evoked),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,ylims=(0.0,1.0));
plot(slow_itpc_plot_norm,fast_itpc_plot_norm,phasereset_itpc_plot_norm,tunedphasereset_itpc_plot_norm,evoked_itpc_plot_norm,layout=(3,2),size=(1500,1500),ylabel="Normalized ITPC",margin=5Plots.mm,plot_title="0.3 noise to stim ratio (Normalized)")



#do a line plot of the normalised bars
slow_itpc_plot_line=plot(stim_rates,stimrate_ITPCs./maximum(stimrate_ITPCs),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="stimrate ITPC",title="4Hz NGNMM 'b'-Stim",size=(800,600),alpha=0.5,marker=:circle);
plot!(slow_itpc_plot_line,stim_rates,fourHz_ITPCs./maximum(stimrate_ITPCs),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,marker=:circle);
fast_itpc_plot_line=plot(stim_rates,stimrate_ITPCs_fastNGNMM./maximum(stimrate_ITPCs_fastNGNMM),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="15Hz NGNMM 'b'-Stim",size=(800,600),alpha=0.5,marker=:circle);
plot!(fast_itpc_plot_line,stim_rates,fourHz_ITPCs_fastNGNMM./maximum(stimrate_ITPCs_fastNGNMM),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,marker=:circle,ylims=(0.0,1.0));
phasereset_itpc_plot_line=plot(stim_rates,stimrate_ITPCs_phasereset./maximum(stimrate_ITPCs_phasereset),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="UnTuned Phase Resetting Model 'b'-Stim",size=(800,600),alpha=0.5,marker=:circle);
plot!(phasereset_itpc_plot_line,stim_rates,fourHz_ITPCs_phasereset./maximum(stimrate_ITPCs_phasereset),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,marker=:circle,ylims=(0.0,1.0));
tunedphasereset_itpc_plot_line=plot(stim_rates,stimrate_ITPCs_tunedphasereset./maximum(stimrate_ITPCs_tunedphasereset),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="Tuned Phase Resetting Model 'b'-Stim",size=(800,600),alpha=0.5,marker=:circle);
plot!(tunedphasereset_itpc_plot_line,stim_rates,fourHz_ITPCs_tunedphasereset./maximum(stimrate_ITPCs_tunedphasereset),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,marker=:circle,ylims=(0.0,1.0));
evoked_itpc_plot_line=plot(stim_rates,stimrate_ITPCs_evoked./maximum(stimrate_ITPCs_evoked),xlabel="Stimulus Rate (Hz)",ylabel="Normalized ITPC",label="Stimulus Rate ITPC",title="Evoked Model 'b'-Stim",size=(800,600),alpha=0.5,marker=:circle);
plot!(evoked_itpc_plot_line,stim_rates,fourHz_ITPCs_evoked./maximum(stimrate_ITPCs_evoked),xlabel="Stimulus Rate (Hz)",label="4Hz ITPC",alpha=0.5,marker=:circle,ylims=(0.0,1.0));
plot(slow_itpc_plot_line,fast_itpc_plot_line,phasereset_itpc_plot_line,tunedphasereset_itpc_plot_line,evoked_itpc_plot_line,layout=(3,2),size=(1500,1500),ylabel="Normalized ITPC",margin=5Plots.mm,plot_title="0.3 noise to stim ratio (Normalized)")

