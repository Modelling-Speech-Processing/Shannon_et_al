using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using ComponentArrays, OrdinaryDiffEq, Plots, Parameters, JLD2
#to vary the stimulus rate I need a new version of each model that scales time in its 
#call of the stimulus interpolation object.
condition="vowel" #which phoneme condition to use for the drive
idx=1 #which test index to use
noisestimratio=0.7 #which noise to stimulus ratio to use
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
drive_amplitude=Π
τ=1.0
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1, τ=τ) 
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,10.0)
saveat=0.0001   

results=Ensemble_NoisyPhoneme_VaryStimRate(prob_func,time_range,p,u0,20,saveat)


#and with faster stimulus rate
τ=4.0
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1, τ=τ) 

results_faststim=Ensemble_NoisyPhoneme_VaryStimRate(prob_func,time_range,p,u0,20,saveat)

data_zsync=abs.(results[1][3,:])
data_zsync_faststim=abs.(results_faststim[1][3,:])
p_norm=plot(range(0.0,10.0,length=length(data_zsync)),data_zsync,label="Normal Stimulus Rate (4Hz) _ scaled",xlabel="Time (s)",ylabel="|Z| Synchrony",legend=:bottomright)
p_fast=plot(range(0.0,10.0,length=length(data_zsync_faststim)),data_zsync_faststim,label="4x Faster Stimulus Rate (16Hz) _ scaled")

τ=0.25
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1, τ=τ)
results_slowstim=Ensemble_NoisyPhoneme_VaryStimRate(prob_func,time_range,p,u0,20,saveat)
data_zsync_slowstim=abs.(results_slowstim[1][3,:])
p_slow=plot(range(0.0,10.0,length=length(data_zsync_slowstim)),data_zsync_slowstim,label="0.25x Slower Stimulus Rate (1Hz) _ scaled")

plot(p_norm,p_fast,p_slow,layout=(3,1),size=(800,900),title="NMM Phoneme Drive with Noise - Varying Stimulus Rate",xlabel="Time (s)",ylabel="|Z| Synchrony")

ITPC_norm=calculate_ITPC(results, 1/saveat, (5.5,10.0),C,vsyn,(5.5,10.0))
ITPC_fast=calculate_ITPC(results_faststim, 1/saveat, (5.5,10.0),C,vsyn,(5.5,10.0))
ITPC_slow=calculate_ITPC(results_slowstim, 1/saveat, (5.5,10.0),C,vsyn,(5.5,10.0))
freqs=ITPC_norm[4]
#[676] is 150Hz
#[19] is 4Hz
plot(freqs[1:676],ITPC_norm[1][1:676],label="Normal Stimulus Rate (4Hz)",xlabel="Frequency (Hz)",ylabel="ITPC",legend=:topright)
plot!(freqs[1:676],ITPC_fast[1][1:676],label="4x Faster Stimulus Rate (16Hz)")
plot!(freqs[1:676],ITPC_slow[1][1:676],label="0.25x Slower Stimulus Rate (1Hz)")
ITPC_norm[2]

bar([ITPC_norm[1][19],ITPC_fast[1][19],ITPC_slow[1][19]],label=["Normal Stimulus Rate (4Hz)" "4x Faster Stimulus Rate (16Hz)" "0.25x Slower Stimulus Rate (1Hz)"],xlabel="Frequency (Hz)",ylabel="ITPC",legend=:topright,ylims=(0.0,1.0))

#what are the 1 Hz and 16Hz idxs?
freqs[19] #4Hz
freqs[6] #1Hz (= 1.1111Hz)
freqs[73] #16Hz

bar([4,16,1.1111],[ITPC_norm[1][19],ITPC_fast[1][73],ITPC_slow[1][6]],label=["Normal Stimulus Rate (4Hz)" "4x Faster Stimulus Rate (16Hz)" "0.25x Slower Stimulus Rate (1Hz)"],xlabel="Frequency (Hz)",ylabel="ITPC at stim rate",legend=:topright,ylims=(0.0,1.0))