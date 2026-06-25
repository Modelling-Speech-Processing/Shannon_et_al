using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using Statistics
using ProgressLogging
using ComponentArrays, OrdinaryDiffEq, Plots, Parameters, JLD2
using LaTeXStrings
using Interpolations, DSP
Plots.default(titlefontsize=16,legendfontsize=12,tickfontsize=9,guidefontsize=11)


condition="b" #which phoneme condition to use for the drive
idx=1 #which test index to use
noisestimratio=0.0#which noise to stimulus ratio to use
noise_filename="NoiseSequences60.csv" #which noise file to use

global NGNMM_NSP_paper_code.interpolators_global = jldopen("./Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]

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
time_range=(0.0,20.0)
saveat=0.0001   

prob_func=vary_noise_and_initial_conditions_NGNMM

#callback to silence the stimulus after 10 seconds (the end of the stimulus period)
#then we can see any ringing effect from the final phoneme. 

function callback_condition(u, t, integrator)
    t==10.0
end

function callback_affect!(integrator)
    integrator.p.drive_amplitude=0.0
end

cb=DiscreteCallback(callback_condition, callback_affect!)

tstops=[10.0]

function Ensemble_NoisyPhoneme_with_callback(prob_func,timerange,p,u0,trajectories,saveat,callback_func,tstops)
    model=NMM_PhonemeDrive_Noisy
    #create ensemble ODE problem using NMM_PhonemeDrive as the model & the given prob_func.
    prob=ODEProblem(model,u0,timerange,p)
    EnsembleProb=EnsembleProblem(prob,prob_func=prob_func)

    #solve for given number of trajectories
    results=solve(EnsembleProb,EnsembleThreads(),trajectories=trajectories,saveat=saveat,callback=callback_func,tstops=tstops)
    return(results)
end
results=Ensemble_NoisyPhoneme_with_callback(prob_func,time_range,p,u0,20,saveat,cb,tstops)

#WONT WORK AS locs is just for the first stimulus stream. 
#getting more results with other phoneme conditions to get better IEPC averages.
global NGNMM_NSP_paper_code.interpolators_global = jldopen("./Phoneme_Drives/drive_interpolators_b_test_2_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
results_b_2=Ensemble_NoisyPhoneme_with_callback(prob_func,time_range,p,u0,20,saveat,cb,tstops)
global NGNMM_NSP_paper_code.interpolators_global = jldopen("./Phoneme_Drives/drive_interpolators_b_test_3_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]
results_b_3=Ensemble_NoisyPhoneme_with_callback(prob_func,time_range,p,u0,20,saveat,cb,tstops)
results=vcat(results,results_b_2,results_b_3); #combine results from the 3 different phoneme conditions, to get more events for the IEPC calculation.


## with 1/f noise on the firing rates rather than abs(Z).
noisy_rates = Vector{Vector{Float64}}(undef, 60)
for (j,res) in enumerate(results)
    for i in 1:20
        rates=get_firing_rate_NMM(res[i],p.C,p.vsyn)[1]
            
        trialwise_seed=1000+i #different seed for each trial.
        noise=seeded_noise(trialwise_seed+i, 1.0, 0.0, length(rates)) #1/f noise with specific seed, different over the 60 trials, but the same over external parameter sets.
            
        noise_power=var(noise)
        rate_power=var(rates)
            
        desired_signal_to_noise_ratio=0.1
        noise_scaling_factor=sqrt(rate_power/(noise_power*desired_signal_to_noise_ratio))
            
        scaled_noise=noise.*noise_scaling_factor
        noisy_rates[(j-1)*20+i]=rates+scaled_noise
    end
end

## now filter the noisy trajectories into frequency bands.
# freq_range=range(0.5,40,length=100) 
freq_range=3:1:5 #actually just the freqs around 4Hz for this.
trajectories=noisy_rates

filtered_trajectories=Array{Float64,3}(undef,length(trajectories),length(freq_range),length(results[1][1][3,:]))
filters=Vector{Any}(undef,length(freq_range))
for (i,freq) in enumerate(freq_range)
   if freq<1
    filters[i]=digitalfilter(Bandpass(0.1,freq+1),Butterworth(2);fs=10000)
   else
    filters[i]=digitalfilter(Bandpass(freq-1,freq+1),Butterworth(2);fs=10000)
   end
end

@progress for i in 1:length(freq_range)
    for j in 1:length(trajectories)
         filtered_trajectories[j,i,:]=filtfilt(filters[i],trajectories[j])
    end
end

hilbert_transforms=Array{ComplexF64,3}(undef,size(filtered_trajectories))
@progress for i in 1:size(filtered_trajectories,2)
    for j in 1:size(filtered_trajectories,1)
        hilbert_transforms[j,i,:]=hilbert(filtered_trajectories[j,i,:])
    end
end

phases=angle.(hilbert_transforms)

##plot IEPC over time for 4Hz band only (it is band pass filtered to +-1Hz)
using Statistics
function moving_mean_circular(data::AbstractVector, window_size::Integer)
    n = length(data)
    smoothed_data = similar(data)
    
    # A centered window requires an odd size to have an equal number of points on both sides
    if iseven(window_size)
        throw(ArgumentError("Window size must be odd for a perfectly centered window."))
    end
    
    half_win = window_size ÷ 2
    
    for i in 1:n
        val_sum = zero(eltype(data))
        
        # Iterate over the relative centered window: [-half_win, ..., 0, ..., half_win]
        for j in -half_win:half_win
            # mod1(x, n) natively handles circular wrapping for 1-based indexing
            # e.g., if n=100 and i+j = 0, idx becomes 100
            # if i+j = 101, idx becomes 1
            idx = mod1(i + j, n)
            val_sum += data[idx]
        end
        
        smoothed_data[i] = val_sum / window_size
    end
    
    @info size(smoothed_data)
    return smoothed_data 
end
Plots.default(titlefontsize=16,legendfontsize=6,tickfontsize=9,guidefontsize=11)

freq_range[2]
four_hz_IEPC=abs.(mean(mean(exp.(im.*phases[:,2:2,1:20*10000]),dims=2)[:,1,:]',dims=2))

#PLOT FIGURE 9_1:

plot(1/10000:1/10000:20.0,four_hz_IEPC,label="4Hz IEPC",lw=1,color=:blue)
#add a line at 7.3 seconds when the stimulus was switched off
plot!([10.0,10.0],[0.0,1.0],color=:red,linewidth=1,linestyle =:dash, label="stim_off")
#and at 5seconds when the stimulus starts:
plot!([5.0,5.0],[0.0,1.0],color=:black,linewidth=1,linestyle =:dot,label="stim_on")
smoothed_four_hz_IEPC=moving_mean_circular(four_hz_IEPC[:,1],12001)
plot!(1/10000:1/10000:20.0,smoothed_four_hz_IEPC,label="4Hz IEPC smoothed",linewidth=2,color=:purple)
plot!(dpi=300,xlabel=L"\textrm{time~(s)}", ylabel=L"\textrm{IEPC~~}4\textrm{Hz}")

#size and save it:
one_column_size=figure_size_tuple(1, aspect_ratio=1.0)
plot!(size=one_column_size,dpi=300,legend=:topright)

savefig("IEPC_4Hz_band.pdf")
savefig("IEPC_4Hz_band.svg")
