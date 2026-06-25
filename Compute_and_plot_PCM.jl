using Pkg; Pkg.activate(".")
using NGNMM_NSP_paper_code
using Arrow, DataFrames, JLD2
using Interpolations, ComponentArrays, Plots, LaTeXStrings

#font sizes:
Plots.default(titlefontsize=16,legendfontsize=12,tickfontsize=9,guidefontsize=11)


frequencies=collect(range(2,15.0,length=60)) #60 frequencies between 2 and 15 Hz
sine_path1="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_1.jld2"
sine_path2="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_2.jld2"
sine_path3="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_3.jld2"

sine_drive_data1=jldopen(sine_path1,"r")
sine_drive_interpolators1=sine_drive_data1["drives"]
sine_drive_data2=jldopen(sine_path2,"r")
sine_drive_interpolators2=sine_drive_data2["drives"]
sine_drive_data3=jldopen(sine_path3,"r")
sine_drive_interpolators3=sine_drive_data3["drives"]
sine_drive_interpolators=vcat(sine_drive_interpolators1,sine_drive_interpolators2,sine_drive_interpolators3)

### using the evoked model:


save_path="./"
save_traj=false

# set parameters
u0=ComponentArray(x1=0.0, x2=0.0)
time_range=(0.0,10.0)
phoneme_sampling_rate=44100
drive_amplitude=20.0
α=1/0.03
p=ComponentArray(α=α, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,noise_case_reference=1)
PCM_range=(5.5,10.0)
freq_range=(2,15.0) #Hz
stimulus_type="peakenvelope"

Evoked_PCMs=get_evoked_model_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(Evoked_PCMs,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="Impulse stim - evoked model",marker_z=frequencies)

#phase-resetting oscillator with frequency tuning:
u0=ComponentArray(θ=0.0, r=1.0)
c=1.0 #will be updated to 0.7*pi*(1/(maximum(interpolators_global[1])*drive_amplitude)) given a particular test condition,.
time_range=(0.0,10.0)
F=4.0
phoneme_sampling_rate=44100
drive_amplitude=20.0
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
# p=ComponentArray(F=F, c=c, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,modulation=phase_modulation,q=q_normalisation,noise_case_reference=1,frequencies=[F for i in 1:60])
PCM_range=(5.5,10.0)

freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

tuned_phase_resetting_oscillator_PCMs=get_coupled_oscillator_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(tuned_phase_resetting_oscillator_PCMs,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="sine+1 stim - tuned coupled oscillator",marker_z=frequencies)

#to inspect trajectory and its response to stimulus:
# traj_path_coupled_oscillator="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/5000sr_oscillator_response_to_sine_stimulus_freq_range_0.5to30.0Hz_peakenvelope_transformed.arrow.lz4"
# trajectories_coupled_oscillator=Matrix(DataFrame(Arrow.Table(traj_path_coupled_oscillator)))
# sine_impulses=get_unitary_peakenv_impulses(sine_drive_interpolators,phoneme_sampling_rate) 
# traj_times=range(5.5,10.0,length=size(trajectories_coupled_oscillator,1))
# stimulus_indexes=traj_times.*phoneme_sampling_rate  
# stimulus_times=stimulus_indexes./phoneme_sampling_rate
# collect(stimulus_times)
# ps_coupled=[]
# for freq_idx in 1:5:50
#     p=plot(traj_times,trajectories_coupled_oscillator[:,freq_idx],xlabel="Time (s)")
#     plot!(p,stimulus_times,sine_impulses[freq_idx](stimulus_indexes),label="Stimulus",xlabel="Time (s)",xlims=(5.0,10.0))
#     #plot sine wave on top:
#     plot!(p,stimulus_times, sine_drive_interpolators[freq_idx](stimulus_indexes))
#     push!(ps_coupled,p)
# end
# plot(ps_coupled...,layout=(5,2),size=(1200,800))






#slow parameter set NMM
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
DAR=12.2 #scaled phoneme DAR -> phoneme magnitude ~0.25, sine is 1. so *0.25.
#Drive 
phoneme_sampling_rate=44100
drive_amplitude=Π
# p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=drive_amplitude,noise_selector=1,noise_case_reference=1) 
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=DAR*η_0,noise_selector=1,noise_case_reference=1)  #with dar amplitude setting
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,10.0)
PCM_range=(5.5,10.0)

freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

slow_NMM_PCMs=get_NGNMM_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(slow_NMM_PCMs,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="DAR set full 1+ sine wave stim - 4Hz NGNMM",marker_z=frequencies)

#fast parameter set NMM
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
DAR=12.2*0.25
#Drive 
phoneme_sampling_rate=44100
drive_amplitude=Π
# p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=drive_amplitude,noise_selector=1,noise_case_reference=1) 
p=ComponentArray(α=α,k=k,C=C,vsyn=vsyn,Δ=Δ, η_0=η_0,α_D=α_D,sampling_rate=phoneme_sampling_rate, drive_amplitude=DAR*η_0,noise_selector=1,noise_case_reference=1)  #with dar amplitude setting
u0=ComponentArray(g_dot=0.0,g=0.82,Z=0.0+im*0.0, A_dot=0.0, A=0.0)
time_range=(0.0,20.0)
PCM_range=(5.5,10.0)

freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

fast_NMM_PCMs=get_NGNMM_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(fast_NMM_PCMs,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="DAR setabs sine wave stim - 15Hz NGNMM",marker_z=frequencies)

#to plot a grey dashed circle of radius 1 to indicate boundary of unit circle
circle_x = cos.(range(0, 2π, length=100))
circle_y = sin.(range(0, 2π, length=100))
#to compute mean resultant of the PCM vectors:
function compute_mean_resultant(PCMs)
    #turn into polar coordinates
    angles=angle.(PCMs)
    magnitudes=abs.(PCMs)
    polarform=magnitudes.*exp.(im.*angles)
    mean_resultant=sum(polarform)/length(PCMs)
    return mean_resultant
end


#plot all four cases on one figure: PPP (evoked impulse stim), coupled oscillator sine stim, NGNMM slow sine stim, NGNMM fast sine stim
one_column_size=figure_size_tuple(1, aspect_ratio=1.0) #square aspect ratio for 4 panel plot, will be adjusted to fit 4 panels in one column width.
p1=scatter(Evoked_PCMs,xlims=(-1.1,1.1),ylims=(-1.1,1.1),marker_z=frequencies,label=nothing, cbar=false,markerstrokewidth=0, markersize=6,xlabel=L"$\textrm{Real}$", ylabel=L"$\textrm{Imaginary}$");
plot!(p1, circle_x, circle_y, linestyle=:dash, linecolor=:gray, label=nothing,size=one_column_size);
quiver!(p1,[0.0],[0.0],quiver=([real(compute_mean_resultant(Evoked_PCMs))],[imag(compute_mean_resultant(Evoked_PCMs))]),arrow=:arrow,linewidth=2,linecolor=:black,label="Mean Resultant");
plot!(dpi=300)

#testing saving one column figure.
# savefig(p1,"Evoked_model_PCMs_one_column_300dpi_72pixperinch.pdf")


p2=scatter(tuned_phase_resetting_oscillator_PCMs,xlims=(-1.1,1.1),ylims=(-1.1,1.1),marker_z=frequencies,label=nothing, cbar=false,markerstrokewidth=0, markersize=6,xlabel=L"\textrm{Real}", ylabel=L"\textrm{Imaginary}");
plot!(p2, circle_x, circle_y, linestyle=:dash, linecolor=:gray, label=nothing,size=one_column_size);
quiver!(p2,[0.0],[0.0],quiver=([real(compute_mean_resultant(tuned_phase_resetting_oscillator_PCMs))],[imag(compute_mean_resultant(tuned_phase_resetting_oscillator_PCMs))]),arrow=:arrow,linewidth=2,linecolor=:black,label="Mean Resultant");
p3=scatter(slow_NMM_PCMs,xlims=(-1.1,1.1),ylims=(-1.1,1.1),marker_z=frequencies,label=nothing, cbar=false,markerstrokewidth=0, markersize=6,xlabel=L"\textrm{Real}", ylabel=L"\textrm{Imaginary}");
plot!(p3, circle_x, circle_y, linestyle=:dash, linecolor=:gray, label=nothing,size=one_column_size);
quiver!(p3,[0.0],[0.0],quiver=([real(compute_mean_resultant(slow_NMM_PCMs))],[imag(compute_mean_resultant(slow_NMM_PCMs))]),arrow=:arrow,linewidth=2,linecolor=:black,label="Mean Resultant");
p4=scatter(fast_NMM_PCMs,xlims=(-1.1,1.1),ylims=(-1.1,1.1),marker_z=frequencies,label=nothing,cbar=false,markerstrokewidth=0, markersize=6,xlabel=L"\textrm{Real}", ylabel=L"\textrm{Imaginary}");
plot!(p4, circle_x, circle_y, linestyle=:dash, linecolor=:gray, label=nothing,size=one_column_size);
quiver!(p4,[0.0],[0.0],quiver=([real(compute_mean_resultant(fast_NMM_PCMs))],[imag(compute_mean_resultant(fast_NMM_PCMs))]),arrow=:arrow,linewidth=2,linecolor=:black,label="Mean Resultant");

clims = extrema(frequencies[1:60])
h2 = scatter([0,0], [0,1], zcolor=[0,1], clims=clims,
                 xlims=(1,1.1), xshowaxis=false, yshowaxis=false, label="", c=:thermal, colorbar_title="Frequency (Hz)", grid=false);
l=@layout [grid(2,2) a{0.01w}]
plot(p1,p2,p3,p4,h2,layout=l,right_margin=3Plots.mm,dpi=300,markersize=4,size=figure_size_tuple(2, aspect_ratio=1.3))
savefig("PCM_comparison_4panel_twocolwidth.pdf")
#format nicely:



















### phase-resetting oscillator case without frequency adjustment:
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
F=4.0 #Hz
p=ComponentArray(F=F, c=c, drive_amplitude=drive_amplitude, noise_selector=1, sampling_rate=phoneme_sampling_rate,modulation=phase_modulation,q=q_normalisation,noise_case_reference=1,frequencies=[F for i in 1:60])
PCM_range=(5.5,10.0)    
save_path="./"
save_traj=true
freq_range=(2.0,15.0) #Hz
stimulus_type="envelope"

PCMs_coupled_oscillator_nofreqadj=get_coupled_oscillator_PCM_across_60freqs_randinitcond(p,u0,time_range,PCM_range,freq_range,save_traj,save_path,stimulus_type)
scatter(PCMs_coupled_oscillator_nofreqadj,xlims=(-1,1),ylims=(-1,1),aspect_ratio=:equal,title="sine wave stim - coupled oscillator no freq adj",marker_z=frequencies)

