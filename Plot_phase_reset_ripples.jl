using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using Statistics
using ComponentArrays, OrdinaryDiffEq, Plots, Parameters, JLD2
using LaTeXStrings


condition="b" #which phoneme condition to use for the drive
idx=1 #which test index to use
noisestimratio=0.0 #which noise to stimulus ratio to use
noise_filename="NoiseSequences60.csv" #which noise file to use

global NGNMM_NSP_paper_code.interpolators_global = jldopen("./Phoneme_Drives/drive_interpolators_$(condition)_test_$(idx)_NSR_$(noisestimratio)_$(split(noise_filename,".")[1]).jld2","r")["drives"]

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

#find the start of each phoneme in the 5 second stimulus window (5.0, 10.0)
#the stimulus is at ~4Hz with some noise, so maybe try for now to just take 5, 5.25 etc.. 
#then plot the mean response in a 200ms window after these times, across all 20 runs and all ~20 phonemes in each run.
plot(abs.(results[2][3,:]),dpi=300,xlims=(25000,30000))
stim_times=collect(range(5.0,9.0,step=0.25))
response_window=0.25
response_sr=Int64(1/saveat) 
response_window_samples=Int64(response_window*response_sr) 

responses=Vector{Vector{Float64}}(undef,340)
for i in 1:20
    for j in 1:17
        responses[(i-1)*17+j]=abs.(results[i][3,Int64((stim_times[j]-response_window)*response_sr):Int64(stim_times[j]*response_sr+response_window_samples)])
    end
end

response_times=range(-0.25,0.25,step=1/response_sr)
plot(response_times,responses[5],label="single trial")
plot!(response_times,mean(responses),label="mean response")
stimulus_time_points=range(-0.25,0.25,step=1/44100)
plot!(stimulus_time_points,interpolators_global[1]((stim_times[2]-0.25)*44100:stim_times[3]*44100),label="example drive")
plot!(legend=:bottomright)
plot!(ylabel=L"\textrm{firing~rate}", xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}")

#size and save it:
one_column_size=figure_size_tuple(1, aspect_ratio=1.7)
plot!(size=one_column_size,legend=:outertopright,dpi=300,xticks=([-0.2,0,0.2],[-0.2,0,0.2]))
savefig("phase_reset_ripples.pdf")


interpolators_global[1]((stim_times[1]-0.25)*44100:stim_times[2]*44100)
for i in 1:20
    for j in 1:17
        println(((i-1)*17)+j)
    end
end


#check state at 5 seconds of the response to make sure it has been randomised...
init_conds=Vector{ComplexF64}(undef,20)
for i in 1:20
    init_conds[i]=results[i][3,Int64(5.0*response_sr)]
    # init_conds[i]=results[i][3,1]#,Int64(5.0*response_sr)]
end
scatter!(init_conds,xlims=(-1,1),ylims=(-1,1))

plot(abs.(results[1][3,:]))
p=plot()
for i in 1:20
    plot!(p,response_times,abs.(results[i][3,Int64((stim_times[1]-response_window)*response_sr):Int64(stim_times[1]*response_sr+response_window_samples)]))
end
plot(p)
plot!(dpi=300,legend=nothing, xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}",ylabel=L"\textrm{firing~rate}")
plot!(size=one_column_size)
savefig("peri_onset_time_responses.pdf")


using ColorSchemes,Colors
Bang_wong_color_palette=[(230,159,0),(0,114,178),(0,158,115),(204,121,167),(86,180,233),(240,228,66),(0,0,0)]
Bang_wong_color_palette_normalised=[(col[1]/255.0,col[2]/255.0,col[3]/255.0) for col in Bang_wong_color_palette]
color_scheme=ColorScheme([Colors.RGB(col...) for col in Bang_wong_color_palette_normalised],"bang wong colorblind friendly")
#font sizes:
Plots.default(titlefontsize=16,legendfontsize=12,tickfontsize=9,guidefontsize=11)

### plot mean response and transparent individual responses over stimulus onset window:
response_plot=plot();
response_times=range(-0.25,0.25,step=1/response_sr)
for i in 1:20
    plot!(response_plot,response_times,abs.(results[i][3,Int64((stim_times[1]-response_window)*response_sr):Int64(stim_times[1]*response_sr+response_window_samples)]),color=color_scheme[2],alpha=0.25,label="",linewidth=1,ylabel=L"\textrm{firing~rate}");
end
#add mean line:
plot!(response_plot,response_times,mean([abs.(results[i][3,Int64((stim_times[1]-response_window)*response_sr):Int64(stim_times[1]*response_sr+response_window_samples)]) for i in 1:20]),color=:blue,alpha=1.0,label=nothing,linewidth=2.0,xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}",ylabel=L"\textrm{firing~rate}");
plot(response_plot)


## and another plot showing the mean phoneme drive here too.
drive_plot=plot();
drive_time_points=range(-0.25,0.25,step=1/44100)
for i in 1:20
    plot!(drive_plot,drive_time_points,NGNMM_NSP_paper_code.interpolators_global[i]((stim_times[1]-0.25)*44100:(stim_times[1]+0.25)*44100),color=color_scheme[1],alpha=0.25,label=nothing,linewidth=1,ylabel=L"\textrm{amplitude}");
end
plot!(drive_plot,drive_time_points,mean([NGNMM_NSP_paper_code.interpolators_global[i]((stim_times[1]-0.25)*44100:(stim_times[1]+0.25)*44100) for i in 1:20]),color=:orange2,alpha=1.0,label=nothing,linewidth=2.0,xlabel=L"\textrm{peri\operatorname{-}onset~time~(s)}",ylabel=L"\textrm{amplitude}");
plot!(drive_plot,yticks=([0.0,0.1,0.2],[0.0,0.1,0.2]))
plot(drive_plot)

# both together
one_column_size_tall=figure_size_tuple(1,aspect_ratio=1.2)

plot((response_plot,drive_plot)...,layout=grid(2,1, heights=[0.75,0.25]), size=one_column_size_tall,dpi=300)

# savefig("slow_NMM_peri_onset_time_responses_transparent_with_mean.pdf")
savefig("fast_NMM_peri_onset_time_responses_transparent_with_mean.pdf")