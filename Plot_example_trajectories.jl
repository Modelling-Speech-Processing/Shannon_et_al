using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code, JLD2
using ComponentArrays, Plots, JSON, WAV, CSV, Statistics, DelimitedFiles, Random, OrdinaryDiffEq, DSP, FFTW
plotlyjs()

raw_audio=wavread("./Stimuli/1.wav")
plot(raw_audio[1][:,1])
raw_audio[2]
global NGNMM_NSP_paper_code.interpolators_global = jldopen("./Phoneme_Drives/drive_interpolators_$("b")_test_$(1)_NSR_$(0.7)_NoiseSequences60.jld2","r")["drives"]
interpolators_global[1]
plot(times,NGNMM_NSP_paper_code.interpolators_global[1][5.0*44100:10*44100])
times=5.0:1/44100:10.0

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
time_range=(0.0,20.0)
ITPCrange=(5.5,10.0)
output=Ensemble_NoisyPhoneme(vary_noise,time_range,p,u0,20,0.0001)

rates=get_firing_rate_NMM(output[1],C,vsyn)
ITPC,fcoeffs,puvs,freqs,rates_fromITPC=NGNMM_NSP_paper_code.calculate_ITPC_1overf_noise(output,Int32(1/0.0001),ITPCrange,C,vsyn,ITPCrange)


individual_size=(1.0*width_px,0.5width_px).*2.5
all_size=(width_px*2.3,width_px*1.1) 
all_bottom_margin=3.0Plots.mm
all_right_margin=1.5Plots.mm
all_legendfontsize=8
all_xtickfontsize=8
all_ytickfontsize=8
all_xlabelfontsize=9
all_ylabelfontsize=9

width_pts = 246.0  # or any other value
inches_per_points = 1.0/72.27
width_inches = width_pts *inches_per_points
width_px= width_inches*100/2  # or  width_inches*DPI.. /2 as two square plots per column approx.

Bang_wong_color_palette=[(230,159,0),(0,114,178),(0,158,115),(204,121,167),(86,180,233),(240,228,66),(0,0,0)]
Bang_wong_color_palette_normalised=[(col[1]/255.0,col[2]/255.0,col[3]/255.0) for col in Bang_wong_color_palette]
using ColorSchemes,Colors
color_scheme=ColorScheme([Colors.RGB(col...) for col in Bang_wong_color_palette_normalised],"bang wong colorblind friendly")


#audiofigure
audioplot=plot(times[1:220313],raw_audio[1],color=color_scheme[1],legend=:topright,xlabel="time (s)",ylabel="amplitude",size=individual_size)
ps=[audioplot for i in 1:4]
plot(ps...,layout=(2,2),size=all_size)
times=5.5:1/44100:10.0
#envelope figure
envfigure=plot(times,NGNMM_NSP_paper_code.interpolators_global[1][5.5*44100:10*44100],color=color_scheme[1],legend=:topright,label="'b' envelope",xlabel="time (s)",ylabel="amplitude",size=individual_size)

##NGNMM rate figure:
normalised_NGNMM_rate_plot=plot(output[1].t[Int32(5.5*(1/0.0001)):Int32(10*(1/0.0001))],rates_fromITPC[1]./maximum(rates_fromITPC[1]),color=color_scheme[2],legend=:topright,label="NGNMM_slow",xlabel="time (s)",ylabel="normalised firing rate",size=individual_size)

### Evoked model rate figure:
times=5.5:1/44100:10.0

α=1/0.030
Π=1.0
#Drive 
phoneme_sampling_rate=44100
drive_amplitude=Π
p=ComponentArray(α=α,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1) 
u0=ComponentArray(x1=0.0,x2=0.0)
time_range=(0.0,20.0)
ITPCrange=(5.5,10.0)
output=Ensemble_EvokedModel(vary_noise,time_range,p,u0,20,0.0001)
ITPC,fcoeffs,puvs,freqs,rates_fromITPC=NGNMM_NSP_paper_code.calculate_ITPC_EvokedModel_noisyrates(output,Int32(1/0.0001),ITPCrange,ITPCrange)


normalised_evoked_rate_plot=plot(output[1].t[Int32(5.5*(1/0.0001)):Int32(10*(1/0.0001))],rates_fromITPC[1]./maximum(rates_fromITPC[1]),color=color_scheme[4],legend=:topright,label="evoked_env",xlabel="time (s)",ylabel="normalised firing rate",size=individual_size)


## phase reset model figure:
times=5.0:1/44100:10.0
phase_modulation=1.0
F=4.0 #hz
Π=1.0
if phase_modulation<=0.5
    q_normalisation=2/((1-phase_modulation)*π) #scaling factor to make absolute area under stimulus modulation curve = 4. 
else
    m=phase_modulation
    q_normalisation=4/(4*m*sqrt(1-((m-1)/m)^2)+2*(m-1)asin((m-1)/m))
end
#Drive 
drive_amplitude=Π
phoneme_sampling_rate=44100

p=ComponentArray(F=F,c=1.0,sampling_rate=phoneme_sampling_rate, drive_amplitude=Π,noise_selector=1,noise_case_reference=1.0,modulation=phase_modulation,q=q_normalisation) 


max_stimulus_amplitude=maximum(maximum.(NGNMM_NSP_paper_code.interpolators_global))
# p.drive_amplitude=1.0/max_stimulus_amplitude

#DAR so that stimulus integral is = 20000 (unmodified first vowel stimulus was 20315, so approximating that.)
# println("mean sum of phoneme envelope", mean(sum.([i[5.0*44100:10.0*44100] for i in interpolators_global])))
p.drive_amplitude=20000/mean(sum.([i[5.0*44100:10.0*44100] for i in NGNMM_NSP_paper_code.interpolators_global]))

# p.c=0.7*pi*(1/(max_stimulus_amplitude*p.drive_amplitude*p.q)) #including q. this is the final maximumum stimulus amplitude after all normalisation. 
#s(t) normalised to 20000, then s(t)*phasemod normalised to 4.0*20000, now set so max phase correction after all this is 0.7π.
#for high c simulation, a max reset that is the full half circle
p.c=1.0*pi*(1/(max_stimulus_amplitude*p.drive_amplitude*p.q)) #including q. this is the final maximumum stimulus amplitude after all normalisation. 
            

u0=ComponentArray(θ=0.0, r=1.0)
time_range=(0.0,20.0)
ITPCrange=(5.5,10.0)
output=Ensemble_CoupledOscillators_modulated(vary_noise,time_range,p,u0,20,0.0001)
ITPC,fcoeffs,puvs,freqs,rates_fromITPC=NGNMM_NSP_paper_code.calculate_ITPC_CoupledOscillators_noisyrates(output,Int32(1/0.0001),ITPCrange,ITPCrange)


normalised_phase_reset_rate_plot=plot(output[1].t[Int32(5.5*(1/0.0001)):Int32(10*(1/0.0001))],rates_fromITPC[1]./maximum(rates_fromITPC[1]),color=color_scheme[3],legend=:topright,label="Phase-reset",xlabel="time (s)",ylabel="normalised firing rate",size=individual_size)



#plot envelope figure, and the three model rate figures together,.
all_plots=[envfigure,normalised_NGNMM_rate_plot,normalised_phase_reset_rate_plot,normalised_evoked_rate_plot]
tosave=plot(all_plots...,layout=(2,2),size=all_size.*2.5,bottom_margin=all_bottom_margin,right_margin=all_right_margin,xtickfontsize=all_xtickfontsize,ytickfontsize=all_ytickfontsize,xlabelfontsize=all_xlabelfontsize,ylabelfontsize=all_ylabelfontsize)
#save it
savefig(tosave,"./phonemepaper_modelactivitycomparison.pdf")