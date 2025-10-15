using Pkg;Pkg.activate(".")
using JSON, DelimitedFiles, Plots, Statistics, MAT
using NGNMM_NSP_paper_code: get_frequencies, get_sorting_indices
using LaTeXStrings
gr()

##
frequencies=get_frequencies(4.5*10000,10000)[1:676]
_4Hz_index=indexin([4.0],frequencies)[1]
driveamplituderatios=[0.1,0.496,2.46,12.2,60.5,100.0]
driveamplituderatiosevoked=[0.1,0.496,1.0,2.46,12.2,60.5,100.0]
noisestimratios=[0.0,0.05,0.1,0.3,0.5,0.7,0.9,1.0]
phonemes = ["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
NGNMM_path="./Results/NMM/"

C0066D025_data=Vector{Any}()
for DAR in driveamplituderatios
    for NSR in noisestimratios
        push!(C0066D025_data,[JSON.parsefile(NGNMM_path*"1overf_Randinitconds_c066D025LowFreqTestallnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end



Aine_data=Vector{Any}()
for DAR in driveamplituderatios
    for NSR in noisestimratios
        push!(Aine_data,[JSON.parsefile(NGNMM_path*"1overf_Randinitconds_AineRerunTestallnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end


#control case data
oscillator_path="./Results/Phase_Reset_Model/"
phase_reset_data=Vector{Any}()
for DAR in driveamplituderatios[1]
    for NSR in noisestimratios
        push!(phase_reset_data,[JSON.parsefile(oscillator_path*"SNR05_high_c_60randinitcond_seeded_noise_normalised_nonrectified_phasemod_$(1.0)_integralnormalised_coupled_oscillator_control_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end

nophasereset_data=Vector{Any}()
for DAR in driveamplituderatios[1]
    for NSR in noisestimratios
        push!(nophasereset_data,[JSON.parsefile(oscillator_path*"SNR05_high_c_60randinitcond_seeded_noise_normalised_nophasereset_phasemod_$(1.0)_integralnormalised_coupled_oscillator_control_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end

#check 4Hz ITPC for b condition of the baseline C0066D025 model for each of the three tests:
baseline_C0066D025_4Hz_ITPCs=Vector{Vector{Float64}}()
for (i,DAR) in enumerate(driveamplituderatios)
    push!(baseline_C0066D025_4Hz_ITPCs,[])
    for j in 1:3
        push!(baseline_C0066D025_4Hz_ITPCs[i],C0066D025_data[8*(i-1)+7][1]["noisestimratio: 0.9"]["ITPCs"]["f"][j][19])
    end
end

#modulated phase-reset data
modulations=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
modulated_nonrectified_data=Vector{Any}()
for DAR in driveamplituderatios[1]
    for modulation in modulations
        push!(modulated_nonrectified_data,[JSON.parsefile(oscillator_path*"SNR05_high_c_60randinitcond_seeded_noise_normalised_nonrectified_modulated_phasemod_$(modulation)_integralnormalised_coupled_oscillator_control_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(0.1)data.json")])
    end
end

Frequency_varied_phasereset_data=Vector{Any}()
Frequencies=[2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0]
for DAR in driveamplituderatios[1]
    for Freq in Frequencies
        push!(Frequency_varied_phasereset_data,[JSON.parsefile(oscillator_path*"FreqVary_freq$(Freq)Hz_SNR05_high_c_60randinitcond_seeded_noise_normalised_nonrectified_modulated_phasemod_1.0_integralnormalised_coupled_oscillator_control_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(0.1)data.json")])
    end
end

#evoked model data
evoked_model_path="./Results/Evoked_Model/"

evoked_derivativestim_data=Vector{Any}()
for DAR in driveamplituderatiosevoked
    for NSR in noisestimratios
        push!(evoked_derivativestim_data,[JSON.parsefile(evoked_model_path*"SNR05_high_c_60randinitcond_seeded_noise_evokedODE_stimulustype_derivative_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end

evoked_envelopestim_data=Vector{Any}()
for DAR in driveamplituderatiosevoked
    for NSR in noisestimratios
        push!(evoked_envelopestim_data,[JSON.parsefile(evoked_model_path*"SNR05_high_c_60randinitcond_seeded_noise_evokedODE_stimulustype_envelope_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end

evoked_peakenvstim_data=Vector{Any}()
for DAR in driveamplituderatiosevoked
    for NSR in noisestimratios
        push!(evoked_peakenvstim_data,[JSON.parsefile(evoked_model_path*"SNR05_high_c_60randinitcond_seeded_noise_evokedODE_stimulustype_peakenvelope_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end

evoked_peakratestim_data=Vector{Any}()
for DAR in driveamplituderatiosevoked
    for NSR in noisestimratios
        push!(evoked_peakratestim_data,[JSON.parsefile(evoked_model_path*"SNR05_high_c_60randinitcond_seeded_noise_evokedODE_stimulustype_peakrate_test_allnoisesDAR$(round(DAR,sigdigits=3))NSR$(NSR)data.json")])
    end
end



#Experimental Data:
oana_data=matread("./ITPC_peakder_condition.mat")
oana_ITPCs=oana_data["meanitc4_cond"]
sorted_oana_ITPCs=oana_ITPCs[get_sorting_indices("./")]

time2pds=Frequency_varied_phasereset_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]


########################
########################
########################
# exploring scatters vs t2pd for each drive setting.
########################
########################
########################
time2pds=Frequency_varied_phasereset_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]
fourHzITPCS=Vector{Array{Float64,3}}()
sorting_indxs=get_sorting_indices("./")

for (i,DAR) in enumerate(driveamplituderatios)
    this_dar_array=Array{Float64,3}(undef,20,8,15) #5 models, 8 NSRs, 15 conditions
    for (j,NSR) in enumerate(noisestimratios)
        #phase reset models and evoked model have constant DAR.
        if NSR==0.05
            #swapped so it uses squared ITPCs before computing means.
            this_dar_array[1,j,:]=reduce(hcat,[mean(x->x.^2,C0066D025_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[2,j,:]=reduce(hcat,[mean(x->x.^2,Aine_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[3,j,:]=reduce(hcat,[mean(x->x.^2,phase_reset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[6,j,:]=reduce(hcat,[mean(x->x.^2,nophasereset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

            if DAR<1.0
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

            else
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

            end
            this_dar_array[11,j,:]=sorted_oana_ITPCs[:]

        else
            this_dar_array[1,j,:]=reduce(hcat,[mean(x->x.^2,C0066D025_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[2,j,:]=reduce(hcat,[mean(x->x.^2,Aine_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[3,j,:]=reduce(hcat,[mean(x->x.^2,phase_reset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[6,j,:]=reduce(hcat,[mean(x->x.^2,nophasereset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

         if DAR<1.0
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

            else
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

            end
            this_dar_array[11,j,:]=sorted_oana_ITPCs[:]
        end
    end
    push!(fourHzITPCS,this_dar_array)
end
scatter_ITPC_vs_t2pd_plots=Array{Any,2}(undef,length(driveamplituderatios),length(noisestimratios))
labels=["NGNMM slow" "NGNMM fast" "Phase Reset" "no Phase Reset" "evoked_derivative" "evoked_envelope" "evoked_peakenv" "evoked_peakrate" "EEG"]
#marker styles for the scatter plots
markerstyles=[:circle :dtriangle :star :diamond :hexagon :star :utriangle :circle :square :utriangle :dtriangle]
colors=[:blue :red :green :orange :purple :brown :cyan :gray :black]
for i in 1:length(driveamplituderatios)
    for j in 1:length(noisestimratios)
        if j==8
            # scatter_ITPC_vs_t2pd_plots[i,j]=scatter(time2pds,fourHzITPCS[i][[1,2,3,4,6,7,8,9,10,11,12,13,14,17,18,19,20],j,:]',label=labels,markerstrokealpha=0.0,markeralpha=1.0,markershape=markerstyles,title="NSR $(noisestimratios[j])",xlabel="Time to peak derivative (s)",ylabel="4Hz ITPC",legend=:outertopright)
            scatter_ITPC_vs_t2pd_plots[i,j]=scatter(time2pds,fourHzITPCS[i][[1,2,3,6,7,8,9,10,11],j,sorting_indxs]',color=colors,label=labels,markerstrokewidth=0.1,markeralpha=1.0,markershape=markerstyles,title=L"\lambda ="*" $(noisestimratios[j])",xlabel="Time to peak derivative (s)",ylabel="4Hz ITPC",legend=:outertopright)
        else
            # scatter_ITPC_vs_t2pd_plots[i,j]=scatter(time2pds,fourHzITPCS[i][[1,2,3,4,6,7,8,9,10,11,12,13,14,17,18,19,20],j,:]',label=nothing,markerstrokealpha=0.0,markeralpha=1.0,markershape=markerstyles,title="NSR $(noisestimratios[j])",xlabel="Time to peak derivative (s)",ylabel="4Hz ITPC",legend=:outertopright)
            scatter_ITPC_vs_t2pd_plots[i,j]=scatter(time2pds,fourHzITPCS[i][[1,2,3,6,7,8,9,10,11],j,sorting_indxs]',color=colors,label=nothing,markerstrokewidth=0.1,markeralpha=1.0,markershape=markerstyles,title=L"\lambda ="*" $(noisestimratios[j])",xlabel="Time to peak derivative (s)",ylabel="4Hz ITPC",legend=:outertopright)
        end
    end
end

#plot the 4HzITPC vs t2pd for each DAR for each NSR
for d in 1:length(driveamplituderatios)
    display(scatter(scatter_ITPC_vs_t2pd_plots[d,:]...,layout=(2,4),ylims=(0.0,1.0),size=(1000,1500),xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",plot_title=L"D="*"$(driveamplituderatios[d])",legend=:topright))
end
#just save DAR=12.2 figure
tosave=plot(scatter_ITPC_vs_t2pd_plots[4,:]...,layout=(3,3),ylims=(0.0,1.0),size=(1200,1200),xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",plot_title="",legend=:topright,markerstrokewidth=0.0,markersize=6,left_margin=5Plots.mm,legendfontsize=12,xlabelfontsize=14,ylabelfontsize=14,ytickfontsize=12,xtickfontsize=12)
savefig(tosave,"./ppaper_ITPC_vs_t2pd_scatter_allNSR_allmodels_DAR12.2.pdf")
########################
########################
########################
# figure that is just the NGNMM results across NSRs each series is a different DAR
########################
########################
########################

time2pds=Frequency_varied_phasereset_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]
fourHzITPCS=Vector{Array{Float64,3}}()
sorting_indxs=get_sorting_indices("./")

for (i,DAR) in enumerate(driveamplituderatios)
    this_dar_array=Array{Float64,3}(undef,11,8,15) #11 models, 8 NSRs, 15 conditions
    for (j,NSR) in enumerate(noisestimratios)
        #phase reset models and evoked model have constant DAR.
        if NSR==0.05
            #swapped so it uses squared ITPCs before computing means.
            this_dar_array[1,j,:]=reduce(hcat,[mean(x->x.^2,C0066D025_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[2,j,:]=reduce(hcat,[mean(x->x.^2,Aine_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[3,j,:]=reduce(hcat,[mean(x->x.^2,phase_reset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[6,j,:]=reduce(hcat,[mean(x->x.^2,nophasereset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            # this_dar_array[2,j,:]=Aine_data[8*(i-1)+j][1]["noisestimratio: $NSR"]["sortedITPCS_t2PD"][1]
            # this_dar_array[3,j,:]=phase_reset_data[j][1]["noisestimratio: $NSR"]["sortedITPCS_t2PD"][1]
            # this_dar_array[6,j,:]=nophasereset_data[j][1]["noisestimratio: $NSR"]["sortedITPCS_t2PD"][1]
            if DAR<1.0
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

            else
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

            end
            this_dar_array[11,j,:]=sorted_oana_ITPCs[:]

        else
            this_dar_array[1,j,:]=reduce(hcat,[mean(x->x.^2,C0066D025_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[2,j,:]=reduce(hcat,[mean(x->x.^2,Aine_data[8*(i-1)+j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[3,j,:]=reduce(hcat,[mean(x->x.^2,phase_reset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
            this_dar_array[6,j,:]=reduce(hcat,[mean(x->x.^2,nophasereset_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])

         if DAR<1.0
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])


            else
                this_dar_array[7,j,:]=reduce(hcat,[mean(x->x.^2,evoked_derivativestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[8,j,:]=reduce(hcat,[mean(x->x.^2,evoked_envelopestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[9,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakenvstim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])
                this_dar_array[10,j,:]=reduce(hcat,[mean(x->x.^2,evoked_peakratestim_data[j][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19] for condition in Condition_keys])


            end
            this_dar_array[11,j,:]=sorted_oana_ITPCs[:]
        end
    end
    push!(fourHzITPCS,this_dar_array)
end

fourHzITPCS_DARfinaldim=stack(fourHzITPCS) #is 20x8x15x6 (20 models, 8 NSRs, 15 conditions, 6 DARs)

#what we are interested in:
NGNMM_EEG_ITPCS=fourHzITPCS_DARfinaldim[[1,2,11],:,:,:] #3x8x15x6 (3 models, 8 NSRs, 15 conditions, 6 DARs)

scatter_ITPC_vs_t2pd_plots=Array{Any,1}(undef,length(noisestimratios)) # 8 noise ratios.
for j in 1:length(noisestimratios)
    scatter_ITPC_vs_t2pd_plots[j]=scatter!()
end

labels=["NGNMM slow" "NGNMM fast" "EEG"]
#marker styles for the scatter plots
markerstyles=[:circle :utriangle :square]# :star :diamond :hexagon :star :square :circle :dtriangle :utriangle :dtriangle]
# for color bars
#marker alpha decreases for decreasing DAR?
driveamplituderatios
log_driveamplituderatios_scaled_to_0_0p1_for_transparency=(log.(driveamplituderatios) .- minimum(log.(driveamplituderatios)))./(maximum(log.(driveamplituderatios))-minimum(log.(driveamplituderatios))).*0.9 .+0.1
markeralphas=log_driveamplituderatios_scaled_to_0_0p1_for_transparency

using ColorSchemes
blues=[(0,0,1.0,markeralphas[1]),(0,0,1.0,markeralphas[2]),(0,0,1.0,markeralphas[3]),(0,0,1.0,markeralphas[4]),(0,0,1.0,markeralphas[5]),(0,0,1.0,markeralphas[6])]
orange=[(1.0,165/255,0,markeralphas[1]),(1.0,165/255,0,markeralphas[2]),(1.0,165/255,0,markeralphas[3]),(1.0,165/255,0,markeralphas[4]),(1.0,165/255,0,markeralphas[5]),(1.0,165/255,0,markeralphas[6])]
using Colors
color_scheme_blue=[Colors.RGBA(col...) for col in blues]
color_scheme_orange=[Colors.RGBA(col...) for col in orange]
color_scheme_black=[Colors.RGBA(0,0,0,1.0) for i in 1:6]


# Make discrete gradients with 6 steps each
blue_grad   = cgrad(:Blues_6, 256)
orange_grad = cgrad(:Oranges_6, 256)

# Black doesn’t need a gradient
black_color = cgrad(color_scheme_black,6, categorical = true)


# Normalise the drive values to [0,1] for mapping into the gradient
function map_drive_to_colour(drive, gradient, min_drive, max_drive)
    t = (log10(drive) - log10(min_drive)) / (log10(max_drive) - log10(min_drive))
    return gradient[t]   # pick colour from gradient
end
min_drive=minimum(driveamplituderatios)
max_drive=maximum(driveamplituderatios)


# colors=[color_scheme_blue color_scheme_orange color_scheme_black]
# for d in 1:length(driveam
plots_for_each_nsr=Array{Any,1}(undef,length(noisestimratios))
for n in 1:length(noisestimratios)
    plots_for_each_nsr[n]=scatter()
    for d in 1:length(driveamplituderatios)
        series_colors=[map_drive_to_colour(driveamplituderatios[d], blue_grad, min_drive, max_drive) map_drive_to_colour(driveamplituderatios[d], orange_grad, min_drive, max_drive) black_color[1]]
        if d==6&&n==8
            # scatter!(plots_for_each_nsr[n],time2pds,NGNMM_EEG_ITPCS[[1,2,3],n,sorting_indxs,d]',label=labels,color=colors,markerstrokewidth=0.5,markersize=6,markeralpha=markeralphas[d],markershape=markerstyles,xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",legend=:bottomright,title=L"\lambda="*" $(noisestimratios[n])")
            scatter!(plots_for_each_nsr[n],time2pds,NGNMM_EEG_ITPCS[[1,2,3],n,sorting_indxs,d]',label=labels,color=series_colors,markerstrokewidth=0.0,markersize=6,markershape=markerstyles,xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",legend=:bottomright,title=L"\lambda="*" $(noisestimratios[n])",cbar=true)
        else
            scatter!(plots_for_each_nsr[n],time2pds,NGNMM_EEG_ITPCS[[1,2,3],n,sorting_indxs,d]',label=nothing,color=series_colors,markerstrokewidth=0.0,markersize=6,markershape=markerstyles,xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",legend=:bottomright,title=L"\lambda="*" $(noisestimratios[n])")
        end
    end
end
display(plot(plots_for_each_nsr...,layout=(3,3),ylims=(0.0,1.0),size=(1000,1500),xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",legend=:topright,markersize=6,markerstrokewidth=0.1,left_margin=5Plots.mm))


# Make a dummy scatter for the colour bar
drivevals = driveamplituderatios
dummy = scatter(drivevals, drivevals;
    zcolor = drivevals,
    c = blue_grad,                # pick one gradient
    colorbar = true,
    label = "",
    markersize = 0)   

display(plot(plots_for_each_nsr...,
    dummy,                        # add dummy subplot with colour bar
    layout = (3,3),               # extra row for colour bar
    size = (1200,1200),
    ylims = (0.0,1.0),
    xlabel = L"L\textrm{~(ms)}",
    ylabel = L"4\textrm{Hz}~~R^2",
    legend = :topright,
    markersize = 6,
    markerstrokewidth = 0.1,
    left_margin = 5Plots.mm
))

tosave=plot(plots_for_each_nsr...,
    dummy,                        # add dummy subplot with colour bar
    layout = (3,3),               # extra row for colour bar
    size = (1200,1200),
    ylims = (0.0,1.0),xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",legend=:topright,markersize=6,markerstrokewidth=0.0,left_margin=5Plots.mm,legendfontsize=12,xlabelfontsize=14,ylabelfontsize=14,ytickfontsize=12,xtickfontsize=12)
savefig(tosave,"./ppaper_NGNMM_EEG_ITPCs_vs_t2pd_across_NSRs_seriesDAR_bluebar.pdf")
savefig(tosave,"./ppaper_NGNMM_EEG_ITPCs_vs_t2pd_across_NSRs_seriesDAR_orangebar.pdf")

########################
#plot the 4HzITPC vs t2pd for each DAR for each NSR
# for d in 1:length(driveamplituderatios)
for d in 1
    display(scatter(scatter_ITPC_vs_t2pd_plots[d,:]...,layout=(2,4),ylims=(0.0,1.0),size=(1200,1200),xlabel=L"L\textrm{~(ms)}",ylabel=L"4\textrm{Hz}~~R^2",plot_title=L"D="*"$(driveamplituderatios[d])",legend=:topright))
end


########################
########################
########################
# getting all correlations, and p values, CI's, and corrected p values for every result. will tabulate in SI.
# using the fourHzITPCS_DARfinaldim data object.
########################
########################
########################
models=["NGNMM slow" "NGNMM fast" "Phase Reset" "no Phase Reset" "evoked_derivative" "evoked_envelope" "evoked_peakenv" "evoked_peakrate"]#"C0066D025_derivativestim" "C0066D025_peakenvelopestim" "C0066D025_peakratestim" "Aine_peakratestim" "pr DerivStim" "pr peakenvstim" "pr peakratestim"]
sorting_indxs=get_sorting_indices("./")
#sorted condition dimension to align with time2pds.
data_for_correlations=fourHzITPCS_DARfinaldim[[1,2,3,6,7,8,9,10,11],:,sorting_indxs,:] # 9 models, 8 NSRs, 15 conditions, 6 DARs
#drop the eeg data
data_for_correlations=data_for_correlations[1:8,:,:,:]

#to correlate to:
time2pds=Frequency_varied_phasereset_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]
#functions to get stats:
using Statistics, Distributions
function get_t_statistic(correlation_coefficient, n)
    # n is the number of samples
    t_statistic = correlation_coefficient * sqrt((n - 2) / (1 - correlation_coefficient^2))
    return t_statistic
end

function get_p_value_two_tailed(t_statistic, n)
    # n is the number of samples
    p_value = 2 * (1 - cdf(TDist(n - 2), abs(t_statistic)))
    return p_value
end


function compute_correlation_confidence_interval(r::Float64, n::Int)
    z =atanh(r)
    se = 1 / sqrt(n - 3)
    z_critical = 1.96 # for 95% confidence
    lower_bound = tanh(z - z_critical * se)
    upper_bound = tanh(z + z_critical * se)
    return (lower_bound, upper_bound)
end

function benjamini_hochberg(p_values::Vector{Float64})
    m = length(p_values)
    
    # Record the original order and sort the p-values
    sorted_indices = sortperm(p_values)
    sorted_p_values = p_values[sorted_indices]
    
    # Calculate the BH adjusted p-values
    # Start from the largest p-value and work backwards
    adjusted_p = similar(sorted_p_values)
    adjusted_p[m] = sorted_p_values[m] # The last one is unadjusted
    
    for i in (m-1):-1:1
        # Calculate the raw adjusted value: (m/i) * p_i
        raw_adjustment = sorted_p_values[i] * (m / i)
        
        # The adjusted p-value must be non-decreasing
        # It's the minimum of its own adjustment and the next highest one
        adjusted_p[i] = min(raw_adjustment, adjusted_p[i+1])
    end
    
    # Ensure values do not exceed 1.0
    clamp!(adjusted_p, 0.0, 1.0)
    
    # Return the adjusted p-values to their original order
    original_order_adjusted_p = similar(p_values)
    original_order_adjusted_p[sorted_indices] = adjusted_p
    
    return original_order_adjusted_p
end

function bonferroni_correction(p_values::Vector{Float64})
    m = length(p_values)
    corrected_p_values = p_values .* m
    corrected_p_values = clamp.(corrected_p_values, 0.0, 1.0) # Ensure values do not exceed 1.0
    return corrected_p_values
end

#check data shape
data_for_correlations[1,1,:,1]
statistics=Array{Float64,4}(undef,length(models),length(noisestimratios),length(driveamplituderatios),6) #6 is for r, CI_lower,CI_upper, p, p_bonferroni, p_FDR
for m in 1:length(models)
    for n in 1:length(noisestimratios)
        for d in 1:length(driveamplituderatios)

            #correlations
            statistics[m,n,d,1]=cor(data_for_correlations[m,n,:,d],time2pds) #r value
            #p value
            statistics[m,n,d,4]=get_p_value_two_tailed(get_t_statistic(statistics[m,n,d,1],length(time2pds)),length(time2pds)) #p value
            #CIs
            CIs=compute_correlation_confidence_interval(statistics[m,n,d,1],length(time2pds))
            #assign lower and upper bounds CIs to array
            statistics[m,n,d,2]=CIs[1] #CI lower
            statistics[m,n,d,3]=CIs[2] #CI upper
            #corrected p values, bonferroni 
            statistics[:,n,d,5]=bonferroni_correction(statistics[:,n,d,4]) #bonferroni
            #and FDR
            statistics[:,n,d,6]=benjamini_hochberg(statistics[:,n,d,4]) #FDR
        end

    end
end

using CSV, DataFrames
df_NGNMM_slow=DataFrame(Model=String[],NSR=Float64[],DAR=Float64[],r=Float64[],CI_lower=Float64[],CI_upper=Float64[],p=Float64[],p_bonferroni=Float64[],p_FDR=Float64[])
for n in 1:length(noisestimratios)
    for d in 1:length(driveamplituderatios)
        push!(df_NGNMM_slow,("NGNMM slow",noisestimratios[n],driveamplituderatios[d],statistics[1,n,d,1],statistics[1,n,d,2],statistics[1,n,d,3],statistics[1,n,d,4],statistics[1,n,d,5],statistics[1,n,d,6]))
    end
end
CSV.write("./ppaper_correlation_statistics_table_NGNMM_slow_SI.csv",df_NGNMM_slow)

df_NGNMM_fast=DataFrame(Model=String[],NSR=Float64[],DAR=Float64[],r=Float64[],CI_lower=Float64[],CI_upper=Float64[],p=Float64[],p_bonferroni=Float64[],p_FDR=Float64[])
for n in 1:length(noisestimratios)
    for d in 1:length(driveamplituderatios)
        push!(df_NGNMM_fast,("NGNMM fast",noisestimratios[n],driveamplituderatios[d],statistics[2,n,d,1],statistics[2,n,d,2],statistics[2,n,d,3],statistics[2,n,d,4],statistics[2,n,d,5],statistics[2,n,d,6]))
    end
end
CSV.write("./ppaper_correlation_statistics_table_NGNMM_fast_SI.csv",df_NGNMM_fast)

for (i,model) in enumerate(models[3:end])
    df_model=DataFrame(Model=String[],NSR=Float64[],r=Float64[],CI_lower=Float64[],CI_upper=Float64[],p=Float64[],p_bonferroni=Float64[],p_FDR=Float64[])
    model_index=i+1
    for n in 1:length(noisestimratios)
        push!(df_model,(model,noisestimratios[n],statistics[model_index,n,1,1],statistics[model_index,n,1,2],statistics[model_index,n,1,3],statistics[model_index,n,1,4],statistics[model_index,n,1,5],statistics[model_index,n,1,6]))
    end
    CSV.write("./ppaper_correlation_statistics_table_$(replace(lowercase(replace(model," "=>"_")),"-" => "_"))_SI.csv",df_model)
end


## 
############################
############################
############################
# Plot the Scatters for the DAR=12.2 and NSR=0.9 case for each model in  a row BUT WITH SQUARED ITPC as in Oanas data.: MAIN RESULT FIGURE.
##############################
############################
#############################
gr()
default(xguidefontsize=11)
default(yguidefontsize=11)
default(xtickfontsize=9)
default(ytickfontsize=9)
default(legendfontsize=9)
default(fontfamily="computer modern")

Bang_wong_color_palette=[(230,159,0),(0,114,178),(0,158,115),(204,121,167),(86,180,233),(240,228,66),(0,0,0)]
Bang_wong_color_palette_normalised=[(col[1]/255.0,col[2]/255.0,col[3]/255.0) for col in Bang_wong_color_palette]
using ColorSchemes,Colors
color_scheme=ColorScheme([Colors.RGB(col...) for col in Bang_wong_color_palette_normalised],"bang wong colorblind friendly")

width_pts = 469*2 #single column.  # or any other value
inches_per_points = 1/72.27 
width_inches = width_pts *inches_per_points
width_px= width_inches*72.27  # or  width_inches*DPI.. /2 as two square plots per column approx.


symbols=[:circle :cross :star :diamond :hexagon :octagon :square :plus :x :utriangle :dtriangle]
linestyles=[:solid :dash :dot :dashdot :dashdotdot :solid :dash :dot :dashdot :dashdotdot]
Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]

time2pds=Frequency_varied_phasereset_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]
fourHzITPCS=Vector{Array{Float64,3}}()
model_data=[C0066D025_data, Aine_data, phase_reset_data, nophasereset_data, evoked_derivativestim_data, evoked_envelopestim_data, evoked_peakenvstim_data, evoked_peakratestim_data];
experimental_plot=plot()
scatter!(experimental_plot, time2pds, sorted_oana_ITPCs[:], label="EEG", xlabel="", ylabel=L"4\textrm{Hz}~~R^2", linewidth=2,color=color_scheme[1],ylims=(0.0,0.31),smooth=true,markerstrokewidth=1,legend=true)
this_dar_array=Array{Float64,2}(undef,10,15) #4 models, 15 conditions
sorting_indxs=get_sorting_indices("./")
ngnmm_lables=["NGNMM slow", "NGNMM fast"]
NGNMM_plot=plot()
for (i,model) in enumerate(model_data[1:2])
    this_dar_array[i,:]=reduce(hcat,[mean(x->x.^2,model[8*3+7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19] for condition in Condition_keys])'
    range_data=reduce(hcat,[reduce(hcat,model[8*3+7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19,:] for condition in Condition_keys])
    @show size(range_data)
    upper_bounds=maximum(x->x^2,range_data,dims=1)[sorting_indxs]
    @show upper_bounds
    lower_bounds=minimum(x->x^2,range_data,dims=1)[sorting_indxs]
    scatter!(NGNMM_plot, time2pds, this_dar_array[i,sorting_indxs], yerror=(this_dar_array[i,sorting_indxs]-lower_bounds,upper_bounds.-this_dar_array[i,sorting_indxs]),label=ngnmm_lables[i], linewidth=2, color=color_scheme[2], smooth=true, markerstrokewidth=1,markershape=symbols[i],linestyle=linestyles[i],ylims=(0.0,0.3),legend=true)
end
display(NGNMM_plot)

phase_reset_plot=plot()
labels_pr=["Phase reset",  "No phase reset"]
idxs=[3,6] #indices of phase reset models in model_data
for (i,model) in enumerate(model_data[[3,4]])
    this_dar_array[idxs[i],:]=reduce(hcat,[mean(x->x.^2,model[7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19] for condition in Condition_keys])'
    range_data=reduce(hcat,[reduce(hcat,model[7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19,:] for condition in Condition_keys])
    @show size(range_data)
    upper_bounds=maximum(x->x^2,range_data,dims=1)[sorting_indxs]
    @show upper_bounds
    lower_bounds=minimum(x->x^2,range_data,dims=1)[sorting_indxs]
    @show lower_bounds

    scatter!(phase_reset_plot, time2pds, this_dar_array[idxs[i],sorting_indxs], yerr=(this_dar_array[idxs[i],sorting_indxs]-lower_bounds,upper_bounds.-this_dar_array[idxs[i],sorting_indxs]),label=labels_pr[i],xlabel="L (s)", linewidth=2, color=color_scheme[3], smooth=true, markerstrokewidth=1,markershape=symbols[i],linestyle=linestyles[i],ylims=(0.0,0.4),legend=true)
end
display(phase_reset_plot)
#saving:

q=plot((experimental_plot, NGNMM_plot, phase_reset_plot, )..., layout=(1,3), size=(width_px,width_px/3),ylims=(0.0,0.3001), xlabel=L"L\textrm{~(ms)}", legend_font_halign=:right,bottom_margin=10.0Plots.mm,right_margin=0.5Plots.mm,left_margin=8Plots.mm)

tosave=q
width, height = tosave.attr[:size]
Plots.prepare_output(tosave)
savefig(tosave,"./phonemepaper_mainresultfig_v2.pdf")


############################
############################
############################
# figure 2: evoked model results"
##############################
############################
#############################
time2pds=Frequency_varied_phasereset_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]
fourHzITPCS=Vector{Array{Float64,3}}()
evoked_model_data=[evoked_envelopestim_data, evoked_derivativestim_data, evoked_peakenvstim_data, evoked_peakratestim_data];
aexperimental_plot=plot()
# scatter!(experimental_plot, time2pds, sorted_oana_ITPCs[:], label="EEG", xlabel="", ylabel="ITPC 4Hz", linewidth=2,color=color_scheme[1],ylims=(0.0,0.31),smooth=true,markerstrokewidth=1,legend=false)
scatter!(experimental_plot, time2pds, sorted_oana_ITPCs[:], label="", xlabel="", ylabel="ITPC 4Hz", linewidth=2,color=color_scheme[1],ylims=(0.0,0.31),smooth=true,markerstrokewidth=1,legend=false)
this_dar_array=Array{Float64,2}(undef,4,15) #4 models, 15 conditions
evoked_model_plots=[]
labels_evkd=["Evoked - envelope", "Evoked - derivative", "Evoked - peakenv", "Evoked - peakrate"]
idxs=[1,2,3,4] #indices of evoked models in model_data
for (i,model) in enumerate(evoked_model_data) 
    p=plot()
    #plot experimental data in the background of each plot.
    if i==1
        scatter!(p, time2pds, sorted_oana_ITPCs[:], label="EEG", xlabel="", ylabel=L"4\textrm{Hz}~~R^2", linewidth=2,color=color_scheme[1],ylims=(0.0,0.31),smooth=true,markerstrokewidth=1,legend=true)
    else
        scatter!(p, time2pds, sorted_oana_ITPCs[:], label="EEG", xlabel="", ylabel="", linewidth=2,color=color_scheme[1],ylims=(0.0,0.31),smooth=true,markerstrokewidth=1,legend=true)
    end
    # scatter!(evoked_model_plot, time2pds, this_dar_array[idxs[i],:].^2, label=labels_evkd[i], linewidth=2, xlabel="Latency to\n peak derivative (s)",color=color_scheme[4], smooth=true, markerstrokewidth=1,markershape=symbols[i],linestyle=linestyles[i],ylims=(0.0,0.55),legend=false)

    this_dar_array[idxs[i],:]=reduce(hcat,[mean(x->x.^2,model[7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19] for condition in Condition_keys])'
    range_data=reduce(hcat,[reduce(hcat,model[7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19,:] for condition in Condition_keys])
    @show size(range_data)
    upper_bounds=maximum(x->x^2,range_data,dims=1)[sorting_indxs]
    @show upper_bounds
    lower_bounds=minimum(x->x^2,range_data,dims=1)[sorting_indxs]
    @show lower_bounds  
    if i==3 #plot both peak rate and peakenv in the same plot
        this_dar_array[4,:]=reduce(hcat,[mean(x->x.^2,evoked_model_data[4][7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19] for condition in Condition_keys])'
        fourth_range_data=reduce(hcat,[reduce(hcat,evoked_model_data[4][7][1]["noisestimratio: 0.9"]["ITPCs"][condition])[19,:] for condition in Condition_keys])
        fourth_upper_bounds=maximum(x->x^2,fourth_range_data,dims=1)[sorting_indxs]
        fourth_lower_bounds=minimum(x->x^2,fourth_range_data,dims=1)[sorting_indxs]


        scatter!(p, time2pds, this_dar_array[i,sorting_indxs], label=labels_evkd[i], yerr=(this_dar_array[i,sorting_indxs]-lower_bounds,upper_bounds.-this_dar_array[i,sorting_indxs]), linewidth=2, xlabel="L (s)",color=color_scheme[4], smooth=true, markerstrokewidth=1,markershape=symbols[3],linestyle=linestyles[1],ylims=(0.0,0.55),
        legend=true)
        println("plotting both peak rate and peakenv in the same plot")
        scatter!(p, time2pds, this_dar_array[4,sorting_indxs], label=labels_evkd[4], yerr=(this_dar_array[4,sorting_indxs]-fourth_lower_bounds,fourth_upper_bounds-this_dar_array[4,sorting_indxs]), linewidth=2, xlabel="L (s)",color=color_scheme[4], smooth=true, markerstrokewidth=1,markershape=symbols[1+3],linestyle=linestyles[2],ylims=(0.0,0.55),legend=true)
        push!(evoked_model_plots,p)
        break
    else
        scatter!(p, time2pds, this_dar_array[i,sorting_indxs], label=labels_evkd[i], yerr=(this_dar_array[i,sorting_indxs]-lower_bounds,upper_bounds.-this_dar_array[i,sorting_indxs]), linewidth=2, xlabel="L (s)",color=color_scheme[4], smooth=true, markerstrokewidth=1,markershape=symbols[3],linestyle=linestyles[1],ylims=(0.0,0.55),
        legend=true)
        push!(evoked_model_plots,p)
    end
    
end
plot(evoked_model_plots..., layout=(1,3), size=(width_px,width_px/3),ylims=(0.0,0.3001), xlabel=L"L\textrm{~(ms)}", legend_font_halign=:right,bottom_margin=10.0Plots.mm,right_margin=0.5Plots.mm,left_margin=8Plots.mm)

#save figure 2:
tosave=plot(evoked_model_plots..., layout=(1,3), size=(width_px,width_px/3),ylims=(0.0,0.3001), xlabel="L (ms)", legend_font_halign=:right,bottom_margin=10.0Plots.mm,right_margin=0.5Plots.mm,left_margin=8Plots.mm)
width, height = tosave.attr[:size]
Plots.prepare_output(tosave)
savefig(tosave,"./ppaper_evoked_model_results_v2.pdf")

#### all six plots together? to get without EEG on evoked, just run the above evoked code commenting out the experimental plot part.
tosave=plot((experimental_plot, NGNMM_plot, phase_reset_plot, evoked_model_plots...)..., layout=(2,3), size=(width_px,width_px/1.5),ylims=(0.0,0.3001), xlabel=L"L\textrm{~(ms)}", legend_font_halign=:right,bottom_margin=10.0Plots.mm,right_margin=0.5Plots.mm,left_margin=8Plots.mm)
width, height = tosave.attr[:size]
Plots.prepare_output(tosave)
savefig(tosave,"./ppaper_combined_main_results.pdf")


##############################
##############################
##############################
# plotting similar, but on x axis is the phonemes, and across in rows are the models, lines connecting each models points so you can see if they overlap or are just a shift of each other.
##############################
##############################
##############################

time2pds=Frequency_varied_phasereset_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]
model_data_select=[sorted_oana_ITPCs, C0066D025_data, phase_reset_data,  evoked_derivativestim_data]
model_lables=["EEG" "Slow NGNMM" "Phase Reset"  "Evoked - derivative stim"]
Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
ITPCs_4Hz_per_model_sorted=Float64.([model_data_select[1] reduce(hcat,[mean(x->x.^2,model_data_select[2][8*(4-1)+7][1]["noisestimratio: 0.9"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,model_data_select[3][8*(1-1)+7][1]["noisestimratio: 0.9"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,model_data_select[4][7][1]["noisestimratio: 0.9"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")]])

scatter(time2pds,ITPCs_4Hz_per_model_sorted[:,1],label="EEG",xlabel="Phoneme",ylabel="4Hz ITPC",legend=:outertopright,ylims=(0.0,0.3),markershape=:circle,markercolor=color_scheme[1],linewidth=2)
for i in 2:4
    scatter!(time2pds,ITPCs_4Hz_per_model_sorted[:,i],label=model_lables[i],markercolor=color_scheme[i],linewidth=2)
end
display(plot!(size=(900,600),xticks=(time2pds,Condition_keys),markerstrokewidth=0.5))
# Add coloured lines joining corresponding points
for i in 1:4
    plot!(time2pds, ITPCs_4Hz_per_model_sorted[:,i], color=color_scheme[i], linewidth=1, label=false)
end
display(plot!())

width_pts = 469*2 #single column.  # or any other value
inches_per_points = 1.0/72.27
width_inches = width_pts *inches_per_points
width_px= width_inches*72.27  # or  width_inches*DPI.. /2 as two square plots per column approx.
label_fontsize=12
legend_fontsize=9
xtick_fontsize=9
ytick_fontsize=9

scatter_symbols=[:circle :hexagon :utriangle :diamond]

scatter(ITPCs_4Hz_per_model_sorted[:,1],label="EEG",xlabel="Phoneme",ylabel="4Hz ITPC",ylims=(0.0,0.3),markershape=scatter_symbols[1],markercolor=color_scheme[1],linewidth=2)
for i in 2:4
    scatter!(ITPCs_4Hz_per_model_sorted[:,i],label=model_lables[i],markercolor=color_scheme[i],linewidth=2,markershape=scatter_symbols[i])
end
display(plot!(size=(width_px,width_px),xticks=(1:15,Condition_keys),markerstrokewidth=0.5))

display(plot!(grid=true))
#re-sort the data so it is grouped by phoneme type:
#group identities: 1:vowel, 2:voiced-stops, 3:unvoiced-stops, 4:nasals, 5:sibilants, 6:liquidsm 7:frivatives
phoneme_groups=[2,4,4,2,1,6,2,3,3,3,6,7,5,7,5]
Condition_keys=Condition_keys
Q=collect(zip(phoneme_groups,1:1:(length(phoneme_groups))))
sortedQ=sort(Q, by=first)
sortingIndexes=[x[2] for x in sortedQ]
grouped_ITPCs_4Hz_per_model_sorted=Array{Float64,2}(undef,15,4)
for i in 1:4
    grouped_ITPCs_4Hz_per_model_sorted[:,i]=ITPCs_4Hz_per_model_sorted[sortingIndexes,i]
end
### without time2pd on x axis, just phonemes
scatter()
plot_colors=[:green,:orange,:yellow,:blue,:pink,:cyan,:purple]
phoneme_colors=[2,4,4,2,1,6,2,3,3,3,6,7,5,7,5]
sorted_phoneme_colors=plot_colors[phoneme_colors]
grouped_phoneme_colors=sorted_phoneme_colors[sortingIndexes]
sorted_conditions=Condition_keys[sortingIndexes]
#### add vspan shaded regions for each phoneme group:
for i in 1:15
    vspan!([i-0.5,i+0.5], color=grouped_phoneme_colors[i], alpha=0.4, label="")
end
scatter!(1:15,grouped_ITPCs_4Hz_per_model_sorted[:,1],label="EEG",xlabel="phoneme",ylabel=L"4\textrm{Hz}~~R^2",ylims=(0.0,0.4),markershape=scatter_symbols[1],markercolor=color_scheme[1],linewidth=2)
for i in 2:4
    scatter!(1:15,grouped_ITPCs_4Hz_per_model_sorted[:,i],label=model_lables[i],markercolor=color_scheme[i],linewidth=2,markershape=scatter_symbols[i])
end
display(plot!(size=(width_px/2,width_px/3),xticks=(1:15,Condition_keys),markerstrokewidth=0.5))

display(plot!(grid=true))
display(plot!())
x=plot!()
tosave=x
width, height = tosave.attr[:size]
Plots.prepare_output(tosave)
PlotlyJS.savefig(Plots.plotlyjs_syncplot(tosave),"./ppaper_ITPC_scatter_phonemex.pdf",width=Int64(ceil(width)), height=Int64(ceil(height)))
savefig(tosave,"./ppaper_ITPC_scatter_phonemex.pdf")


##############
##############
##############
##############
#final modulation and frequency vary plots. just want to plot the mean ITPC across conditions as we vary each parameter. 
##############
##############
##############
symbols_varyfreq=[:circle :cross :diamond :hexagon :octagon :square :x :utriangle :dtriangle :star]
Frequency_varied_phasereset_sqrITPCs=Array{Float64,2}(undef,length(Frequencies[1:2:end]),15) #frequencies, 15 conditions
raw3_freq_varied_phasereset_sqrITPCs=Array{Float64,3}(undef,length(Frequencies[1:2:end]),15,3) #frequencies, 15 conditions, 3 raw ITPCs.
for (fidx,Freq) in enumerate(Frequencies[1:2:end])
    Frequency_varied_phasereset_sqrITPCs[fidx,:]=reduce(hcat,[mean(x->x.^2,Frequency_varied_phasereset_data[1+2*(fidx-1)][1]["noisestimratio: $(0.1)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")]
    raw3_freq_varied_phasereset_sqrITPCs[fidx,:,1]=reduce(hcat,[reduce(hcat,Frequency_varied_phasereset_data[1+2*(fidx-1)][1]["noisestimratio: $(0.1)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2
    raw3_freq_varied_phasereset_sqrITPCs[fidx,:,2]=reduce(hcat,[reduce(hcat,Frequency_varied_phasereset_data[1+2*(fidx-1)][1]["noisestimratio: $(0.1)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2
    raw3_freq_varied_phasereset_sqrITPCs[fidx,:,3]=reduce(hcat,[reduce(hcat,Frequency_varied_phasereset_data[1+2*(fidx-1)][1]["noisestimratio: $(0.1)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2
end


modulated_nonrectified_ITPCs=Array{Float64,2}(undef,11,15) #11 modulations, 15 conditions
raw3_modulated_nonrectified_ITPCs=Array{Float64,3}(undef,11,15,3) #11 modulations, 15 conditions, 3 raw ITPCs.
for (midx,modulation) in enumerate(modulations)
    modulated_nonrectified_ITPCs[midx,:]=reduce(hcat,[mean(x->x.^2,modulated_nonrectified_data[midx][1]["noisestimratio: $(0.1)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")]
    raw3_modulated_nonrectified_ITPCs[midx,:,1]=reduce(hcat,[reduce(hcat,modulated_nonrectified_data[midx][1]["noisestimratio: $(0.1)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2
    raw3_modulated_nonrectified_ITPCs[midx,:,2]=reduce(hcat,[reduce(hcat,modulated_nonrectified_data[midx][1]["noisestimratio: $(0.1)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2
    raw3_modulated_nonrectified_ITPCs[midx,:,3]=reduce(hcat,[reduce(hcat,modulated_nonrectified_data[midx][1]["noisestimratio: $(0.1)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2
end
time2pds=modulated_nonrectified_data[1][1]["noisestimratio: 0.1"]["sortedITPCS_t2PD"][4]
# labels=reshape(string.(modulations),1,11)
#add "m =" to labels
# labels=["m = "*i for i in labels]
# vary_mod_plot=scatter(time2pds,modulated_nonrectified_ITPCs',size=(width_px,width_px/2),label=labels,xlabel="latency to peak derivative (s)",ylabel="ITPC",markershape=symbols_varyfreq,color_palette=color_scheme,legend=:outertopright,ylims=(0.0,1.0))

raw3_modulated_nonrectified_ITPCs[:,:,1]
#plot mean ITPC across conditions.
mean_modulated_nonrectified_ITPCs=mean(modulated_nonrectified_ITPCs,dims=2)
#raw 3 test 1 scatter
modulationplot=scatter(modulations,mean(raw3_modulated_nonrectified_ITPCs[:,:,1],dims=2),
    xlabel="phase resetting modulation",
    ylabel="Mean ITPC",
    legend=:false,
    xticks=(modulations, modulations),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols_varyfreq,
    markersize=3,
    markeralpha=0.5,
    color=color_scheme[3],
    size=(width_px/2,width_px/3))

scatter!(modulations,mean(raw3_modulated_nonrectified_ITPCs[:,:,2],dims=2),
    xlabel="phase resetting modulation",
    ylabel="Mean ITPC",
    legend=:false,
    xticks=(modulations, modulations),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols_varyfreq,
    markersize=3,
    markeralpha=0.5,
    color=color_scheme[3])

scatter!(modulations,mean(raw3_modulated_nonrectified_ITPCs[:,:,3],dims=2),
    xlabel="phase resetting modulation",
    ylabel="Mean ITPC",
    legend=:false,
    xticks=(modulations, modulations),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols_varyfreq,
    markersize=3,
    markeralpha=0.5,
    color=color_scheme[3])


plot!(modulationplot,modulations,mean_modulated_nonrectified_ITPCs,
    xlabel=L"m",
    ylabel=L"\textrm{mean}~~4\textrm{Hz}~~R^2",
    legend=:false,
    xticks=(modulations, modulations),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols_varyfreq,
    markersize=3,
    color=color_scheme[3])

#and next to this (below it) plot the effect of frequency.
raw3_freq_varied_phasereset_sqrITPCs
frequencyplot=scatter(Frequencies[1:2:end],mean(raw3_freq_varied_phasereset_sqrITPCs[:,:,1],dims=2),
    xlabel="Frequency (Hz)",
    ylabel="Mean ITPC",
    legend=:false,
    xticks=(Frequencies[1:2:end], Frequencies[1:2:end]),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols[1],
    markersize=3,
    markeralpha=0.5,
    color=color_scheme[3],
    size=(width_px/2,width_px/3))

scatter!(frequencyplot,Frequencies[1:2:end],mean(raw3_freq_varied_phasereset_sqrITPCs[:,:,2],dims=2),
    xlabel="Frequency (Hz)",
    ylabel="Mean ITPC",
    legend=:false,
    xticks=(Frequencies[1:2:end], Frequencies[1:2:end]),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols[1],
    markersize=3,
    markeralpha=0.5,
    color=color_scheme[3])

scatter!(frequencyplot,Frequencies[1:2:end],mean(raw3_freq_varied_phasereset_sqrITPCs[:,:,3],dims=2),
    xlabel="Frequency (Hz)",
    ylabel="Mean ITPC",
    legend=:false,
    xticks=(Frequencies[1:2:end], Frequencies[1:2:end]),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols[1],
    markersize=3,
    markeralpha=0.5,
    color=color_scheme[3])      

plot!(frequencyplot,Frequencies[1:2:end],mean(Frequency_varied_phasereset_sqrITPCs,dims=2),
    xlabel="frequency (Hz)",
    ylabel=L"\textrm{mean}~~4\textrm{Hz}~~R^2",
    legend=:false,
    xticks=(Frequencies[1:2:end], Frequencies[1:2:end]),
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols[1],
    markersize=3,
    color=color_scheme[3])

plot(frequencyplot,modulationplot,
    layout=(2,1),
    size=(width_px/2,width_px/3))

#savefig

tosave=plot(frequencyplot, modulationplot,
    layout=(2,1),
    size=(width_px/2,width_px/3))
width, height = tosave.attr[:size]
Plots.prepare_output(tosave)
PlotlyJS.savefig(Plots.plotlyjs_syncplot(tosave),"./ppaper_modulation_frequency_vary.pdf",width=Int64(ceil(width)), height=Int64(ceil(height)))
savefig(tosave,"./ppaper_modulation_frequency_vary.pdf")

##########################
##########################
##########################
#evaluating mean 4Hz ITPC over phonemes for each model across λ (Noise-stimulus-ratios).
##########################
##########################
##########################
sorting_indxs=get_sorting_indices("./")
Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
gr()
default(xguidefontsize=11)
default(yguidefontsize=11)
default(xtickfontsize=9)
default(ytickfontsize=9)
default(legendfontsize=9)
default(fontfamily="computer modern")

Bang_wong_color_palette=[(230,159,0),(0,114,178),(0,158,115),(204,121,167),(86,180,233),(240,228,66),(0,0,0)]
Bang_wong_color_palette_normalised=[(col[1]/255.0,col[2]/255.0,col[3]/255.0) for col in Bang_wong_color_palette]
using ColorSchemes,Colors
color_scheme=ColorScheme([Colors.RGB(col...) for col in Bang_wong_color_palette_normalised],"bang wong colorblind friendly")

symbols=[:circle  :hexagon  :utriangle :diamond :plus :x :rect :cross :star]
linestyles=[:solid :dash :dot :dashdot :dashdotdot :solid :solid :solid :solid] # for the 9 models.

width_pts = 469*2 #single column.  # or any other value
inches_per_points = 1/72.27 
width_inches = width_pts *inches_per_points
width_px= width_inches*72.27  # or  width_inches*DPI.. /2 as two square plots per column approx.

NSRplot_model_data=[sorted_oana_ITPCs, Aine_data, C0066D025_data, phase_reset_data, nophasereset_data,  evoked_envelopestim_data,evoked_derivativestim_data,evoked_peakenvstim_data,evoked_peakratestim_data];
Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
sqr_ITPCs_4Hz_per_model_sorted=Array{Float64,3}(undef, 15,9,length(noisestimratios))
max_ITPCs_4Hz_per_model_sorted=Array{Float64,3}(undef, 15,9,length(noisestimratios)) # to store the max ITPCs for each model and NSR.
min_ITPCs_4Hz_per_model_sorted=Array{Float64,3}(undef, 15,9,length(noisestimratios)) # to store the min ITPCs for each model and NSR.
raw_3_ITPCs_4Hz_per_model_sorted=Array{Float64,4}(undef, 15,9,length(noisestimratios),3) # to store the raw ITPCs for each model and NSR.
for (i,NSR) in enumerate(noisestimratios)
    max_ITPCs_4Hz_per_model_sorted[:,:,i]=Float64.([NSRplot_model_data[1] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[2][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[3][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[4][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[5][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[6][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[7][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[8][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] maximum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[9][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs]])
    min_ITPCs_4Hz_per_model_sorted[:,:,i]=Float64.([NSRplot_model_data[1] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[2][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[3][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[4][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[5][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[6][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[7][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[8][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs] minimum(x->x^2,reduce(hcat,[reduce(hcat,NSRplot_model_data[9][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[19,:] for condition in Condition_keys]),dims=1)[sorting_indxs]])
    sqr_ITPCs_4Hz_per_model_sorted[:,:,i]=Float64.([NSRplot_model_data[1] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[2][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[3][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[4][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[5][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[6][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[7][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[8][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")] reduce(hcat,[mean(x->x.^2,NSRplot_model_data[9][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")]])
    raw_3_ITPCs_4Hz_per_model_sorted[:,:,i,1]=Float64.([NSRplot_model_data[1] reduce(hcat,[reduce(hcat,NSRplot_model_data[2][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[3][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[4][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[5][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[6][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[7][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[8][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[9][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2])
    raw_3_ITPCs_4Hz_per_model_sorted[:,:,i,2]=Float64.([NSRplot_model_data[1] reduce(hcat,[reduce(hcat,NSRplot_model_data[2][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[3][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[4][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[5][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[6][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[7][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[8][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[9][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2])
    raw_3_ITPCs_4Hz_per_model_sorted[:,:,i,3]=Float64.([NSRplot_model_data[1] reduce(hcat,[reduce(hcat,NSRplot_model_data[2][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[3][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[4][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[5][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[6][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[7][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[8][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[9][i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2])
end

raw_3_ITPCs_4Hz_per_model_sorted[:,:,7,1] # this is the raw ITPCs for the 4Hz band, for the 7th NSR.
raw_3_ITPCs_4Hz_per_model_sorted[:,:,7,2] # this is the raw ITPCs for the 4Hz band, for the 7th NSR.

sqr_ITPCs_4Hz_per_model_sorted[:,:,7]
max_ITPCs_4Hz_per_model_sorted[:,:,7]
min_ITPCs_4Hz_per_model_sorted[:,:,7]

mean4Hz_ITPC_oana_to_models=Array{Float64,2}(undef,9,length(noisestimratios)) #9 model, 8 NSRs.
max_data_mean4Hz_ITPC_oana_to_models=Array{Float64,2}(undef,9,length(noisestimratios)) #9 model, 8 NSRs.
min_data_mean4Hz_ITPC_oana_to_models=Array{Float64,2}(undef,9,length(noisestimratios)) #9 model, 8 NSRs.
mean4Hz_ITPC_oana_to_models_test1=Array{Float64,2}(undef,9,length(noisestimratios)) #9 model, 8 NSRs.
mean4Hz_ITPC_oana_to_models_test2=Array{Float64,2}(undef,9,length(noisestimratios)) #9 model, 8 NSRs.
mean4Hz_ITPC_oana_to_models_test3=Array{Float64,2}(undef,9,length(noisestimratios)) #9 model, 8 NSRs.

for (i,NSR) in enumerate(noisestimratios)
    for j in 1:9
        mean4Hz_ITPC_oana_to_models[j,i]=mean(sqr_ITPCs_4Hz_per_model_sorted[:,j,i])
        max_data_mean4Hz_ITPC_oana_to_models[j,i]=mean(max_ITPCs_4Hz_per_model_sorted[:,j,i])
        min_data_mean4Hz_ITPC_oana_to_models[j,i]=mean(min_ITPCs_4Hz_per_model_sorted[:,j,i])

        mean4Hz_ITPC_oana_to_models_test1[j,i]=mean(raw_3_ITPCs_4Hz_per_model_sorted[:,j,i,1])
        mean4Hz_ITPC_oana_to_models_test2[j,i]=mean(raw_3_ITPCs_4Hz_per_model_sorted[:,j,i,2])
        mean4Hz_ITPC_oana_to_models_test3[j,i]=mean(raw_3_ITPCs_4Hz_per_model_sorted[:,j,i,3])
    end
end

mean4Hz_ITPC_oana_to_models #columns are NSRs, rows are models.
max_data_mean4Hz_ITPC_oana_to_models
min_data_mean4Hz_ITPC_oana_to_models
@show (mean4Hz_ITPC_oana_to_models.-min_data_mean4Hz_ITPC_oana_to_models,max_data_mean4Hz_ITPC_oana_to_models.-mean4Hz_ITPC_oana_to_models)
mean4Hz_ITPC_oana_to_models_plot=scatter(noisestimratios,mean4Hz_ITPC_oana_to_models_test1', markeralpha=0.5,label=["EEG" "NGNMM fast" "NGNMM slow" "Phase Reset" "No Phase Reset" "Evoked - env stim" "Evoked - deriv stim" "Evoked - peak env stim" "Evoked - peak rate stim"],
    xlabel=L"\lambda",
    ylabel=L"\textrm{mean}\;\;4\textrm{Hz}\;\; R^2",
    legend=:false,
    ylims=(-0.0,1.0),
    linewidth=2,
    markershape=symbols,
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2] color_scheme[3] color_scheme[3] color_scheme[4] color_scheme[4] color_scheme[4] color_scheme[4]],size=(width_px/2,width_px/3))
scatter!(mean4Hz_ITPC_oana_to_models_plot,noisestimratios,mean4Hz_ITPC_oana_to_models_test2', markeralpha=0.5,label=["EEG" "NGNMM fast" "NGNMM slow" "Phase Reset" "No Phase Reset" "Evoked - env stim" "Evoked - deriv stim" "Evoked - peak env stim" "Evoked - peak rate stim"],
    legend=:false,
    linewidth=2,
    markershape=symbols,
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2] color_scheme[3] color_scheme[3] color_scheme[4] color_scheme[4] color_scheme[4] color_scheme[4]],size=(width_px/2,width_px/3))
scatter!(mean4Hz_ITPC_oana_to_models_plot,noisestimratios,mean4Hz_ITPC_oana_to_models_test3', markeralpha=0.5,label=["EEG" "NGNMM fast" "NGNMM slow" "Phase Reset" "No Phase Reset" "Evoked - env stim" "Evoked - deriv stim" "Evoked - peak env stim" "Evoked - peak rate stim"],
    legend=:false,
    linewidth=2,
    markershape=symbols,
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2] color_scheme[3] color_scheme[3] color_scheme[4] color_scheme[4] color_scheme[4] color_scheme[4]],size=(width_px/2,width_px/3))
plot!(mean4Hz_ITPC_oana_to_models_plot,noisestimratios,mean4Hz_ITPC_oana_to_models', markeralpha=0,label=["EEG" "NGNMM fast" "NGNMM slow" "Phase Reset" "No Phase Reset" "Evoked - env stim" "Evoked - deriv stim" "Evoked - peak env stim" "Evoked - peak rate stim"],
    legend=:false,
    linewidth=2,
    markershape=symbols,
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2] color_scheme[3] color_scheme[3] color_scheme[4] color_scheme[4] color_scheme[4] color_scheme[4]],size=(width_px/2,width_px/3))


mock_plot_for_legend_1=plot(abs.(mean4Hz_ITPC_oana_to_models[1:5,:])'.*0.0, label=["EEG" "NGNMM fast" "NGNMM slow" "Phase Reset" "No Phase Reset"],
    xlabel="",
    grid=false,
    axis=false,
    yaxis=false,
    ylabel="",
    legend=:top,
    ylims=(0.9,1.0),
    linewidth=2,
    markershape=permutedims(symbols[1:5]),
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2] color_scheme[3] color_scheme[3]],size=(width_px/2,width_px/3))   
mock_plot_for_legend_2=plot(abs.(mean4Hz_ITPC_oana_to_models[6:9,:])'.*0.0, label=["Evoked - env stim" "Evoked - deriv stim" "Evoked - peak env stim" "Evoked - peak rate stim"],
    xlabel="",
    grid=false,
    axis=false,
    yaxis=false,
    ylabel="",
    legend=:top,
    ylims=(0.9,1.0),
    linewidth=2,
    markershape=permutedims(symbols[1,6:9]),
    markersize=3,
    color=[color_scheme[4] color_scheme[4] color_scheme[4] color_scheme[4]],size=(width_px/2,width_px/3))  



l = @layout [
    a{0.7h}; [ b{0.2h} c{0.2h} ]
]   

plot(mean4Hz_ITPC_oana_to_models_plot, mock_plot_for_legend_1,mock_plot_for_legend_2,
    layout=l,size=(width_px/2,width_px/3),
)

#savefig
tosave=plot(mean4Hz_ITPC_oana_to_models_plot, mock_plot_for_legend_1,mock_plot_for_legend_2,
    layout=l,size=(width_px/2,width_px/3),
)
width, height = tosave.attr[:size]
Plots.prepare_output(tosave)
PlotlyJS.savefig(Plots.plotlyjs_syncplot(tosave),"./ppaper_varynsr_mean4Hz_ITPC_v2.pdf",width=Int64(ceil(width)), height=Int64(ceil(height)))
savefig(tosave,"./ppaper_varynsr_mean_ITPC_0_1_yrange_v2.pdf")


###############################
###############################
###############################
#now plotting the effect of drive strength on the two NGNMM models' mean 4Hz ITPCs across phonemes. Same as we did for the noisestimratios but for the driveamplituderatios.
###############################
###############################
###############################

drive_strength_model_data=[sorted_oana_ITPCs, Aine_data, C0066D025_data]
drive_strength_plot_model_ITPCs=Array{Float64,3}(undef,15,length(driveamplituderatios),length(drive_strength_model_data)) # 3 models, 8 driveamplituderatios

driveamplituderatios
noisestimratios
Aine_data
#NGNMM data is stored having looped per DAR, and for each DAR it adds data for each NSR. there are 6 DARs and 8 NSRs, so 48 datasets. and i want the 0.9 NSR (second to last NSR) 
raw_3_drive_strength_model_ITPCs=Array{Float64,4}(undef,15,length(driveamplituderatios),length(drive_strength_model_data),3) # 3 models, 8 driveamplituderatios, 3 ITPCs.
for (i,driveamplituderatio) in enumerate(driveamplituderatios)
    drive_strength_plot_model_ITPCs[:,i,1]=Float64.(drive_strength_model_data[1])
    drive_strength_plot_model_ITPCs[:,i,2]=Float64.(reduce(hcat,[mean(x->x.^2,drive_strength_model_data[2][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")])
    drive_strength_plot_model_ITPCs[:,i,3]=Float64.(reduce(hcat,[mean(x->x.^2,drive_strength_model_data[3][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition]) for condition in Condition_keys])'[:,19][get_sorting_indices("./")])
    raw_3_drive_strength_model_ITPCs[:,i,1,1]=Float64.(drive_strength_model_data[1])
    raw_3_drive_strength_model_ITPCs[:,i,1,2]=Float64.(drive_strength_model_data[1])
    raw_3_drive_strength_model_ITPCs[:,i,1,3]=Float64.(drive_strength_model_data[1])
        # raw_3_ITPCs_4Hz_per_model_sorted[:,:,i,1]=Float64.([NSRplot_model_data[1] reduce(hcat,[reduce(hcat,NSRplot_model_data[2][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[3][8*(4-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[4][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[5][8*(1-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[6][8*(3-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[7][8*(3-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[8][8*(3-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2 reduce(hcat,[reduce(hcat,NSRplot_model_data[9][8*(3-1)+i][1]["noisestimratio: $(NSR)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2])

    raw_3_drive_strength_model_ITPCs[:,i,2,1]=Float64.(reduce(hcat,[reduce(hcat,drive_strength_model_data[2][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2)
    raw_3_drive_strength_model_ITPCs[:,i,2,2]=Float64.(reduce(hcat,[reduce(hcat,drive_strength_model_data[2][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2)
    raw_3_drive_strength_model_ITPCs[:,i,2,3]=Float64.(reduce(hcat,[reduce(hcat,drive_strength_model_data[2][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2)
    raw_3_drive_strength_model_ITPCs[:,i,3,1]=Float64.(reduce(hcat,[reduce(hcat,drive_strength_model_data[3][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition])[:,1][19] for condition in Condition_keys])[sorting_indxs].^2)
    raw_3_drive_strength_model_ITPCs[:,i,3,2]=Float64.(reduce(hcat,[reduce(hcat,drive_strength_model_data[3][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition])[:,2][19] for condition in Condition_keys])[sorting_indxs].^2)
    raw_3_drive_strength_model_ITPCs[:,i,3,3]=Float64.(reduce(hcat,[reduce(hcat,drive_strength_model_data[3][7+8*(i-1)][1]["noisestimratio: $(0.9)"]["ITPCs"][condition])[:,3][19] for condition in Condition_keys])[sorting_indxs].^2)
end


drive_strength_plot_model_ITPCs
drive_strength_plot_model_ITPCs[:,1,:]
#plot the ITPCs for each model and drive strength.
drive_strength_plot=Array{Any,1}() # 8 driveamplituderatios
for (i,driveamplituderatio) in enumerate(driveamplituderatios)
    pp=scatter(time2pds,drive_strength_plot_model_ITPCs[:,i,:],
        label=permutedims(driveamplituderatios),
        xlabel="latency to peak derivative (s)",
        ylabel="ITPC",
        markershape=symbols_varyfreq,
        color_palette=color_scheme,
        ylims=(0.0,1.0),
        markersize=3,
        grid=true,
        legend=nothing)
        push!(drive_strength_plot,pp)

end

plot(drive_strength_plot..., layout=(1,6),size=(width_px*3,width_px))


mean4Hz_ITPC_drive_strength=Array{Float64,2}(undef,length(driveamplituderatios),3) # 3 models, 8 driveamplituderatios
#mean4Hz_ITPC for each test
mean4Hz_ITPC_drive_strength_test1=Array{Float64,2}(undef,length(driveamplituderatios),3) # 3 models, 8 driveamplituderatios
mean4Hz_ITPC_drive_strength_test2=Array{Float64,2}(undef,length(driveamplituderatios),3) # 3 models, 8 driveamplituderatios
mean4Hz_ITPC_drive_strength_test3=Array{Float64,2}(undef,length(driveamplituderatios),3) # 3 models, 8 driveamplituderatios
for (i,driveamplituderatio) in enumerate(driveamplituderatios)
    for j in 1:3
        mean4Hz_ITPC_drive_strength[i,j]=mean(drive_strength_plot_model_ITPCs[:,i,j])
        mean4Hz_ITPC_drive_strength_test1[i,j]=mean(raw_3_drive_strength_model_ITPCs[:,i,j,1])
        mean4Hz_ITPC_drive_strength_test2[i,j]=mean(raw_3_drive_strength_model_ITPCs[:,i,j,2])
        mean4Hz_ITPC_drive_strength_test3[i,j]=mean(raw_3_drive_strength_model_ITPCs[:,i,j,3])
    end
end
DAR_mean4Hz_ITPCplot=scatter(driveamplituderatios,mean4Hz_ITPC_drive_strength_test1, markeralpha=0.5,
    xlabel=L"D",
    ylabel=L"\textrm{mean}\;\;4\textrm{Hz}\;\;R^2",
    legend=:false,
    label="",
    ylims=(0.0,1.0),
    linewidth=2,
    markershape=symbols,
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2]],size=(width_px/2,width_px/3))
scatter!(DAR_mean4Hz_ITPCplot,driveamplituderatios,mean4Hz_ITPC_drive_strength_test2, markeralpha=0.5,
    legend=:false,
    label="",
    linewidth=2,
    markershape=symbols,
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2]],size=(width_px/2,width_px/3))
scatter!(DAR_mean4Hz_ITPCplot,driveamplituderatios,mean4Hz_ITPC_drive_strength_test3, markeralpha=0.5,
    legend=:false,
    label="",
    linewidth=2,
    markershape=symbols,
    markersize=3,
    color=[color_scheme[1] color_scheme[2] color_scheme[2]],size=(width_px/2,width_px/3))

plot!(DAR_mean4Hz_ITPCplot,driveamplituderatios,mean4Hz_ITPC_drive_strength,
    legend=:topright,
    xscale=:log10,
    xticks=(driveamplituderatios, driveamplituderatios),
    linewidth=2,
    markersize=3,
    markershape=[symbols[1] symbols[2] symbols[3]], # for the 3 models.
    label=["EEG" "NGNMM fast" "NGNMM slow"],
    color=[color_scheme[1] color_scheme[2] color_scheme[2]],
    linestyle=[linestyles[1]  linestyles[1] linestyles[2]], # for the 3 models.
    size=(width_px/2,width_px/3))

#save figure:
tosave=DAR_mean4Hz_ITPCplot
width, height = tosave.attr[:size]
Plots.prepare_output(tosave)
PlotlyJS.savefig(Plots.plotlyjs_syncplot(tosave),"./ppaper_drivestrength_mean4Hz_ITPC.pdf",width=Int64(ceil(width)), height=Int64(ceil(height)))
savefig(tosave,"./ppaper_drivestrength_mean_ITPC_fixedlabels.pdf")

