## initialise and load packages
using Pkg;Pkg.activate(".")
using JSON, DelimitedFiles, Plots, Statistics, MAT
using NGNMM_NSP_paper_code: get_frequencies, get_sorting_indices, get_concordance_correlation_coefficient, figure_size_tuple
using LaTeXStrings, Plots
gr()
Plots.default(titlefontsize=16,legendfontsize=12,tickfontsize=9,guidefontsize=11)

##
#get the data:

frequencies=get_frequencies(4.5*10000,10000)[1:676]
_4Hz_index=indexin([4.0],frequencies)[1]
driveamplituderatios=[0.1,0.496,2.46,12.2,60.5,100.0]
driveamplituderatiosevoked=[12.2]
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

## here I want to plot the concordance coefficient between the condition ITPCs and the EEG itpcs, for each lambda and D and slow and fast model, instead of the ITPC per condition per model per lambda per D.
#this will be in two heatmaps. One for each model.

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

fourHzITPCS_DARfinaldim=stack(fourHzITPCS) #is 11x8x15x6 (9 models, 8 NSRs, 15 conditions, 6 DARs)  #indices 4:5,:,:,: are undef.
model_labels=["NMM slow" "NMM fast" "Phase Reset" "no Phase Reset" "evoked_derivative" "evoked_envelope" "evoked_peakenv" "evoked_peakrate" "EEG"]

CCCs=zeros(11,8,6) #still with 4:5 model indices as undefined.

for model_idx in 1:11 #loop over models.
    for nsr_idx in 1:8 #loop over NSRs
        for dar_idx in 1:6 #loop over DARs
            if model_idx==11 #do not need to re sort the sorted oana experiemental itpcs at index 11
                CCCs[model_idx,nsr_idx,dar_idx]=get_concordance_correlation_coefficient(fourHzITPCS_DARfinaldim[model_idx,nsr_idx,:,dar_idx],sorted_oana_ITPCs[:])
            else
                CCCs[model_idx,nsr_idx,dar_idx]=get_concordance_correlation_coefficient(fourHzITPCS_DARfinaldim[model_idx,nsr_idx,get_sorting_indices("./"),dar_idx],sorted_oana_ITPCs[:])
            end
        end
    end
end
CCCs


heatmaps=[]
j=1
clims=(0,1)
for i in vcat([1,2,3],[6,7,8,9,10,11])
    hm=heatmap(CCCs[i,:,:],title=model_labels[j],ylabel=L"$\lambda$",xlabel=L"D",clims=clims)
    plot!(hm,xticks=([1,2,3,4,5,6],driveamplituderatios),yticks=([1,2,3,4,5,6,7,8],noisestimratios))
    push!(heatmaps,hm)

    j+=1
end
two_column_size=figure_size_tuple(2, aspect_ratio=1)

plot(heatmaps...,layout=(3,3),size=two_column_size.*2,dpi=300)
savefig("all_concordance_corr_coeffs.pdf")



### #plot a heatmap of all the models on x axis vs the noise level on the y axis:

clims=(0,1)
noise_heatmap=heatmap([1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],CCCs[vcat([1,2,3],[6,7,8,9,10]),:,4]',clims=clims,xlabel="",ylabel=L"$\lambda$",cbar_title=L"$\rho$")
plot!(noise_heatmap,xticks=([1,2,3,4,5,6,7,8],model_labels[1:end-1]),xrotation=60)
plot!(noise_heatmap,yticks=([1,2,3,4,5,6,7,8],noisestimratios))
#retrieve maximum CCC for evoked models:
max_evoked_CCCs=maximum(CCCs[vcat([3]),:,4]',dims=1)
max_evoked_CCCs=(CCCs[vcat([3]),:,4]')

#plot just the NMM heatmaps against noise and drive strength:

NMM_heatmaps=[]
j=1
clims=(0,1)
for i in [1,2]
    # hm=heatmap(CCCs[i,:,:],title=model_labels[j],ylabel=L"$\lambda$",xlabel=L"D",clims=clims)
    hm=heatmap(CCCs[i,:,:],ylabel=L"$\lambda$",xlabel=L"D",clims=clims)
    plot!(hm,xticks=([1,2,3,4,5,6],driveamplituderatios),yticks=([1,2,3,4,5,6,7,8],noisestimratios))
    push!(NMM_heatmaps,hm)

    j+=1
end
plot(NMM_heatmaps...)
two_column_size=figure_size_tuple(2, aspect_ratio=3)

#extract maximum CCC for each NMM:
max_CCCs=[maximum(CCCs[i,:,:]) for i in [1,2]]

plot(vcat(NMM_heatmaps,noise_heatmap)..., layout=(1,3),size=two_column_size.*2, dpi=300,bottom_margin=20Plots.mm,cbar_title=L"$\rho$")

savefig("NMM_CCC_and_noiseheatmap.pdf")


