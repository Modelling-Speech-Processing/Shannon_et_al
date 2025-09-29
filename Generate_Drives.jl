using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using ComponentArrays, Plots, JSON, WAV, CSV, Statistics, DelimitedFiles, Random, DataFrames, Arrow
using ProgressLogging, JLD2



phoneme_sequence_envelopes=JSON.parsefile("PhonemeEnvelopes_allconditions.json")

noise_names=["NoiseSequences60.csv","NoiseSequences60_3.csv","NoiseSequences60_old.csv"]


noise_stimulus_ratios=[0.0,0.05,0.1,0.3,0.7,0.9,1.0]
Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]

#create the drive interpolators vector for each phoneme for each noise ratio and save it to a file for later use.

# @progress for noise_stimulus_ratio in noise_stimulus_ratios
for noise_name in noise_names
    @progress for phoneme in Condition_keys
        for test in 1:3
        stimulus_envelope=phoneme_sequence_envelopes[phoneme][test]
        drives=generate_drive_interpolators_specify_noise2stimulus_ratio_forsaving(stimulus_envelope,noise_stimulus_ratios[1],noise_name)
        noise_filename_short=split(noise_name,".")[1]
        jldsave("./Phoneme_Drives/drive_interpolators_$(phoneme)_test_$(test)_NSR_$(noise_stimulus_ratios[1])_$(noise_filename_short).jld2";drives)
        end
    end
end
# end

for noise_filename in noise_names
for noise_stimulus_ratio in noise_stimulus_ratios
    for phoneme in Condition_keys
        for test in 1:3
            noise_filename_short=split(noise_names[1],".")[1]
            if isfile("./Phoneme_Drives/drive_interpolators_$(phoneme)_test_$(test)_NSR_$(noise_stimulus_ratio)_$(noise_filename_short).jld2")
                println("File exists")
            else
                println("File does not exist")
            end
        end
    end
end
end
#all created. 125GB.

#check can open
xx=jldopen("./Phoneme_Drives/drive_interpolators_vowel_test_1_NSR_0.0_NoiseSequences60_new.jld2","r")
plot(xx["drives"][1])