using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using ComponentArrays, Plots, JSON, WAV, CSV, Statistics, DelimitedFiles, Random, DataFrames, Arrow
using ProgressLogging, JLD2

# phoneme_sequence_envelopes=JSON.parsefile("PhonemeEnvelopes_allconditions.json")
# noise_names=["NoiseSequences60.csv","NoiseSequences60_3.csv","NoiseSequences60_old.csv"]
# noise_stimulus_ratios=[0.0,0.05,0.1,0.3,0.7,0.9,1.0]

frequency_range=(2.0,15.0) #Hz

sampling_rate=44100.0 #Hz


#create the drive interpolators vector for each phoneme for each noise ratio and save it to a file for later use.
#check for directory to store them, make one if not present:
if !isdir("./Sine_Drives/")
    mkdir("./Sine_Drives/")
end

# @progress for noise_stimulus_ratio in noise_stimulus_ratios
for test in 1:3
    drives=generate_sine_drive_interpolators_for_saving(sampling_rate, frequency_range,test)
    # jldsave("./Sine_Drives/drive_interpolators_$(test).jld2";drives)
    jldsave("./Sine_Drives/$(sampling_rate)sr_wide_drive_interpolators_$(test).jld2";drives)
end

for test in 1:3
    if isfile("./Sine_Drives/drive_interpolators_$(test).jld2")
        println("File exists")
    else
        println("File does not exist")
    end
end

#check it looks ok
xx=jldopen("./Sine_Drives/drive_interpolators_1.jld2","r")
plot(range(0.0,10.0,length=Int(10*sampling_rate)), xx["drives"][1])
xx["drives"][1]
length(range(0.0,10.0,length=Int(10*sampling_rate)))