using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using ComponentArrays, Plots, JSON, WAV, CSV, Statistics, DelimitedFiles, Random

#load all except 95.wav (it is shorter than the rest)
stimuli_filenames=readdir("Stimuli")[1:end-1] #ignoring 95.wav
individual_phoneme_envelopes=Vector{Vector}(undef, 20*length(stimuli_filenames)) #20 phonemes per sequence.
for (idx,filename) in enumerate(stimuli_filenames)
    X,fs,_,_ = wavread("./Stimuli/"*filename)
    envelope=sum(abs.(narrowband_envelopes(X,fs,80,8000,32)),dims=2)[:] #sum of hilbert transform of the signal passed through 32 band pass filters between 80Hz and 8000Hz (A Cochlear filter?)
    for i in 1:20
        phoneme_start_index=Int64(1+((i-1)*fs*5/20))
        phoneme_end_index=Int64(phoneme_start_index+(fs*5/20))
        if i==20
            individual_phoneme_envelopes[i+(idx-1)*20]=envelope[phoneme_start_index:end]
        else
            individual_phoneme_envelopes[i+(idx-1)*20]=envelope[phoneme_start_index:phoneme_end_index]
        end
    end
end

#cut any phonemes with non standard length
median(length.(individual_phoneme_envelopes))
#most common size, 11026
individual_phoneme_envelopes=filter(x->length(x)==11026,individual_phoneme_envelopes)
writedlm("IndividualPhonemes.csv",reduce(hcat,individual_phoneme_envelopes))

#draw from the individual phonemes data, and construct Noise sequences by stretching and squeezing the phoneme envelopes
#creating 5 second sequences at non-4Hz frequency.
individual_phoneme_envelopes=readdlm("IndividualPhonemes.csv",Float64)
# noise_rng=MersenneTwister(1234)
noise_rng=MersenneTwister(1234+3) #for noise set 3. (1 is 'old', 2 is 'new', 3 is '3')
NoiseDataSet=Vector{Vector}(undef,60)
for i in 1:60
    NoiseDataSet[i]=select_and_modify_20_random_phonemes([individual_phoneme_envelopes[:,i] for i in 1:size(individual_phoneme_envelopes,2)],noise_rng)
end
#save the 60 noise sequences. 
# writedlm("NoiseSequences60.csv",NoiseDataSet) #each row is a sequence.
writedlm("NoiseSequences60_3.csv",NoiseDataSet) #each row is a sequence.

#test loading 
# loaded_noise_sequences=readdlm("NoiseSequences60.csv",Float64)
#plot the first sequence
# plot(loaded_noise_sequences[1,:],legend=false)
