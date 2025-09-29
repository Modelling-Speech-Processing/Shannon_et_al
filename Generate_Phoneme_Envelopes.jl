using Pkg; Pkg.activate("$(pwd())")
using NGNMM_NSP_paper_code
using ComponentArrays, Plots, JSON

Condition_keys=["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"]
stim_paths=readdir("./StimuliNorm/")
path_condition_pairs=[(Condition_keys[i],stim_paths[1+(3*(i-1)):3+(3*(i-1))]) for i in eachindex(Condition_keys)]
stim_path_dict=Dict{String,Vector{String}}(path_condition_pairs)
Phoneme_Envelopes=Dict{String,Vector{Vector{Float64}}}()
for condition in Condition_keys
    store_envelopes!(Phoneme_Envelopes,condition,stim_path_dict)
end

envelopedatastring=json(Phoneme_Envelopes)

open("PhonemeEnvelopes_allconditions.json", "w") do f
    write(f, envelopedatastring)
end

#check saved data
# Phoneme_Envelopes=JSON.parsefile("PhonemeEnvelopes_allconditions.json")
# plot(Phoneme_Envelopes["vowel"][1])
