# Modelling-sharpness-specific-entrainment-to-speech-data-and-code
Code and data for our work modelling sharpness-specific tuning of neural entrainment to speech.

## Julia version and packages
This code uses Julia 1.10.3, and the following packages: Arrow v2.8.0, CSV v0.10.15, ColorSchemes v3.29.0, Colors v0.13.0, ComponentArrays v0.15.25, DSP v0.8.2, DataFrames v1.7.0, DelimitedFiles v1.9.1, DiffEqCallbacks v4.4.1, Distributions v0.25.118, FFTW v1.8.1, Interpolations v0.15.1, JLD2 v0.5.11, JSON v0.21.4, LaTeXStrings v1.4.0, MAT v0.10.7, OrdinaryDiffEq v6.92.0, Parameters v0.12.3, Plots v1.40.10, PowerLawNoise v0.0.1, WAV v1.2.0, LinearAlgebra, Random, Statistics v1.10.0. These are defined in the Project.toml file so you can create the julia project/environment using them.

## File structure:
All of the functions used to pre-process the stimuli, and run the tests are defined in ./src/NGNMM_NSP_paper_code.jl. This contains the module NGNMM_NSP_paper_code that is imported in the other scripts where necessary.

## Input Data
The raw audio stimuli are in ./Stimuli/ from which the noise is generated. The raw audio stimuli that are the 45 test streams (3 streams per consonant-vowel condition) are in ./StimuliNorm .

## Code
There are several "Generate_x.jl" scripts.  Generate_Noise.jl creates the 60 noise streams constructed by squeezing and stretching random phonemes from ./Stimuli . Generate_Phoneme_Envelopes.jl turns the raw audio streasm in. ./StimuliNorm into envelopes. Generate_Drives.jl creates the drives using the envelopes, noise, and chosen noise-stimulus ratios.

The "Compute_x.jl" scripts run the simulations for the various models tested here. "ControlCase" refers to the Phase-resetting oscillator model. "Compute_ITPCs_Across_Phonemes_randinitconds_1of.jl" runs the tests with the neural mass model. And "EvokedODEModel" is  the Evoked model. These scripts were setup to run on an HPC using SLURM, which is the reason for all the preamble that helped avoid precompilation errors when being run over many nodes at once.

There are two plotting/statistics scripts. "Plot_example_trajectories.jl" will produce supplementary figure S1, showing an example envelope stream, and resulting response from each model. "Plot_paper_results_figures_and_stats.jl" produces all of the other results figures in the paper, inclduing csv's of the statistics tables presented in the supplementary material. It uses the data stored in the ./Results/ folder of the main directory.

## Results Data
The results are divided into directories for each model in ./Results/. They are in JSON format with the following architecture: 
```
{"noisestimratio: NSR":{
                        "correlations":[cor4Hz, cor8Hz, cor12Hz], #from un-squared ITPC's
                        "ITPCs":{#the un-squared ITPCs (i.e. = R)
                                "f":[[ITPCs from f condition test 1], 
                                     [ITPCs from f condition test 2],
                                     [ITPCs from f condition test 3]], #ITPCs from 0-150Hz.
                                "other conditions....":[[],[],[]],...
                        },
                        "sortedITPCS_t2PD": [mean unsquared 4Hz ITPCs in order of increasing latency to peak derivative],
                                            [mean unsquared 8Hz ITPCs in order of increasing latency to peak derivative],
                                            [mean unsquared 12Hz ITPCs in order of increasing latency to peak derivative],[sorted latency to peak derivative],
                        "sqr_correlations": N/A - dont use.
                        }

}
```
The result filenames correspond to the simulation settings used. 
./Results/Statistics/ contains csv's for each of the statistics tables presented in the supplementary material.

## other files
mean_times_to_peak_deriv.csv contains the unsorted latencies to peak derivatives for the conditions in the following order: ["vowel","b","d","g","k","p","t","m","n","s","z","l","r","f","v"].
ITPC_peakder_condition.mat contains the original ITPCs and time to peak derivative data from Cucu et al's paper available at https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.826105/full.

NoiseSequences60_old.csv, NoiseSequences60.csv, and NoiseSequences60_3.csv contain 60 noise streams each generated with a different random seed. Each of these streams is actually only 5 seconds long, so the Generate_Drives.jl script, calling "NGNMM_NSP_paper_code.generate_drive_interpolators_specify_noise2stimulus_ratio_forsaving()", will combine 3, 5-second streams together for each test stream (5 seconds of noise before the stimulus, 5 seconds during, and 5 seconds after (not used)), so each NoiseSequences csv file provides 20 noise streams. The simulation running functions then run the simulations across each set of 20 drives, before computing the ITPC over the 60 trials per test.

IndividualPhonemes.csv contains the envelopes of all the individual phonemes present in ./Stimuli . It is created by Generate_Noise.jl. PhonemeEnvelopes_allconditions.json contains the stimuli envelope streams, constructed from the audio files in ./StimuliNorm, by Generate_Phoneme_Envelopes.jl.