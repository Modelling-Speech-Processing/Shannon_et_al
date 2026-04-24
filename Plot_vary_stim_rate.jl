## initialise and load packages
using Pkg;Pkg.activate(".")
using JSON, DelimitedFiles, Plots, Statistics, MAT, DataFrames, CSV
using LaTeXStrings
using NGNMM_NSP_paper_code
gr()

## This is how the data was saved:

# save_as_df_with_headings(tuned_pro_stimrate_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/tuned_pro_stimrate_sqrITPCs.csv")
# save_as_df_with_headings(evoked_stimrate_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/evoked_stimrate_sqrITPCs.csv")
# save_as_df_with_headings(stimrate_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/slow_NGNMM_stimrate_sqrITPCs.csv")
# save_as_df_with_headings(fast_NGNMM_stimrate_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/fast_NGNMM_stimrate_sqrITPCs.csv")

# #4 Hz ITPCs too:
# save_as_df_with_headings(untuned_pro_fourHz_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/untuned_pro_fourHz_sqrITPCs.csv")
# save_as_df_with_headings(tuned_pro_fourHz_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/tuned_pro_fourHz_sqrITPCs.csv")  
# save_as_df_with_headings(evoked_fourHz_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/evoked_fourHz_sqrITPCs.csv")
# save_as_df_with_headings(fourHz_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/slow_NGNMM_fourHz_sqrITPCs.csv")
# save_as_df_with_headings(fast_NGNMM_fourHz_ITPCs, Condition_keys, stim_rates, "./Results/Speech_rate_test/fast_NGNMM_fourHz_sqrITPCs.csv")


## I want to plot the mean ITPC over conditions (the column headings of the dataframes) for each model vs stimulus rate (the first column of the data frames).
#with the stdev over conditions as a ribbon. Using plot (i.e a line graph, with all models on the same plot)

##load the dataframes
name_extension="15secondstimulus_peaksearch"
tuned_pro_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/tuned_pro_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
untuned_pro_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/untuned_pro_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
evoked_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/evoked_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
slow_NGNMM_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/slow_NGNMM_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
fast_NGNMM_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/fast_NGNMM_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)

tuned_pro_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/tuned_pro_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
untuned_pro_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/untuned_pro_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
evoked_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/evoked_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
slow_NGNMM_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/slow_NGNMM_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
fast_NGNMM_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/fast_NGNMM_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)

#plot the data on a plot for stimrate and a plot for 4Hz itpcs, each model is a series.
Plots.default(titlefontsize=16,legendfontsize=8,tickfontsize=9,guidefontsize=11)

p_stimrateITPCs=plot(tuned_pro_stimrate_ITPCs_df[:,1], mean.(eachrow(tuned_pro_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(tuned_pro_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{tuned~PR}", xlabel=L"\textrm{syllabic~rate~(Hz)}", ylabel=L"\textrm{mean}~R^2", legend=:topleft)
plot!(p_stimrateITPCs, untuned_pro_stimrate_ITPCs_df[:,1], mean.(eachrow(untuned_pro_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(untuned_pro_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{untuned~PR}")
plot!(p_stimrateITPCs, evoked_stimrate_ITPCs_df[:,1], mean.(eachrow(evoked_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(evoked_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{evoked}")
plot!(p_stimrateITPCs, slow_NGNMM_stimrate_ITPCs_df[:,1], mean.(eachrow(slow_NGNMM_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(slow_NGNMM_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{NMM~slow}")
plot!(p_stimrateITPCs, fast_NGNMM_stimrate_ITPCs_df[:,1], mean.(eachrow(fast_NGNMM_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(fast_NGNMM_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{NMM~fast}")   
plot!(p_stimrateITPCs,dpi=300, legend=:outertopright,xticks=(vcat(tuned_pro_stimrate_ITPCs_df[1:4:end,1],tuned_pro_stimrate_ITPCs_df[end,1]),vcat(tuned_pro_stimrate_ITPCs_df[1:4:end,1],tuned_pro_stimrate_ITPCs_df[end,1])))
tuned_pro_stimrate_ITPCs_df[:,1]

#plot in one column width figure:
one_column_size=figure_size_tuple(1, aspect_ratio=1.5)
plot!(p_stimrateITPCs,size=one_column_size)
savefig("15second_vary_syllabic_rate_sr_ITPCs_peaksearch.pdf")
#format nicely:










# #plot the 4Hz itpcs: DOES NOT WORK ANYMORE, WITH LONGER TIMESPAN< THE 4Hz freq idx changed.
# #so the data from the 15secondstimulus run has ITPCs at about 1.25Hz here, so all very low.
# p_4HzITPCs=plot(tuned_pro_4Hz_ITPCs_df[:,1], mean.(eachrow(tuned_pro_4Hz_ITPCs_df[:,2:end])), ribbon=std.(eachrow(tuned_pro_4Hz_ITPCs_df[:,2:end])), label="Tuned Pro", xlabel="Stimulus Rate (Hz)", ylabel="ITPC", title="4Hz ITPC vs Stimulus Rate", legend=:topleft)
# plot!(p_4HzITPCs, untuned_pro_4Hz_ITPCs_df[:,1], mean.(eachrow(untuned_pro_4Hz_ITPCs_df[:,2:end])), ribbon=std.(eachrow(untuned_pro_4Hz_ITPCs_df[:,2:end])), label="Untuned Pro")
# plot!(p_4HzITPCs, evoked_4Hz_ITPCs_df[:,1], mean.(eachrow(evoked_4Hz_ITPCs_df[:,2:end])), ribbon=std.(eachrow(evoked_4Hz_ITPCs_df[:,2:end])), label="Evoked")
# plot!(p_4HzITPCs, slow_NGNMM_4Hz_ITPCs_df[:,1], mean.(eachrow(slow_NGNMM_4Hz_ITPCs_df[:,2:end])), ribbon=std.(eachrow(slow_NGNMM_4Hz_ITPCs_df[:,2:end])), label="Slow NGNMM")
# plot!(p_4HzITPCs, fast_NGNMM_4Hz_ITPCs_df[:,1], mean.(eachrow(fast_NGNMM_4Hz_ITPCs_df[:,2:end])), ribbon=std.(eachrow(fast_NGNMM_4Hz_ITPCs_df[:,2:end])), label="Fast NGNMM")   
# plot!(dpi=300, legend=:outertopright)


# #bar chart of 4Hz itpcs per condition for the fast NGNMM
# bar(fast_NGNMM_4Hz_ITPCs_df[:,1],Matrix(fast_NGNMM_4Hz_ITPCs_df[:,2:end]), label="Fast NGNMM", xlabel="Stimulus Rate (Hz)", ylabel="4Hz ITPC", title="4Hz ITPC vs Stimulus Rate for Fast NGNMM", legend=false, dpi=300,ylims=(0.0,1.0))