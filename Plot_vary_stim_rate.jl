## initialise and load packages
using Pkg;Pkg.activate(".")
using JSON, DelimitedFiles, Plots, Statistics, MAT, DataFrames, CSV
using LaTeXStrings
using NGNMM_NSP_paper_code
gr()

##load the dataframes
name_extension="15secondstimulus_peaksearch"
tuned_pro_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/tuned_pro_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
untuned_pro_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/untuned_pro_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
evoked_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/evoked_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
slow_NGNMM_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/slow_NGNMM_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)
fast_NGNMM_stimrate_ITPCs_df=CSV.read("./Results/Speech_rate_test/fast_NGNMM_stimrate_sqrITPCs_$(name_extension).csv",DataFrame)

#if you want to access just the 4Hz ITPCs:
# name_extension="15secondstimulus"
# tuned_pro_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/tuned_pro_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
# untuned_pro_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/untuned_pro_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
# evoked_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/evoked_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
# slow_NGNMM_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/slow_NGNMM_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)
# fast_NGNMM_4Hz_ITPCs_df=CSV.read("./Results/Speech_rate_test/fast_NGNMM_fourHz_sqrITPCs_$(name_extension).csv",DataFrame)

#plot figure 6: 
Plots.default(titlefontsize=16,legendfontsize=8,tickfontsize=9,guidefontsize=11)
line_types=[:solid, :dash, :dot, :dashdot, :dashdotdot]
p_stimrateITPCs=plot(tuned_pro_stimrate_ITPCs_df[:,1], mean.(eachrow(tuned_pro_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(tuned_pro_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{tuned~PR}", xlabel=L"\textrm{syllabic~rate~(Hz)}", ylabel=L"\textrm{mean}~R^2", legend=:topleft,linestyle=line_types[1])
plot!(p_stimrateITPCs, untuned_pro_stimrate_ITPCs_df[:,1], mean.(eachrow(untuned_pro_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(untuned_pro_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{untuned~PR}", linestyle=line_types[5])
plot!(p_stimrateITPCs, evoked_stimrate_ITPCs_df[:,1], mean.(eachrow(evoked_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(evoked_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{evoked}", linestyle=line_types[3])
plot!(p_stimrateITPCs, slow_NGNMM_stimrate_ITPCs_df[:,1], mean.(eachrow(slow_NGNMM_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(slow_NGNMM_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{NMM~slow}", linestyle=line_types[4])
plot!(p_stimrateITPCs, fast_NGNMM_stimrate_ITPCs_df[:,1], mean.(eachrow(fast_NGNMM_stimrate_ITPCs_df[:,2:end])), ribbon=std.(eachrow(fast_NGNMM_stimrate_ITPCs_df[:,2:end])), label=L"\textrm{NMM~fast}", linestyle=line_types[2])   
plot!(p_stimrateITPCs,dpi=300, legend=:outertopright,xticks=(vcat(tuned_pro_stimrate_ITPCs_df[1:4:end,1],tuned_pro_stimrate_ITPCs_df[end,1]),vcat(tuned_pro_stimrate_ITPCs_df[1:4:end,1],tuned_pro_stimrate_ITPCs_df[end,1])))
tuned_pro_stimrate_ITPCs_df[:,1]

#plot in one column width figure:
one_column_size=figure_size_tuple(1, aspect_ratio=1.5)
plot!(p_stimrateITPCs,size=one_column_size,dpi=300)
savefig("15second_vary_syllabic_rate_sr_ITPCs_peaksearch.pdf")


