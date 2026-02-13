
using Pkg; Pkg.activate(".")
using NGNMM_NSP_paper_code
using Arrow, DataFrames, JLD2, Arrow
using Interpolations, ComponentArrays, Plots, LaTeXStrings, CSV


frequencies=collect(range(2,15.0,length=60)) #60 frequencies between 2 and 15 Hz
sine_path1="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_1.jld2"
sine_path2="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_2.jld2"
sine_path3="/Users/as15635/Documents/Projects/Shannon_et_al/Sine_Drives/drive_interpolators_3.jld2"

sine_drive_data1=jldopen(sine_path1,"r")
sine_drive_interpolators1=sine_drive_data1["drives"]
sine_drive_data2=jldopen(sine_path2,"r")
sine_drive_interpolators2=sine_drive_data2["drives"]
sine_drive_data3=jldopen(sine_path3,"r")
sine_drive_interpolators3=sine_drive_data3["drives"]
sine_drive_interpolators=vcat(sine_drive_interpolators1,sine_drive_interpolators2,sine_drive_interpolators3)



#load /Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/fast_NGNMM_response_to_abs_sine_stimulus_freq_range_2.0to15.0Hz_envelope_transformed.arrow.lz4
savepath="/Users/as15635/Documents/Projects/Shannon_et_al/Trajectories/" 
freq_range=(2.0,15.0) 
stimulus_type="envelope"

fourHz_NGNNM_abs_data=Matrix(DataFrame(Arrow.Table(savepath*"slow_NGNMM_response_to_abs_sine_stimulus_freq_range_2.0to15.0Hz_envelope_transformed.arrow.lz4")))
sine_sampling_rate=44100
time_range=(5.5,10.0)
plot(abs.(sine_drive_interpolators[1][time_range[1]*sine_sampling_rate:time_range[2]*sine_sampling_rate]))
plot!(fourHz_NGNNM_abs_data[1:end,1]./maximum(fourHz_NGNNM_abs_data[:,1]))
plot!(ylims=(0,10))

frequencies[38]
#check gaussian filtering:
filtered_4Hz=gaussian_filter_and_hilbert_for_PCM(fourHz_NGNNM_abs_data[:,1],44100,4.0)
plot!(angle.(filtered_4Hz))
plot!(ylims=(-1,2),xlims=(0,100000))

filtered_rect_sin=gaussian_filter_and_hilbert_for_PCM(abs.(sine_drive_interpolators[1][time_range[1]*sine_sampling_rate:time_range[2]*sine_sampling_rate]),44100,4.0) 

plot(angle.(filtered_4Hz))
plot!(angle.(filtered_rect_sin)) 
plot!(ylims=(-1,2),xlims=(0,100000))
phase_diffs=angle.(filtered_4Hz).-angle.(filtered_rect_sin)
scatter(exp.(im.*phase_diffs),xlims=(-1,1),ylims=(-1,1))

mean(exp.(im.*phase_diffs))