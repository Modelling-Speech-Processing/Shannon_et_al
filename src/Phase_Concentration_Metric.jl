


### reformed drive interpolators generator to create sine waves of varying frequency for PCM tests.
#Will use the same global interpolators object type as before.
global sine_interpolators_global::Vector{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Linear{Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Linear{Throw{OnGrid}}}, Tuple{StepRange{Int64, Int64}}}, BSpline{Linear{Throw{OnGrid}}}, Line{Nothing}}}

"""
Will make an "interpolators" array that contains interpolations of 20 sine waves of different frequency and put it in global variable Interpolators.
This can then be used in the ensemble simulation via the interpolation selector variable in a prob_func

before the stimulus is applied, noise is applied of amplitude scale 1.0 (which will be multiplied by the drive_amplitude parameter of the ODE model later)
the noisestimratio parameter sets the ratio of noise to stimulus amplitude to use when the stimulus is applied. for instance, if a ratio of 0.4 noise to 0.6 stimulus is required
just set noisestimratio to 0.4. This results in a constant amplitude drive being applied to the model of format (1.0*noise, (0.4*noise + 0.6*stimulus), 1.0*noise).
"""
function generate_drive_interpolators_specify_noise2stimulus_ratio(stimulus_envelope_sum,noisestimratio,noise_filename)
    stimlength=length(stimulus_envelope_sum)
    noise=readdlm("./$(noise_filename)",',')
    if typeof(noise[1])==SubString{String}
        noise=readdlm("./$(noise_filename)")
    end
    all_noise_vector=[noise[i,:] for i in 1:size(noise,1)]
    noise_length=5*44100
    out_of_stim_noise_scale=1.0
    during_stim_noise_scale=noisestimratio #scales size of amplitude for the during-stimulus noise.

    for trial in 1:20   
        #will go through the noise samples in triplets: 1 is the before stimulus noise, 1 is during stimulus and 1 is after stimulus
        #will loop over this, and create a single interpolator for each set of noise samples, with the stimulus added to the middle and appropriate scale applied
        #drive is [noise,stimulus,noise,zeroes] (in case it runs over, making extrapolation to be 0 only)
        trial_drive_input=[all_noise_vector[1+(trial-1)*3]*out_of_stim_noise_scale;(stimulus_envelope_sum.*(1-during_stim_noise_scale)).+(all_noise_vector[2+(trial-1)*3][noise_length-stimlength+1:end]*during_stim_noise_scale);all_noise_vector[3+(trial-1)*3]*out_of_stim_noise_scale;zeros(length(all_noise_vector[1]))]
        
        xs=1:1:length(trial_drive_input)
       
        interpolators[][trial]=linear_interpolation(xs,Vector(trial_drive_input),extrapolation_bc=Line())
        #interpolators[trial]=interpolate(Vector(trial_drive_input),BSpline(Linear()))#,extrapolation_bc=Line())
    end
    return interpolators[] #look into typed global
end
