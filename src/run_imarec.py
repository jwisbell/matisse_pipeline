#------- Script to run the esorex MATISSE Image Reconstruction ------------------
#------- Author: Jacob Isbell (MPIA) --------------------------------------------
#------- v0.0.1 -- 27 Sept 2018 -------------------------------------------------


from subprocess import call
import argparse
import set_of_frames

#--------------------------------------------------------------------------------
#------- initialize the default filepaths for settings and recipes --------------
#--------------------------------------------------------------------------------
config_loc = '/home/isbell/.esorex/esorex.rc'
recipe_dir = '/usr/local/misc/matisse/install/lib64/esopipes-plugins'
data_dir = '.'
#need to specify the input and output directories once I have data


#--------------------------------------------------------------------------------
#------------ Set up the argument parser for command line arguments -------------
#---- This really shouldn't take many options besides the config file -----------
#---- and the input and output directories for data -----------------------------
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Wrapper to Run the MATISSE DRL and keep your sanity. \n\n Developed by Jacob Isbell (MPIA) --- v0.0.1 -- 27 Sept 2018')
parser.add_argument('-c','--config_file',dest='config_file', metavar='CONFIG_FILEPATH', type=str, default=config_loc, help='The filepath containing your esorex.rc file.')
#parser.add_argument('-r','--recipe_dir', dest='recipe_dir', metavar='RD', type=str, default=recipe_loc, help='Enter the filepath containing the directory containing your recipes.')
parser.add_argument('-d', '--data_dir', dest='data_dir', metavar='DATA_DIRECTORY', type=str, default=data_dir, help='The directory containing your raw data frames.')
args = parser.parse_args()


#Load the initial set of frames from a directory
init_sof = set_of_frames.SOF(args.data_dir)



#--------------------------------------------------------------------------------
#----------------------------- mat_cal_imarec -----------------------------------
#---- Reconstructs an image based on calibrated interferometric measurements ----
#---- REQUIRES: TARGET_CAL_INT
#---- OUTPUT: TARGET_REC
#--------------------------------------------------------------------------------

    fov             = 40.0          # Field of view for the reconstructed image in [mas]. [40.0]
    npix            = 256           # Size of the reconstructed image in pixels. Powers of 2 should
#                                       be used (speeds up the FFT), but this is not enforced. [256]
    nbresult        = 0             # Number of reconstructions written to the result file. The
#                                       best result is always stored. If 0 is given, all created
#                                       reconstructions are stored in the result file. [0]
    lambda_from     = 1.0           # Shortest wavelength for the input data in [um]. [1.0]
    lambda_to       = 15.0          # Longest wavelength for the input data in [um]. [15.0]
    lambda_list     = 'none'        # A list of lambda ranges (pairs of lower and upper
#                                       wavelength). It is a sequence of comma separated floating
#                                       point numbers or 'none'. This list overwrites the
#                                       --lambda_from and --lambda_to parameters. [none]
    engine          = 1             # Specifies the optimization engine used for the image
#                                       reconstruction. 1 = ASA-CG, 2 = L-BFGS-B. [1]
    algo_mode       = 1             # Specifies if bispectrum and/or complex visibilities are used
#                                       for reconstruction. 1 = use bispectrum, 2 = use complex
#                                       visibilities, 3 = use bispectrum and complex visibilities.[1]
    calc_t3amp      = 1             # Flag if the T3 amplitude and error should be calculated. [1]
    calc_vis2f0     = 1             # Flag if an artificial squared visibility and error for f=0 should be calculated. [1]
    calc_visamp     = 1             # Flag if the VIS amplitude and error should be calculated. [1]
    calc_visf0      = 1             # Flag if an artificial complex visibility and error for f=0 should be calculated. [1]

    start_mode      = 0             # The mode for reading/creating the start image. 0 = read from
#                                       file, 1 = point source, 2 = gaussian disc, 3 = uniform disc,
#                                       4 = fully darkened disc, 5 = Lorentz disc. [0]
    start_param     = 0.0           # Additional parameter for the start image creation (mode=0 ->
#                                       scale [mas/px], mode=2 -> FWHM [mas], mode=3 -> diameter
#                                       [mas], mode=4 -> diameter [mas], mode=5 -> FWHM [mas]). [0.0]
    start_select    = 4             # Mode selection for the start image. [4]
    prior_mode      = 0             # The mode for reading/creating the prior image. 0 = read from
#                                       file, 1 = point source, 2 = gaussian disc, 3 = uniform disc,
#                                       4 = fully darkened disc, 5 = Lorentz disc. [0]
    prior_param     = 0.0           # Additional parameter for the prior image creation (mode=0 ->
#                                       scale [mas/px], mode=2 -> FWHM [mas], mode=3 -> diameter
#                                       [mas], mode=4 -> diameter [mas], mode=5 -> FWHM [mas]). [0.0]
#                           
    prior_select    = 4             # Mode selection for the prior image. [4]
    model_scale     = 0.1           # Pixel scale of the optional model image [mas/px]. [0.1]
    weight_power    = 1.0           # Weight power (uv weight calculation). [1.0]
    om_start        = 1.0           # Start radius of the object mask [mas]. [1.0]
    om_step         = 1.0           # Step size for the object mask radius scan [mas]. [1.0]
    om_count        = 6             # Number of object mask radius scans. [6]
    om_scale        = 0.1           # Pixel scale of the optional object mask image [mas/px]. [0.1]
    mu_start        = 1.0           # Start value for the regularization parameter mu. [1.0]
    mu_factor       = 0.1           # Factor between two consecutive regularization parameter values. [0.1]
    mu_count        = 6             # Number of regularization parameter scans. [6]
    reg_func        = 0             # Regularisation function (0 = no regularization). [0]
    reg_eps         = 0.0           # Epsilon for regularisation function 4 (edge preserving). [0.0]
    grad_tol        = 1e-11         # Tolerance value for ASA_CG. [1e-11]
    conv_scale      = 1.0           # Scale factor for the convolution (1.0 means max baseline is
#                                       used, negative value -> gaussian PSF, positive value -> tent
#                                       PSF). [1.0]
    cost_func       = 1             # Cost function (1 or 2). [1]
    cost_weight     = 0.0           # Weight for the cost function 2 (weight between closure phase nd modulus term). [0.0]
    ncorr           = 5             # Number of  corrections for L-BFGS-B. [5]
    factr           = 10.0          # L-BFGS-B tolerance for termination test. [10.0]
    pg_tol          = 1e-6          # L-BFGS-B projected gradient tolerance for termination test. [1e-06]
    asa_count       = 1             # Number of ASA-CG iterations per reconstruction. [1]
    cc_threshold    = 0.05          # Threshold for cross correlation. [0.05]
    mjd_tol         = 0.0001        # Maximum allowed MJD difference for finding a VIS2 element for  a T3 element [d]. [0.0001]
    bl_tol          = 0.05          # Maximum allowed baseline difference for finding a VIS2  element for a T3 element [m]. [0.05]
    wl_tol          = 0.0           # Maximum allowed wavelength difference wavelength filter [um] [0.0].
    precision       = 0             # Number of digits after decimal point for gradient and cost
#                                       value (precision < 0 : round relative, precision == 0 : no
#                                       round, precision > 0 : round absolute). [0]
    wiener_filter   = 0             # Flag if a Wiener filter is applied to the gradient. [0]
    filter_fwhm     = 0.0           # Start value for a gaussian gradient filter(FWHM) in [mas] [0.0].                       
    filter_factor   = 0.99          # Factor between two consecutive gradient filter sizes. [0.99]
    guess           = 0             # Flag if only a model fit is requested. [0]
    fit_fwhm        = 2.0           # Start FWHM for the model fit [mas]. [2.0]
    noise_seed      = 42            # Seed value for the noise random generator. [42]
    noise_factor    = 0.0           # Noise factor (noise_sigma = error*factor) for the noise
#                                       random generator. [0.0]
    info_flags      = 'param'       # Flags controlling the information printed during
#                                       reconstruction. [param]
    vis2_name       = ''            # ASCII file for measured and reconstructed squared
#                                       visibilities. []
    cp_name         = ''            # ASCII file for measured and reconstructed closure phases. []
    vis_name        = ''            # ASCII file for measured and reconstructed complex
#                                       visibilities. []

#make the set of frames
init_sof.mat_cal_imarec_sof()
#mci_sof = init_sof._write_sof_file(init_sof.mat_cal_imarec_dict, 'mat_cal_imarec.sof')
mci_sof = 'mat_cal_imarec.sof'


call('esorex --config=\"{}\" --fov={} --npix={} --nbresult={} --lambda_from={} --lambda_to={} --lambda_list={} --engine={} --algo_mode={} --calc_t3amp={} --calc_vis2f0={} \
    --calc_visamp={} --calc_visf0={} --start_mode={} --start_param={} --start_select={} --prior_mode={} --prior_param={} --prior_select={} --model_scale={} --weight_power={} \
    --om_start={} --om_step={} --om_count={} --om_scale={} --mu_start={} --mu_factor={} --mu_count={} --reg_func={} --reg_eps={} --grad_tol={} --conv_scale={} --cost_func={} \
    --cost_weight={} --ncorr={} --factr={} --pg_tol={} --asa_count={} --cc_threshold={} --mjd_tol={} --bl_tol={} --wl_tol={} --precision={} --wiener_filter={} filter_fwhm={}\
    --filter_factor={} --guess={} --fit_fwhm={} --noise_seed={} --noise_factor={} --info_flags=\"{}\" --vis2_name={} --cp_name={} --vis_name={} {}'.format(\
        fov, npix, nbresult, lambda_from, lambda_to, lambda_list, engine, algo_mode, calc_t3amp, calc_visf0,\
        calc_visamp, calc_visf0, start_mode, start_param, start_select, prior_mode, prior_param, prior_select, model_scale, weight_power, \
        om_start, om_count, om_scale, mu_start, mu_factor, mu_count, reg_func, reg_eps, grad_tol, conv_scale, cost_func,\
        cost_weight, ncorr, factr, pg_tol, asa_count, cc_threshold, mjd_tol, bl_tol, wl_tol, precision, wiener_filter, filter_fwhm,\
        filter_factor, guess, fit_fwhm, noise_seed, noise_factor, info_flags, vis2_name, cp_name,vis_name, imarec_sof), shell=True)