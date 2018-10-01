#------- Script to run the esorex MATISSE Data Reduction Pipeline (DPL) ----------
#------- Author: Jacob Isbell (MPIA) --------------------------------------------
#------- v0.0.1 -- 27 Sept 2018 -------------------------------------------------


from subprocess import call
import argparse
import set_of_frames

#--------------------------------------------------------------------------------
#----------------- initialize the default filepaths for settings and recipes ----
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
parser = argparse.ArgumentParser(description='Wrapper to Run the MATISSE DRL and keep your sanity. Developed by Jacob Isbell (MPIA) --- v0.0.1 -- 27 Sept 2018')
parser.add_argument('-c','--config_file',dest='config_file', metavar='CONFIG_FILEPATH', type=str, default=config_loc, help='The filepath containing your esorex.rc file.')
#parser.add_argument('-r','--recipe_dir', dest='recipe_dir', metavar='RD', type=str, default=recipe_loc, help='Enter the filepath containing the directory containing your recipes.')
parser.add_argument('-d', '--data_dir', dest='data_dir', metavar='DATA_DIRECTORY', type=str, default=data_dir, help='The directory containing your raw data frames.')
args = parser.parse_args()


#Load the initial set of frames from a directory
init_sof = set_of_frames.SOF('.')


#test to make sure it is working
call('esorex --config=\"{}\" --recipes'.format(args.config_file), shell=True)


#--------------------------------------------------------------------------------
#----------------------start the actual calibration -----------------------------
#----- call the commands in the order specified by p28 of ----------------------- 
#----- the DRL Design Doc (Berio+2017) ------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------- The recipes are called with the syntax: ------------------------
#---- esorex [esorex-options] [myrecipe [recipe options] [sof]] -----------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#----------- the set of frames file(s) (sof(s)) are of the form -----------------
#---- /full/path/to/file.fits CLASSIFICATION ------------------------------------
#---- classification keywords are e.g. SCIENCE, MASTER_BIAS, DARK, ... ----------
#--------------------------------------------------------------------------------

"""
CURRENTLY SECTIONING OFF CODE FOR DEBUGGING PURPOSES
"""
do_mcd = True 
do_mef = False 
do_mes = False 
do_mek = False
do_mre_targ = False 
do_mre_cal = False
do_mco = False 


#--------------------------------------------------------------------------------
#----------------------------- mat_cal_det --------------------------------------
#---- Compute the FLATFIELD, BADPIX, and NONLINEARITY maps for the detector -----
#---- REQUIRES: FLAT, DARK  
#---- OUTPUT: FLATFIELD, BADPIX, NONLINEARITY
#---- FILES: mat_master_flat.fits, mat_master_bpm.fits, mat_master_nonlinearity.fits
#--------------------------------------------------------------------------------

if do_mcd:
	gain= 0.0               	# Default conversion gain in [e-/DU]. [0.0]
	nditskip = 0            	# number of skipped frames. [0]
	cosmics = 0            		# flag if cosmics should be detected. [0]
	darklimit = 100.0           # Absolute limit used for darks: good = [-limit ... +limit]. [100.0]
	flatlimit = 0.2          	# Relative limit used for flats: good = [ref*(1 - limit) ... ref*(1 + limit)]. [0.2]
	min_linear_range = 100.0    # Threshold for flatfield and non-linearity map calculation: LOW < threshold <= LINEAR. [100.0]
	max_linear_range = 1.5e4   	# Threshold for flatfield and non-linearity map calculation: LINEAR < threshold <= NON_LINEAR. [1.5e+04]
	max_nonlinear_range = 3e4 	# Threshold for flatfield and non-linearity map calculation: NON_LINEAR < threshold <= SATURATED. [3e+04]
	max_rel_deviation = 0.01  	# Maximum relative deviation for the nonlinearity fit. [0.01]
	max_abs_deviation =10.0  	# Maximum absolute deviation for the nonlinearity fit [e-]. [10.0]
	compensate = '[pb,ct]'      # Defines which kind of compensation should be applied (none = no compensation at all, all = all compensations possible, dd = detector specific defaults, pb = subtract pixel bias,
	# 								gb = subtract global bias, cb = subtract detector channel bias, rb = subtract row bias, ct = subtract crosstalk). [pb,ct]
	nlstart = 1024            	# first nonlinearity mapping sampling point. [1024]
	nlstep = 64             	# distance between two nonlinearity mapping sampling points. [64]
	poi = '[0,0]'               # pixel of interest: <x>,<y>. [0,0]
	expert = 0             		# expert flag. [0]
	nt = 'TRUE'                 # Flag if a new method should be used for detector calibration. [TRUE]
	#init_sof.mat_cal_det_sof()
	#mcd_sof = init_sof._write_sof_file(init_sof.mat_cal_det_dict, 'mat_cal_det.sof')
	mcd_sof = 'mat_cal_det.sof'

	call('esorex --config=\"{}\" mat_cal_det --gain={} --nditskip={} --cosmics={} --darklimit={} --flatlimit={} --min_linear_range={} --max_linear_range={}\
		 --max_nonlinear_range={} --max_rel_deviation={} --max_abs_deviation={} --compensate={} --nlstart={} --nlstep={} --poi={} --expert={} --nt={} {}'\
		.format(args.config_file, gain, nditskip, cosmics, darklimit, flatlimit, min_linear_range, max_linear_range, max_nonlinear_range, max_rel_deviation,\
			max_abs_deviation, compensate, nlstart, nlstep, poi, expert, nt, mcd_sof), shell=True)

	#init_sof.add_file('','FLATFIELD')
	#init_sof.add_file('','BADPIX')
	#init_sof.add_file('', 'NONLINEARITY')


#--------------------------------------------------------------------------------
#----------------------------- mat_est_flat -------------------------------------
#---- Compute the observed flatfield for later calibration ----------------------
#---- REQUIRES: OBSFLAT, OBSDARK, MASTERFLAT, MASTERBPM, MASTERNL  
#---- OUTPUT: OBS_FLATFIELD
#--------------------------------------------------------------------------------
if do_mef:
	obsflat_type = 'det'    # Defines which kind of OBS_FLATFIELD is created 
	#							(const = 1.0 for all pixels, det = detector flat, instr = instrument flat). [det]
	recalc_flat  = 0     	# Flag if the detector flatfield should be recalculated. [0]
	gain         = 0.0    	# Default conversion gain in [e-/DU]. [0.0]
	#init_sof.mat_est_flat_sof()
	#mef_sof = init_sof._write_sof_file(init_sof.mat_est_flat_dict, 'mat_est_flat.sof')
	mef_sof = 'mat_est_flat.sof'

	call('esorex --config=\"{}\" mat_est_flat --obsflat_type={} --recalc_flat={} --gain={} {}'\
		.format(args.config_file, obsflat_type, recalc_flat, gain, mef_sof), shell=True)

	#init_sof.add_file('', 'OBS_FLATFIELD')

#--------------------------------------------------------------------------------
#----------------------------- mat_est_shift ------------------------------------
#---- Compute the SHIFT MAP of the observation for later calibration ------------
#---- REQUIRES: SPECTRA_HOTDARK, SPECTRA_IMAGES, DISTOR_HOTDARK, DISTOR_IMAGES --
#---- OUTPUT: SHIFT_MAP
#--------------------------------------------------------------------------------
if do_mes:
	debug = 'FALSE'  #This parameter allows to printing debugging information. [FALSE]
	#init_sof.mat_est_shift_sof()
	#mes_sof = init_sof._write_sof_file(init_sof.mat_est_shift_dict, 'mat_est_shift.sof')
	mes_sof = 'mat_est_shift.sof'

	call('esorex --config=\"{}\" mat_est_shift --debug={} {}'\
		.format(args.config_file, debug, mes_sof), shell=True)

	#init_sof.add_file('', 'SHIFT_MAP')

#--------------------------------------------------------------------------------
#----------------------------- mat_est_kappa ------------------------------------
#---- Compute the KAPPA_MATRIX for the instrument -------------------------------
#---- REQUIRES: KAPPA_HOTDARK, KAPPA_SRC 
#---- OUTPUT: KAPPA_MATRIX
#--------------------------------------------------------------------------------
if do_mek:
	#no options
	init_sof.mat_est_kappa_sof()
	#mek_sof = init_sof._write_sof_file(init_sof.mat_est_kappa_dict, 'mat_est_kappa.sof')
	mek_sof = 'mat_est_kappa.sof'
	call('esorex --config=\"{}\" mat_est_shift {}'\
		.format(args.config_file, mek_sof), shell=True)

	#init_sof.add_file('','KAPPA_MATRIX')

#--------------------------------------------------------------------------------
#----------------------------- mat_raw_estimates --------------------------------
#---- Compute the raw intensities from the science frame ------------------------
#---- REQUIRES: TARGET_RAW
#---- OUTPUT: TARGET_RAW_INT
#--------------------------------------------------------------------------------
if do_mre_targ:
	compensate        = '[pb,cb,rb,nl,if,bp,od]' 	# Defines which kind of compensation should be applied (none = no compensation at all, all = all compensations possible, dd = detector specific defaults, pb = subtract pixel bias, gb =
	# 													subtract global bias, cb = subtract detector channel bias, rb = subtract row bias, ct = subtract crosstalk, nl = nonlinearity compensation, if = divide by instrument flat, df = divide by
	# 													detector flat, bp = bad pixel interpolation, el = convert to electrons, od = remove optical distortion). [pb,cb,rb,nl,if,bp,od]
	gain              = 0.0 						# Default conversion gain in [e-/DU]. [0.0]
	reduce_flag       = 'TRUE'						# Flag if the reference sub-windows should be removed from the result. [TRUE]
	ioi               = '[0,0]'						# images of interest: <first>,<count>. [0,0]
	useAvgSky         = 'FALSE'						# useAvgSky option. [FALSE]
	useKappaMatrix    = 'TRUE'						# useKappaMatrix option. [TRUE]
	replaceTel        = 0							# Replace Photometry of one telescope by the mean of the 3 others. (0: none, 1: AT1/UT1, 2: AT2/UT2, 3: AT3/UT3, 4: AT4/UT4). [0]
	useOpdMod         = 'FALSE'						# useOpdMod option. [FALSE]
	coherentIntegTime = 0.0							# Specify a coherent integration time (in s). [0.0]
	ChromaticOpdFit   = 'FALSE'						# Chromatic OPD Fit Option. [FALSE]
	corrFlux          = 'FALSE'						# corrFlux option. [FALSE]
	catalog           = 1							# calibrator catalog. [1]
	diamStar          = 1.0							# calibrator angular diameter. [1.0]
	diamErr           = 0.01						# calibrator angular diameter error. [0.01]
	cumulBlock        = 'FALSE'						# cumul all blocks of an OB. [FALSE]

	init_sof.mat_raw_estimates_targ_sof()
	#mre_targ_sof = init_sof._write_sof_file(init_sof.mat_raw_esimates_targ_dict, 'mat_raw_estimates_targ.sof')
	mre_targ_sof = 'mat_raw_estimates_targ.sof'

	call('esorex --config=\"{}\" mat_raw_estimates --gain={} --reduce_flag={} --ioi={} --useAvgSky={} --useKappaMatrix={} --replaceTel={} \
		--useOpdMod={} --coherentIntegTime={} --ChromaticOpdFit={} --corrFlux={} --catalog={} --diamStar={} --diamErr={} --cumulBlock={} {}'\
		.format(args.config_file, gain, reduce_flag, ioi, useAvgSky, useKappaMatrix, replaceTel,\
		 useOpdMod, coherentIntegTime, ChromaticOpdFit, corrFlux, catalog, diamStar, diamErr, cumulBlock, mre_targ_sof), shell=True)

	init_sof.add_file('', 'TARGET_RAW_INT')


#--------------------------------------------------------------------------------
#----------------------------- mat_raw_estimates --------------------------------
#---- Compute the raw intensities from the calib frame --------------------------
#---- REQUIRES: CALIB_RAW, HOT_DARK, CALIB_SRC_RAW
#---- OUTPUT: CALIB_RAW_INT
#--------------------------------------------------------------------------------
if do_mre_cal:
	compensate        = '[pb,cb,rb,nl,if,bp,od]' 	# Defines which kind of compensation should be applied (none = no compensation at all, all = all compensations possible, dd = detector specific defaults, pb = subtract pixel bias, gb =
	# 													subtract global bias, cb = subtract detector channel bias, rb = subtract row bias, ct = subtract crosstalk, nl = nonlinearity compensation, if = divide by instrument flat, df = divide by
	# 													detector flat, bp = bad pixel interpolation, el = convert to electrons, od = remove optical distortion). [pb,cb,rb,nl,if,bp,od]
	gain              = 0.0 						# Default conversion gain in [e-/DU]. [0.0]
	reduce_flag       = 'TRUE'						# Flag if the reference sub-windows should be removed from the result. [TRUE]
	ioi               = '[0,0]'						# images of interest: <first>,<count>. [0,0]
	useAvgSky         = 'FALSE'						# useAvgSky option. [FALSE]
	useKappaMatrix    = 'TRUE'						# useKappaMatrix option. [TRUE]
	replaceTel        = 0							# Replace Photometry of one telescope by the mean of the 3 others. (0: none, 1: AT1/UT1, 2: AT2/UT2, 3: AT3/UT3, 4: AT4/UT4). [0]
	useOpdMod         = 'FALSE'						# useOpdMod option. [FALSE]
	coherentIntegTime = 0.0							# Specify a coherent integration time (in s). [0.0]
	ChromaticOpdFit   = 'FALSE'						# Chromatic OPD Fit Option. [FALSE]
	corrFlux          = 'FALSE'						# corrFlux option. [FALSE]
	catalog           = 1							# calibrator catalog. [1]
	diamStar          = 1.0							# calibrator angular diameter. [1.0]
	diamErr           = 0.01						# calibrator angular diameter error. [0.01]
	cumulBlock        = 'FALSE'						# cumul all blocks of an OB. [FALSE]

	init_sof.mat_raw_estimates_cal_sof()
	#mre_targ_sof = init_sof._write_sof_file(init_sof.mat_raw_esimates_cal_dict, 'mat_raw_estimates_cal.sof')
	mre_cal_sof = 'mat_raw_estimates_cal.sof'

	call('esorex --config=\"{}\" mat_raw_estimates --gain={} --reduce_flag={} --ioi={} --useAvgSky={} --useKappaMatrix={} --replaceTel={} \
		--useOpdMod={} --coherentIntegTime={} --ChromaticOpdFit={} --corrFlux={} --catalog={} --diamStar={} --diamErr={} --cumulBlock={} {}'\
		.format(args.config_file, gain, reduce_flag, ioi, useAvgSky, useKappaMatrix, replaceTel,\
		 useOpdMod, coherentIntegTime, ChromaticOpdFit, corrFlux, catalog, diamStar, diamErr, cumulBlock, mre_targ_sof), shell=True)

	init_sof.add_file('', 'CALIB_RAW_INT')

#--------------------------------------------------------------------------------
#----------------------------- mat_cal_oifits -----------------------------------
#---- All calibrate steps are applied to the SCIENCE frame ----------------------
#---- REQUIRES: TARGET_RAW_INT, CALIB_RAW_INT
#---- OUTPUT: TARGET_CAL_INT, INTERP_TF2
#--------------------------------------------------------------------------------
if do_mco:
	tfKeep      = 0      		# store interpolate function. [0]
	tfInterp    = 0      		# transfer function interpolation method (0: average, 2:linear function). [0]
	cumulBlock  = 'FALSE'   	# cumul all blocks of an OB. [FALSE]

	init_sof.mat_cal_oifits_sof()
	#mco_sof = init_sof._write_sof_file(init_sof.mat_cal_oifits_dict, 'mat_cal_oifits.sof' )
	mco_sof = 'mat_cal_oifits.sof'
	call('esorex --config=\"{}\" mat_cal_oifits --tfKeep={} --tfInterp={} --cumulBlock={} {}'\
		.format(args.config_file, tfKeep, tfInterp, cumulBlock, mco_sof), shell=True)

	#init_sof.add_file('','TARGET_CAL_INT')
	#init_sof.add_file('','INTERP_TF2')


print "Calibration has been completed for the Visibilities and Phase information.\
		 To construct an image of the observation, please see 'run_imarec.py'"




#this will be a better casapy ;)








