from astropy.io import fits 
import numpy as np 
import glob
import sys

class SOF():
	frame_dict = {}
	def __init__(self, folder):
		#for each .fits file in the folder, open it and get the CLASSIFICATION info from the header
		#save a dictionary with keys as CLASSIFICATION and entries as file paths
		files =  glob.glob(folder+'/*.fit')+ glob.glob(folder+'/*.fits') + glob.glob(folder+'/*.ofits') 
		self.parent_dir = folder
		for f in files:
			try:
				hdu = fits.open(f)
				classification = hdu[2].data['target'][0]
				if classification in self.frame_dict.keys():
					self.frame_dict[classification].append(f)
				else:
					self.frame_dict[classification] = [f]
			except KeyError as e:
				print "Keyword {!r} not found in extension {} of {}.".format('target', 2, f)
				print "Continuing, but SOF may be incomplete... \n"

	def add_file(self, fp, classification):
		#------------------------------------------------------------
		#---- Add a filepath to the dictionary key: classification --
		#---- If the key already exists, append to the list of fps --
		#---- otherwise, add the key: value pair --------------------
		#------------------------------------------------------------
		if classification in self.frame_dict.keys():
			self.frame_dict[classification].append(fp)
		else:
			self.frame_dict[classification] = [fp]

	def _write_sof_file(self, subdict, outname):
		#-----------------------------------------------------------
		#---- Write an sof file to the location outname ------------
		#---- It will contain all entries with the classifications -
		#---- within the list keys ---------------------------------
		#-----------------------------------------------------------
		try:
			f = open(outname, 'w')
			for k in subdict.keys:
				for entry in self.subdict[k]:
					f.write("{} \t {} \n".format(entry, k))
			f.close()
			print "File saved to {} successfully".format(outname)
			return outname
		except IOError:
			print "No such file or directory: {}".format(outname)
			print "Cannot continue. Exiting with code 1..."
			sys.exit(1)

	def mat_cal_det_sof(self):
		#---------------------------------------------------------
		#---- mat_cal_det requires the FLATFIELD, BADPIX, and ----
		#---- NONLINEARITY maps for the detector. ----------------
		#-----Create a subdictionary of only these files ---------
		#---------------------------------------------------------
		keys = ['FLAT', 'DARK']
		try:
			self.mat_cal_det_dict = {key:self.frame_dict[key] for key in keys}
		except KeyError as e:
			print "Key {} not in master dict. Did you enter the correct directory ({})?".format(e, self.parent_dir)

	def mat_est_flat_sof(self):
		#---------------------------------------------------------
		#---- mat_est_flat requires the OBSFLAT, OBSDARK,  -------
		#---- FLATFIELD, BADPIX, NONLINEARITY (from mat_cal_det)--
		#---------------------------------------------------------
		#-----Create a subdictionary of only these files ---------
		#---------------------------------------------------------
		keys = ['OBSFLAT', 'OBSDARK','FLATFIELD', 'BADPIX', 'NONLINEARITY']
		try:
			self.mat_est_flat_dict = {key:self.frame_dict[key] for key in keys}
		except KeyError as e:
			print "Key {} not in master dict. Did you enter the correct directory ({})?".format(e, self.parent_dir)

	def mat_est_shift_sof(self):
		#---------------------------------------------------------
		#---- mat_est_shift requires the SPECTRA_HOTDARK, --------
		#---- SPECTRA_IMAGES, DISTOR_HOTDARK, DISTOR_IMAGES, -----
		#---- BADPIX, NONLINEARITY, OBS_FLATFIELD (from previous -
		#---- calibration steps) ---------------------------------
		#---------------------------------------------------------
		#-----Create a subdictionary of only these files ---------
		#---------------------------------------------------------
		keys = ['SPECTRA_HOTDARK', 'SPECTRA_IMAGES', 'DISTOR_HOTDARK', 'DISTOR_IMAGES', 'BADPIX', 'NONLINEARITY', 'OBS_FLATFIELD' ]
		try:
			self.mat_est_shift_dict = {key:self.frame_dict[key] for key in keys}
		except KeyError as e:
			print "Key {} not in master dict. Did you enter the correct directory ({})?".format(e, self.parent_dir)

	def mat_est_kappa_sof(self):
		#------------------------------------------------------------------------
		#---- mat_est_kappa requires the KAPPA_HOTDARK, KAPPA_SRC (both new),---- 
		#---- BADPIX, NONLINEARITY, OBS_FLATFIELD, SHIFT_MAP (from prev. steps)--
		#------------------------------------------------------------------------
		#-----Create a subdictionary of only these files ------------------------
		#------------------------------------------------------------------------
		keys = ['KAPPA_HOTDARK', 'KAPPA_SRC', 'BADPIX', 'NONLINEARITY', 'OBS_FLATFIELD', 'SHIFT_MAP']
		try:
			self.mat_est_kappa_dict = {key:self.frame_dict[key] for key in keys}
		except KeyError as e:
			print "Key {} not in master dict. Did you enter the correct directory ({})?".format(e, self.parent_dir)

	def mat_raw_estimates_targ_sof(self):
		#------------------------------------------------------------------------
		#---- mat_raw_estimates for the TARGET_RAW ------------------------------
		#---- Needs also BADPIX, NONLINEARITY, OBS_FLATFIELD, SHIFT_MAP, --------
		#---- KAPPA_MATRIX from previous steps ----------------------------------
		#------------------------------------------------------------------------
		#-----Create a subdictionary of only these files ------------------------
		#------------------------------------------------------------------------
		keys = ['TARGET_RAW', 'BADPIX', 'NONLINEARITY', 'OBS_FLATFIELD', 'SHIFT_MAP', 'KAPPA_MATRIX']
		try:
			self.mat_raw_estimates_targ_dict = {key:self.frame_dict[key] for key in keys}
		except KeyError as e:
			print "Key {} not in master dict. Did you enter the correct directory ({})?".format(e, self.parent_dir)

	def mat_raw_estimates_cal_sof(self):
		#------------------------------------------------------------------------
		#---- mat_raw_estimates for the CALIB_RAW, HOT_DARK, CALIB_SRC_RAW ------
		#---- Needs also BADPIX, NONLINEARITY, OBS_FLATFIELD, SHIFT_MAP, --------
		#---- KAPPA_MATRIX from previous steps ----------------------------------
		#------------------------------------------------------------------------
		#-----Create a subdictionary of only these files ------------------------
		#------------------------------------------------------------------------
		keys = ['CALIB_RAW','HOT_DARK','CALIB_SRC_RAW', 'BADPIX', 'NONLINEARITY', 'OBS_FLATFIELD', 'SHIFT_MAP', 'KAPPA_MATRIX']
		try:
			self.mat_raw_estimates_cal_dict = {key:self.frame_dict[key] for key in keys}
		except KeyError as e:
			print "Key {} not in master dict. Did you enter the correct directory ({})?".format(e, self.parent_dir)


	def mat_cal_oifits_sof(self):
		#------------------------------------------------------------------------
		#---- mat_cal_oifits takes all the prevous calib steps to make the ------
		#---- TARGET_CAL_INT and INTERP_TF2 -------------------------------------
		#---- Needs TARGET_RAW_INT, CALIB_RAW_INT -------------------------------
		#------------------------------------------------------------------------
		#-----Create a subdictionary of only these files ------------------------
		#------------------------------------------------------------------------

		keys = ['TARGET_RAW_INT','CALIB_RAW_INT']
		try:
			self.mat_cal_oifits_dict = {key:self.frame_dict[key] for key in keys}
		except KeyError as e:
			print "Key {} not in master dict. Did you enter the correct directory ({})?".format(e, self.parent_dir)


if __name__ =="__main__":
	sof = SOF('/home/isbell/matisse_pl/') 


