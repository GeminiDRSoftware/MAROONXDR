Parameter defaults and options
------------------------------
::

   suffix               '_reduced'           Filename suffix
   target_name          None                 Target name to downselect files
   simbad_target_name   None                 SIMBAD resolvable target name
   use_coords           False                Use telescope pointing coordinates instead of target name
   zp_pc                0.0                  Zeropoint for counts_pc channel. Determined from data if not provided.
   zp_frd               0.0                  Zeropoint for counts_frd channel. Determined from data if not provided.
   start_time           'filename'           Time to consider to compute exposure start
      Allowed values:
      	mjd_start	Telescope MJD written at start of exposure
      	mjd_end	Telescope MJD written at end of readout
      	filename	UTC from filename
      
   report               True                 Generate diagnostic PDF of exposure meter time series.
