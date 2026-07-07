Parameter defaults and options
------------------------------
::

   suffix               '_structureStandardized' Filename suffix
   bad_wcs              'exit'               Method for WCS handling
      Allowed values:
      	exit	Exit reduction if discrepant WCS found
      	fix	Attempt to fix discrepant WCS using offsets
      	new	Create new WCS from target coordinates and offsets
      	ignore	Do not check or fix the WCS
      
