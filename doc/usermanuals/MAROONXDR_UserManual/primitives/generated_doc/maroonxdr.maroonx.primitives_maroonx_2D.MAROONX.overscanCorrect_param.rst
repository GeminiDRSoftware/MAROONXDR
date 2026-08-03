Parameter defaults and options
------------------------------
::

   suffix               '_overscanCorrected' Filename suffix
   function             'none'               Fitting function
      Allowed values:
      	none	Row-by-row values
      	spline3	Cubic spline
      	chebyshev	Chebyshev polynomial
      
   order                None                 Order of fitting function
      	Valid Range = [0,inf)
   lsigma               3.0                  Low rejection in sigma of fit
      	Valid Range = [0,inf)
   hsigma               3.0                  High rejection in sigma of fit
      	Valid Range = [0,inf)
   niter                0                    Maximum number of rejection iterations
      	Valid Range = [0,inf)
   grow                 0                    Rejection growing radius
      	Valid Range = [0,inf)
   nbiascontam          0                    Number of columns to exclude from averaging
      	Valid Range = [0,inf)
   bias_type            None                 Overscan type
      Allowed values:
      	serial	Serial
      	parallel	Parallel
      	None	Field is optional
      
