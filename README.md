PolySolvy-Version-x---Bivariate.py contains all the individual files used by Working_Example.py to interpolate and graph several examples of bivariate functions. 

Contents: 
  - vmond.py contains the Vandermonde() function along with the Gauss() function it calls; these are used to calculate the coefficient vector used by bivariate.py
  - bivariate.py combines the terms from a coefficient vector produced by Vandermonde() with their corresponding (x^n)(y^n) terms (this is dependent on the size of the data set)
  - graph_bivariate.py contains the graphing functionality of the project. It uses bivariate.py to generate appropriate z-values for a generated range of data points (including those which are being interpolated) and then graphs them using plot_trisurf() from Axes3D. It overlays the interpolating function with the original points that are being interpolated.
  - multivar_interp_generic.py is a file that was used to visualize the power matrix and terms of the bivariate function. It has no functional use in EXAMPLE.py
  - EXAMPLE.py the working example which contains all of the code files compiled along with a series of examples. The intention is that these be run in a jupyter notebook using the IPython kernel (this allows for the interactive graphing).
  
  How To Use EXAMPLE.py:
    - Copy all the code in EXAMPLE.py
    - Open web-browser tab and navigate to https://jupyter.org/try
    - Once there, click on "Try Classic Notebook" and then wait for the page to load
    - Once the page has loaded click on "insert" tab at the top menue, from there click on "Insert Cell Below".
    - Paste the entirety of the copied code from EXAMPLE.py into the "In [ ]" cell that appears
    - Hit "Run" on the top menue and enjoy.
    
    *NOTE: to try other examples enter them at the bottom of the cell (near the other examples) and then press "Run" again. Input must follow the same patterning as the other examples and contain a "square" amount of data points (i.e. 1, 4, 9, etc.). The limitation of using jupyter notebook is that it can't handle running more than 9 points currently. Further development for sets of 16+ is needed.
  
