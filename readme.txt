This folder contains the codes I wrote for covariance matrix estimation, in the framework of the collaboration with A. Musolas.

The main files are :
	- run_opti_bezier_quotient.m, 
	- run_opti_bezier_section.m,
	- run_opti_piecewise_bilinear.m

Those three files load the data, and perform the tests in order to recover the results from the paper. 

--------------------------------------------------------------------------------------------------
IMPORTANT WARNING : the data are not contained in this git repo, you should download them first. 
Here is how to proceed for the wind field data:
- create a subfolder data_points/data_surface_extracted
- download the data, and store them in this folder, in the following form. Each data is stored in a separate .mat file. The file should contain directly the Y factorization of the data, i.e., if the true covariance matrix is C = YY', we store directly the Y factor. This Y factor is expected to be a matrix of size 3024 x 3024, whose columns have decreasing magnitude, so that to truncate the rank to r, we just have to keep the r first columns. 
- rename the mat file as Yuv_xyz.mat. For example, the file containing the data associated with a heading = 1.5pi/32 and a magnitude 1 will be named Y15_010.mat, while the one associated with heading = 4 pi /32 and magnitude 13 will be Y40_130.mat (so, basically, we keep one number to characterize the decimal part of the parameter)
--------------------------------------------------------------------------------------------------


--------------------------------------------------------------------------------------------------
SECOND IMPORTANT WARNING : Notice also that the three above mentioned codes take some time to run (between 5 and 8 hours on a Windows 7 platform with 8 cores at 3.60 GHz and 16 GB ram). Not all parts of the code are required to solve the optimization problem on the surface, and some parts of the code can be commented in order to reduce further computation time. See the documentation inside the three above mentioned codes for more information.
--------------------------------------------------------------------------------------------------

Author: E. Massart

Last update : October 24th, 2018




