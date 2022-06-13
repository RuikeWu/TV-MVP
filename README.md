# TV-MVP
We provide detailed step by step implementation procedure for the paper "Time-varying Minimum Variance Portfolio" (Fan et al. 2022).  We first introduce some main functions including estimating time-varying factor loading, time-varying covariance matrix(also sparse residual covariance matrix estimation), rolling window procedure, and the two-step hypothesis test in the main paper.	Further, we introduce some auxiliary functions, such as the tools for selecting values of various tuning parameters and some functions used in listed examples. Finally, we provide a large example to show a whole picture of our codes. All example codes could be found in the file "main_process.m" for users' convenience. The first version of paper is available at http://dx.doi.org/10.2139/ssrn.3956956

# Preliminary setting
Before running provided code, some preliminary setting need to be done for adapting to the local environment. We recommend that you put the provided folder 'TV-MVP' under the folder 'bin' in MATLAB.<br>
(1)Install cvx in MATLAB(we provide the CVX toolbox named "cvx.rar", the user only need to unzip the file, and run "cvx_setup.m" in MATLAB command window, more instruction of CVX toolbox can be referred to http://cvxr.com/cvx/)<br>
(2)Install packages 'spcov', 'PDSCE' and 'R.matlab' locally in R.<br>
(3)Open 'Rspcov.bat' in text form, edit the first path to the location of local "Rscript.exe", e.g. F:/R-3.6.1/bin/Rscript.exe, and the user need to replace F:/R-3.6.1/bin/ by their local path. And further edit the second path to the location of  'TV-MVP' folder, e.g. F:/matlab/bin/TV-MVP/spcov_test.R, and the user need to replace F:/matlab/bin/TV-MVP by their local path. We use 'bat' to run the code written by R and 'bat' file called by MATLAB.<br> 
(4)Open MATLAB, and change the working path to current folder, e.g. F:/matlab/bin/TV-MVP/

# Example
For contained functions and some examples, please refer to the provided document 'Manual.pdf'.
