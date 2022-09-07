# TV-MVP
We provide a detailed step-by-step implementation procedure for the paper "Time-varying Minimum Variance Portfolio" (Fan et al. 2022). We first introduce some main functions, including estimating time-varying factor loading, time-varying covariance matrix(also sparse residual covariance matrix estimation), rolling window procedure, and the two-step hypothesis test in the main paper. Further, we introduce some auxiliary functions, such as the tools for selecting values of various tuning parameters and functions used in the listed examples for each main function. Finally, we provide a comprehensive example to show the whole picture of our codes. All example codes could be found in the file "main_process.m" for users' convenience. The early version of the paper is available at http://dx.doi.org/10.2139/ssrn.3956956

# Preliminary setting
Before running provided code, some preliminary settings must be done to adapt to the local environment. We recommend the user to put the provided folder 'TV-MVP' under the folder 'bin' in MATLAB.
(1) Install CVX toolbox in MATLAB (the user needs first to download CVX toolbox named ``cvx.rar'', and then unzip the file, and run cvxsetup.m in MATLAB command window, more instructions of CVX toolbox can be referred to http://cvxr.com/cvx/)
(2) Install packages 'spcov', 'PDSCE' and 'R.matlab' locally in R.
(3) Open 'Rspcov.bat' in text form, edit the first path to the location of local "Rscript.exe", e.g. F:/R-3.6.1/bin/Rscript.exe, and the user need to replace F:/R-3.6.1/bin/ by their local path. Then, edit the second path to the location of 'TV-MVP' folder, e.g., F:/matlab/bin/TV-MVP/spcov_test.R, and the user needs to replace F:/matlab/bin/TV-MVP with their local path. We use 'bat' to run the code written by R and 'bat' file called by MATLAB.
(4) Open 'spcov\_test.R', edit the variable loc\_path2 in sixth row to your local path, e.g. F:/matlab/bin/TV-MVP/
(5) Open MATLAB and change the working path to the current folder, e.g., F:/matlab/bin/TV-MVP/

# Examples
Please refer to the provided document 'Manual.pdf' for contained functions and some examples.
