================================================================
Probability Density Estimator, Version 1.0, April 2019
Jenny Farmer jfarmer@carolina.rr.com
Donald Jacobs djacobs1@uncc.edu
University of North Carolina at Charlotte


================================================================
GENERAL INFORMATION
================================================================

The Probability Density Estimator is a non-parametric estimator based on the principle of Maximum Entropy.  Given a sample of independent and identically distributed continuous random variables, this program estimates the probability density function (PDF). For details of the method, please see the following publication:

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0196937

Please cite this publication if you use this code for your research.


=================================================================
INSTALLATION FOR MATLAB (2018 or greater)
=================================================================

The PDF Estimator consists of the following files:

EstimatePDF.cpp; EstimatePDF.h 
callPDF.cpp;  callPDF.h
ChebyShev.cpp; ChebyShev.h
InputData.cpp; InputData.h
InputParameters.cpp; InputParameters.h
MinimizeScore.cpp; MinimizeScore.h
Partition.cpp; Partition.h
Score.cpp; Score.h
WriteResults.cpp; WriteResults.h
OutputControl.cpp; OutputControl.h
example.m; CompilePDF.m; README.txt


Installation Steps


1. Prior to installing the PDF Estimator, the MingGW C/C++ compiler for Windows must be installed as a MATLAB Add-on.  To install, select [Add-Ons/Get Add-Ons] from the HOME menu within MATLAB and search for ‘MinGW’.  Select and install MinGW-w64.



2. Compile the source code to produce a matlab executable (mex) file by executing the CompilePDF.m script in MATLAB.


=================================================================
USAGE
=================================================================

Invoke the PDF Estimator in a MATLAB script, function, or from the MATLAB command prompt, as follows:

[failed, y, pdf, cdf, sqr, lagrange, score, confidence, SURD] =... 	EstimatePDF(data, parameters)


The 'data' input parameter is required, and should consist of a single vector of random sample data.  For example, to estimate the PDF from a Gaussian distribution sample with default parameters, execute the following command:

[failed, x, pdf, cdf, sqr, lagrange, score, confidence, SURD] =... 	EstimatePDF(randn(1000))

For additional information about the 'parameters' input argument, see the OPTIONAL INPUT PARAMETERS section.  An example MATLAB script, example.m, is also included with this package, demonstrating the usage of non-default parameters and suggestions for plotting results.


===============================================================
OUTPUT PARAMETERS
===============================================================

failed          non-zero if a solution was not found
y               range of values in PDF (independent variable)	
pdf             Probability Density Function (PDF)
cdf             Cummulative Denstiy Function
sqr             Scaled Quantile Residual
lagrange        Lagrange coefficients
score           Value returned by the score-type selected
confidence      SURD threshold achieved 
SURD            Sample Uniform Random Data            


For additional information see NOTES section.


=================================================================
OPTIONAL INPUT PARAMETERS
=================================================================

To customize the output and access advanced features, any of the following optional parameters can be included in structure and passed to the EstimatePDF function in the second argument (see USAGE).  Default values are in brackets. Only set the parameter values you wish to overide. 


parameters.SURDtarget          [40]      
parameters.SURDmin             [5]
parameters.SURDmax             [100] 
parameters.LagrangeMin         [1]
parameters.LagrangeMax         [200]
parameters.lowBound            [calculated]
parameters.highBound           [calculated]      
parameters.integrationPoints   [calculated]
parameters.debug               [false]
parameters.partition           [1025]         
parameters.fuzz                [false]     
parameters.scoreType           ['QZ']
parameters.minVariance         [true]
parameters.outlierCutoff       [7]
parameters.adaptiveDx          [false]

For additional information see NOTES section.


=================================================================
NOTES
=================================================================

The following section includes a few brief notes concerning more advanced input and output options available, and how they may affect performance of the estimation.  For a greater understanding of the methodology used, please see the publication referenced in the GENERAL INFORMATION section.


1. SURD Threshold

Sample Uniform Random Data (SURD) loosely correlates with the strength of the solution, with higher thresholds indicating more probably solutions for the PDF.  Note that the solution will not 
be considered a failure if it reaches at least a 5% SURD thereshold, which may be lower than the target specified with the parameters.SURDtarget option.  


2. Scaled Quantile Residual

The equation for Scaled Quantile Residual (SQR) is given by SQR = sqrt(N+2)*(u - uniform-u) where N is the number of data samples. SQR plots are very useful as a diagnostic measure because they are sample size invariant and have universal characteristics independent of the true PDF.


3. Lagrange Coefficients

Each Lagrange multiplier returned as output is an expansion coefficient in the series of orthogonal functions within an exponential. The more complex the shape of the distribution, the more Lagrange multipliers are required to accurately define the PDF. 


4. Greater accuracy in numerical integration can be controlled

Increasing the number of integration points will improve the resolution of the PDF, but could increase runtime.  Decreasing the integration points is not recommended, as it may produce poor solutions.  


5. Failed solutions 

Two safety measures are implemented to prevent the program from continuing an unreasonably long time without finding a solution. 

i) If progress stalls and the score is not improving significantly after many attempts, or 
ii) If the maximum number of Lagrange multipliers has been reached. 

If the maximum number of Lagrange multipliers is reached, this indicates that the solution is likely not yet converged. The user can increase the maximum. However, the default maximum of 200 is set to prevent cases that may never converge. 


6. Parametric maximum entropy method can be used with this program.

If a user desires an exact number of Lagrange multipliers, the  minLagrange and maxLagrange parameter options can be set to equal values. For example, if the user knows the distribution is a Gaussian, then the user could set both the minimum and maximum Lagrange mutlipliers to 3.  In this case, the output will be equivalent to a parametric maximum entropy method, where the number of Lagrange multipliers is known in advance. 


7. Verbose outputs for debugging. 

For more details on the progress of the program and explanations of possible warnings and outcomes, set the debug parameter option to true.


