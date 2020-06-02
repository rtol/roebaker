# roebaker
Matlab scripts and functions to fit the Roe-Baker distribution to estimates of the climate sensitivity.

Main script is FitRoeBaker.m

For each of the papers surveyed by Knutti, this script fits either a Normal or a Roe-Baker distribution. The Normal ones are fitted directly in the script by the method of moments. The Roe-Baker distributions are fitted by calibrating the central tendency and lower and upper bound. As Knutti reports statistics in nine different formats, there are nine auxiliary functions.

There are also auxiliary functions to return the characteristics of the Roe-Baker distribution and its likelihood. For symmetry, the same functions are there for the Normal distribution.

Necessary data are in Knutti.mat. If you want to trace that data, look at the Excel file.
