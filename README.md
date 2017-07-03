# b_blocks
A simple Bayesian block decomposition utility.  
The algorithm is given in Scargle, Jeffrey D.; Norris, Jay P.; Jackson, Brad; Chiang, James
Studies in Astronomical Time Series Analysis. VI. Bayesian Block Representations,
http://adsabs.harvard.edu/abs/2013ApJ...764..167S and 
Jackson et al., An algorithm for optimal partitioning of data on an interval,
http://adsabs.harvard.edu/abs/2005ISPL...12..105J  

The code is partially based on battblocks utility
https://heasarc.gsfc.nasa.gov/ftools/caldb/help/battblocks.html

# Usage
b_blocks --PriorRatio=9.0 --FileName=input.thc --bbFileName=output.thi --LcType=Poisson

PriorRatio:
The higher the PriorRatio the less blocks are produced. 

## Poisson counts (LcType=Poisson)
input.thc format:  
Ti Tf Counts  
Ti and Tf - a time bin begin and end  
Count - number of counts in the bin  

input.thc example:
 2.919   5.863   488.010  
 5.863   8.807   451.890  
 8.807  11.751   453.890  
11.751  14.695   473.960  

## Gaussian distributed quantity (LcType=Gauss)
input.thi format:  
Ti Tf y yErr  
Ti and Tf - a time bin begin and end  
y - a quantity in the bin  
yErr - standard deviation of the quantity  

input.thi example:  
0.000 0.064  10.000  1.000  
0.064 0.128  10.000  1.000  
0.128 0.192  10.000  1.000  

## Output
output.thi (containing Bayesian blocks) format:  
Ti Tf CountRate CountRateErr  
Ti and Tf - a block begin and end  
CountRate - count rate in the block (for Poisson data) or averaged "x" (for Gaussian data)
CountRateErr - standard deviation of the count rate in the block, for Poisson data CountRateErr = sqrt(Rate/(Tf-Ti))

output.thi example:  
    2.416 14919.662  156.425    0.102  
14919.662 14931.438  170.200    3.802  
14931.438 56659.690  156.925    0.061  

# Command line arguments
-h, --help         print this help and exit  
-p, --PriorRatio   prior ratio  
-f, --FileName     input THC (LcType=Poisson) or THI (LcType=Gauss) file name  
-o, --bbFileName   output THI file name  
-t, --LcType       Gauss|Poisson  

# Installation
In the source code folder run  
cmake .  
make  

*gsl* and *gslcblas* are required.
