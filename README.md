# IGRF12-SGP4-Routine

This is a routine compiled from IGRF12 website and SGP4 routines from STK and various sources. 

## Routine to calculate magnetic field of a celestial object

This C++ files are able to calcualte magnetic field in the X, Y and Z directions with a given Lat, Long, Alt and Date. 
The error is within 20nT per year as given by the documentations.

## Routine to calculate position and velocity vectors of a celestial object

It also contains predictive model to calculate position and velocity vectors of a celestial object (satellite) with a given Two Line Element. 
The preditive model is within error of 150nT for estimation of 10 days with exponential increase of error day by day.

## Routine to calculate the sun vector of a celestial object

With the given TLE, Julian date can be easily calculated with the SGP4 routine. The Julian date can be passed to Sun vector routine in the code to calculate the Sun vector of the celestial object. 

## Conclusion

With the Sun vector and Magnetic Field vector, the Rotation matrix for attitude control and determination can be done for a certain celestial object. 

## Usage

Download the files and load it with an IDE with C++ compiler. I used it with Eclipse IDE. Load the files and read the comments on the routine functions inside the code to specific your need. 
