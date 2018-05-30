#!/bin/bash

base="/Users/sharon/code"
exe_dir="$base/solid"

make 

inc=90 # Observer's inclination (measured from spin axis)
mass="1.4" # in Solar Masses
radius="12" # Equatorial radius in km
temp=0.35 # in keV (Surface temperature)
distance=6.1713552e18 # 200 pc in meters

#Flags:
# -n = Number of azimuthal bins the circumference is cut up into
# -q = Shape model: q=3 means spherical; q=1 means oblate


#while test "$inc" -lt 91
#    do

# Spherical star; Quad = 0
    "$exe_dir/solid" -q 3 -m "$mass" -r "$radius" -f 1 -n 16 -i "$inc"  -o "sphere-$mass-$radius.txt" -T "$temp" -D "$distance" 



#    let inc=inc+5
#done


  
#
# -m mass in solar mass units
# -r radius (at equator) in km
# -f frequency in Hz
# -n number of bins
# -i inclination angle (degrees)
# -o output file
# -d data file

