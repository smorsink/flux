# flux
Computes the flux from a uniformly emitting rotating neutron star 

Compile: Type make at the command line. This will create an executable "flux"

To Run, edit the bash script "go". You will need to change the name of your "base" directory.

The script go runs as its default a star with the following parameters:
M = 1.4 Msun; R=12 km; Spin Frequency = 1 Hz (Essentially non-rotating)
T = 0.35 keV (in the star's rest frame)

The output goes to the file specified in the bash script. The output looks like:

#M = 1.4 Msun #Req = 12 km
#spin      Flux(1keV)      BolFlux/Fs      v_{eq}/c             Area             A/4piR^2   SolidAng         SolidAng/Ss
1          17.253          1.0016          8.46615e-05          1809.56          1          1.81543e-29          1.0016

Flux(1keV) measured at 1keV (units to be discussed later)
BolFlux/Fs = The bolometric flux divided by the bolometric flux for a spherical star. Should be 1 for a non rotating star.
v_{eq}/c = velocity at the equator divided by the speed of light
Area = Surface area of the star (in km^2)
A/4piR^2 = Surface area divided by surface area of a spherical star with the same radius. Should be 1 for a non-rotating star.
SolidAng = Solid angle subtended by the star, in steradians
SoldAng/Ss = Solid angle divided by the value for a spherical star. Should be 1 for a non-rotating star.
