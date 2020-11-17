# Generation of non stationary accelerograms
Function for the generation of non-stationary accelerograms adapted from original fortran program of Sabetta F. and Pugliese A. (1995)

## Input Parameters
Mw = moment magnitude

Re = epicentral distance (km)

isito = site conditions (0=rock; 1=shallow all.; 2=deep alluvium)

isig = st. dev. of GMPE (0=median value,1=84th percentile)

dt = time step of output accelerogram

nacc = number of accelerograms to be generated

tot_dur = total duration of accelerogram

scale = scale factor (1 for cm/s/s)

## Output Parameters
t = time vector (s)

acc = acceleration time history (in cm/s/s)

vel = velocity time history (in cm/s)

disp = displacement time history (in cm)

## Citation
Sabetta, Fabio, and Antonio Pugliese. "Estimation of response spectra and simulation of nonstationary earthquake ground motions." Bulletin of the Seismological Society of America 86.2 (1996): 337-352.
