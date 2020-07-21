# tau_AS

Scripts for calculating tau_AS (raw, normalized, downsampled, and normalized + downsampled) from rMATS data for male-biased, female-biased, and unbiased genes.
Each script is run as follows:

> python Tau_##.py /path/to/SE.MATS.JC.txt /path/to/MXE.MATS.JC.txt readcoverage /path/to/coordinates

Coordinates should be in tab separated format like so:
'ENSAPLG00000021936 MB'
