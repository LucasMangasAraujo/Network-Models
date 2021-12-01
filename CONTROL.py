## CONTROL CODE TO SIMULATE THE BEHAVIOR OF A GIVEN RUBBER ELATICITY MODEL
""""Incompressebility is assumed in this code"""
#===================================================================
## IMPORT BUILT IN MODULES
from math import *
import numpy as np
import sys, os
import seaborn as sns
import pandas as p

#===================================================================
## IMPORT USER DEFINED FUNCTIONS
from UTILS import*
from INTERFACE import MATISU

#===================================================================
## VARIABLES DECLARATIONS

# Constant Paramters
R0=float(0.0); R1=float(1.0);

# Strings
MODEL = str(''); LOAD = str('') ; QUADSCHEME = str('')

# Integers
I = int(0); J = int(0); NINCR = int(R0); NPROP = int(R0); NSTRES = int(R0);

# Scalars
LAMBDA = R0; STRETCHINCR = R0;


# Second Order Tensors
T = np.zeros((3,3),float); TR = np.zeros((3,3),float);


# Unkowns a priori
RPROP = []; ETOT = []; STRETCH = [];

# Unkwowns Dictionaries
POLY = {}

#===================================================================
## GET MATERIAL PROPERTIES

# Set INPUT file name
INPUT = 'PARAMETERS.txt'

# Call MATPROP to assemble Material Model and loading conditions
(RPROPS,MODEL,LOADTYPE,NINCR,STRETCHINCR, QUADSCHEME, POLY) = MATPROP(INPUT)



#===================================================================
## Assemble loading histoty

# Initialize
STRETCH = np.zeros((NINCR+1,),float);

# Assemble
for I in range(0,NINCR+1):
	if (I == 0):
		STRETCH[I] = R1;
	else:
		STRETCH[I] = STRETCH[I-1]+STRETCHINCR;


#===================================================================
## INITIALIZE THE OUTPUT FILES

# Dummy output file, equal for all cases
DUMMYFILE = open('out.txt','w+')

# Based on the Input File, create the specific output file
EXT = '.txt'; # File Extension

if (LOADTYPE == 'UNIAXIAL'):
	if ('POLY' in MODEL):
		DIS = [*POLY][0] #Ditribution name
		OUT = 'UNI' + MODEL + '_' + DIS + EXT;
	else:
		OUT = 'UNI' + MODEL + EXT;
elif (LOADTYPE == 'BIAXIAL'):
	if 'POLY' in MODEL:
		OUT = 'BI' + MODEL + '_' + DIS + EXT;
	else:
		OUT = 'BI' + MODEL + EXT;
elif (LOADTYPE == 'SHEAR'):
	if 'POLY' in MODEL:
		OUT = 'SHEAR' + MODEL + EXT;
	else:
		OUT = 'SHEAR' + MODEL + '_' + DIS + EXT;
elif (LOADTYPE == 'COMPRESSUNI'):
	if 'POLY' in MODEL:
		DIS = [*POLY][0] #Ditribution name
		OUT = 'COMPRESSUNI' + MODEL + '_' + DIS + EXT;
	else:
		OUT = 'COMPRESSUNI' + MODEL + EXT;
elif (LOADTYPE == 'COMPSHEAR'):
	if 'POLY' in MODEL:
		DIS = [*POLY][0] #Ditribution name
		OUT = 'COMPSHEAR' + MODEL + '_' + DIS + EXT;
	else:
		OUT = 'COMPSHEAR' + MODEL + EXT;
		
#===================================================================
## STATE UPDATE

# Print the initial point
DUMMYFILE.write('{:2.4f} {:3.4f} {:3.4f}'.format(STRETCH[0],R0,R0))
DUMMYFILE.write('\n')

# Run Loop
for I in range(1,NINCR+1):
	# Call State Update functions
	if 'FULL' in MODEL:
		if 'POLY' in MODEL:
			TR,T = MATISU(MODEL,LOADTYPE,RPROPS,STRETCH[I],QUADSCHEME,POLY);
		else:
			TR,T = MATISU(MODEL,LOADTYPE,RPROPS,STRETCH[I],QUADSCHEME);
			
	else:
		TR,T = MATISU(MODEL,LOADTYPE,RPROPS,STRETCH[I]);
		
	# Print on output file
	DUMMYFILE.write('{:2.4f} {:3.4f} {:3.4f}'.format(STRETCH[I],TR[0,0],T[0,0]))
	DUMMYFILE.write('\n')
	#ENDFOR
# ENDIF		    

#===================================================================
## COPY SDUMMY FILE TO SPECIFIC OUTPUT FILE 

# Close file
DUMMYFILE.close()

# Copy and Move file to proper folder
os.system('copy out.txt %s' %OUT) # Copy file
if (MODEL == 'GAUSS'):
	os.system('move %s Gaussian' %OUT);
elif (MODEL == '3CHAIN'):
	os.system('move %s 3Chain' %OUT);
elif (MODEL == '8CHAIN'):
	os.system('move %s 8Chain' %OUT);
elif 'FULL' in MODEL:
	if 'POLY' in MODEL:
		os.system('move %s Polydisperse_Full_Network' %OUT);
	else:
		os.system('move %s Full_Network' %OUT);
		


#===================================================================
## Post Processing 

import matlab.engine
