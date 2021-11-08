## CONTROL CODE TO SIMULATE THE BEHAVIOR OF A GIVEN RUBBER ELATICITY MODEL
""""Incompressebility is assumed in this code"""
#===================================================================
## IMPORT BUILT IN MODULES
import numpy as np
import matplotlib.pylab as plt

#===================================================================
## IMPORT USER DEFINED FUNCTIONS
from PREFUN import*
from INTERFACE import MATISU


#===================================================================
## OPEN INPUT FILES
PROPFILE=open('PROP.txt','r')


#===================================================================
## VARIABLES DECLARATIONS

# Constant Paramters
R0=float(0.0); R1=float(1.0);

# Strings
MODEL=str(''); LOAD=str('') 

# Integers
I=int(0); J=int(0); NINCR=int(R0); NPROP=int(R0); NSTRES=int(R0);

# Scalars
LAMBDA=R0; 

# Vectors


# Second Order Tensors
DEFGRAD=np.zeros((3,1),float);


# Unkowns a priori
RPROP=[]; ETOT=[]; STRETCH=[];

#===================================================================
## GET MATERIAL PROPERTIES

# Ask which model will be uses
MODEL=input('Which Network Model?\n')
#MODEL='GAUSS' # Testing Code

# Check Model: Call MATPROP
RPROPS=MATPROP(PROPFILE,MODEL)

# Close File
PROPFILE.close()

#===================================================================
## ASSEMBLE LOADING CONDITIONS

# Ask Loading Conditions if Needed
LOADTYPE=input('Which Loading Condition?\n ')
#LOADTYPE='UNIAXIAL'; # Testing the code

# Ask Numbe of Increments
NINCR=int(eval(input('Numer of Increments:\n')))
#NINCR=100; # Default Value for Testing Code

# Targe Deformation
LAMBDA=int(eval(input('Target Stretch:\n ')));

# Call Loading Fuction

ETOT=LOADING(LOADTYPE, MODEL,NINCR,LAMBDA);

#===================================================================
## STATE UPDATE

# Assemble Nominal Stress array

if(LOADTYPE=='UNIAXIAL'):
    NSTRES=1; #
    NOMSTRES=R0;
    if(MODEL=='GAUSS'):
					OUT=open('UNIGAUSS.txt','w+')
    elif(MODEL=='3CHAIN'):
					OUT=open('UNI3CHAIN.txt','w+')
					
elif(LOADTYPE=='BIAXIAL'):
    NOMSTRESS=R0;
    if(MODEL=='GAUSS'):
					OUT=open('BIGAUSS.txt','w+')
elif(LOADTYPE=='SHEAR'):
    NOMSTRES=R0;
    if(MODEL=='GAUSS'):
					OUT=open('SHEARGAUSS.txt','w+')
					
    
    
# Run Loop
# Print the initial point
OUT.write('{:2.4f} {:3.4f}'.format(ETOT[0],R0))
OUT.write('\n')
STRETCH=R0
for I in range(1,NINCR+1):
	STRETCH=ETOT[I];
	NOMSTRES=MATISU(LOADTYPE,MODEL,RPROPS,STRETCH,NSTRES);
		
	# Print on output file
	OUT.write('{:2.4f} {:3.4f}'.format(STRETCH,NOMSTRES))
	OUT.write('\n')
	#ENDFOR
# ENDIF		    
#===================================================================
## Close File
OUT.close()

		

		



