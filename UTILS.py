## FUNCTIONS USED TO ASSEMBLE ARRAYS, LOADING PATHS AND ETC (PRE-PROCESSING IF YOU WILL)

#-----------------------------------------------------------------------------------
# Import relavant modules
import numpy as np
from math import *
import sys, os

#======================================================================================
''' MATPROP: THIS FUNCTION ASSEMBLE ARRAY OR SCALAR CARRYING THE MATERIAL PROPERTIES

Inputs: 
 		INFILE: TXT FILE CONTAINING THE MATERIAL PROPERTIES		

Outputs:
		RPROP: Array or scalar containing the material properties
		MODEL: STRING CONTANING THE MODEL IDENTIFIER
		LOADTYPE: String with the loading condition
		NINCR: Numer of increments
		STRETCHINCR: Delta stretch
		QUADSCHEME: Particular integration scheme to be used
		
		
'''
def MATPROP(INFILE):
	#-----------------------------------------------------------------------------------
	# Declating Local Variables:
	# Constants
	R0 = float(0.0); R1 = float(1.0)
	
	# Strings:
	DUMMY = str(); MODEL = str(); LOADTYPE = str(); QUADSCHEME = str();
	
	# Integers
	I = int()
	
	# Unkowns
	RPROP = []; DATA = [];
	
	#-----------------------------------------------------------------------------------
	# Start Reading Input File
	FILE = open(INFILE,'r')
	FILE.readline();
	
	#-----------------------------------------------------------------------------------
	# Model Name
	FILE.readline();
	MODEL = FILE.readline().strip('\n');
	
	# Loading type
	FILE.readline();
	LOADTYPE = FILE.readline().strip('\n');
	
	#-----------------------------------------------------------------------------------
	# Basde on the previous lines, assemble material properties array
	
	FILE.readline();
	DATA = FILE.readline().split();
	
	if(MODEL == 'GAUSS'): # Gaussian Model
		RPROP = R0; # Only one material par (shear modulus)
		RPROP = float(DATA[0]);
		
	elif(MODEL == '3CHAIN'): # 3-chain model
		RPROP = np.zeros((3,),float)
		for I in range(0,len(DATA)):
			RPROP[I] = float(DATA[I]);
		
	elif(MODEL == '8CHAIN'): # 8-chain model
		RPROP = np.zeros((len(DATA),),float)
		for I in range(0,len(DATA)):
			RPROP[I] = float(DATA[I]);
	elif(MODEL == 'FULL'): # Full-Network Model
		RPROP = np.zeros((len(DATA),),float)
		for I in range(0,len(DATA)):
			RPROP[I] = float(DATA[I]);
			
	#-----------------------------------------------------------------------------------
	# Number of increments and Target Stretch
	FILE.readline();
	DATA = FILE.readline().split();
	NINCR = int(DATA[0]); MAXSTRETCH = float(DATA[1]);
	
	# Calculate 
	STRETCHINCR = (MAXSTRETCH - R1)/NINCR;
	
	#-----------------------------------------------------------------------------------
	# For the Full Network model we nedd the string of the Integration Scheme
	if (MODEL == 'FULL'):
		FILE.readline().strip('\n');
		QUADSCHEME = FILE.readline();
	
	#-----------------------------------------------------------------------------------
	# Close File
	FILE.close()
		
	return RPROP, MODEL, LOADTYPE, NINCR, STRETCHINCR, QUADSCHEME
	
#****************************************************************************************	
	

#======================================================================================
""" DEFGRAD: ASSEMBLE THE DEFORMATION GRADIENT BASED ON THE TYPE LOADING

Inputs:
		LOATYPE: String Informting loading condition;
		STRETCH: Stretch
		
Outputs:
		F: Deformation Gradient Tensor

"""

def DEFGRAD(LOADTYPE,STRETCH):
	#-----------------------------------------------------------------------------------
	# Declare Local Variables
	
	# Constants
	R0 = float(0.0); R1 = float(1.0); R2 = float(2.0);
	
	# Second order tensors
	F =  np.zeros((3,3),float);
	
	#-----------------------------------------------------------------------------------
	# Assemble deformation gradient
	
	F[0,0]=STRETCH; # Equal in all cases
	
	if (LOADTYPE == 'UNIAXIAL'):
		F[1,1] = R1/sqrt(F[0,0]);
		F[2,2] = F[1,1];
	elif (LOADTYPE == 'BIAXIAL'):
		F[1,1] = F[0,0];
		F[2,2] = R1/(F[0,0]**R2)
	elif (LOADTYPE == 'SHEAR'):
		F[1,1] = R1;
		F[2,2] = R1/F[0,0];
		
	return F
	
#****************************************************************************************

#======================================================================================
""" toPIOLA: Converte Cauchy Stress tensor to Piola  

Inputs:
		F: Deformation Gradient;
		T: Cauchy stress tensor
		JACOB: Volumetic ration or jacobian
Outputs:
		TR: Piola Stress tensor (First)

"""

def toPIOLA(F,T,JACOB):
	#-----------------------------------------------------------------------------------
	# Import relavant modules
	from numpy.linalg import inv
	
	#-----------------------------------------------------------------------------------
	# Declare Local Variables
	
	# Constants
	R0 = float(0.0); R1 = float(1.0); R2 = float(2.0);
	
	# Integers
	I = int(R0); J = int(R0);  K = int(R0); 
	
	# Second order tensors
	TR =  np.zeros((3,3),float);
	
	#-----------------------------------------------------------------------------------
	# Compute Piola Stress
	
	# Invert deformation gradient
	FINVT = inv(F)
	FINVT = FINVT.T
	
	for I in range(0,3):
		for J in range(0,3):
			for K in range(0,3):
				TR[I,J]=TR[I,J]+(T[I,K]*FINVT[K,J])
	
	TR=JACOB*TR
	
	return TR
	
#****************************************************************************************

#======================================================================================
""" INVLANGEVIN:   

Inputs:
		X: Point where the Langevin Inverse must be evaluated
		APPROX: Type of approximation desired
Outputs:
		INVL: Inverse Langevin 

"""

def INVLANGEVIN(X,APPROX = 'PADE'):
	#-----------------------------------------------------------------------------------
	# Declare Local Variables
	
	# Constants
	R0 = float(0.0); R1 = float(1.0); R2 = float(2.0); R3 = float(3.0);
	
	# Scalars
	INVL = R0;
	
	#-----------------------------------------------------------------------------------
	# Compute L^{-1}
	
	if(APPROX == 'PADE'): # Best Pade Aproximation x*(3.-(x**2))/(1.-(x**2))
		INVL = (X*( R3 - (X**R2) ) )/( R1 - (X**R2) );
	
	return INVL
	
#****************************************************************************************

#======================================================================================
""" QUAD:   

Inputs:
		NAME: String with the Quadrature Rule Desired
		obs: The input f
Outputs:
		Q: Number of quadrature points;
		DIRQUAD : Dictionary containinf quadrature nodes and weitghs

"""

def QUAD(NAME):
	#-----------------------------------------------------------------------------------
	# Declare Local Variables
	
	
	# Constants
	R0 = float(0.0); R1 = float(1.0); R2 = float(2.0); R3 = float(3.0);
	
	# Integers 
	I = int(R0); Q = int(R0);
	
	# Scalars
	
	# Logicals 
	FLAG = True
	
	# Dictionary
	DIRQUAD ={};
	
	# Unkowns
	DATA = [];
	
	#-----------------------------------------------------------------------------------
	# Open File 

	
	FILE = open(NAME,'r');
	
	#-----------------------------------------------------------------------------------
	# Start Reading the File
	
	# Read First Line of the file wich contain if the scheme is half-sphere based
	KEY = FILE.readline().strip('\n')
	
	if 'Half' in KEY:
		HALF = True
	
	while (FLAG):
		KEY = FILE.readline().strip('\n')
		if 'Node' in KEY:
			FLAG = False
			 
	
	# Reset Flag to true
	FLAG = True
	
	while (FLAG):
		KEY = FILE.readline().strip('\n') # Read line after line
		if 'END' in KEY:
			FLAG = False;
		else:
			DATA = KEY.split();
			DIRQUAD[int(DATA[0])] = [float(DATA[1]),float(DATA[2]),float(DATA[3]),float(DATA[4])]
		
	# Close File	
	FILE.close();
	
	# Number of Nodes
	Q = len(DIRQUAD);
	


	return Q, DIRQUAD, HALF
