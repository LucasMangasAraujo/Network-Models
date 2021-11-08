## FUNCTIONS USED TO ASSEMBLE ARRAYS, LOADING PATHS AND ETC (PRE-PROCESSING IF YOU WILL)


#======================================================================================
## MATPROP: THIS FUNCTION ASSEMBLE ARRAY OR SCALAR CARRYING THE MATERIAL PROPERTIES

# Inputs: 
# 		INFILE: TXT FILE CONTAINING THE MATERIAL PROPERTIES
#		MODEL: STRING CONTANING THE MODEL IDENTIFIER

# Outputs
#		RPROP: Array or scalar containing the material properties
def MATPROP(INFILE,MODEL):
	#-----------------------------------------------------------------------------------
	# Declating Local Variables:
	# Constants
	R0=float(0.0)
	
	# Strings:
	DUMMY=str()
	
	# Integers
	I=int()
	
	# Unkowns
	RPROP=[]
	
	#-----------------------------------------------------------------------------------
	# Carry out assembly
	if(MODEL=='GAUSS'): # Gaussian Model
		RPROP=R0;
	# Only one material parameter
		for line in INFILE:
			DUMMY=line.split()
			RPROP=eval(DUMMY[0])
	#elif(MODEL=='3Chain'): # Non-Gaussian 3-chain model
	
	


	return RPROP
	
#****************************************************************************************	
	

#======================================================================================
""" LOADING: THIS FUNCTION ASSEMBLES THE DEF HISTORY

OBS: The inputs are stretches

Inputs: 
 		LOADTYPE: String informing the laoding condition
		MODEL: Strinf informing the model 
		NINCR: Number of increments used to discretize the load
                DEF1:  Axial deformation in the first principal direction

Outputs
		RPROP: Array or scalar containing the material properties
"""
def LOADING(LOADTYPE,MODEL,NINCR,DEF1):
	#-----------------------------------------------------------------------------------
	# Import Relevant Modules
	import numpy as np
	#-----------------------------------------------------------------------------------
	# Declaring Local Variables:
	# Constants
	R0=float(0.0); R1=float(1);
	
	# Integers
	I=int(R0); J=int(R0);
	
	# Arrays
	
	#-----------------------------------------------------------------------------------
	# Assemble Load
	
	if(LOADTYPE=='UNIAXIAL' or LOADTYPE=='BIAXIAL' or LOADTYPE=='SHEAR'): 
		ETOT=np.zeros((NINCR+1,),float);
		ETOT[0]=R1
		EINCR=(DEF1-R1)/NINCR;
		for I in range(1,NINCR+1):
			ETOT[I]=ETOT[I-1]+EINCR;
		
		
	




	return ETOT