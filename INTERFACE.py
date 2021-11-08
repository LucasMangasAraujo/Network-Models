## INTERFACE: FUNCTIONS WICH WILL BE USED TO ASSEMBLE THE MECHANICAL RESPONSE


#**************************************************************************************


#======================================================================================
## MATISU: FUNCTION WHICH SELECTS THE FUNCTION USED TO COMPUTE THE NOMIMAL STRESS

# Inputs:
# 			LOADTYPE:
#			MODEL:
#           RPROPS: Material Properties Array
#			STRETCH: Chain stretch (SCALAR OR TENSOR) 
#			NSTRES: Dimensions of the stress


#******************

# Outputs:
# 			STRES: NOMINAL STRESS
#  			
 
#******************

def MATISU(LOADTYPE,MODEL,RPROPS,STRETCH,NSTRES):
	#-----------------------------------------------------------------------------------
	# Import Relevant Modules
	import numpy as np
	import STATEUPDATE as SU
	
	#-----------------------------------------------------------------------------------
	# Declaring Local Variables
	
	# Constants
	R0=float(0.0);
	
	# Integers
	I=int(R0); J=int(R0);
	
	# Arrays
	
	#-----------------------------------------------------------------------------------
	# Select Model to be used
	# Gauss
	if(MODEL=='GAUSS'):
		NOMSTRES=R0;
		NOMSTRES=SU.GAUSS(LOADTYPE,RPROPS,STRETCH,NSTRES)
	
	
	return NOMSTRES
	
	
#**************************************************************************************


 
 