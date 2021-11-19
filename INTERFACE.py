## INTERFACE: FUNCTIONS WICH WILL BE USED TO ASSEMBLE THE MECHANICAL RESPONSE


#-----------------------------------------------------------------------------------
# Import Relevant Modules
import numpy as np
import STATEUPDATE as SU

#**************************************************************************************


#======================================================================================
''' MATISU: FUNCTION WHICH SELECTS THE FUNCTION USED TO COMPUTE THE NOMIMAL STRESS

 Inputs:
 		LOADTYPE:
		MODEL:
        RPROPS: Material Properties Array
		STRETCH: Chain stretch (SCALAR OR TENSOR) 
		NSTRES: Dimensions of the stress


******************

 Outputs:
 		STRES: NOMINAL STRESS
  			
 
******************
'''

def MATISU(MODEL,LOADTYPE,RPROPS,STRETCH):	
	#-----------------------------------------------------------------------------------
	# Declaring Local Variables
	
	# Constants
	R0=float(0.0);
	
	# Integers
	I=int(R0); J=int(R0);
	
	# Second Order Tensors
	T = np.zeros((3,3),float); TR = np.zeros((3,3),float); 
	
	#-----------------------------------------------------------------------------------
	# Select Model to be used
	# Gauss
	if(MODEL == 'GAUSS'):
		TR,T = SU.GAUSS(LOADTYPE,RPROPS,STRETCH);
	elif(MODEL == '3CHAIN'):
		# For Default, The initial end-to-end distance is equal to N^0.5*B
		TR,T = SU.THREECHAIN(LOADTYPE,RPROPS,STRETCH); # 
	elif(MODEL == '8CHAIN'):
		# For Default, The initial end-to-end distance is equal to N^0.5*B
		TR,T = SU.EIGHTCHAIN(LOADTYPE,RPROPS,STRETCH);
	
	
	return TR,T
	
	
#**************************************************************************************


 
 