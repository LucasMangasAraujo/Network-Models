## STATEUPDATE: COLLECTION OF NETWORK MODELS

#**************************************************************************************


#======================================================================================
""" SUGAUSS: THIS FUNCTION COMPUTES THE NOMINAL STRESS FOR A GIVEN LOAD CONDITION CONDERING
GAUSS MODEL

Inputs: 
		LOADTYPE: String with the laoding partition
		RPROPS: Real array or sclar with material properties
		STRETCH: Stretch
		NSTRES: Integer with nunber of stress components
 			
******************

 Outputs:	
		NOMSTRES: First entry of the Firt Piola Stress tensor
		
 
******************
"""

def GAUSS(LOADTYPE,RPROPS,STRETCH,NSTRES):
  	#-----------------------------------------------------------------------------------
	## Import Relevant Modules
	import numpy as np
	from numpy.linalg import inv
	
	#-----------------------------------------------------------------------------------
	## Declaring Local Variables
	
	# Constants
	R0=float(0.0); RP5=float(0.5); R1=float(1.0); R2=float(2.0); R3=float(3.)
	
	# Integers
	I=int(R0); J=int(R0); K=int(R0);
	
	# Scalar
	R1D3=R0; GMODU=R0; DFREEDI1=R0; P=R0; I1=R0;
	
	# Vectors
	
	# Second Order Tensors
	DEFGRAD=np.zeros((3,3),float); B=np.zeros((3,3),float);
	T=np.zeros((3,3),float);       TR=np.zeros((3,3),float);
	DEVB=np.zeros((3,3),float);    DEFGRADINVT=np.zeros((3,3),float); 
	
	# Unknoes
	NOMSTRES=[];
	
	#-----------------------------------------------------------------------------------
	## Setting Useful Parameters
	R1D3=R1/R3
	
	
	#-----------------------------------------------------------------------------------
	## Getting Material Parameters
	GMODU=np.copy(RPROPS);
	
	
	#-----------------------------------------------------------------------------------
	## Define form of the constitutive relatio
	
	DFREEDI1=RP5*GMODU # Free-Energy Derivative
	
	# Kinematics and related variables depending on the loading
	if(LOADTYPE=='UNIAXIAL'): # Uniaxial Loading
		# Assemble Deformatinon Gradient
		DEFGRAD[0,0]=STRETCH; DEFGRAD[1,1]=R1/np.sqrt(STRETCH);
		DEFGRAD[2,2]=R1/np.sqrt(STRETCH);
		
		for I in range(0,3):
			# Arrange Left Cauchy Green Tenso
			for J in range(0,3):
				for K in range(0,3):
					B[I,J]=B[I,J]+DEFGRAD[I,K]*DEFGRAD[J,K]
		
		# Assemble Deviatoric Part of B
		I1=B[0,0]+B[1,1]+B[2,2] # First Invariant of B
		
		DEVB[0,0]=B[0,0]-(R1D3*I1); 
		DEVB[1,1]=B[1,1]-(R1D3*I1);
		DEVB[2,2]=B[2,2]-(R1D3*I1);
		
		# Compute Indeterminate pressure
		P=-R2*DFREEDI1*DEVB[1,1]
		
	elif(LOADTYPE=='BIAXIAL'): # EQUIBIAXIAL loading
		DEFGRAD[0,0]=STRETCH; DEFGRAD[1,1]=STRETCH;
		DEFGRAD[2,2]=R1/(STRETCH**R2);
		
		# Assemble left Cauchy-Green Tensor
		for I in range(0,3):
			for J in range(0,3):
				for K in range(0,3):
					B[I,J]=B[I,J]+DEFGRAD[I,K]*DEFGRAD[J,K]
		
		# Compute First Invariant 
		I1=B[0,0]+B[1,1]+B[2,2]
		
		# Extract dev(B)
		DEVB[0,0]=B[0,0]-(R1D3*I1); 
		DEVB[1,1]=B[1,1]-(R1D3*I1);
		DEVB[2,2]=B[2,2]-(R1D3*I1);
		
		# Pressure Stress From Boundary Conditions
		P=-R2*DFREEDI1*DEVB[2,2]
	
	elif(LOADTYPE=='SHEAR'):
		DEFGRAD[0,0]=STRETCH; DEFGRAD[1,1]=R1;
		DEFGRAD[2,2]=R1/(STRETCH);
		
		# Assemble left Cauchy-Green Tensor
		for I in range(0,3):
			for J in range(0,3):
				for K in range(0,3):
					B[I,J]=B[I,J]+DEFGRAD[I,K]*DEFGRAD[J,K]
		
		# Compute First Invariant 
		I1=B[0,0]+B[1,1]+B[2,2]
		
		# Extract dev(B)
		DEVB[0,0]=B[0,0]-(R1D3*I1); 
		DEVB[1,1]=B[1,1]-(R1D3*I1);
		DEVB[2,2]=B[2,2]-(R1D3*I1);
		
		# Pressure Stress From Boundary Conditions
		P=-R2*DFREEDI1*DEVB[2,2]		
		
	
	#----------------------------------------------------------------------------------------	
	## Compute Stress Tensors
	
	# Cauchy Stress Tensor
	for I in range(0,3):
		for J in range(0,3):
			if(I==J):
				T[I,J]=R2*DFREEDI1*DEVB[I,J]+P
			else:
				T[I,J]=R2*DFREEDI1*DEVB[I,J]
	
	# Compute First Piola stress tensor
	DEFGRADINVT=DEFGRADINVT.T
	DEFGRADINVT=inv(DEFGRAD)
	
	for I in range(0,3):
		for J in range(0,3):
			for K in range(0,3):
				TR[I,J]=TR[I,J]+(T[I,K]*DEFGRADINVT[K,J])
			
	# First entry of the Piola Stress Tensor	
	NOMSTRES=TR[0,0]
  
  
	return NOMSTRES
  
 #======================================================================================
""" 3CHAIN: THIS FUNCTION COMPUTES THE NOMINAL STRESS FOR A GIVEN LOAD CONDITION CONDERING
3 CHAIN MODEL

Inputs: 
		LOADTYPE: String with the laoding partition
		RPROPS: Real array or sclar with material properties
		STRETCH: Stretch
		NSTRES: Integer with nunber of stress components
 			
******************

 Outputs:	
		NOMSTRES: First entry of the Firt Piola Stress tensor
		
 
******************
"""
 