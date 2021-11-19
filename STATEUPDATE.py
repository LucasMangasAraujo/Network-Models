## STATEUPDATE: COLLECTION OF NETWORK MODELS

#-----------------------------------------------------------------------------------
## Import Relevant Modules
import numpy as np
from math import *
from UTILS import *

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
		TR: Piola Stress Tensor
		T: Cauchy stress tensor
		
 
******************
"""

def GAUSS(LOADTYPE,RPROPS,STRETCH):
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
	F=np.zeros((3,3),float); B=np.zeros((3,3),float);
	T=np.zeros((3,3),float);       TR=np.zeros((3,3),float);
	DEVB=np.zeros((3,3),float);    FINVT=np.zeros((3,3),float); 
	
	
	#-----------------------------------------------------------------------------------
	## Setting Useful Parameters
	R1D3=R1/R3
	
	
	#-----------------------------------------------------------------------------------
	## Getting Material Parameters
	GMODU=np.copy(RPROPS);
	
	
	#-----------------------------------------------------------------------------------
	## Define form of the constitutive relatio
	
	DFREEDI1=RP5*GMODU # Free-Energy Derivative
	
	# Assable Deformatinon Gradient
	F = DEFGRAD(LOADTYPE,STRETCH)
	
	# Assemble Left Cauchy Green Tensor
	for I in range(0,3):
		for J in range(0,3):
			for K in range(0,3):
				B[I,J]=B[I,J]+F[I,K]*F[J,K]
	
	# Extract first invariant of B and its deviatoric partition
	I1 = B[0,0]+B[1,1]+B[2,2] # First Invariant of B
	
	DEVB[0,0]=B[0,0]-(R1D3*I1); 
	DEVB[1,1]=B[1,1]-(R1D3*I1);
	DEVB[2,2]=B[2,2]-(R1D3*I1);
	
	# Hydrostatic stress assuing 3 -directions is a stress-free
	P = -R2*DFREEDI1*DEVB[2,2]
	
	#----------------------------------------------------------------------------------------	
	## Compute Stress Tensors
	
	# Cauchy Stress Tensor
	for I in range(0,3):
		for J in range(0,3):
			if(I == J):
				T[I,J]=R2*DFREEDI1*DEVB[I,J]+P
			else:
				T[I,J]=R2*DFREEDI1*DEVB[I,J]
	
	TR = toPIOLA(F,T,R1) 
  
  
	return TR,T
  
#======================================================================================


""" THREECHAIN: THIS FUNCTION COMPUTES THE NOMINAL STRESS FOR A GIVEN LOAD CONDITION CONDERING
3 CHAIN MODEL

Inputs: 
		LOADTYPE: String with the laoding partition
		RPROPS: Real array or sclar with material properties
		STRETCH: Stretch
 			
******************

 Outputs:	
		NOMSTRES: First entry of the Firt Piola Stress tensor
		
 
******************
"""

def THREECHAIN(LOADTYPE,RPROPS,STRETCH,RMSFLAG = True):
	#-----------------------------------------------------------------------------------
	## Declaring Local Variables
	
	# Constants
	R0=float(0.0); RP5=float(0.5); R1=float(1.0); R2=float(2.0); R3=float(3.)
	
	# Integers
	I = int(R0); J = int(R0); K = int(R0); NKUHN = int(R0);
	
	# Scalar
	R1D3 = R0; GMODU = R0; P = R0; I1 = R0; BKUHN = R1;
	
	
	# Vectors
	DFREEDL = np.zeros((3,),float);
	
	# Second Order Tensors
	F = np.zeros((3,3),float); B=np.zeros((3,3),float);
	T = np.zeros((3,3),float);       TR=np.zeros((3,3),float);
	DEVB = np.zeros((3,3),float);    FINVT = np.zeros((3,3),float);
	DFREEDBI = np.zeros((3,3),float); DFREEDBI_BI = np.zeros((3,3),float);
	SOID = np.zeros((3,3),float);     SHAPE = np.zeros((3,3),float);
	
	#-----------------------------------------------------------------------------------
	## Setting Useful Parameters
	R1D3=R1/R3
	
	#-----------------------------------------------------------------------------------
	## Getting material properties of the model
	GMODU = RPROPS[0]; NKUHN = RPROPS[1]; BKUHN = RPROPS[2];
	
	# Reference configuration end-to-end distance
	if (RMSFLAG):
		R_REF = sqrt(NKUHN)*BKUHN;
	else:
		R_REF = RPROPS[3];
		
	#-----------------------------------------------------------------------------------
	## Define form of the constitutive relation
	
	# Assemble Deformatinon Gradient
	F = DEFGRAD(LOADTYPE,STRETCH)
	
	# Assemble Left Cauchy Green Tensor
	for I in range(0,3):
		for J in range(0,3):
			for K in range(0,3):
				B[I,J]=B[I,J]+F[I,K]*F[J,K]
	

	#%%%%%%%%%%%%%%%%%%%%%%%%%%
	# Free energy derivative
	
	# 
	for I in range(0,3):
		DFREEDBI[I,I] = ( R1/(R2*F[I,I]) )*( GMODU*R1D3*(NKUHN*BKUHN/R_REF)*INVLANGEVIN(F[I,I]*R_REF/(NKUHN*BKUHN)) )
	
	for I in range(0,3):
		for J in range(0,3):
			for K in range(0,3):
				DFREEDBI_BI[I,J] = DFREEDBI_BI[I,J] + DFREEDBI[I,K]*B[K,J]
	
	DC_DFREEDBI_BI = DFREEDBI[0,0]*B[0,0] + DFREEDBI[1,1]*B[1,1] + DFREEDBI[2,2]*B[2,2]  
				
	# Assebmle shape tensor
	for I in range(0,3):
		for J in range(0,3):
			SHAPE[I,J] = R2*DFREEDBI_BI[I,J] - (R2*R1D3*DC_DFREEDBI_BI*SOID[I,J]);
	
	#%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	# Pressure from the boundary conditions
	P = -R1*SHAPE[2,2];
	
	#----------------------------------------------------------------------------------------	
	## Compute Stress Tensors
	
	# Cauchy Stress Tensor
	for I in range(0,3):
		for J in range(0,3):
			if(I == J):
				T[I,J] = SHAPE[I,J] + P
			else:
				T[I,J] = SHAPE[I,J]
	
	TR = toPIOLA(F,T,R1)

	return TR,T
	
#======================================================================================
""" EIGHTCHAIN: THIS FUNCTION COMPUTES THE NOMINAL STRESS FOR A GIVEN LOAD CONDITION CONDERING
3 CHAIN MODEL

Inputs: 
		LOADTYPE: String with the laoding partition
		RPROPS: Real array or sclar with material properties
		STRETCH: Stretch
 			
******************

 Outputs:	
		TR: Piola Stres Tensor
		T:  Cauchy Stress Tensor
		
 
******************
"""

def EIGHTCHAIN(LOADTYPE,RPROPS,STRETCH,RMSFLAG = True):
	#----------------------------------------------------------------------------------------
	## Declaring Local Variables
	
	# Constants
	R0=float(0.0); RP5=float(0.5); R1=float(1.0); R2=float(2.0); R3=float(3.)
	
	# Integers
	I = int(R0); J = int(R0); K = int(R0);
	
	# Scalar
	R1D3 = R0; GMODU = R0; DFREEDI1 = R0; P = R0; I1 = R0; LAMBDA = R0;
	
	# Vectors
	
	# Second Order Tensors
	F=np.zeros((3,3),float); B=np.zeros((3,3),float);
	T=np.zeros((3,3),float);       TR=np.zeros((3,3),float);
	DEVB=np.zeros((3,3),float);     
	
	
	#-----------------------------------------------------------------------------------
	## Setting Useful Parameters
	R1D3=R1/R3; 
	
	#-----------------------------------------------------------------------------------
	## Getting material properties of the model
	GMODU = RPROPS[0]; NKUHN = RPROPS[1]; BKUHN = RPROPS[2];
	
	# Reference configuration end-to-end distance
	if (RMSFLAG):
		R_REF = sqrt(NKUHN)*BKUHN;
	else:
		R_REF = RPROPS[3];
		
	#-----------------------------------------------------------------------------------
	## Define form of the constitutive relation
	
	# Assable Deformatinon Gradient
	F = DEFGRAD(LOADTYPE,STRETCH)
	
	# Assemble Left Cauchy Green Tensor
	for I in range(0,3):
		for J in range(0,3):
			for K in range(0,3):
				B[I,J]=B[I,J]+F[I,K]*F[J,K]
	
	# Extract first invariant of B and its deviatoric partition
	I1 = B[0,0] + B[1,1] + B[2,2] # First Invariant of B
	
	DEVB[0,0]=B[0,0] - (R1D3*I1); 
	DEVB[1,1]=B[1,1] - (R1D3*I1);
	DEVB[2,2]=B[2,2] - (R1D3*I1);
		
	# Assebmle Chain Stretch
	LAMBDA = sqrt(I1/R3);
	
	# Free-Energy Derivative
	DFREEDI1=RP5*(GMODU*R1D3)*(NKUHN*BKUHN/R_REF)*(R1/LAMBDA)*INVLANGEVIN( R_REF*LAMBDA/(NKUHN*BKUHN) ); # Free-Energy Derivative
	
	# Hydrostatic stress assuing 3 -directions is a stress-free
	P = -R2*DFREEDI1*DEVB[2,2]
	
	#-----------------------------------------------------------------------------------
	## Compute Stress Tensors
	
	# Cauchy Stress Tensor
	for I in range(0,3):
		for J in range(0,3):
			if(I == J):
				T[I,J]=R2*DFREEDI1*DEVB[I,J]+P
			else:
				T[I,J]=R2*DFREEDI1*DEVB[I,J]
	
	TR = toPIOLA(F,T,R1) 


	return TR,T