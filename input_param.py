import numpy as np
from constants import * 

#INPUT:
#--------------------------------------------------------------------------------------
#Number of electrodes
Number_lead=2
#Number of phonon baths
Number_thermal_bath=1
#Number of Pade poles for fermi function
Number_pole_Fermi=20 
#Number of Pade poles for Bose function
Number_pole_Bose=2
#Fermionic hierarchical turncation tier
Truncation_tier_Fermi=2
#Bosonic hierarchical turncation tier
Truncation_tier_Bose=1
N_index=Number_thermal_bath*(Number_pole_Bose+1)
# sigma +/- (creation/annihilation)
Sign_pm=2
#Number of electronic levels
Number_el_state=1
dim_fockspace=2**Number_el_state

#chemical potential [mu_L,mu_R] in units of eV
chem_potentials=np.array([2.0,-2.0])
#temperature of electrodes in units of K
temperature_Fermi=300
#temperature of phonon baths in units of K
temperature_Bose=300
High_temperature_approximation=False
Ishizaki_trunc_on=False                                                      

#ADO filtering threshold
threshold_on=False
if threshold_on:    threshold_value=1.e-4
# True: wide band limit, False: lorentzian broadening
WBL_on=True
if not WBL_on:  broadenings=np.array([5.0,5.0])  # in units of eV


#reorganization energy in units of eV
lambda_T=np.array([0.05])            
# charateristic time of the bath in units of fs
tau_c= 50
# charateristic frequency of the bath (=1/tau_c)    
wc= np.array([1./tau_c*(hbar_eVs*1e15)])                   

# Time step in units of fs
dt=2.e-2
# Output at every dt/dt_out points, dt_out in units of fs
dt_out=2.0e-1       
# Time steps
nsteps=5000
# Propagation from t=0 or not
Init_propagation=True
if not Init_propagation:    time_intermediate=0.

# nuclear mass
mass=amu
# Parameters of Morse potential surface
# width (in units of m^{-1})
a=1.028/Bohr 
# Dissociation energy in units of eV
D_e=2.38 

mass=mass/au2mass
Morse_fre=a*(au2m)*np.sqrt(2.0*D_e/au2eV/(mass))*au2eV
chi=Morse_fre/D_e

# for sparsity check
rho_0_el=np.zeros((dim_fockspace,dim_fockspace),float)
rho_0_el[0,0]=1.
energies_el=np.array([1.,1.])
el_couplings_el=np.zeros((2,2),float)

#Neutral and charged state potential energy surfaces
energies=[]
energies.append("2.38*(np.exp(-1.028*(x-1.78e-10)/"+str(Bohr)+")-1.0)**2")
energies.append("2.38*(np.exp(-1.028*(x-1.98e-10)/"+str(Bohr)+")-1.0)**2+1.0")
# electronic couplings
el_couplings=[["0"]*dim_fockspace]*dim_fockspace

# coupling of the electronic level to left and right lead
el_lead_coupling_left=np.sqrt(0.05/(2.0*np.pi) )
el_lead_coupling_right=np.sqrt(0.05/(2.0*np.pi) )
max_el_lead_coupling=max(el_lead_coupling_left,el_lead_coupling_right)
el_lead_couplings=[]
el_lead_couplings.append(["(0.475*(1.0-np.tanh(2*(x-4.0e-10)/"+str(Bohr)+"))+0.05)*"+str(el_lead_coupling_left)])
el_lead_couplings.append(["(0.475*(1.0-np.tanh(2*(x-4.0e-10)/"+str(Bohr)+"))+0.05)*"+str(el_lead_coupling_right)])

# coupling to the phonon bath
el_bath_couplings=np.array([["(x-1.78e-10)/"+str(Bohr)+"*np.exp(-1.0*"+str(chi)+"*((x-1.78e-10)/"+str(Bohr)+")**2)","(x-1.98e-10)/"+str(Bohr)+"*np.exp(-1.0*"+str(chi)+"*((x-1.98e-10)/"+str(Bohr)+")**2)"]])

# complex absorbing potential 
CAP_thresholds=[4.0e-10,4.0e-10]
CAP_pots=["0.*((x-"+str(CAP_thresholds[0])+")/1e-10)**4*np.heaviside(x-"+str(CAP_thresholds[0])+",0)","0.*((x-"+str(CAP_thresholds[1])+")/1e-10)**4*np.heaviside(x-"+str(CAP_thresholds[1])+",0)"]

#DVR grids  # in units of meter
x_0=1.3e-10
x_N=5.e-10
# N_DVR must be odd
N_DVR=63
x_infties=1.0
