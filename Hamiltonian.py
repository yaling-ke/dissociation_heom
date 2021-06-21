import numpy as np
from constants import *

class nuclear_DVR():
    def __init__(self,number_el_state,mass,N_DVR,x_0,x_N,x_infties,CAP_thresholds,CAP_pots,energies,el_couplings,CAP_button='on'):
        self.number_el_state=number_el_state
        self.mass=mass
        self.N_DVR=N_DVR
        self.x_0=x_0
        self.x_N=x_N
        self.x_infties=x_infties
        self.CAP_button=CAP_button
        self.CAP_thresholds=CAP_thresholds
        self.CAP_pots=CAP_pots
        self.energies=energies
        self.el_couplings=el_couplings
        self.L=self.x_N-self.x_0
        self.delta_x=self.L/(self.N_DVR+1)
        self.alpha_list=[alpha for alpha in range(1,self.N_DVR+1)]
        self.gridpoints=[self.x_0+alpha*self.delta_x for alpha in self.alpha_list]

        Hamil=self.H_DVR()
        eigenvals,eigenvecs=np.linalg.eigh(Hamil[0:self.N_DVR,0:self.N_DVR])
        self.initial_nuclear_distribution_vector=eigenvecs[:,0].copy()
        if self.CAP_button=='on':
            self.alpha_list.append(self.alpha_list[-1]+1)
            self.gridpoints.append(x_infties)
            self.N_DVR +=1
            self.initial_nuclear_distribution_vector=np.zeros(self.N_DVR,dtype=complex)
            self.initial_nuclear_distribution_vector[0:self.N_DVR-1]=eigenvecs[:,0].copy()


    # Creation and Annihilation Operators for the whole DVR basis:
    def Operators_Creation_Annihilation(self):
        operators=[]
        d_dagger=[]
        d=[]
        for i in range(1,self.number_el_state+1):
            op=np.zeros((2**self.number_el_state,2**self.number_el_state),float)
            for j in range (2**self.number_el_state):
                if len(str(bin(j)))-2 <i or str(bin(j))[-i]=='0':
                    pref=1.
                    for k in range(1,i):
                        if str(bin(j))[-k]=='1':
                            pref=pref*(-1)
                    op[j+2**(i-1)][j]=pref
            d_dagger.append(op)
            d.append(op.transpose())
        operators.append(d_dagger)
        operators.append(d)
        return np.array(operators)
        
        
    def H_DVR(self):
        #Define second derivative in DVR representation
        def D_2_DVR(alpha,beta):
            if alpha==beta:
                return -(np.pi/(self.delta_x/au2m))**2*(1./3+1./(6*(self.N_DVR+1)**2)-1./(2.*(self.N_DVR+1)**2*(np.sin(alpha*np.pi/(self.N_DVR+1))**2)))
            else:
                a=alpha*np.pi/(self.N_DVR+1)
                b=beta*np.pi/(self.N_DVR+1)
                return -(np.pi/(self.delta_x/au2m))**2*(float(2*(-1)**(alpha-beta))/((self.N_DVR+1)**2)*np.sin(a)*np.sin(b)/(np.cos(a)-np.cos(b))**2)

        h=np.zeros((2**self.number_el_state*self.N_DVR,2**self.number_el_state*self.N_DVR),complex)
        for fock_state in range(2**self.number_el_state):
            for alpha,gridpoint in zip(self.alpha_list,self.gridpoints):
                x=gridpoint
                h[fock_state*self.N_DVR+(alpha-1),fock_state*self.N_DVR+(alpha-1)]+=eval(self.energies[fock_state])
                for beta in self.alpha_list:
                    h[fock_state*self.N_DVR+(alpha-1),fock_state*self.N_DVR+(beta-1)]+= -1./(2.*self.mass)*D_2_DVR(alpha,beta)*au2eV

        for fock_state1 in range(2**self.number_el_state):
            for fock_state2 in range(fock_state1,2**self.number_el_state):
                for alpha, gridpoint in zip(self.alpha_list,self.gridpoints):
                    x=gridpoint
                    h[fock_state1*self.N_DVR+(alpha-1),fock_state2*self.N_DVR+(alpha-1)]+=eval(self.el_couplings[fock_state1][fock_state2])
                    h[fock_state2*self.N_DVR+(alpha-1),fock_state1*self.N_DVR+(alpha-1)]+=eval(self.el_couplings[fock_state1][fock_state2])
        return h


    def H_DVR_Lindblad(self):

        #Define second derivative in DVR representation
        def D_2_DVR(alpha,beta):
            if alpha==beta:
                return -(np.pi/(self.delta_x/au2m))**2*(1./3+1./(6*(self.N_DVR)**2)-1./(2.*(self.N_DVR)**2*(np.sin(alpha*np.pi/(self.N_DVR))**2)))
            else:
                a=alpha*np.pi/(self.N_DVR)
                b=beta*np.pi/(self.N_DVR)
                return -(np.pi/(self.delta_x/au2m))**2*(float(2*(-1)**(alpha-beta))/((self.N_DVR)**2)*np.sin(a)*np.sin(b)/(np.cos(a)-np.cos(b))**2)

        h=np.zeros((2**self.number_el_state*self.N_DVR,2**self.number_el_state*self.N_DVR),complex)
        h_CAP=np.zeros((2**self.number_el_state*self.N_DVR,2**self.number_el_state*self.N_DVR),complex)
        for fock_state in range(2**self.number_el_state):
            for alpha,gridpoint in zip(self.alpha_list[:-1],self.gridpoints[:-1]):
                x=gridpoint
                h[fock_state*self.N_DVR+(alpha-1),fock_state*self.N_DVR+(alpha-1)]+=eval(self.energies[fock_state])
                if self.CAP_button=="on":
                    if gridpoint>self.CAP_thresholds[fock_state]:   
                        h_CAP[fock_state*self.N_DVR+(alpha-1),fock_state*self.N_DVR+(alpha-1)]+= -1j*eval(self.CAP_pots[fock_state])
                for beta in self.alpha_list[:-1]:
                    h[fock_state*self.N_DVR+(alpha-1),fock_state*self.N_DVR+(beta-1)]+= -1./(2.*self.mass)*D_2_DVR(alpha,beta)*au2eV

            alphaa=self.alpha_list[-1]
            x=self.gridpoints[-1]
            h[fock_state*self.N_DVR+alphaa-1,fock_state*self.N_DVR+alphaa-1]+=eval(self.energies[fock_state])

        for fock_state1 in range(2**self.number_el_state):
            for fock_state2 in range(fock_state1,2**self.number_el_state):
                for alpha, gridpoint in zip(self.alpha_list,self.gridpoints):
                    h[fock_state1*self.N_DVR+(alpha-1),fock_state2*self.N_DVR+(alpha-1)]+=eval(self.el_couplings[fock_state1][fock_state2])
                    h[fock_state2*self.N_DVR+(alpha-1),fock_state1*self.N_DVR+(alpha-1)]+=eval(self.el_couplings[fock_state1][fock_state2])

                alphaa=self.alpha_list[-1]
                x_1=self.gridpoints[-1]
                h[fock_state1*self.N_DVR+(alphaa-1),fock_state2*self.N_DVR+(alphaa-1)]+=eval(self.el_couplings[fock_state1][fock_state2])
                h[fock_state2*self.N_DVR+(alphaa-1),fock_state1*self.N_DVR+(alphaa-1)]+=eval(self.el_couplings[fock_state1][fock_state2])
        return h,h_CAP

    def get_electronic_operators(self,energies_el,el_couplings_el):
        hamil_el=np.zeros((2**self.number_el_state,2**self.number_el_state),float)
        for fock_state in range(2**self.number_el_state):
            hamil_el[fock_state,fock_state]=energies_el[fock_state]
        for fock_state1 in range(2**self.number_el_state):
            for fock_state2 in range(fock_state1+1,2**self.number_el_state):
                hamil_el[fock_state1,fock_state2]=el_couplings_el[fock_state1][fock_state2]
                hamil_el[fock_state2,fock_state1]=el_couplings_el[fock_state1][fock_state2]
        return self.Operators_Creation_Annihilation(), hamil_el

    def get_operators(self):
        if self.CAP_button=='on':
            H,H_CAP=self.H_DVR_Lindblad()
        else:
            H,H_CAP=self.H_DVR()
        return H,H_CAP

    def get_rho_0(self,rho_0):
        self.rho_0_DVR=np.zeros((2**self.number_el_state*self.N_DVR,2**self.number_el_state*self.N_DVR),complex)
        for fock_state in range(2**self.number_el_state):
            for alpha in range(len(self.initial_nuclear_distribution_vector)):
                for beta in range(len(self.initial_nuclear_distribution_vector)):
                    self.rho_0_DVR[fock_state*self.N_DVR+alpha,fock_state*self.N_DVR+beta]=complex(rho_0[fock_state,fock_state]*(self.initial_nuclear_distribution_vector[alpha].conjugate()*self.initial_nuclear_distribution_vector[beta]))
        if abs(abs(np.trace(self.rho_0_DVR))-1.0)>1e-10:
            print("Initial density matrix is not normalized.")
            exit()
        else:
            return self.rho_0_DVR

    def get_CAP(self):
        cap=np.zeros((2**self.number_el_state,self.N_DVR),float)
        if self.CAP_button =='on':
            for alpha,gridpoint in zip(self.alpha_list,self.gridpoints):
                x=gridpoint
                for fock_state in range(2**self.number_el_state):
                    if gridpoint >= self.CAP_thresholds[fock_state]:
                        cap[fock_state,(alpha-1)]+=eval(self.CAP_pots[fock_state])
        return cap

    def get_V_DVR(self,el_lead_couplings,number_lead):
        Vs=np.zeros((number_lead,self.number_el_state,self.N_DVR),float)

        for lead in range(number_lead):
            for el_state in range(self.number_el_state):
                for alpha,gridpoint in zip(self.alpha_list,self.gridpoints):
                    x=gridpoint
                    Vs[lead,el_state,(alpha-1)]=eval(el_lead_couplings[lead][el_state])
        return Vs

    def get_V_DVR_thermal(self,el_bath_couplings,number_bath):
        Vs_T=np.zeros((number_bath,2**self.number_el_state*self.N_DVR),float)

        for bath_state in range(number_bath):
            for fock_state in range(2**self.number_el_state):
                for alpha,gridpoint in zip(self.alpha_list,self.gridpoints):
                    x=gridpoint
                    Vs_T[bath_state,fock_state*self.N_DVR+(alpha-1)]=eval(el_bath_couplings[bath_state][fock_state])
        return Vs_T

    def get_phi_DVR(self,j,x): #DVR basisfunctions
        if self.x_0 <= x <=self.x_N:
            return np.sqrt(2./self.L)*np.sin(j*np.pi*(x-self.x_0)/self.L)
        else:
            return 0

    def get_U_DVR(self, j,alpha): # transformation matrix between DVR basisfunctions and the eigenfunctions of the position operator
        return np.sqrt(2./(self.N_DVR+1))*np.sin(np.pi*j*alpha/(self.N_DVR+1))

    def get_chi_DVR(self,alpha,x): #eigenfunctions of the position operator in the DVR basis
        s=0.
        for j in range(1,self.N_DVR+1):
            s+=self.get_phi_DVR(j,x)*self.get_U_DVR(j,alpha)
        return s

    def reconstruct_eigenstate(self,i,eigenvecs,x):
        s=0.
        for j in range(self.N_DVR):
            s+=eigenvecs[j,i]*self.get_chi_DVR(j+1,x)
        return s

