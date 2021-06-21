import numpy as np
from constants import *

class Pade_decomposition():

    def __init__(self,number_Pade_poles,Spec):
        self.number_Pade_poles=number_Pade_poles
        self.Spec=Spec
        if 'Fermi' in self.Spec or 'fermi' in self.Spec:
            self.Generate_Fermi_Pade_Parameters()
        elif 'Bose' in self.Spec or 'bose' in self.Spec:
            self.Generate_Bose_Pade_Parameters()


    def Generate_Fermi_Pade_Parameters(self):
        """
        Function generating the Pade poles and their respective weights. The routine used here is the [N-1,N] PSD scheme proposed in JCP,134,244106(2011).
        """
        def Kronecker_Delta(m,n):
            if m==n:
                return 1.
            else:
                return 0.
        def b_Fermi(m):
            return 2.*m-1.
    
        M=2*self.number_Pade_poles
        Delta       =np.zeros((M,M),float)
        Delta_tilde =np.zeros((M-1,M-1),float) 
    
        for m in range(M):
            for n in range(M):
                Delta[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Fermi(m+1)*b_Fermi(n+1))
    
        for m in range(M-1):
            for n in range(M-1):
                Delta_tilde[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Fermi(m+2)*b_Fermi(n+2))
    
        roots_Delta         =np.linalg.eigvalsh(Delta)
        roots_Delta_tilde   =np.linalg.eigvalsh(Delta_tilde)
    
        self.Xis_Fermi=[]
        for root in roots_Delta:
            if root>0:
                self.Xis_Fermi.append(2./abs(root))
        self.Xis_Fermi.sort()
        self.Chis_Fermi=[]
        for root in roots_Delta_tilde:
            if abs(root)>=1.e-13:
                if root>0:
                    self.Chis_Fermi.append(2./abs(root))
        self.Chis_Fermi.sort()
    
        self.kappas_Fermi=[]
        for j in range(self.number_Pade_poles):
            Quotient=1.
            for k in range(self.number_Pade_poles):
                if k!=(self.number_Pade_poles-1) and k!=j:
                    Quotient*=(self.Chis_Fermi[k]**2-self.Xis_Fermi[j]**2)/(self.Xis_Fermi[k]**2-self.Xis_Fermi[j]**2)
                elif k==(self.number_Pade_poles-1) and k!=j:
                    Quotient*=1./(self.Xis_Fermi[k]**2-self.Xis_Fermi[j]**2)
                elif k!=(self.number_Pade_poles-1) and k==j:
                    Quotient*=(self.Chis_Fermi[k]**2-self.Xis_Fermi[j]**2)
            Quotient*=0.5*self.number_Pade_poles*b_Fermi(self.number_Pade_poles+1)
            self.kappas_Fermi.append(Quotient)
    
    def Generate_Bose_Pade_Parameters(self):
        """
        Function generating the Pade poles and their respective weights. The routine used here is the [N-1,N] PSD scheme proposed in JCP,134,244106(2011).
        """
        def Kronecker_Delta(m,n):
            if m==n:
                return 1.
            else:
                return 0.
        
        def b_Bose(m):
            return 2.*m+1.
    
        M=2*self.number_Pade_poles
        Delta       =np.zeros((M,M),float)
        Delta_tilde =np.zeros((M-1,M-1),float) 
        for m in range(M):
            for n in range(M):
                Delta[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Bose(m+1)*b_Bose(n+1))
    
        for m in range(M-1):
            for n in range(M-1):
                Delta_tilde[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Bose(m+2)*b_Bose(n+2))
    
        roots_Delta         =np.linalg.eigvalsh(Delta)
        roots_Delta_tilde   =np.linalg.eigvalsh(Delta_tilde)
    
        self.Xis_Bose=[]
        for root in roots_Delta:
            if root >0:
                self.Xis_Bose.append(2./root)
        self.Xis_Bose.sort()
        self.Chis_Bose=[]
        for root in roots_Delta_tilde:
            if abs(root)>=1.e-13:
                if root >0:
                    self.Chis_Bose.append(2./abs(root))
        self.Chis_Bose.sort()
    
        self.kappas_Bose=[]
        for j in range(self.number_Pade_poles):
            Quotient=1.
            for k in range(self.number_Pade_poles):
                if k!=(self.number_Pade_poles-1) and k!=j:
                    Quotient*=(self.Chis_Bose[k]**2-self.Xis_Bose[j]**2)/(self.Xis_Bose[k]**2-self.Xis_Bose[j]**2)
                elif k==(self.number_Pade_poles-1) and k!=j:
                    Quotient*=1./(self.Xis_Bose[k]**2-self.Xis_Bose[j]**2)
                elif k!=(self.number_Pade_poles-1) and k==j:
                    Quotient*=(self.Chis_Bose[k]**2-self.Xis_Bose[j]**2)
            Quotient*=0.5*self.number_Pade_poles*b_Bose(self.number_Pade_poles+1)
            self.kappas_Bose.append(Quotient)
    
    
    def Fermi_Pade_Approx(self, energy=0.,mu=0.,temperature=20.):
        """
        Approximation of the Fermi-Dirac distribution function using the Pade poles 'etas' and 'Xis'.
        Input: energy, chemical potential mu, and temperatue are in units of eV, eV, and K, respectively. Pade poles, etas and Xis, are dimensionless.
        Output: approximate Fermi-Dirac distribution function.
        """
        f=0.5
        x=(energy-mu)/(k_B*temperature)
        for j in range(len(self.kappas_Fermi)):
            f-=2.*self.kappas_Fermi[j]*x/(x**2+self.Xis_Fermi[j]**2)
        return f
    
    def Fermi_V(self,energy=0.,mu=0.,temperature=20.):
        """
        Input: energy, the chemical potential mu, and the temperature are in units of eV, eV, and K, respectively.
        Output: Fermi-Dirac Distribution Function.
        """
        try:
            return 1./( 1.+np.exp( (energy-mu)/(k_B*temperature) ) )
        except:
            return 0.
    
    def Bose_Pade_Approx(self,energy=0.,temperature=10.):
        """
        Approximation of the Bose-Einstein distribution function using the Pade poles 'etas' and 'Xis'.
        Input: energy, chemical potential mu, and temperatue are in units of eV, eV, and K, respectively. Pade poles, etas and Xis, are dimensionless.
        Output: approximate Bose-Einstein distribution function.
        """
        x=energy/(k_B*temperature)
        f=0.5+1./x
        for j in range(len(self.kappas_Bose)):
            f+=2.*self.kappas_Bose[j]*x/(x**2+self.Xis_Bose[j]**2)
        return f

    def Bose_V(self,energy=0.1,temperature=20.):
        """
        Input: energy and the temperature are in units of eV and K, respectively.
        Output: Bose-Einstein Distribution Function.
        """
        try:
            return 1./( 1.-np.exp( -energy/(k_B*temperature) ) )
        except:
            return 0.

    def get_Pade_parameters(self):
        if 'Fermi' in self.Spec or 'fermi' in self.Spec:
            return self.kappas_Fermi,self.Xis_Fermi
        elif 'Bose' in self.Spec or 'bose' in self.Spec:
            return self.kappas_Bose,self.Xis_Bose

    def get_etas_gammas_Fermi(self, chem_potentials, temperature, broadenings):
        etas=[]
        gammas=[]  
        for i in range(len(chem_potentials)):
            etas_sub=[np.pi * broadenings[i] * self.Fermi_Pade_Approx(1j*broadenings[i],0.,temperature)]
            for j in range(len(self.Xis_Fermi)):
                etas_sub.append(-1j * 2*np.pi * k_B*temperature * self.kappas_Fermi[j] * broadenings[i]**2 /( broadenings[i]**2 - (k_B*temperature*self.Xis_Fermi[j])**2 ) )
            etas.append(etas_sub)

            gammas_sub=[[broadenings[i]-1j*chem_potentials[i],broadenings[i]+1j*chem_potentials[i]]]
            for j in range(len(self.Xis_Fermi)):
                gammas_sub.append([self.Xis_Fermi[j]*k_B*temperature-1j*chem_potentials[i],self.Xis_Fermi[j]*k_B*temperature+1j*chem_potentials[i]])
            gammas.append(gammas_sub)
        return np.array(etas), np.array(gammas)

    def get_etas_gammas_Fermi_wbl(self, chem_potentials, temperature):
        etas=[]
        gammas=[]  
        for i in range(len(chem_potentials)):
            etas_sub=[]
            for j in range(len(self.Xis_Fermi)):
                etas_sub.append(-1j * 2*np.pi * k_B*temperature * self.kappas_Fermi[j])
            etas.append(etas_sub)

            gammas_sub=[]
            for j in range(len(self.Xis_Fermi)):
                gammas_sub.append([self.Xis_Fermi[j]*k_B*temperature-1j*chem_potentials[i],self.Xis_Fermi[j]*k_B*temperature+1j*chem_potentials[i]])
            gammas.append(gammas_sub)
        return np.array(etas), np.array(gammas)

    def get_etas_gammas_Bose(self, couplings, broadenings, temperature):
        etas=[]
        gammas=[]  
        K_trunc=[]  
        for i in range(len(couplings)):
            etas.append( couplings[i] * broadenings[i] * (1./np.tan(broadenings[i]/(2.*k_B*temperature))-1j))
            gammas.append(broadenings[i])
            K_trunc_sub = 2. * couplings[i] * k_B*temperature / broadenings[i] - couplings[i] * 1./np.tan(broadenings[i]/(2.*k_B*temperature))
            for j in range(len(self.Xis_Bose)):
                etas.append(4.* couplings[i] * self.Xis_Bose[j] * self.kappas_Bose[j] * broadenings[i] /( (self.Xis_Bose[j])**2 -(broadenings[i]/(k_B*temperature))**2 ) )
                gammas.append(self.Xis_Bose[j]*k_B*temperature)
                K_trunc_sub -= 4.* couplings[i] * self.kappas_Bose[j] * broadenings[i] /(k_B*temperature)/( (self.Xis_Bose[j])**2 -(broadenings[i]/(k_B*temperature))**2 ) 
            K_trunc.append(K_trunc_sub)
        return np.array(etas), np.array(gammas),np.array(K_trunc)

    def get_etas_gammas_Bose_high_temperature(self, couplings, broadenings, temperature):
        etas=[]
        gammas=[]  
        K_trunc_T=np.zeros(len(couplings),float)
        for i in range(len(couplings)):
            etas.append( couplings[i] * broadenings[i] * (2.*k_B*temperature/broadenings[i]-1j))
            gammas.append(broadenings[i])
        return np.array(etas), np.array(gammas),np.array(K_trunc)


if __name__=='__main__':
    number_pade_modes=5
    spec='Bose'
    pade=Pade_decomposition(number_pade_modes,spec)
    KAPPA,XI=pade.get_Pade_parameters()
    for i in range(len(KAPPA)):
        print(i,KAPPA[i],XI[i])
    if 'Fermi' in spec or 'fermi' in spec:
        print(pade.Fermi_V(energy=1.e-2,temperature=50.))
        print(pade.Fermi_Pade_Approx(energy=1.e-2,temperature=50.))
    elif 'Bose' in spec or 'bose' in spec:
        print(pade.Bose_V(energy=1.e-2,temperature=50.))
        print(pade.Bose_Pade_Approx(energy=1.e-2,temperature=50.))



