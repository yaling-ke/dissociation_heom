import numpy as np
import itertools
import sys

class Hierarchy_index():

    def __init__(self,number_lead,number_el_state,number_pole,sign_pm,truncation_tier,WBL_on):
        """
        Generating a index list indices[],indices_plus[],indices_minus[] 
        """
        self.number_lead=number_lead
        self.number_el_state=number_el_state
        self.number_pole=number_pole
        self.sign_pm=sign_pm
        self.truncation_tier=truncation_tier
        self.WBL_on=WBL_on
        self.index_mode()
        self.indices_sym_generation()

    def index_mode(self):
        """
        Generating the index_mode: mode in {[0:number_lead,0:number_el_state,0:number_pole+1,0:sign_pm]}
        i.e., number_lead=2(left:0, right:1),    number_el_state=1,  number_pole=1,  sign_pm=2(0:+, 1:-)
        modes: 0:[0,0,0,0]     1:[0,0,0,1]     2:[0,0,1,0]     3:[0,0,1,1]
               4:[1,0,0,0]     5:[1,0,0,1]     6:[1,0,1,0]     7:[1,0,1,1]
        number_mode=8
        """
        leads=np.arange(self.number_lead)
        el_states=np.arange(self.number_el_state)
        if self.WBL_on:
            poles=np.arange(self.number_pole)
        else:
            poles=np.arange(self.number_pole+1)
        sign=np.arange(self.sign_pm)
    
        self.modes=[]
        for elem in itertools.product(leads,el_states,poles,sign):
            self.modes.append(list(elem))
        self.number_mode=len(self.modes)
    
    def indices_sym_generation(self):
        """
        Generating the index list:
        i.e., number_lead=2(left:0, right:1),    number_el_state=1,  number_pole=1,  sign_pm=2(0:+, 1:-)
        modes:0,1,2,3,4,5,6,7
        Truncation_tier: L=4
        0th_tier:indices[0]=[-1,-1,-1,-1]  !!!len(indices[0])=L
                 nth_tier[0]=0
        1th_tier:indices[1]=[0,-1,-1,-1]     indices[2]=[2,-1,-1,-1]     indices[3]=[4,-1,-1,-1]     indices[4]=[3,-1,-1,-1]
                 nth_tier[1:4]=1
        2th_tier:indices[5]=[0,1,-1,-1]      indices[6]=[0,2,-1,-1]      ......                      indices[11]=[0,7,-1,-1]
                 indices[12]=[2,3,-1,-1]     indices[13]=[2,4,-1,-1]     .....                       indices[16]=[2,7,-1,-1]
                 ...                         ...                         ...                         ...
                 indices[20]=[6,7,-1,-1]
                 nth_tier[5:20]=2
        3th_tier:indices[21]=[0,1,2,-1]      indices[22]=[0,1,3,-1]      ......                      indices[26]=[0,1,7,-1]
                 indices[27]=[0,2,3,-1]      ...                         ....                        indices[31]=[0,2,7,-1]
                 ...                         ...                         ...                         ...
                 indices[60]=[5,6,7,-1]
                 nth_tier[21:60]=3
        4th_tier:indices[61]=[0,1,2,3]       indices[62]=[0,1,2,4]       ......                      indices[65]=[0,1,2,7]
                 indices[66]=[0,1,3,4]       ...                         ....                        indices[69]=[0,1,3,7]
                 ...                         ...                         ...                         ...
                 indices[122]=[4,5,6,7]
                 nth_tier[61:122]=4

        for example, [1,4,-1,-1] is not included, because rho_{[1,4,-1,-1]}=(-1)^{2/2} * rho_{[0,5,-1,-1]

        --------------------------------------------------------------------------------------------------------------------------
        Finding the relation between nth_tier and (n+1)th_tier ADOs, and it is stored in index_plus[] 

        dot(rho^{(n)}_{j_1,...,j_n})=[[nth_tier]]+[[(n-1)th_tier]]+[[-i\sum_{j (!= j_1,...,j_n)} A_{bar{j}}*rho^{(n+1)}_{j_1,...,j_n,j}]]_{(n+1)th_tier}

        index_plus[s=[j1,j2,j3,-1]]=[sublist[0], sublist[1], sublist[2], sublist[3],..., sublist[number_mode-1]]  with 
                                              sublist[0]=[index_of_[0,j_1,j_2,j_3],permutation times, sym_turn_on=0 or 1]
                                              sublist[1]=[index_of_[1,j_1,j_2,j_3],permutation times, sym_turn_on=0 or 1]
                                              sublist[2]=[index_of_[j_1,2,j_2,j_3],permutation times, sym_turn_on=0 or 1]
                                              ...
                                              sublist[7]=[index_of_[j_1,j_2,j_3,7],permutation times, sym_turn_on=0 or 1]
          
        for example indices[50]=[2,4,6,-1]
                    sublist[0]=[81,3,0], indices[81]=[0,2,4,6] permutation times=3 from [2,4,6,0], sym_turn_on=0
                    sublist[1]=[90,3,1], indices[90]=[0,3,5,7] permutation times=3 from [2,4,6,1], which is (-1)^2*hermite_conjugate_of(rho^4_[1,2,4,6]), so sym_turn_on=1
                    sublist[2]=[-1,-1,-1]
                    sublist[3]=[109,2,0], indices[109]=[2,3,4,6] permutation times=2 from [2,4,6,3] sym_turn_on=0
                    sublist[4]=[-1,-1,-1]
                    sublist[5]=[114,1,0], indices[114]=[2,4,5,6] permutation times=1 from [2,4,6,5] sym_turn_on=0
                    sublist[6]=[-1,-1,-1] 
                    sublist[7]=[116,0,0], indices[116]=[2,4,6,7] permutation times=0 from [2,4,6,7] sym_turn_on=0
        so, index_plus[50]=[sublist[0], sublist[1], sublist[2], sublist[3],sublist[4], sublist[5], sublist[6], sublist[7]]

        --------------------------------------------------------------------------------------------------------------------------
        Finding the relation between nth_tier and (n-1)th_tier ADOs, and it is stored in index_minus[] 

        dot(rho^{(n)}_{j_1,...,j_n})=[[nth_tier]]+[[-i* \sum_{r=1}^{n} (-1)^{n-r}*C_{j_r}*rho^{(n-1)}_{j_1,...,j_{r-1},j_{r+1},...,j_n} ]](n-1)th_tier+[[(n+1)th_tier]]
                            index_minus[s=[j1,j2,j3,j4]]=[sublist[0], sublist[1], sublist[2], sublist[3]]  with 
                                              sublist[0]=[j_1, index_of_[j_2,j_3,j_4,-1], permutation_times, sym_turn_on=0 or 1]
                                              sublist[1]=[j_2, index_of_[j_1,j_3,j_4,-1], permutation_times, sym_turn_on=0 or 1]
                                              sublist[2]=[j_3, index_of_[j_1,j_2,j_4,-1], permutation_times, sym_turn_on=0 or 1]
                                              sublist[3]=[j_4, index_of_[j_1,j_2,j_3,-1], permutation_times, sym_turn_on=0 or 1]
          
        for example indices[87]=[0,3,4,6]
                    sublist[0]=[0,53,0, 1], indices[53]=[2,5,7,-1] is -1*hermite_conjugate_of(rho^3_[3,4,6,-1]), so sym_turn_on=1
                    sublist[1]=[3,37,0, 0], indices[37]=[0,4,6,-1], so sym_turn_on=0
                    sublist[2]=[4,34,0, 0], indices[34]=[0,3,6,-1], so sym_turn_on=0
                    sublist[3]=[6,32,0, 0], indices[32]=[0,3,4,-1], so sym_turn_on=0
        so, index_minus[87]=[sublist[0], sublist[1], sublist[2], sublist[3]]
        """
    
        self.indices_sym = []
        self.nth_tier_sym = []
        self.index_plus_sym = []
        self.index_minus_sym = []
        tier_range_sym = [0]

        i_tier = 0
        appendix = [-1] * self.truncation_tier
        self.indices_sym.append(appendix)
        self.nth_tier_sym.append(i_tier)
        tier_range_sym.append(1)
        sublist_minus=[]
        for i in range(self.truncation_tier):
            sublist_minus.append([-1]*4)
        self.index_minus_sym.append(sublist_minus)
    
        i_tier = 1
        appendix = [-1] * (self.truncation_tier - i_tier)
        i_count = 1
        sublist_plus = []
        for i_mode in range(self.number_mode):
            if i_mode%2 == 0:
                i_count += 1
                self.indices_sym.append([i_mode] + appendix)
                self.nth_tier_sym.append(i_tier)
                sublist_plus.append([ i_mode//2 + 1, 0, 0])
            else:
                sublist_plus.append([ (i_mode+1)//2, 0, 1])
        self.index_plus_sym.append(sublist_plus)
        tier_range_sym.append(i_count)
        for i_index in range(tier_range_sym[i_tier], tier_range_sym[i_tier+1]):
            sublist_minus=[]
            sublist_minus.append([self.indices_sym[i_index][0],0,0,0])
            for i in range(i_tier,self.truncation_tier):
                sublist_minus.append([-1]*4)
            self.index_minus_sym.append(sublist_minus)

        while ( i_tier < self.truncation_tier):
            i_tier += 1
            appendix = [-1] * (self.truncation_tier - i_tier)
            for i_index in range(tier_range_sym[i_tier-1], tier_range_sym[i_tier]):
                sublist_plus = []
                for i_mode in range(self.number_mode):
                    elem = self.indices_sym[i_index][i_tier-2]
                    if (i_mode > elem):
                        sublist_plus.append([i_count,0,0])
                        i_count += 1
                        self.indices_sym.append( self.indices_sym[i_index][:i_tier-1] + [i_mode] +  appendix)
                        self.nth_tier_sym.append(i_tier)
                    else:
                        if ( self.elem_in_index(i_mode, self.indices_sym[i_index][:i_tier-1])):
                            sublist_plus.append([-1, -1, -1])
                        else:
                            temp=self.indices_sym[i_index][:i_tier-1].copy() + [i_mode]
                            temp.sort()
                            permutation_number = i_tier-1-temp.index(i_mode)
                            sym_on=0
                            if (permutation_number == (i_tier-1) ) and (i_mode%2 == 1):
                                temp_sym=self.find_sym_index(temp)
                                temp=sorted(temp_sym)
                                permutation_number += ( np.sum(np.array(temp_sym) != np.array(temp)) // 2 )
                                sym_on=1

                            for i_pos in range(i_tier-1):
                                if (i_pos==0): 
                                    elem_1 = temp[0]//2+1
                                else:
                                    elem_1=hierarchy_plus
                                elem_2=temp[i_pos+1]
                                hierarchy_plus=self.index_plus_sym[elem_1][elem_2][0]
                            sublist_plus.append([hierarchy_plus, permutation_number, sym_on])
                self.index_plus_sym.append(sublist_plus)
            tier_range_sym.append(i_count)

            for i_index in range(tier_range_sym[i_tier], tier_range_sym[i_tier+1]):
                sublist_minus=[]
                for i in range(i_tier):
                    elem=self.indices_sym[i_index][i]
                    temp=self.indices_sym[i_index][:i_tier].copy()
                    temp.remove(elem)
                    sym_on = 0
                    permutation_number = 0
                    if ( temp[0]%2 != 0 ):
                        temp_sym = self.find_sym_index(temp)
                        temp = sorted(temp_sym)
                        sym_on = 1
                        permutation_number = ( np.sum(np.array(temp_sym) != np.array(temp)) // 2 )
                    hierarchy_minus = temp[0]//2+1
                    for i_pos in range(i_tier-2):
                        elem_1=hierarchy_minus
                        elem_2=temp[i_pos+1]
                        hierarchy_minus=self.index_plus_sym[elem_1][elem_2][0]
                    sublist_minus.append([elem, hierarchy_minus, permutation_number, sym_on])
                for i in range(i_tier,self.truncation_tier):
                    sublist_minus.append([-1]*4)
                self.index_minus_sym.append(sublist_minus)
    
    def filtering(self, index_in, tier, etas, gammas, max_el_lead_coupling):
        factor1=1.0
        for k in range(tier):
            i_mode=index_in[k]
            i_lead=self.modes[i_mode][0]
            i_pole=self.modes[i_mode][2]
            i_sign=self.modes[i_mode][3]
            factor1=factor1*etas[i_lead][i_pole]/gammas[i_lead][i_pole][i_sign].real 
        factor2=1.0
        for k in range(tier-1):
            factor2=factor2*max_el_lead_coupling**2
            fac_sum=0.0
            for m in range(k+1):
                i_mode=index_in[m]
                i_lead=self.modes[i_mode][0]
                i_pole=self.modes[i_mode][2]
                i_sign=self.modes[i_mode][3]
                fac_sum=fac_sum+gammas[i_lead][i_pole][i_sign].real
            factor2=factor2/fac_sum
        return abs(factor1*factor2)

    def threshold_reduce_indices_sym(self,threshold_value,etas,gammas,max_el_lead_coupling):
        index_length = len(self.indices_sym)
        threshold=np.full(index_length,True)
        correspondence=[]

        i_hierarchy=0
        for i in range(index_length):
            tier=self.nth_tier_sym[i]
            if tier < 2:
                correspondence.append(i_hierarchy)
                i_hierarchy += 1
            else:
                factor=self.filtering(self.indices_sym[i], tier, etas, gammas, max_el_lead_coupling) 
                if factor < threshold_value:
                    threshold[i]=False
                    correspondence.append(-1)
                else:
                    correspondence.append(i_hierarchy)
                    i_hierarchy += 1

        self.indices_sym_filtered=[]
        self.nth_tier_sym_filtered=[]
        self.index_minus_sym_filtered=[]
        self.index_plus_sym_filtered=[]
        for i in range(index_length):
            tier=self.nth_tier_sym[i]
            if threshold[i]==True:
                self.indices_sym_filtered.append( self.indices_sym[i] )
                self.nth_tier_sym_filtered.append( self.nth_tier_sym[i] )

                sub_index_minus_sym = self.index_minus_sym[i].copy()
                for i_tier in range(tier):
                    i_minus=sub_index_minus[i_tier][1]
                    if (i_minus != -1): sub_index_minus[i_tier][1]=correspondence[i_minus]
                self.index_minus_sym_filtered.append( sub_index_minus )

                if (tier < self.truncation_tier):
                    sub_index_plus = self.index_plus_sym[i].copy()
                    for i_mode in range(self.number_mode):
                        i_plus=sub_index_plus[i_mode][0]
                        if (i_plus != -1): sub_index_plus[i_mode][0]=correspondence[i_plus]
                    self.index_plus_sym_filtered.append( sub_index_plus )


    def elem_in_index(self,elem,index_original):
        """
        To check whether (elem) equals to an element in [index_original] or not; return True if it is in there
        """
        _in = False
        for i in range(len(index_original)):
            if index_original[i] == elem:
                _in = True
                break
        return _in

    def find_sym_index(self,index_original):
        """
        We can reduce the number of ADOs by utilizing the symmetry relation: 
                       rho_{j_1,j_2,....,j_n}=(-1)^{floor[n/2]} * rho_{bar{j_1},bar{j_2},...,bar{j_n}}.
        Because the latter is dependent on the former, we only need to include the first ADO in the propagation.

        It is observed that, if j_k is even, bar{j_k}=j_k+1, elif j_k is odd, bar{j_k}=j_k-1
        For example, index_original=[1,4,-1,-1] return index_sym=[0,5,-1,-1]
        """
        index_sym=[]
        for i in range(len(index_original)):
            if index_original[i] != -1:
                if index_original[i]%2==0:
                    index_sym.append(index_original[i]+1)
                else:
                    index_sym.append(index_original[i]-1)
            else:
                index_sym.append(-1)
        return index_sym

    def get_modes(self):
        return self.number_mode, np.array(self.modes)

    def get_indices_sym(self):
        return np.array(self.nth_tier_sym), np.array(self.indices_sym), np.array(self.index_plus_sym), np.array(self.index_minus_sym)

    def get_filtered_indices_sym(self,threshold_value,etas,gammas,max_el_lead_coupling):
        self.threshold_reduce_indices_sym(threshold_value,etas,gammas,max_el_lead_coupling) 
        return self.nth_tier_sym_filtered, self.indices_sym_filtered, self.index_plus_sym_filtered, self.index_minus_sym_filtered


if __name__=="__main__":
    Number_lead=2
    Number_pole=3 
    Number_el_state=1
    Sign_pm=2
    Truncation_depth=2

    hierarchy=Hierarchy_index(Number_lead,Number_el_state,Number_pole,Sign_pm,Truncation_depth,False)

    number_mode,modes=hierarchy.get_modes()
    output_modes=open('modes.txt','w')
    for i in range(number_mode):
        output_modes.write('index_modes['+str(i)+']\t=\t'+str(modes[i])+'\n')
    output_modes.close()        
    
    nth_tier, indices, index_plus, index_minus=hierarchy.get_indices_sym()
    output_indices=open('indices.txt','w')
    for i in range(len(indices)):
        output_indices.write('nth_tier['+str(i)+']\t= '+str(nth_tier[i])+'\tindices['+str(i)+']\t= '+str(indices[i])+'\n')
    output_indices.close()

    output_index_minus=open('index_minus.txt','w')
    for i in range(len(index_minus)):
        output_index_minus.write('nth_tier['+str(i)+']\t= '+str(nth_tier[i])+'\tindex_minus['+str(i)+']\t= '+str(index_minus[i])+'\n')
    output_index_minus.close()

    output_index_plus=open('index_plus.txt','w')
    for i in range(len(index_plus)):
        output_index_plus.write('nth_tier['+str(i)+']\t= '+str(nth_tier[i])+'\tindex_plus['+str(i)+']\t= '+str(index_plus[i])+'\n')
    output_index_plus.close()
