import sys,os
import numpy as np
from constants import * 
from input_param import *
import Pade 
import hierarchy_fermi
import hierarchy_bose
import Hamiltonian
import sparsity
from time import perf_counter, asctime

t_start_program = asctime()
output_info = open( 'simulation_info.txt', 'w')
output_info.write( 'The simulation starts at '+str(t_start_program)+'\n')
t_end = perf_counter()

if High_temperature_approximation:
    etas_Bose,gammas_Bose,K_trunc_T=Pade_Bose.get_etas_gammas_Bose_high_temperature(lambda_T, wc, temperature_Bose)
else:
    Pade_Bose=Pade.Pade_decomposition(Number_pole_Bose,'Bose')
    KAPPA_Bose,XI_Bose=Pade_Bose.get_Pade_parameters()
    etas_Bose,gammas_Bose,K_trunc_T=Pade_Bose.get_etas_gammas_Bose(lambda_T, wc, temperature_Bose)
N_hierarchy=hierarchy_bose.aw_number(Truncation_tier_Bose,N_index)
indices_Bose,index_plus_Bose,index_minus_Bose=hierarchy_bose.set_aw(N_index,Truncation_tier_Bose,N_hierarchy)

Pade_Fermi=Pade.Pade_decomposition(Number_pole_Fermi,'Fermi')
KAPPA_Fermi,XI_Fermi=Pade_Fermi.get_Pade_parameters()
if (WBL_on):
    etas_Fermi,gammas_Fermi=Pade_Fermi.get_etas_gammas_Fermi_wbl(chem_potentials, temperature_Fermi)
else:
    etas_Fermi,gammas_Fermi=Pade_Fermi.get_etas_gammas_Fermi(chem_potentials, temperature_Fermi, broadenings)
index_fermi=hierarchy_fermi.Hierarchy_index(Number_lead,Number_el_state,Number_pole_Fermi,Sign_pm,Truncation_tier_Fermi,WBL_on)

number_mode,modes=index_fermi.get_modes()
if  threshold_on :
    nth_tier_Fermi,indices_Fermi,index_plus_Fermi,index_minus_Fermi=index_fermi.get_filtered_indices_sym(threshold_value,etas_Fermi,gammas_Fermi,max_el_lead_coupling)
else:
    nth_tier_Fermi,indices_Fermi,index_plus_Fermi,index_minus_Fermi=index_fermi.get_indices_sym()

t_start = t_end
t_end = perf_counter()
output_info.write("Elapsed time of indices generation: " + str(t_end-t_start) +'\n')

systems=Hamiltonian.nuclear_DVR(Number_el_state,mass,N_DVR,x_0,x_N,x_infties,CAP_thresholds,CAP_pots,energies,el_couplings,CAP_button='on')
dim_DVR=systems.N_DVR
hamil,hamil_CAP=systems.get_operators()
V_DVR=systems.get_V_DVR(el_lead_couplings,Number_lead)
Vd_B=systems.get_V_DVR_thermal(el_bath_couplings,Number_thermal_bath)
for i in range(Number_thermal_bath):
    hamil+= np.diag(lambda_T*Vd_B[i,:]**2)
rho_0_DVR=systems.get_rho_0(rho_0_el)
cap=systems.get_CAP()

d_operators_el,hamil_el=systems.get_electronic_operators(energies_el,el_couplings_el)
d_operators_logical=np.array(d_operators_el,dtype=bool)
hamil_logical=np.array(hamil_el,dtype=bool)
rho_0_logical=np.array(rho_0_el,dtype=bool)

rho_sparsity,number_nonzeros=sparsity.sparsity(modes, nth_tier_Fermi, index_minus_Fermi, index_plus_Fermi, d_operators_logical, hamil_logical, rho_0_logical, dim_fockspace, len(indices_Fermi), len(index_plus_Fermi), Truncation_tier_Fermi, number_mode, Number_el_state)

if (WBL_on):
    rho_nonzeros,number_pair,number_pair_hamil_l,number_pair_hamil_r,number_pair_wbl=sparsity.sparse_matrix_elements_a_wbl(modes, nth_tier_Fermi, index_minus_Fermi, index_plus_Fermi, d_operators_el, hamil_el, rho_sparsity, number_nonzeros, dim_fockspace, len(indices_Fermi), len(index_plus_Fermi), Truncation_tier_Fermi, number_mode, Number_el_state)
    matrix_pair,matrix_pair_hamil_l,matrix_pair_hamil_r,matrix_pair_wbl,matrix_values_pair,matrix_gammas=sparsity.sparse_matrix_elements_b_wbl(modes, nth_tier_Fermi, indices_Fermi, index_minus_Fermi, index_plus_Fermi, gammas_Fermi, etas_Fermi, d_operators_el, hamil_el, rho_sparsity, rho_nonzeros, number_pair, number_pair_hamil_l, number_pair_hamil_r, number_pair_wbl, number_nonzeros, dim_fockspace, len(indices_Fermi), len(index_plus_Fermi), Truncation_tier_Fermi, number_mode, Number_lead, Number_el_state, Number_pole_Fermi, Sign_pm)
else:
    rho_nonzeros,number_pair,number_pair_hamil_l,number_pair_hamil_r=sparsity.sparse_matrix_elements_a(modes, nth_tier_Fermi, index_minus_Fermi, index_plus_Fermi, d_operators_el, hamil_el, rho_sparsity, number_nonzeros, dim_fockspace, len(indices_Fermi), len(index_plus_Fermi), Truncation_tier_Fermi, number_mode, Number_el_state)
    matrix_pair,matrix_pair_hamil_l,matrix_pair_hamil_r,matrix_values_pair,matrix_gammas=sparsity.sparse_matrix_elements_b(modes, nth_tier_Fermi, indices_Fermi, index_minus_Fermi, index_plus_Fermi, gammas_Fermi, etas_Fermi, d_operators_el, hamil_el, rho_sparsity, rho_nonzeros, number_pair, number_pair_hamil_l,number_pair_hamil_r, number_nonzeros, dim_fockspace, len(indices_Fermi), len(index_plus_Fermi), Truncation_tier_Fermi, number_mode, Number_lead, Number_el_state, Number_pole_Fermi, Sign_pm)

t_start = t_end
t_end = perf_counter()
output_info.write("Elapsed time sparsity: " + str(t_end-t_start) +'\n')

with open ('index_info.txt','w') as f:
    f.write('i\tN_lead(0=L,1=R)\tN_elec\tN_pole\t+-\n')
    for i in range(number_mode):
        f.write(str(i)+'\t')
        for j in range(4):
            f.write(str(modes[i][j])+"\t")
        f.write('\n')

    f.write('i\tnth_tier\t')
    for j in range(Truncation_tier_Fermi):
        f.write('index[i]['+str(j)+']'+'\t')
    f.write('\n')
    for i in range(len(nth_tier_Fermi)):
        f.write(str(i)+"\t"+str(nth_tier_Fermi[i])+"\t")
        for j in range(Truncation_tier_Fermi):
            f.write(str(indices_Fermi[i][j])+'\t')
        f.write('\n')

    f.write("thermal bath index information\n")
    f.write("i\t( j_1 ... j_n ... j_N )\t for n = 1, N_index\n")
    for i in range(N_hierarchy+1):
        f.write(str(i)+"\t")
        for j in range(N_index):
            f.write(str(indices_Bose[j][i])+' ')
        f.write('\n')

    f.write('The first-tier indices for the heat flow\n')
    for i in range(N_hierarchy+1):
        if sum(indices_Bose[:,i]) == 1:
            f.write(str(i)+' '+str(indices_Bose[:,i].argmax()+1))
            f.write('\n')

count_plus_T=0
count_minus_T=0
with open('matrix_thermal.txt','w') as f:
	f.write('i_n\t'+'n_index\t'+'i_n^+\n')
	for i in range(N_hierarchy+1):
		for j in range(N_index):
			if ( index_plus_Bose[j][i] != -1 ):
				count_plus_T  += 1
				f.write(str(i)+'\t'+str(index_plus_Bose[j][i])+'\t'+str(j+1)+'\n')
	f.write('i_n\t'+'n_index\t'+'i_n^-\n')
	for i in range(N_hierarchy+1):
		for j in range(N_index):
			if ( index_minus_Bose[j][i] != -1 ):
				count_minus_T += 1
				f.write(str(i)+'\t'+str(index_minus_Bose[j][i])+'\t'+str(j+1)+'\n')

with open('rho_nonzeros.txt','w') as f:
    f.write('i\ti_hierarchy\trow\tcol\n')
    for j in range(number_nonzeros):
        f.write(str(rho_nonzeros[j,0])+'\t'+str(rho_nonzeros[j,1])+'\t'+str(rho_nonzeros[j,2])+'\n')
with open('matrix.txt','w') as f:
    f.write('i_new\t i_old\t i_sym(0,1)\t i_leftright(0=left,1=right)\t i_lead\t i_elec\n')
    for i in range(number_pair):
        f.write(str(matrix_pair[i,0])+'\t'+str(matrix_pair[i,1])+'\t'+str(matrix_pair[i,2])+'\t'+str(matrix_pair[i,3])+'\t'+str(matrix_pair[i,4])+'\t'+str(matrix_pair[i,5])+'\n')
with open('matrix_hamil.txt','w') as f:
    f.write('left: i_new\t i_old\t hamil_el(i_row,i_column)\n')
    for i in range(number_pair_hamil_l):
        f.write(str(matrix_pair_hamil_l[i,0])+'\t'+str(matrix_pair_hamil_l[i,1])+'\t'+str(matrix_pair_hamil_l[i,2])+'\t'+str(matrix_pair_hamil_l[i,3])+'\n')
    f.write('right: i_new\t i_old\t hamil_el(i_row,i_column)\n')
    for i in range(number_pair_hamil_r):
        f.write(str(matrix_pair_hamil_r[i,0])+'\t'+str(matrix_pair_hamil_r[i,1])+'\t'+str(matrix_pair_hamil_r[i,2])+'\t'+str(matrix_pair_hamil_r[i,3])+'\n')
if (WBL_on):
    with open('matrix_wbl.txt','w') as f:
        f.write('i_new\t i_old\t i_elec\t i_tier\t i_sign (-1,1)\n')
        for i in range(number_pair_wbl):
            f.write(str(matrix_pair_wbl[i,0])+'\t'+str(matrix_pair_wbl[i,1])+'\t'+str(matrix_pair_wbl[i,2])+'\t'+str(matrix_pair_wbl[i,3])+'\t'+str(matrix_pair_wbl[i,4])+'\n')
with open('matrix_value.txt','w') as f:
    f.write('matrix_gammas\n')
    for i in range(number_nonzeros):
        f.write('('+str(matrix_gammas[i].real)+','+str(matrix_gammas[i].imag)+')\t')
    f.write('\nmatrix_values_pair\n')
    for i in range(number_pair):
        f.write('('+str(matrix_values_pair[i].real)+','+str(matrix_values_pair[i].imag)+')\t')
    f.write('\n')

with open('system_info.txt','w') as f:
    f.write('0:bath-1\t0:pole\tk-th\tgammas\tetas\n')
    k = 0                                           
    for i in range(Number_thermal_bath):
        for j in range(Number_pole_Bose+1):
            f.write(str(i)+'\t'+str(j)+'\t'+str(k)+'\t('+str(gammas_Bose[k].real)+','+str(gammas_Bose[k].imag)+')\t')
            f.write('('+str(etas_Bose[k].real)+','+str(etas_Bose[k].imag)+')\n')
            k += 1                        

    f.write('0:bath-1\tK_trunc_T\n')
    for i in range(Number_thermal_bath):
        f.write(str(i)+'\t'+str(K_trunc_T[i])+'\n')
    for i in range(Number_el_state):
        f.write('d_dagger for electronic DOF:\t'+str(i+1)+'\n')
        for j1 in range(dim_fockspace):
            for j2 in range(dim_fockspace):
                f.write(str(d_operators_el[0,i,j1,j2])+'\t')
            f.write('\n')
        f.write('d for electronic DOF:\t'+str(i+1)+'\n')
        for j1 in range(dim_fockspace):
            for j2 in range(dim_fockspace):
                f.write(str(d_operators_el[1,i,j1,j2])+'\t')
            f.write('\n')

    f.write('Hamiltonian:\n')
    for j1 in range(dim_fockspace*dim_DVR):
        for j2 in range(dim_fockspace*dim_DVR):
            f.write(str(hamil[j1,j2].real)+'\t')
            #f.write("("+str(hamil[j1,j2].real)+','+str(hamil[j1,j2].imag)+')\t')
        f.write('\n')

    f.write('rho_0_DVR:\n')
    for j1 in range(dim_fockspace*dim_DVR):
        for j2 in range(dim_fockspace*dim_DVR):
            f.write("("+str(rho_0_DVR[j1,j2].real)+','+str(rho_0_DVR[j1,j2].imag)+')\t')
        f.write('\n')

    f.write('cap potential:\n')
    for fock_state in range(dim_fockspace):
        for i in range(dim_DVR):
            f.write(str(cap[fock_state,i])+'\t')
        f.write('\n')

    for i in range(Number_lead):
        for j in range(Number_el_state):
            f.write('V_DVR('+str(i+1)+","+str(j+1)+',:)\n')
            for k in range(dim_DVR):
                f.write(str(V_DVR[i,j,k])+'\t')
            f.write('\n')

    for i in range(Number_thermal_bath):
        f.write('Vd_B('+str(i+1)+',:,:)\n')
        for k1 in range(dim_fockspace*dim_DVR):
            f.write(str(Vd_B[i,k1])+'\t')
        f.write('\n')

output_info.write("number of leads = "+str(Number_lead)+"\n")
output_info.write("number of thermal baths = "+str(Number_thermal_bath)+"\n")
output_info.write("number of electronic states = "+str(Number_el_state)+"\n")
output_info.write("number of Fermi Pade poles = "+str(Number_pole_Fermi)+"\n")
output_info.write("number of Bose Pade poles = "+str(Number_pole_Bose)+"\n")
output_info.write("number of Bose indices = "+str(N_index)+"\n")
output_info.write("Sign pm = "+str(Sign_pm)+"\n")
output_info.write("dimension of fockspace 2^(number_el_state) = "+str(dim_fockspace)+"\n")
output_info.write("Fermi truncation tier = "+str(Truncation_tier_Fermi)+"\n")
output_info.write("Bose truncation tier = "+str(Truncation_tier_Bose)+"\n")
output_info.write("Bose hierarchy = "+str(N_hierarchy)+"\n")
output_info.write("ishizaki truncation = "+str(Ishizaki_trunc_on)+'\n')
output_info.write("wide band limit = "+str(WBL_on)+"\n")
output_info.write("Gamma_L [eV] = %.3f\n"%(el_lead_coupling_left**2*2.0*np.pi))
output_info.write("Gamma_R [eV] = %.3f\n"%(el_lead_coupling_right**2*2.0*np.pi))
output_info.write("time step (dt,[fs])= "+str(dt)+"\n")
output_info.write("time outout step (dt_out,[fs])= "+str(dt_out)+"\n")
output_info.write("Initial propagation = "+str(Init_propagation)+'\n')
if not Init_propagation:
	output_info.write("time intermediate (dt,[fs])= "+str(time_intermediate)+"\n")
output_info.write("number of time steps = "+str(nsteps)+" \n")
output_info.write("number of DVR grid points = "+str(dim_DVR)+"\n")
output_info.write("**********hierarchy parameters **********\n")
output_info.write('number of mode = '+str(number_mode)+'\n')
output_info.write('number of hierarchy = '+str(len(indices_Fermi))+'\n')
tiers_string=["1st tier: ", "2nd tier: ","3rd tier: ","4th tier: ", "5th tier: ","6th tier: ","7th tier: ", "8th tier: "]
for i in range(Truncation_tier_Fermi):
    output_info.write("**{ "+tiers_string[i]+str(np.sum(np.array(nth_tier_Fermi)==(i+1)))+" }**\n")
output_info.write('number of nonzero entries = '+str(number_nonzeros)+'\n')
output_info.write("number_pair_hierarchy = "+str(number_pair)+'\n')
output_info.write("number_pair_hamil_l = "+str(number_pair_hamil_l)+'\n')
output_info.write("number_pair_hamil_r = "+str(number_pair_hamil_r)+'\n')
output_info.write("count_plus_Bose = "+str(count_plus_T)+'\n')
output_info.write("count_minus_Bose = "+str(count_minus_T)+'\n')
if (WBL_on):
    output_info.write("number_pair_wbl = "+str(number_pair_wbl)+'\n')
memory_1=number_pair*16/1024
memory_2=(number_nonzeros*dim_DVR**2*3)*16/1024
output_info.write('memory of coefficients : %.3f [KB] (%.3f [MB]) \n'%(memory_1,memory_1/1024))
output_info.write("*****************************************\n")
output_info.write("memory of ADOs in RK4 : %.3f [KB] (%.3f [MB]) \n"%(memory_2,memory_2/1024))

t_end = asctime()
output_info.write('The simulation ends at '+str(t_end)+'\n')
output_info.close()
