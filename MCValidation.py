# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 11:02:21 2022

I go, factor by factor, calculating the transition probability P_ij, corres-
ponding to the epidemic progression spanning from the chain binomial system of 
Adams et. al. 

@author: Alejandro Alarcón and Guillermo Pérez
"""

import numpy as np
from scipy.stats import binom
import os
os.chdir("C:/Users/Alejandro/Documents/Research/Antonio's internship")
from functions import *


##########Parameters
h = 1/24
GAMMA = 1/1.373
p = np.repeat(0.5,10) 
THETA = 1/2.108
PSI = 1.00
PHI_1 = 0.75 
OMEGA = 1.00       
DELTA_1 = 1/4.176
DELTA_2 = 1 
DELTA_3 = 1    
DELTA_4 = 1
TAU_1 = 1
TAU_2 = 1    
SIGMA_asym = 1/6.284

#Population size
N = 100                                      
#Number of age classes
K = 10                                         

#Number of individuals per compartment at source state
m = np.stack((source_state[0,], source_state[1,], source_state[2,], source_state[3,], source_state[4,], source_state[5,], source_state[6,], source_state[7,], source_state[8,], source_state[9,]) , axis = 0)

#Number of individuals per compartment at destination state
n = np.stack((j_state[0,], j_state[1,], j_state[2,], j_state[3,], j_state[4,], j_state[5,], j_state[6,], j_state[7,], j_state[8,], j_state[9,]) , axis = 0)


#########################A_1################################

#Defining p_star
p_star = np.zeros(10) 

#Calling the contact matrices
C_asym = C_asym()
C_sym = C_sym()
q_sym=1
q_asym=1


#This holds true if broadcasting holds true
for k in range(K):
    p_star[k] = np.sum(np.multiply(q_asym*C_asym[k,], I_PreSym + I_Asym) + \
                       np.multiply(q_sym*C_sym[k,],I_Mild + I_Sev))

p_star = 1-np.exp(-h*p_star)


#The next one can be calculated for values of n_1 running from zero to m_1
#Use binom.pmf(k,n,p) to get probability of k successes out of n trials
#I am relying on broadcasting again
a_1 = binom.pmf(n[0],m[0],1-p_star)

#Computing the product of probabilities
a_1 = np.prod(a_1)



#########################A_2################################

#Defining p_2
p_2 = 1-np.exp(-h*GAMMA)

assert(n[0]+n[1]-m[0] >= 0 and n[0]+n[1]-m[0] <= m[1])
a_2 = binom.pmf(n[0]+n[1]-m[0],m[1],p_2)
a_2 = np.prod(a_2)


#########################A_3################################

#Defining p_3
p_3 = 1-np.exp(-h*p*THETA)
p_3_comp = 1-np.exp(-h*(1-p)*THETA)

#First check point
#If 0 ≤ k1 ≤ m3:
k_1 = m[0]+m[1]+m[2]-n[0]-n[1]-n[2]
print(k_1)

a_3 = np.zeros(10)
#Actually, I am making a sum per age group
for k in range(K):
    if k_1[k] >= 0 and k_1[k] <= m[2,k]:
        for i in range(k_1[k]+1):  #I need to sum one because of python indexing 
            a_3[k] += binom.pmf(i,m[2,k],p_3[k])*binom.pmf(k_1[k]-i,m[2,k]-i,p_3_comp[k])

a_3 = np.prod(a_3)


#########################A_4################################

#Defining p_4
p_4 = 1-np.exp(-h*DELTA_1)

#To be used in case A_4-i: If k1 ≥ n4 − m4 ≥ 0 and 0 ≤ k1 ≤ m3
k_2 = np.zeros(10)
for k in range(K):    
    k_2[k] = min(k_1[k]-(n[3,k]-m[3,k]), m[3,k])
    
#To be used in case A_4-ii: If n4 − m4 < 0 
k_3 = np.zeros(K)
for k in range(K):    
    k_3[k] = min(k_1[k], n[3,k])    
    

#Integrating the two possibilities in one loop
a_4 = np.zeros(K)
for k in range(K):
    #Case i) If k1 ≥ n4 − m4 ≥ 0 and 0 ≤ k1 ≤ m3
    if n[3,k]-m[3,k] >= 0 and n[3,k]-m[3,k] <= k_1[k]:
        for i in range(int(k_2[k])+1):
            a_4[k] += binom.pmf(n[3,k]-m[3,k]+i,m[2,k],p_3[k])*binom.pmf(i,m[3,k],p_4)
    #Case ii) If n4 − m4 < 0 
    if n[3,k]-m[3,k] < 0:
        for i in range(int(k_3[k])+1):
            a_4[k] += binom.pmf(i,m[2,k],p_3[k])*binom.pmf(m[3,k]-n[3,k]+i,m[3,k],p_4)
    
a_4 = np.prod(a_4) 


#########################A_5################################

#Defining p_5
p_5 = 1-np.exp(-h*PSI)
p_6 = 1-np.exp(-h*DELTA_2)



#This is good
a_5 = np.zeros(K)    
for k in range(K):
    #Case A_4-i: 0\leq n_4-m_4\leq k_1        
    if n[3,k]-m[3,k] >= 0 and n[3,k]-m[3,k] <= k_1[k]:
        for i in range(n[3,k]-m[3,k], int(n[3,k]-m[3,k]+k_2[k]+1)): 
            print(i)
            p = m[4,k]-n[4,k]+k_1[k]-i
            #We care for values of p such that 0<=p<=m_5. Otherwise, we get zero 
            #probabiliy, which is consistent with a_5 = np.zeros(K)    
            if p>=0 and p<=m[4,k]:
                print("hello"+str(i))
                #Computing the inner sum. Adding a one, cause it's open interval
                for j in range(p+1):
                    print(j)
                    a_5[k] += binom.pmf(j,m[4,k],p_5)*binom.pmf(p-j,m[4,k]-j,p_6)
    else:
    #Case A_4-ii: n_4-m_4 < 0        
        assert(n[3,k]-m[3,k] < 0)
        for i in range(int(k_3[k])+1):
            p = m[4,k]-n[4,k]+k_1[k]-i
            #We care for values of p such that 0<=p<=m_5. Otherwise, we get zero 
            #probabiliy, which is consistent with a_5 = np.zeros(K)    
            if p>=0 and p<=m[4,k]:
                #Computing the inner sum
                for j in range(p+1):
                    print(j)
                    a_5[k] += binom.pmf(j,m[4,k],p_5)*binom.pmf(p-j,m[4,k]-j,p_6)
a_5 = np.prod(a_5)
    
                        
#########################A_6################################

#Defining p_5
p_5 = 1-np.exp(-h*PSI)
p_6 = 1-np.exp(-h*DELTA_2)
p_7 = 1-np.exp(-h*PHI_1*OMEGA)
p_7_comp = 1-np.exp(-h*(1-PHI_1)*OMEGA)


#Fourth check point
#Assesing which case it belongs to
print(n[5]-m[5])


#Case i) If k3 ≥ n6 − m6 ≥ 0 and 0\leq n_5-m_5\leq k_1    
k_5 = np.zeros(K)
for k in range(K):    
    k_5[k] = min(k_3[k]-(n[5,k]-m[5,k]), m[5,k])


a_6 = np.zeros(K)    
for k in range(K):
    for i in range(int(k_5[k])+1):
        print("hello"+str(i))
        #First, calculating in a the new individuals in I_5 and I_6
        a = np.zeros(K)
        for j in range(i+1):
            a[k] += binom.pmf(j,m[5,k],p_7)*binom.pmf(i-j,m[5,k]-j,p_7_comp)
        a_6[k] += binom.pmf(n[5,k]-m[5,k]+i,m[4,k],p_5)*a[k]   

a_6 = np.prod(a_6)


#This following part has not been checked
#Case iii) If n_6 − m_6 < 0:

k_6 = min(m_5,n_6)
limit = m_6-n_6+i
    
for i in range(k_6):
    a = 0
    for j in range(limit):
        a += binom.pmf(j,m_6,p_7)*binom.pmf(m_6-n_6+i-j,m_6-j,p_7_comp)
    A_6 += binom.pmf(i,m_5,p_5)*a
    
    
#########################A_7################################
p_7 = 1-np.exp(-h*PHI_1*OMEGA)    
p_8 = 1-np.exp(-h*DELTA_3)
p_9 = 1-np.exp(-h*TAU_1)    


#Fifth check point
#Assesing which case it belongs to
print(n[6]-m[6])


#Case i) If k_5 ≥ n7 − m7 ≥ 0 and 0\leq n_6-m_6\leq k_3:   
k_7 = np.zeros(K)
for k in range(K):    
    k_7[k] = min(k_5[k]-(n[6,k]-m[6,k]), m[6,k])


a_7 = np.zeros(K)    
for k in range(K):
    for i in range(int(k_7[k])+1):
        print("hello"+str(i))
        #First, calculating in a the new individuals in D_1' and R_3'
        a = np.zeros(K)
        for j in range(i+1):
            a[k] += binom.pmf(j,m[6,k],p_8)*binom.pmf(i-j,m[6,k]-j,p_9)
        a_7[k] += binom.pmf(n[6,k]-m[6,k]+i,m[5,k],p_7)*a[k]   

a_7 = np.prod(a_7)


#This following part has not been checked
#Case iii) If n7 − m7 < 0

k_8 = min(m_6,n_7)
    
for i in range(k_8):
    a = 0
    for j in range(i):
        a += binom.pmf(j,m_7,p_8)*binom.pmf(m_7-n_7+i-j,m_7-j,p_9)
    A_7 += binom.pmf(i,m_6,p_7)*a   


#########################A_8################################
p_7 = 1-np.exp(-h*PHI_1*OMEGA)    
p_7_comp = 1-np.exp(-h*(1-PHI_1)*OMEGA)
p_8 = 1-np.exp(-h*DELTA_3)
p_9 = 1-np.exp(-h*TAU_1)    
p_10 = 1-np.exp(-h*DELTA_4)
p_11 = 1-np.exp(-h*TAU_2)    

#Sixth check point
#Assesing which case it belongs to
print(n[7]-m[7])


#Case i) If k_5 ≥ n8 − m8 ≥ 0 and 0\leq n_6-m_6\leq k_3:   
k_9 = np.zeros(K)
for k in range(K):    
    k_9[k] = min(k_5[k]-(n[7,k]-m[7,k]), m[7,k])



a_8 = np.zeros(K)    
for k in range(K):
    for i in range(int(k_5[k])+1):
        print(i)
        for j in range(int(k_9[k])+1):
            print(j)
            #We care for ordinary factorials only 
            if n[7,k]-m[7,k]+j <= m[5,k]-i:
                print("hello"+str(i))
                #First, calculating in a the new individuals in R_4 and D_2
                a = np.zeros(K)
                for l in range(j+1):
                    a[k] += binom.pmf(l,m[7,k],p_10)*binom.pmf(j-l,m[7,k]-l,p_11)
                a_8[k] += binom.pmf(n[7,k]-m[7,k]+j,m[5,k]-i,p_7_comp)*a[k]   
                
a_8 = np.prod(a_8)
 

#This following part has yet to be checked    
#Case ii) If n8 − m8 < 0:    

k_10 = min(m_6-I_5_prima,n_8)
limit = m_8-n_8+i

    
for i in range(k_10):
    a = 0
    limit = m_8-n_8+i        #I cannot recycle i, because it will be used later on 
    for j in range(limit):
        a += binom.pmf(j,m_8,p_10)*binom.pmf(m_8-n_8+i-j,m_8-j,p_11)
    A_8 += binom.pmf(i,m_6-I_5_prima,p_7_comp)*a   
    


#########################A_9################################




#########################A_10################################
    
    
