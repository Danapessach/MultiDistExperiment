
# coding: utf-8

# In[ ]:

my_id=2
print("ID:" + str(my_id))


# In[ ]:

DELAY=1

#this means that real numbers will be truncated to this number of places after decimal point
#this must be the same value for all parties
decimals=4


# In[ ]:

import numpy as np
import pandas as pd
import time
import shutil
import subprocess
import os
import ast
import sys

from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype

import time


# In[ ]:

os.getcwd()


# In[ ]:

public_path=os.getcwd()+"/MyDistExperiment/public"

#number of parties
with open(public_path+"/constants/parties.ini", 'r') as f:
    parties=int(f.read())
print("parties:",parties)
private_path=os.getcwd()+"/MyDistExperiment/private/" + str(parties) + "p"

player_file_path=private_path+'/party'+str(my_id)+'/player-' + str(my_id) + '.ini'

python2_path='C:\Anaconda2\envs\py_27\python.exe'


# In[ ]:

#read all public constants

with open(public_path+"/constants/repairs.ini", 'r') as f:
    repairs=ast.literal_eval(f.read())
print("repairs:",repairs)

with open(public_path+"/constants/dataset_name.ini", 'r') as f:
    dataset_name=f.read()
print("dataset_name:",dataset_name)

#lambda
with open(public_path+"/constants/lam.ini", 'r') as f:
    lams=ast.literal_eval(f.read())
    print("lambdas:",lams)

#non-sensitive features list
with open(public_path+"/constants/features.ini", 'r') as f:
    features=ast.literal_eval(f.read())
print("features:",features)
print("feature 1:",features[0])
print("#features:",len(features))

#ignore features list
with open(public_path+"/constants/ignore_features.ini", 'r') as f:
    ignore_features=ast.literal_eval(f.read())
print("ignore_features:",ignore_features)
print("ignore_features 1:",ignore_features[0])
print("#ignore_features:",len(ignore_features))

#sensitive feature and privileged value
with open(public_path+"/constants/sensitive_feature.ini", 'r') as f:
    sensitive_feature=(f.read().split(','))
print("sensitive_feature:",sensitive_feature)
sensitive_attr=sensitive_feature[0]
priv_value=sensitive_feature[1]
print("sensitive_feature:",sensitive_feature[0])
print("privileged value:",sensitive_feature[1])

#target feature
with open(public_path+"/constants/target_feature.ini", 'r') as f:
    target_feature=(f.read())
print("target_feature:",target_feature)

#alphas list
with open(public_path+"/constants/alphas.ini", 'r') as f:
    alphas=ast.literal_eval(f.read())
    print("alphas:",alphas)
    
#betas list
with open(public_path+"/constants/betas.ini", 'r') as f:
    betas=ast.literal_eval(f.read())
    print("betas:",betas)
    
#number of bins list
with open(public_path+"/constants/bin_nums.ini", 'r') as f:
    bin_nums=ast.literal_eval(f.read())
    print("bin_nums:",bin_nums)
    
#size lists
with open(public_path+"/constants/priv_sizes.ini", 'r') as f:
    priv_sizes=ast.literal_eval(f.read())
    print("priv_sizes:",priv_sizes)
    
with open(public_path+"/constants/unpriv_sizes.ini", 'r') as f:
    unpriv_sizes=ast.literal_eval(f.read())
    print("unpriv_sizes:",unpriv_sizes)

if (sum(priv_sizes)<max(bin_nums)) | (sum(unpriv_sizes)<max(bin_nums)):
    print("priv_size:",sum(priv_sizes))
    print("unpriv_size:",sum(unpriv_sizes))
    print("max_bin_nums:",max(bin_nums))
    raise RuntimeError("Dataset size of at least one group is smaller than number of bins")


# In[ ]:

#Load private data
df=pd.read_table(private_path+"/party"+str(my_id)+"/"+ dataset_name + "_partial.csv",delimiter=",")
df


# In[ ]:




# In[ ]:

#perform calcs
#Each party ranks its elements in ascending order

X={'priv':[], 'unpriv':[]} 

#if is_numeric_dtype(df[sensitive_attr]):
if df[sensitive_attr].dtype.kind in 'i':
    priv_value=int(priv_value)
if df[sensitive_attr].dtype.kind in 'f':
    priv_value=float(priv_value)
#print(priv_value)
    
for attr in features:
    X["priv"].append(np.sort(round(df[attr][df[sensitive_attr] == priv_value]*10**(decimals))))
    X["unpriv"].append(np.sort(round(df[attr][df[sensitive_attr] != priv_value]*10**(decimals))))

print(X["priv"])
print(X["unpriv"])


# In[ ]:




# In[ ]:

n = {'priv':np.sum(priv_sizes), 'unpriv':np.sum(unpriv_sizes)} 


# In[ ]:

n


# In[ ]:

#Initialize public variables of K-lists 
#using public constants - n and bin numbers

import os
from pathlib import Path


K={'priv':[], 'unpriv':[]} 
kdics={'priv':{}, 'unpriv':{}}

q={'priv':[], 'unpriv':[]} 
r={'priv':[], 'unpriv':[]} 
for group in ["priv","unpriv"]:
    kdic={}
    for i in range(0,len(features)):
        q[group].append(n[group]//bin_nums[i]) #quotient
        r[group].append(n[group]-bin_nums[i]*q[group][i]) #remainder
        K1=np.arange(1,r[group][i]*(q[group][i]+1)+1,q[group][i]+1)
        K2=np.arange(r[group][i]*(q[group][i]+1)+1,n[group]+1,q[group][i])
        K[group]=[*K1,*K2,n[group]]
        kdic[features[i]]=K[group]
    kdics[group]=kdic
    
print(kdics)




# In[ ]:

kdics_med={'priv':{}, 'unpriv':{}}

for group in ["priv","unpriv"]:
    kdic_med={}
    for attr in features:
        k_min_list=kdics[group][attr]
        k_med_list=[]
        for i in range(0,len(k_min_list)-1):
            k_med_list.append(int(np.ceil((k_min_list[i]+k_min_list[i+1]-1)/2)))
        kdic_med[attr]=k_med_list
    kdics_med[group]=kdic_med

print(kdics_med[group])
print(kdics_med)


# In[ ]:




# In[ ]:

#Compute m-list using "Find-Ranked-Element-Multiparty"
#Initialize the current range [a,b] to [α,β]
#Use our secure application *my_sum_compare.py* app to compute secure computation
#based on viff and TUeVIFFie packages
#we also added an option to use real numbers with truncation


#if ("III_min" in repairs) | ("I_min" in repairs):
mdics={'priv':{}, 'unpriv':{}}

if my_id==0:    
    with open(os.getcwd() + '/runtimes_' + dataset_name + "_" + sensitive_feature[0]+ "_p" + str(parties) + ".csv", 'a') as f:
        f.write("dataset_name,sensitive_feature,parties,repair_type,step,feature,group,max_bin,lam,runtime[sec]\n")

for group in ["priv","unpriv"]:
    mdic={}
    for i in range(0,len(features)):
        
        start_time_mlist_per_feature_per_group = time.time()
        
        mlist=[]
        Klist=kdics[group][features[i]]
        print("Read Klist", Klist)
        for k in Klist:
        #for k in [1]:
            print("#################################################################################")
            print("Begin searching for the *" + str(k) + "*-ranked element for feature " + features[i] + " and group " + group)
            print("#################################################################################")
            a=int(alphas[i]*(10**(decimals)))
            b=int(betas[i]*(10**(decimals)))
            print("Read a", alphas[i])
            print("Read b", betas[i])
            iteration = 0
            done=False
            while not done:
                m_guess = np.ceil((a+b)/2)
                print("a:", a/(10**(decimals)))
                print("b:", b/(10**(decimals)))                
                print("Calculate m_guess", m_guess/(10**(decimals)))
                l=np.sum(X[group][i]<m_guess)
                g=np.sum(X[group][i]>m_guess)
                print("Calculate l")
                print(l)
                print("Calculate g")
                print(g)
                my_l = l
                print("my_l:",my_l)
                print("k-1:",k-1)
                print("Waiting for other parties to respond...")
                python2_command = python2_path + ' ' + os.getcwd() + '\\MyDistExperiment\\my_secure_apps\\my_sum_compare_' + str(parties) + 'p.py --no-ssl ' + player_file_path + ' ' + str(my_l) + ' ' + str(k-1) + ' le'
                process = subprocess.Popen(python2_command.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()

                print(output)
                output_string=str(output)
                start_string=output_string.find("comp: {")
                end_string=output_string[start_string:].find("}")
                output_comp1=float(output_string[(start_string+len("comp: {")):(start_string+end_string)])
                print("sigma(l)<=k-1:" ,output_comp1)

                my_g = g
                print("my_g:",my_g)
                my_comp_value = n[group]-k
                print("n-k:",my_comp_value)
                print("Waiting for other parties to respond...")
                python2_command = python2_path + ' ' + os.getcwd() + '\\MyDistExperiment\\my_secure_apps\\my_sum_compare_' + str(parties) + 'p.py --no-ssl ' + player_file_path + ' ' + str(my_g) + ' ' + str(my_comp_value) + ' le'
                process = subprocess.Popen(python2_command.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()

                print(output)
                output_string=str(output)
                start_string=output_string.find("comp: {")
                end_string=output_string[start_string:].find("}")
                output_comp2=float(output_string[(start_string+len("comp: {")):(start_string+end_string)])
                print("sigma(g)<=n-k: " , output_comp2)

                if (output_comp1==1) & (output_comp2==1):
                    done=True
                    print("#################################################################################")
                    print("The *" + str(k) + "*-ranked element for feature " + features[i] + " and group " + group + " is " + str(m_guess/(10**(decimals))))
                    print("#################################################################################")
                    mlist.append(m_guess/(10**(decimals)))
                else:
                    if (output_comp1==0):
                        b=m_guess-1
                    if (output_comp2==0):
                        a = m_guess+1                
                print("I am party " + str(my_id) + ", I ended round " + str(iteration))
                iteration=iteration+1

        mdic[features[i]]=mlist
        
        runtime_mlist_per_feature_per_group = time.time()-start_time_mlist_per_feature_per_group
        if my_id==0:                    
            valList=[dataset_name,sensitive_feature[0],parties,"III_min","mlist",features[i],group,max(bin_nums),"-",runtime_mlist_per_feature_per_group]
            with open(os.getcwd()+'/runtimes_' + dataset_name + "_" + sensitive_feature[0]+ "_p" + str(parties) + ".csv", 'a') as f:
                for item in valList:
                    f.write('"%s",' % item)
                f.write('\n')
        
        print(mdic)
    mdics[group]=mdic    
    print(mdics)

#after finishing this step - calculate the repair


# In[ ]:

#Compute m-list using "Find-Ranked-Element-Multiparty"
#Initialize the current range [a,b] to [α,β]
#Use our secure application *my_sum_compare.py* app to compute secure computation
#based on viff and TUeVIFFie packages
#we also added option to use real numbers with truncation

if "I_med" in repairs:
    mdics_med={'priv':{}, 'unpriv':{}}

    for group in ["priv","unpriv"]:
        mdic_med={}
        for i in range(0,len(features)):
            mlist_med=[]
            Klist_med=kdics_med[group][features[i]]
            print("Read Klist", Klist_med)
            for k in Klist_med:
            #for k in [1]:
                print("#################################################################################")
                print("Begin searching for the *" + str(k) + "*-ranked element for feature " + features[i] + " and group " + group)
                print("#################################################################################")
                a=int(alphas[i]*(10**(decimals)))
                b=int(betas[i]*(10**(decimals)))
                print("Read a", alphas[i])
                print("Read b", betas[i])
                iteration = 0
                done=False
                while not done:
                    m_guess = np.ceil((a+b)/2)
                    print("a:", a/(10**(decimals)))
                    print("b:", b/(10**(decimals)))                
                    print("Calculate m_guess", m_guess/(10**(decimals)))
                    l=np.sum(X[group][i]<m_guess)
                    g=np.sum(X[group][i]>m_guess)
                    print("Calculate l")
                    print(l)
                    print("Calculate g")
                    print(g)
                    my_l = l
                    print("my_l:",my_l)
                    print("k-1:",k-1)
                    print("Waiting for other parties to respond...")
                    python2_command = python2_path + ' ' + os.getcwd() + '\\MyDistExperiment\\my_secure_apps\\my_sum_compare_' + str(parties) + 'p.py --no-ssl ' + player_file_path + ' ' + str(my_l) + ' ' + str(k-1) + ' le'
                    process = subprocess.Popen(python2_command.split(), stdout=subprocess.PIPE)
                    output, error = process.communicate()

                    print(output)
                    output_string=str(output)
                    start_string=output_string.find("comp: {")
                    end_string=output_string[start_string:].find("}")
                    output_comp1=float(output_string[(start_string+len("comp: {")):(start_string+end_string)])
                    print("sigma(l)<=k-1:" ,output_comp1)

                    my_g = g
                    print("my_g:",my_g)
                    my_comp_value = n[group]-k
                    print("n-k:",my_comp_value)
                    print("Waiting for other parties to respond...")
                    python2_command = python2_path + ' ' + os.getcwd() + '\\MyDistExperiment\\my_secure_apps\my_sum_compare_' + str(parties) + 'p.py --no-ssl ' + player_file_path + ' ' + str(my_g) + ' ' + str(my_comp_value) + ' le'
                    process = subprocess.Popen(python2_command.split(), stdout=subprocess.PIPE)
                    output, error = process.communicate()

                    print(output)
                    output_string=str(output)
                    start_string=output_string.find("comp: {")
                    end_string=output_string[start_string:].find("}")
                    output_comp2=float(output_string[(start_string+len("comp: {")):(start_string+end_string)])
                    print("sigma(g)<=n-k: " , output_comp2)

                    if (output_comp1==1) & (output_comp2==1):
                        done=True
                        print("#################################################################################")
                        print("The *" + str(k) + "*-ranked element for feature " + features[i] + " and group " + group + " is " + str(m_guess/(10**(decimals))))
                        print("#################################################################################")
                        mlist_med.append(m_guess/(10**(decimals)))
                    else:
                        if (output_comp1==0):
                            b=m_guess-1
                        if (output_comp2==0):
                            a = m_guess+1                
                    print("I am party " + str(my_id) + ", I ended round " + str(iteration))
                    iteration=iteration+1

            mdic_med[features[i]]=mlist_med
            print(mdic_med)
        mdics_med[group]=mdic_med    
        print(mdics_med)

    #after finishing this step - calculate the repair


# In[ ]:




# In[ ]:

for repair_type in repairs:
    if repair_type in ['III_min', 'I_min']:
        mdics_rep=mdics
        mdics_find_bin=mdics
    if repair_type in ['I_med']:
        mdics_rep=mdics_med
        mdics_find_bin=mdics
    for lam in lams:
        df_priv_repaired = pd.DataFrame({}) 

        for attr in features:

            start_time_repair_per_feature_per_group = time.time()
            
            ####################################
            #if x is privileged - repair
            #only repairing privileged values
            ####################################
            x=df[attr][df[sensitive_attr] == priv_value]
            x

            bins_priv=mdics_find_bin["priv"][attr]
            print(bins_priv)

            bins_unpriv=mdics_find_bin["unpriv"][attr]
            print(bins_unpriv)


            ####################################
            # Find the relevant bins of 𝑥, where 𝑥 is between the bounds of the bin
            ####################################

            df_tmp = pd.DataFrame({}) 

            for i in range(0, len(bins_priv)-1):
                if i <(len(bins_priv)-2):
                    print(i,",",bins_priv[i],",",bins_priv[i+1])
                    df_tmp.insert(len(df_tmp.columns),str(i),(bins_priv[i] <= x) & (x < bins_priv[i+1]) )
                else:
                    print(i,",",bins_priv[i],",",bins_priv[i+1])
                    df_tmp.insert(len(df_tmp.columns),str(i),(bins_priv[i] <= x) & (x <= bins_priv[i+1]))

            print(df_tmp)


            ####################################
            #Out of the relevant bins of 𝑥, select the lower one for the repair
            ####################################
            #This is a search in each instance row - for the smallest bin with "True"
            #This will be the bin for repair.
            binlist=[]
            for i in range(0,len(df_tmp.index)):
                binlist.append(int(df_tmp[i:i+1].idxmax(axis=1)))

            print(binlist)


            bins_priv=mdics_rep["priv"][attr]
            print(bins_priv)

            bins_unpriv=mdics_rep["unpriv"][attr]
            print(bins_unpriv)
            
            
            ####################################
            #Perform repair
            ####################################
            if repair_type in ["III_min"]:
                x1=(np.array(bins_unpriv)[binlist]+(np.array(x)-np.array(bins_priv)[binlist])/(np.array(bins_priv)[np.array(binlist)+1]-np.array(bins_priv)[binlist]+sys.float_info.epsilon)*(np.array(bins_unpriv)[np.array(binlist)+1]-np.array(bins_unpriv)[binlist]))
            if repair_type in ["I_min","I_med"]:
                x1=(np.array(bins_unpriv)[binlist])
            #print(x1)
            print("lambda:",lam)
            print("original x:",x)
            x_repaired = (1-lam)*x+lam*x1
            print("x_repaired:",x_repaired)


            df_priv_repaired.insert(len(df_priv_repaired.columns),str(attr),x_repaired)

            runtime_repair_per_feature_per_group = time.time()-start_time_repair_per_feature_per_group
            
            if my_id==0:                    
                valList=[dataset_name,sensitive_feature[0],parties,repair_type,"repair",attr,"-",max(bin_nums),lam,runtime_repair_per_feature_per_group]
                with open(os.getcwd()+'/runtimes_' + dataset_name + "_" + sensitive_feature[0]+ "_p" + str(parties) + ".csv", 'a') as f:
                    for item in valList:
                        f.write('"%s",' % item)
                    f.write('\n')                        
            
        df_priv_repaired.insert(0,sensitive_attr,priv_value)
        df_priv_repaired = pd.concat([df_priv_repaired, df[ignore_features][df[sensitive_attr] == priv_value]], axis=1)

        df_priv_repaired.insert(len(df_priv_repaired.columns),target_feature,df[target_feature][df[sensitive_attr] == priv_value])

        print(df_priv_repaired)

        df_unpriv=df[np.array([sensitive_attr,*features,*ignore_features,target_feature])][df[sensitive_attr] != priv_value]
        print(df_unpriv)

        df_repaired=pd.concat([df_priv_repaired,df_unpriv])
        print("################################################")
        print("The repaired dataframe for lambda=" + str(lam) + " is:")
        print("################################################")
        print(df_repaired)
        print("-----------------------------------------------------------------------")

        ################################################
        ###save repaired csv with lambda in file name
        ################################################    
        df_repaired.to_csv(private_path+"/party"+str(my_id)+"/" + dataset_name + "_sens-" + str(sensitive_attr) + "_rep_lam"+str(lam)+"_bins" + str(bin_nums[0]) + "_rep_type" + repair_type + ".csv",index=False)    


# In[ ]:




# In[ ]:




# In[ ]:



