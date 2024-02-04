#!/usr/bin/env python3

import msprime
import tskit
import numpy as np
import itertools
import scipy
import matplotlib.pyplot as plt
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-z",dest="Z_Ne", type=int, required=True, help="")
    parser.add_argument("-m",dest="Mig_rate", type=float, required=True, help="")
    parser.add_argument("-abxd",dest="ABX_div", type=int, required=True, help="")
    parser.add_argument("-abd",dest="AB_div", type=int, required=True, help="")
    parser.add_argument("-ab",dest="AB_Ne", type=int, required=True, help="")
    parser.add_argument("-x",dest="X_Ne", type=int, required=True, help="")
    parser.add_argument("-abx",dest="ABX_Ne", type=int, required=True, help="")

args = parser.parse_args()

Z_Ne=args.Z_Ne
Mig_rate=args.Mig_rate
ABX_div=args.ABX_div
AB_div=args.AB_div
AB_Ne=args.AB_Ne
X_Ne=args.X_Ne
ABX_Ne=args.ABX_Ne

###define demography

demography2 = msprime.Demography()
demography2.add_population(name="A", initial_size=5_000)
demography2.add_population(name="B", initial_size=5_000)
demography2.add_population(name="ABX", initial_size=ABX_Ne)
demography2.add_population(name="O", initial_size=5_000)
demography2.add_population(name="ABXZ", initial_size=5_000)
demography2.add_population(name="X", initial_size=X_Ne)
demography2.add_population(name="AB", initial_size=AB_Ne)
demography2.add_population(name="Z", initial_size=Z_Ne)
demography2.add_population(name="ABXZO", initial_size=5_000)
demography2.set_migration_rate("A","Z", Mig_rate)
demography2.set_migration_rate("B","Z", Mig_rate)
demography2.add_population_split(time=AB_div, derived=["A", "B"], ancestral="AB")
demography2.add_population_split(time=ABX_div, derived=["X","AB"], ancestral="ABX")
demography2.add_population_split(time=180, derived=["ABX","Z"], ancestral="ABXZ")
demography2.add_population_split(time=200, derived=["ABXZ","O"], ancestral="ABXZO")


###simulate and sample

def sim_replicates(num_replicates,d):
    ancestry_reps = msprime.sim_ancestry(
    samples={"A": 25, "X": 25,
             "B": 25, "O": 25, "Z":25},
    demography=d,
    ploidy=2,
    sequence_length=100_000_000,
    recombination_rate=1.1e-8,
    num_replicates=num_replicates,
    )
    for ts in ancestry_reps:
        mutated_ts = msprime.sim_mutations(ts, rate=1e-8)
        yield mutated_ts
        

###calculate f4###

def f4sim(rep,d):

    
    f4s=[]

    for replicate_index, ts in enumerate(sim_replicates(rep,d)):

        samplesX_gf=ts.samples(population=5)
        samplesA_gf=ts.samples(population=0)
        samplesB_gf=ts.samples(population=1)
        samples_D_gf=ts.samples(population=3)

        f4=ts.f4([samplesA_gf,samples_D_gf,samplesX_gf,samplesB_gf])

        f4s.append(f4)



    
    return f4s



f4sim=f4sim(100,demography2) #replicates

f4res=[np.array(f4sim),np.array([Z_Ne]*100),np.array([Mig_rate]*100),
np.array([ABX_div]*100),np.array([AB_div]*100),np.array([AB_Ne]*100),
np.array([X_Ne]*100),np.array([ABX_Ne]*100)]

f4res_t=pd.DataFrame(f4res).T
f4res_t2 = f4res_t.rename(columns={0: 'f4', 1: 'Z_Ne',2: 'Mig_rate', 3: 'ABX_div',
    4: 'AB_div', 5: 'AB_Ne',6: 'X_Ne', 7: 'ABX_Ne'})

path = ["f4sim_","Z",str(Z_Ne),"_","M",str(Mig_rate),"_",
    "ABXd",str(ABX_div),"_","ABd",str(AB_div),"_",
    "AB",str(AB_Ne),"_","X",str(X_Ne),"_","ABX",str(ABX_Ne),".csv"]

f4res_path= "".join(path)

f4res_t2.to_csv(f4res_path,index=False)



###calculate f3###

def f3sim_ind(rep,d):


    for replicate_index, ts in enumerate(sim_replicates(rep,d)):

        samplesX_gf=ts.samples(population=5)
        samplesA_gf=ts.samples(population=0)
        samplesB_gf=ts.samples(population=1)
        samples_D_gf=ts.samples(population=3)
        samples_Z_gf=ts.samples(population=7)

        #to get diploid individuals, merge haploid chromosomes per individual

        samples_A_dip=[]
        for i in range(0,49,2):
            samples_A_dip.append([samplesA_gf[i],samplesA_gf[i+1]])

        samples_X_dip=[]
        for i in range(0,49,2):
            samples_X_dip.append([samplesX_gf[i],samplesX_gf[i+1]])

        samples_B_dip=[]
        for i in range(0,49,2):
            samples_B_dip.append([samplesB_gf[i],samplesB_gf[i+1]])

        samples_Z_dip=[]
        for i in range(0,49,2):
            samples_Z_dip.append([samples_Z_gf[i],samples_Z_gf[i+1]])

        pop_pairs=samples_A_dip+samples_X_dip+samples_B_dip+samples_Z_dip

        #individual combinations for pairwise f3
        
        samples_A2_dip_pairs=itertools.combinations(pop_pairs, 2)
        samples_A2_dip_pairs_list=[list(x) for x in samples_A2_dip_pairs]

        f3_pairs=[]
        for i in samples_A2_dip_pairs_list:
            f3=ts.f3([samples_D_gf,np.array(i[0]),np.array(i[1])])
            f3_pairs.append(f3)


    allp_f3s= [np.array(samples_A2_dip_pairs_list),np.array(f3_pairs)]

    return allp_f3s



allp_f3s_res=f3sim_ind(1,demography2)


#optional export results in four different files

pairs=allp_f3s_res[0]
pair1=[x[0] for x in pairs]
pair2=[x[1] for x in pairs]


allp_f3s_f3val=pd.DataFrame(allp_f3s_res[1])
allp_f3s_pair1=pd.DataFrame(pair1)
allp_f3s_pair2=pd.DataFrame(pair2)

f3params=[np.array([Z_Ne]),np.array([Mig_rate]),
np.array([ABX_div]),np.array([AB_div]),np.array([AB_Ne]),
np.array([X_Ne]),np.array([ABX_Ne])]
f3params_df = pd.DataFrame(f3params)


path_f3val = ["f3sim_","Z",str(Z_Ne),"_","M",str(Mig_rate),"_",
    "ABXd",str(ABX_div),"_","ABd",str(AB_div),"_",
    "AB",str(AB_Ne),"_","X",str(X_Ne),"_","ABX",str(ABX_Ne),"_f3val.csv"]
    
path_pair1 = ["f3sim_","Z",str(Z_Ne),"_","M",str(Mig_rate),"_",
    "ABXd",str(ABX_div),"_","ABd",str(AB_div),"_",
    "AB",str(AB_Ne),"_","X",str(X_Ne),"_","ABX",str(ABX_Ne),"_pair1.csv"]
    
path_pair2 = ["f3sim_","Z",str(Z_Ne),"_","M",str(Mig_rate),"_",
    "ABXd",str(ABX_div),"_","ABd",str(AB_div),"_",
    "AB",str(AB_Ne),"_","X",str(X_Ne),"_","ABX",str(ABX_Ne),"_pair2.csv"]
    
path_parameters = ["f3sim_","Z",str(Z_Ne),"_","M",str(Mig_rate),"_",
    "ABXd",str(ABX_div),"_","ABd",str(AB_div),"_",
    "AB",str(AB_Ne),"_","X",str(X_Ne),"_","ABX",str(ABX_Ne),"_parameters.csv"]


path_f3val_fin= "".join(path_f3val)
path_pair1_fin= "".join(path_pair1)
path_pair2_fin= "".join(path_pair2)
path_parameters_fin= "".join(path_parameters)


allp_f3s_f3val.to_csv(path_f3val_fin,index=False)
allp_f3s_pair1.to_csv(path_pair1_fin,index=False)
allp_f3s_pair2.to_csv(path_pair2_fin,index=False)
f3params_df.to_csv(path_parameters_fin,index=False)
