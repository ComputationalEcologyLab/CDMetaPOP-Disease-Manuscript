#nohup bash -c 'python CDMetaPOPFit.py &'

import numpy as np
import matplotlib.pyplot as plt
import os
import gc
from scipy.optimize import minimize, approx_fprime
from scipy.integrate import odeint
import re
from pathlib import Path
import pandas as pd
from scipy.stats import poisson
import emcee
import corner

def load_fixed_simulation_data(base_path, folder_prefix, compartments, env=False):
    column = "States_AfterAnyAddedInds"
    all_data = {c: [] for c in compartments}
    for i in range(num_sims):
        folder_name = f"{folder_prefix}{i}species0"
        folder_path = os.path.join(base_path, folder_name)
        file_path = os.path.join(folder_path, "summary_popAllTime_DiseaseStates.csv")

        if os.path.isfile(file_path):
            df = pd.read_csv(file_path)
            temp = {c: [] for c in compartments}
            for j in range(runtime):
                try:
                    if "SEIR" in base_path:
                        values = df[column][j].split('|')[0].split(';')
                        for k, c in enumerate(compartments):
                            temp[c].append(float(values[k]))
                    elif "SIDP" in base_path:
                        values = df[column][j].split('|')[0].split(';')
                        temp["S"].append(float(values[0]))
                        temp["I"].append(float(values[1]))
                        temp["R"].append(float(values[2]))
                        temp["P"].append(float(df["EnvResivoir_GetMetrics"][j].split('|')[0]))
                    else:
                        values = df[column][j].split('|')[0].split(';')
                        for k, c in enumerate(compartments):
                            temp[c].append(float(values[k]))
                except:
                    continue
            if all(len(temp[c]) == runtime for c in compartments):
                for c in compartments:
                    all_data[c].append(temp[c])

    means = {c: np.mean(np.array(all_data[c]), axis=0) if all_data[c] else np.zeros(runtime) for c in compartments}
    return means, all_data


############# CI
runtime = 50
num_sims = 100
time = np.linspace(0, runtime, runtime)

folder_path = Path(r'../NewFiles/OnePatch_SIRS/output/')

# List only the directories (folders)
folders = [f.name for f in folder_path.iterdir() if f.is_dir()]

S_trajectories = np.zeros((len(folders), runtime))
I_trajectories = np.zeros((len(folders), runtime))
R_trajectories = np.zeros((len(folders), runtime))
for i, folder_name in enumerate(folders):

    print(i)

    base_path = os.path.join(folder_path,folder_name)
    prefix = 'run0batch0mc'
    comps = ['S', 'I', 'R']
    means, all_data = load_fixed_simulation_data(base_path, prefix, comps, env=True)

    S_trajectories[i, :] = means['S']
    I_trajectories[i, :] = means['I']
    R_trajectories[i, :] = means['R']

np.save("S_trajectories.npy", S_trajectories)
np.save("I_trajectories.npy", I_trajectories)
np.save("R_trajectories.npy", R_trajectories)

#Ulower_bound = np.percentile(S_trajectories, 2.5, axis=0)
#Uupper_bound = np.percentile(S_trajectories, 97.5, axis=0)
#Umedian = np.percentile(S_trajectories, 50, axis=0)

#Ilower_bound = np.percentile(I_trajectories, 2.5, axis=0)
#Iupper_bound = np.percentile(I_trajectories, 97.5, axis=0)
#Imedian = np.percentile(I_trajectories, 50, axis=0)

#Rlower_bound = np.percentile(R_trajectories, 2.5, axis=0)
#Rupper_bound = np.percentile(R_trajectories, 97.5, axis=0)
#Rmedian = np.percentile(R_trajectories, 50, axis=0)


## 3. Create the Plot
#plt.figure(figsize=(10, 10))

#plt.fill_between(time, Ulower_bound, Uupper_bound, color='#56B4E9', alpha=0.2, label='95% Credible Interval (MCMC)')
#plt.plot(time, Umedian, color='#56B4E9', linestyle='--', label='MCMC Median')

#plt.fill_between(time, Ilower_bound, Iupper_bound, color='#D55E00', alpha=0.2, label='95% Credible Interval (MCMC)')
#plt.plot(time, Imedian, color='#D55E00', linestyle='--', label='MCMC Median')

#plt.fill_between(time, Rlower_bound, Rupper_bound, color='#009E73', alpha=0.2, label='95% Credible Interval (MCMC)')
#plt.plot(time, Rmedian, color='#009E73', linestyle='--', label='MCMC Median')


#T_ave = np.load(r"../mcmc/bestfit/ODE_Data/T_ave.npy")
#I_ave = np.load(r"../mcmc/bestfit/ODE_Data/I_ave.npy")
#D_ave = np.load(r"../mcmc/bestfit/ODE_Data/D_ave.npy")
#plt.plot(time, T_ave, color='#56B4E9', label='Observed Data')
#plt.plot(time, I_ave, color='#D55E00', label='Observed Data')
#plt.plot(time, D_ave, color='#009E73', label='Observed Data')

#plt.xlabel('Year')
#plt.ylabel('Population Size')
#plt.legend()
#plt.savefig("images/CI")











