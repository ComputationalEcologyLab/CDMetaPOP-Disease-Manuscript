import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.lines import Line2D


# Plotting
model_configs = {
    "SIR": {
        "base_path": r"../Figure_Reproduction/Figure1/Figure_1_summary_data/SIR",
        "prefix": "run0batch0mc",

        "comps": ['S', 'I', 'R']
    },
    "SIRS": {
        "base_path": r"../Figure_Reproduction/Figure1/Figure_1_summary_data/SIRS",
        "prefix": "run0batch0mc",
        "comps": ['S', 'I', 'R']
    },
    "SEIR": {
        "base_path": r"../Figure_Reproduction/Figure1/Figure_1_summary_data/SEIR",
        "prefix": "run0batch0mc",
        "comps": ['S', 'E', 'I', 'R']
    },
    "SIRP": {
        "base_path": r"../Figure_Reproduction/Figure1/Figure_1_summary_data/SIRP",
        "prefix": "run0batch0mc",
        "comps": ['S', 'I', 'R', 'P']
    }
}

timesteps = range(51)

for model, cfg in model_configs.items():
    
    file_path = f"{cfg['base_path']}/{model}_results_"
    means = np.load(file_path+'means.npy', allow_pickle=True).item()
    all_data = np.load(file_path+'all_data.npy', allow_pickle=True).item()
    t_ode = np.load(file_path+'t_ode.npy', allow_pickle=True)
    sol = np.load(file_path+'sol.npy', allow_pickle=True)
    comps = list(np.load(file_path+'comps.npy', allow_pickle=True))


    for comp in cfg["comps"]:
        idx = comps.index(comp)
        observed = sol[:, idx][::20]
        simulated = means[comp][:-1]
        def calculate_kge(observed, simulated):
            # 1. Correlation coefficient
            r = np.corrcoef(observed, simulated)[0, 1]
            
            # 2. Variability ratio (Alpha)
            alpha = np.std(simulated) / np.std(observed)
            
            # 3. Bias ratio (Beta)
            beta = np.mean(simulated) / np.mean(observed)
            
            # Calculate KGE
            kge = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
            
            return kge, r, alpha, beta

        # Usage
        score, corr, var, bias = calculate_kge(observed, simulated)

        print(f"Total KGE Score: {score:.3f}")
        print(f"  - Correlation: {corr:.3f}")
        print(f"  - Variability Match: {var:.3f} (1.0 is perfect)")
        print(f"  - Bias Match: {bias:.3f} (1.0 is perfect)")
        print("\n")





