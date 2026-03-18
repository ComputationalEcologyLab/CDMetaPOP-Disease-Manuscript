import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import pandas as pd
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

def log_likelihood(params, t, S_ave, I_ave, R_ave, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput):
    """
    Calculates how well the model fits the data.
    """
    try:
        output = CDMetaPOP(params, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput)
        model_s = output['S'] # Adjust 'S' to match your target data column
        model_i = output['I']        
        model_r = output['R']

        sigma2 = 100
        ssr = np.sum((S_ave - model_s)**2) + np.sum((I_ave - model_i)**2) + np.sum((R_ave - model_r)**2)
        print(ssr)
        
        return -0.5 * (ssr / sigma2)
    except Exception as e:
        return -np.inf

def log_prior(params):
    """
    Strictly enforces your (0, 1) boundaries.
    """
    p1, p2 = params
    if 0.0 < p1 < 1.0 and 0.0 < p2 < 1.0:
        return 0.0
    return -np.inf

def log_probability(params, t, S_ave, I_ave, R_ave, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, t, S_ave, I_ave, R_ave, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput)


# Comparison ###############################################

sub_label = ['a', 'b', 'c', 'd']

plt.rcParams['axes.labelsize'] = 35
plt.rcParams['xtick.labelsize'] = 21
plt.rcParams['ytick.labelsize'] = 21
plt.rcParams['legend.fontsize'] = 26
label_size = 35

runtime = 50
num_sims = 100
time = np.linspace(0, runtime, runtime)
base_path = r"Figure_2_summary_data/bestfit/Bestfit_Data/"
prefix = 'run0batch0mc'
comps = ['S', 'I', 'R']
means, all_data = load_fixed_simulation_data(base_path, prefix, comps, env=True)

S_ave = np.load(r"Figure_2_summary_data/bestfit/ODE_Data/S_ave.npy")
I_ave = np.load(r"Figure_2_summary_data/bestfit/ODE_Data/I_ave.npy")
R_ave = np.load(r"Figure_2_summary_data/bestfit/ODE_Data/R_ave.npy")

fig = plt.figure(figsize=(10, 10))

plt.plot(time, S_ave, "-", color='#56B4E9', label="S: ODE")
plt.plot(time, I_ave, "-", color='#D55E00', label="I: ODE")
plt.plot(time, R_ave, "-", color='#009E73', label="R: ODE")

plt.plot(time, means[comps[0]], color='#56B4E9', linestyle='--', label="S: Best Fit")
plt.plot(time, means[comps[1]], color='#D55E00', linestyle='--', label="I: Best Fit")
plt.plot(time, means[comps[2]], color='#009E73', linestyle='--', label="R: Best Fit")

for line in all_data[comps[0]]:
    plt.plot(time, line, color='#56B4E9', linestyle='--', alpha=0.1, linewidth=0.5)
for line in all_data[comps[1]]:
    plt.plot(time, line, color='#D55E00', linestyle='--', alpha=0.1, linewidth=0.5)
for line in all_data[comps[2]]:
    plt.plot(time, line, color='#009E73', linestyle='--', alpha=0.1, linewidth=0.5)

plt.text(0.02, 0.94, f"({sub_label[0]})", transform=fig.transFigure, fontsize=label_size, fontweight='bold')
plt.ylim(0,100)
plt.grid(True)
plt.xlabel("Year")
plt.ylabel("Population size")
plt.legend()
plt.tight_layout()
plt.savefig("Figure_outputs/bestfit.png")
plt.close()

############# CI
S_trajectories = np.load("Figure_2_summary_data/S_trajectories.npy")
I_trajectories = np.load("Figure_2_summary_data/I_trajectories.npy")
R_trajectories = np.load("Figure_2_summary_data/R_trajectories.npy")

Ulower_bound = np.percentile(S_trajectories, 2.5, axis=0)
Uupper_bound = np.percentile(S_trajectories, 97.5, axis=0)
Umedian = np.percentile(S_trajectories, 50, axis=0)

Ilower_bound = np.percentile(I_trajectories, 2.5, axis=0)
Iupper_bound = np.percentile(I_trajectories, 97.5, axis=0)
Imedian = np.percentile(I_trajectories, 50, axis=0)

Rlower_bound = np.percentile(R_trajectories, 2.5, axis=0)
Rupper_bound = np.percentile(R_trajectories, 97.5, axis=0)
Rmedian = np.percentile(R_trajectories, 50, axis=0)


# 3. Create the Plot
plt.figure(figsize=(10, 10))

plt.fill_between(time, Ulower_bound, Uupper_bound, color='#56B4E9', alpha=0.2)
plt.plot(time, Umedian, color='#56B4E9', linestyle='--', label='S: MCMC')

plt.fill_between(time, Ilower_bound, Iupper_bound, color='#D55E00', alpha=0.2)
plt.plot(time, Imedian, color='#D55E00', linestyle='--', label='I: MCMC')

plt.fill_between(time, Rlower_bound, Rupper_bound, color='#009E73', alpha=0.2)
plt.plot(time, Rmedian, color='#009E73', linestyle='--', label='R: MCMC')


S_ave = np.load(r"Figure_2_summary_data/bestfit/ODE_Data/S_ave.npy")
I_ave = np.load(r"Figure_2_summary_data/bestfit/ODE_Data/I_ave.npy")
R_ave = np.load(r"Figure_2_summary_data/bestfit/ODE_Data/R_ave.npy")
plt.plot(time, S_ave, color='#56B4E9', label='S: ODE')
plt.plot(time, I_ave, color='#D55E00', label='I: ODE')
plt.plot(time, R_ave, color='#009E73', label='R: ODE')

plt.text(0.02, 0.94, f"({sub_label[3]})", transform=fig.transFigure, fontsize=label_size, fontweight='bold')
plt.ylim(0,100)
plt.grid(True)
plt.xlabel('Year')
plt.ylabel('Population Size')
plt.legend()
plt.tight_layout()
plt.savefig("Figure_outputs/CI")
plt.close()

######################################################3
sampler = np.load(r"Figure_2_summary_data/mcmc_sampler.npy", allow_pickle=True).item()
flat_samples = sampler.get_chain(discard=100, thin=10, flat=True)

labels = [r"$\beta$", r"$\gamma$"]
ndim = len(labels)

fig = corner.corner(
    flat_samples, 
    labels=labels, 
    truths=[0.51457746, 0.20896594],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    plot_datapoints=True,
    smooth=0.75,
    title_kwargs={"fontsize": 26},
    fig=plt.figure(figsize=(10, 10))
)
axes = np.array(fig.axes).reshape((ndim, ndim))
org_truth = [0.5, 0.2]
for i in range(ndim):
    for j in range(i + 1):
        ax = axes[i, j]
        ax.axvline(org_truth[j], color="red", linestyle="--", lw=1.5)
        if i > j:
            ax.axhline(org_truth[i], color="red", linestyle="--", lw=1.5)
            ax.plot(org_truth[j], org_truth[i], marker="s", color="red")
plt.legend(
    [plt.Line2D([0], [0]), 
     plt.Line2D([0], [0], color='red', linestyle='--'),
     plt.Line2D([0], [0], color='black', linestyle='--')],
    ['Best Fit', 'ODE Value', 'Quantiles'],
    bbox_to_anchor=(1, ndim), loc='upper right'
)

axes[0,0].text(0.02, 0.94, f"({sub_label[1]})", transform=fig.transFigure, fontsize=label_size, fontweight='bold')
plt.savefig("Figure_outputs/corner.png")
plt.close()




#####################################################################
samples = sampler.get_chain()
fig, axs = plt.subplots(2, figsize=(10, 10), sharex=True)

axs[0].plot(samples[:, :, 0], "k", alpha=0.3)
axs[0].set_ylabel(labels[0])

axs[1].plot(samples[:, :, 1], "k", alpha=0.3)
axs[1].set_ylabel(labels[1])

axs[0].text(0.02, 0.94, f"({sub_label[2]})", transform=fig.transFigure, fontsize=label_size, fontweight='bold')
plt.xlabel("Number of Iterations")
plt.tight_layout()
plt.savefig("Figure_outputs/walker.png")
plt.close()


#############################################
import matplotlib.image as mpimg

fig, axs = plt.subplots(2, 2, figsize=(20, 20))

axs[0, 0].imshow(mpimg.imread("Figure_outputs/bestfit.png"))
axs[0, 0].axis('off')
axs[0, 1].imshow(mpimg.imread("Figure_outputs/corner.png"))
axs[0, 1].axis('off')
axs[1, 0].imshow(mpimg.imread("Figure_outputs/walker.png"))
axs[1, 0].axis('off')
axs[1, 1].imshow(mpimg.imread("Figure_outputs/CI.png"))
axs[1, 1].axis('off')

plt.tight_layout()
plt.savefig("Figure_outputs/Figure_2.png")
plt.close()











