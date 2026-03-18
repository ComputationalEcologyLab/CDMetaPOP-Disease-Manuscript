#nohup bash -c 'python CDMetaPOPFit_SIRS.py &'

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
import shutil

def order_of_magnitude(value):
    if value == 0:
        return 0
    return int(np.floor(np.log10(abs(value))))

def get_recent_folder(path, phrase):
    folders = [f.name for f in path.iterdir() if f.is_dir() and phrase in f.name]

    def get_number(folder_name):
        match = re.search(rf"{phrase}(\d+)", folder_name)
        return int(match.group(1)) if match else 0

    folders.sort(key=get_number)
    
    return folders[-1]

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

def CDMetaPOP(params0, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput):

    print(params0)
    beta = params0[0]
    gamma = params0[1]
    phi = params0[2]

    Transition_Matrix = [
        [0.0,  0.0,  phi],
        [beta, 0.0,  0.0],
        [0.0,  gamma, 0.0],
    ]

    TransitionMatrix_SIDP_SSR = pd.DataFrame(Transition_Matrix)
    TransitionMatrix_SIDP_SSR.to_csv(LocOfTM, index=False, header=False)   

    # Run CDMetaPOP
    os.system("cd "+ LocOfsrc +"&& python CDmetaPOP.py "+ LocOfRunVars +" RunVars.csv output/output")
    gc.collect()

    # Get new data folder
    path = Path(LocOfOutput)
    phrase = "output"
    recent_folder = get_recent_folder(path, phrase)

    # Make copy of transitions matrix
    TMFile = LocOfTM.split("/")
    TransitionMatrix_SIDP_SSR.to_csv(LocOfOutput + recent_folder + '/' + TMFile[-1], index=False, header=False) 

    # Create data from CDMetaPOP runs
    base_path = LocOfOutput+recent_folder
    prefix = 'run0batch0mc'
    comps = ['S', 'I', 'R']
    means, all_data = load_fixed_simulation_data(base_path, prefix, comps, env=True)

    gc.collect()
    return means

def SSR(params0, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput):

    output = CDMetaPOP(params0, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput)

    SSR = 0
    SSR += np.sum((T_ave - output['S'])**2)
    SSR += np.sum((I_ave - output['I'])**2)
    SSR += np.sum((D_ave - output['R'])**2)

    print(params0)
    print(SSR)
    gc.collect()
    return SSR

def create_bestfit(params0, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput):

    dir_path = Path(r"Figure_3_from_source_data/bestfit/Bestfit_Data")
    if dir_path.exists() and dir_path.is_dir():
        shutil.rmtree(dir_path)
    else:
        "Nothing"

    output = CDMetaPOP(params0, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput)

    path = Path(LocOfOutput)
    phrase = "output"
    recent_folder = get_recent_folder(path, phrase)

    source = Path(LocOfOutput+recent_folder)
    destination = Path(r"Figure_3_from_source_data/bestfit/")

    if source.exists() and source.is_dir():
        shutil.move(str(source), str(destination))
    else:
        print("Source directory not found.")
        exit(0)

    old_folder = Path(r"Figure_3_from_source_data/bestfit/"+recent_folder)
    new_folder = Path(r"Figure_3_from_source_data/bestfit/"+"Bestfit_Data")
    old_folder.rename(new_folder)

    return

LocOfTM = r'CDMetaPOP_inputs/otherfiles/disease/TransitionMatrix_SIRS_SSR.csv'

LocOfsrc = r'CDMetaPOP/'

LocOfRunVars = r'CDMetaPOP_inputs/SIRS/'


LocOfOutput = r'CDMetaPOP_inputs/SIRS/output/'


filename = LocOfRunVars+'RunVars.csv'
RunVars = pd.read_csv(filename)

runtime = RunVars.at[0, 'runtime']
num_sims = RunVars.at[0, 'mcruns']

time = np.linspace(0, runtime, runtime)
dt = time[1] - time[0]
N = 100

beta = 0.5
gamma = 0.2
phi = 0.2
vital = 0.0
params0 = [N, beta, gamma, phi]
y0 = [N - 10, 10, 0]

endtime = runtime
tau = dt
n_simulations = 100
Time = np.zeros([n_simulations, 1],dtype=tuple)
Ux = np.zeros([n_simulations, 1],dtype=tuple)
Ix = np.zeros([n_simulations, 1],dtype=tuple)
Dx = np.zeros([n_simulations, 1],dtype=tuple)

for i in range(n_simulations):
    Time[i][0] = [0.0]
    Ux[i][0] = [y0[0]]
    Ix[i][0] = [y0[1]]
    Dx[i][0] = [y0[2]]

Rate = np.zeros((5))

for i in range(n_simulations):

    timetime = 0
    while timetime < endtime:

        Rate[0] = beta * (Ix[i][-1][-1] / N) * Ux[i][-1][-1]
        Rate[1] = gamma * Ix[i][-1][-1]
        Rate[2] = phi * Dx[i][-1][-1]

        leap = tau

        uT  = np.random.random()
        uI  = np.random.random()
        uR  = np.random.random()
        
        NT  = poisson.ppf(uT, Rate[0]*leap)
        NI  = poisson.ppf(uI, Rate[1]*leap)
        NR  = poisson.ppf(uR, Rate[2]*tau)

        NT  = min(NT, Ux[i][-1][-1])
        NI  = min(NI, Ix[i][-1][-1])
        NR  = min(NR, Dx[i][-1][-1])

        Ux[i][0] = np.append(Ux[i][0], Ux[i][-1][-1] - NT + NR)
        Ix[i][0] = np.append(Ix[i][0], Ix[i][-1][-1] + NT - NI)
        Dx[i][0] = np.append(Dx[i][0], Dx[i][-1][-1] + NI - NR)

        Time[i][0] = np.append(Time[i][0], Time[i][-1][-1] + leap)
        timetime = Time[i][-1][-1]

Time_ave = Time[0][0]
T_ave = Ux.mean(axis=0)[0]
I_ave = Ix.mean(axis=0)[0]
D_ave = Dx.mean(axis=0)[0]
np.save("Figure_3_from_source_data/bestfit/ODE_Data/T_ave.npy", T_ave)
np.save("Figure_3_from_source_data/bestfit/ODE_Data/I_ave.npy", I_ave)
np.save("Figure_3_from_source_data/bestfit/ODE_Data/D_ave.npy", D_ave)

# Fitting #####################################
test_params = np.random.rand(3)

simplex = np.array([test_params, [test_params[0]*2.0, test_params[1], test_params[2]], [test_params[0], test_params[1]*2.0, test_params[2]], [test_params[0], test_params[1], test_params[2]*2.0]]) 
res = minimize(SSR, x0=test_params, args=(LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput), method="Nelder-Mead", options={'initial_simplex': simplex, 'xatol': 1e-3})
np.save("res.npy", res)

print("ODE values")
print(params0)
print(res.x)
with open(r"Figure_3_from_source_data/bestfit/bestfit_values.txt", w) as outfile:
    print(params0,file=outfile)
    print(res.x,file=outfile)
print("\n")
print("DONE")


# MCMC

def log_likelihood(params, t, T_ave, I_ave, D_ave, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput):
    """
    Calculates how well the model fits the data.
    """
    try:
        output = CDMetaPOP(params, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput)
        model_s = output['S']
        model_i = output['I']        
        model_r = output['R']

        sigma2 = 100
        ssr = np.sum((T_ave - model_s)**2) + np.sum((I_ave - model_i)**2) + np.sum((D_ave - model_r)**2)
        print(ssr)
        
        return -0.5 * (ssr / sigma2)
    except Exception as e:
        return -np.inf

def log_prior(params):
    """
    Strictly enforces your (0, 1) boundaries.
    """
    p1, p2, p3 = params
    if 0.0 < p1 < 1.0 and 0.0 < p2 < 1.0 and 0.0 < p3 < 1.0:
        return 0.0
    return -np.inf

def log_probability(params, t, T_ave, I_ave, D_ave, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, t, T_ave, I_ave, D_ave, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput)


nwalkers = 16
nsteps = 200
ndim = 3

initial_guess = res.x 
pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)

print("Starting MCMC")

sampler = emcee.EnsembleSampler(
    nwalkers, 
    ndim, 
    log_probability, 
    args=(time, T_ave, I_ave, D_ave, LocOfTM, LocOfsrc, LocOfRunVars, LocOfOutput)
)

# Run the MCMC
sampler.run_mcmc(pos, nsteps, progress=True)
np.save("mcmc_sampler.npy", sampler)

# Discard 100 steps for burn-in
flat_samples = sampler.get_chain(discard=tau, thin=10, flat=True)
np.save("flat_samples.npy", flat_samples)

print("\n--- Parameter Estimation Results ---")
labels = ["p1", "p2", "p3"]
with open(r"Figure_3_from_source_data/uncertainty_values.txt", w) as outfile:
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(f"{labels[i]}: {mcmc[1]:.4f} (+{q[1]:.4f} / -{q[0]:.4f})", file=outfile)

# CI ###########################################
runtime = 50
num_sims = 100
time = np.linspace(0, runtime, runtime)

folder_path = Path(r'CDMetaPOP_inputs/SIRS/output/')

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

np.save("Figure_3_from_source_data/S_trajectories.npy", S_trajectories)
np.save("Figure_3_from_source_data/I_trajectories.npy", I_trajectories)
np.save("Figure_3_from_source_data/R_trajectories.npy", R_trajectories)


# Plot ###############################################################

output_dir = Path(r"Figure_outputs/from_source")
output_dir.mkdir(parents=True, exist_ok=True)

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
base_path = r"Figure_3_from_source_data/bestfit/Bestfit_Data/"
prefix = 'run0batch0mc'
comps = ['S', 'I', 'R']
means, all_data = load_fixed_simulation_data(base_path, prefix, comps, env=True)

T_ave = np.load(r"Figure_3_from_source_data/bestfit/ODE_Data/T_ave.npy")
I_ave = np.load(r"Figure_3_from_source_data/bestfit/ODE_Data/I_ave.npy")
D_ave = np.load(r"Figure_3_from_source_data/bestfit/ODE_Data/D_ave.npy")

fig = plt.figure(figsize=(10, 10))

plt.plot(time, T_ave, "-", color='#56B4E9', label="S: ODE")
plt.plot(time, I_ave, "-", color='#D55E00', label="I: ODE")
plt.plot(time, D_ave, "-", color='#009E73', label="R: ODE")

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
plt.savefig("Figure_outputs/from_source/bestfit.png")
plt.close()

############# CI
S_trajectories = np.load("Figure_3_from_source_data/S_trajectories.npy")
I_trajectories = np.load("Figure_3_from_source_data/I_trajectories.npy")
R_trajectories = np.load("Figure_3_from_source_data/R_trajectories.npy")

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


T_ave = np.load(r"Figure_3_from_source_data/bestfit/ODE_Data/T_ave.npy")
I_ave = np.load(r"Figure_3_from_source_data/bestfit/ODE_Data/I_ave.npy")
D_ave = np.load(r"Figure_3_from_source_data/bestfit/ODE_Data/D_ave.npy")
plt.plot(time, T_ave, color='#56B4E9', label='S: ODE')
plt.plot(time, I_ave, color='#D55E00', label='I: ODE')
plt.plot(time, D_ave, color='#009E73', label='R: ODE')

plt.text(0.02, 0.94, f"({sub_label[3]})", transform=fig.transFigure, fontsize=label_size, fontweight='bold')
plt.ylim(0,100)
plt.grid(True)
plt.xlabel('Year')
plt.ylabel('Population Size')
plt.legend()
plt.tight_layout()
plt.savefig("Figure_outputs/from_source/CI")
plt.close()

######################################################3
sampler = np.load(r"Figure_3_from_source_data/mcmc_sampler.npy", allow_pickle=True).item()
flat_samples = sampler.get_chain(discard=100, thin=10, flat=True)

labels = [r"$\beta$", r"$\gamma$", r"$\phi$"]
ndim = len(labels)

fig = corner.corner(
    flat_samples, 
    labels=labels, 
    truths=[0.54985047, 0.22113357, 0.21202528],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    plot_datapoints=True,
    smooth=0.75,
    title_kwargs={"fontsize": 26},
    fig=plt.figure(figsize=(10, 10))
)
axes = np.array(fig.axes).reshape((ndim, ndim))
org_truth = [0.5, 0.2, 0.2]
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
plt.savefig("Figure_outputs/from_source/corner.png")
plt.close()




#####################################################################
samples = sampler.get_chain()
fig, axs = plt.subplots(3, figsize=(10, 10), sharex=True)

axs[0].plot(samples[:, :, 0], "k", alpha=0.3)
axs[0].set_ylabel(labels[0])

axs[1].plot(samples[:, :, 1], "k", alpha=0.3)
axs[1].set_ylabel(labels[1])

axs[2].plot(samples[:, :, 2], "k", alpha=0.3)
axs[2].set_ylabel(labels[2])

axs[0].text(0.02, 0.94, f"({sub_label[2]})", transform=fig.transFigure, fontsize=label_size, fontweight='bold')
plt.xlabel("Number of Iterations")
plt.tight_layout()
plt.savefig("Figure_outputs/from_source/walker.png")
plt.close()


#############################################
import matplotlib.image as mpimg

fig, axs = plt.subplots(2, 2, figsize=(20, 20))

axs[0, 0].imshow(mpimg.imread("Figure_outputs/from_source/bestfit.png"))
axs[0, 0].axis('off')
axs[0, 1].imshow(mpimg.imread("Figure_outputs/from_source/corner.png"))
axs[0, 1].axis('off')
axs[1, 0].imshow(mpimg.imread("Figure_outputs/from_source/walker.png"))
axs[1, 0].axis('off')
axs[1, 1].imshow(mpimg.imread("Figure_outputs/from_source/CI.png"))
axs[1, 1].axis('off')

plt.tight_layout()
plt.savefig("Figure_outputs/from_source/Figure_3.png")
plt.close()














