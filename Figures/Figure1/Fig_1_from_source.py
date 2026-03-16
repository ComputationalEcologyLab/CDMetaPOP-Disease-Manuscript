import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import poisson

def tauleap_SIR(propensities, y0, tau, runtime, N, beta, gamma):

    n_simulations = 100
    Time = np.zeros([n_simulations, 1],dtype=tuple)
    S = np.zeros([n_simulations, 1],dtype=tuple)
    I = np.zeros([n_simulations, 1],dtype=tuple)
    R = np.zeros([n_simulations, 1],dtype=tuple)

    for i in range(n_simulations):
        Time[i][0] = [0.0]
        S[i][0] = [y0[0]]
        I[i][0] = [y0[1]]
        R[i][0] = [y0[2]]

    for i in range(n_simulations):

        timetime = 0
        while timetime < runtime:

            Rate = propensities(S[i][-1][-1], I[i][-1][-1], N, beta, gamma)

            uT  = np.random.random()
            uI  = np.random.random()
            
            NT  = poisson.ppf(uT, Rate[0]*tau)
            NI  = poisson.ppf(uI, Rate[1]*tau)

            NT  = min(NT, S[i][-1][-1])
            NI  = min(NI, I[i][-1][-1])

            S[i][0] = np.append(S[i][0], S[i][-1][-1] - NT)
            I[i][0] = np.append(I[i][0], I[i][-1][-1] + NT - NI)
            R[i][0] = np.append(R[i][0], R[i][-1][-1] + NI)

            Time[i][0] = np.append(Time[i][0], Time[i][-1][-1] + tau)
            timetime = Time[i][-1][-1]

    Time_ave = Time[0][0]
    S_ave = S.mean(axis=0)[0]
    I_ave = I.mean(axis=0)[0]
    R_ave = R.mean(axis=0)[0]
    output = np.column_stack((S_ave, I_ave, R_ave))

    return output

def tauleap_SIRS(propensities, y0, tau, runtime, N, beta, gamma, phi):

    n_simulations = 100
    Time = np.zeros([n_simulations, 1],dtype=tuple)
    S = np.zeros([n_simulations, 1],dtype=tuple)
    I = np.zeros([n_simulations, 1],dtype=tuple)
    R = np.zeros([n_simulations, 1],dtype=tuple)

    for i in range(n_simulations):
        Time[i][0] = [0.0]
        S[i][0] = [y0[0]]
        I[i][0] = [y0[1]]
        R[i][0] = [y0[2]]

    for i in range(n_simulations):

        timetime = 0
        while timetime < runtime:

            Rate = propensities(S[i][-1][-1], I[i][-1][-1], R[i][-1][-1], N, beta, gamma, phi)

            uT  = np.random.random()
            uI  = np.random.random()
            uR  = np.random.random()
            
            NT  = poisson.ppf(uT, Rate[0]*tau)
            NI  = poisson.ppf(uI, Rate[1]*tau)
            NR  = poisson.ppf(uR, Rate[2]*tau)

            NT  = min(NT, S[i][-1][-1])
            NI  = min(NI, I[i][-1][-1])
            NR  = min(NR, R[i][-1][-1])

            S[i][0] = np.append(S[i][0], S[i][-1][-1] - NT + NR)
            I[i][0] = np.append(I[i][0], I[i][-1][-1] + NT - NI)
            R[i][0] = np.append(R[i][0], R[i][-1][-1] + NI - NR)

            Time[i][0] = np.append(Time[i][0], Time[i][-1][-1] + tau)
            timetime = Time[i][-1][-1]

    Time_ave = Time[0][0]
    S_ave = S.mean(axis=0)[0]
    I_ave = I.mean(axis=0)[0]
    R_ave = R.mean(axis=0)[0]
    output = np.column_stack((S_ave, I_ave, R_ave))

    return output

def tauleap_SEIR(propensities, y0, tau, runtime, N, beta, epsilon, gamma):

    n_simulations = 100
    Time = np.zeros([n_simulations, 1],dtype=tuple)
    S = np.zeros([n_simulations, 1],dtype=tuple)
    E = np.zeros([n_simulations, 1],dtype=tuple)
    I = np.zeros([n_simulations, 1],dtype=tuple)
    R = np.zeros([n_simulations, 1],dtype=tuple)

    for i in range(n_simulations):
        Time[i][0] = [0.0]
        S[i][0] = [y0[0]]
        E[i][0] = [y0[1]]
        I[i][0] = [y0[2]]
        R[i][0] = [y0[3]]

    for i in range(n_simulations):

        timetime = 0
        while timetime < runtime:

            Rate = propensities(S[i][-1][-1], E[i][-1][-1], I[i][-1][-1], N, beta, epsilon, gamma)

            uT  = np.random.random()
            uE  = np.random.random()
            uI  = np.random.random()
            
            NT  = poisson.ppf(uT, Rate[0]*tau)
            NE  = poisson.ppf(uE, Rate[1]*tau)
            NI  = poisson.ppf(uI, Rate[2]*tau)

            NT  = min(NT, S[i][-1][-1])
            NE  = min(NE, E[i][-1][-1])
            NI  = min(NI, I[i][-1][-1])

            S[i][0] = np.append(S[i][0], S[i][-1][-1] - NT)
            E[i][0] = np.append(E[i][0], E[i][-1][-1] + NT - NE)
            I[i][0] = np.append(I[i][0], I[i][-1][-1] + NE - NI)
            R[i][0] = np.append(R[i][0], R[i][-1][-1] + NI)

            Time[i][0] = np.append(Time[i][0], Time[i][-1][-1] + tau)
            timetime = Time[i][-1][-1]

    Time_ave = Time[0][0]
    S_ave = S.mean(axis=0)[0]
    E_ave = E.mean(axis=0)[0]
    I_ave = I.mean(axis=0)[0]
    R_ave = R.mean(axis=0)[0]
    output = np.column_stack((S_ave, E_ave, I_ave, R_ave))

    return output

def tauleap_SIRP(propensities, y0, tau, runtime, N, beta_d, beta_i, gamma, alpha, mu):

    n_simulations = 100
    Time = np.zeros([n_simulations, 1],dtype=tuple)
    S = np.zeros([n_simulations, 1],dtype=tuple)
    I = np.zeros([n_simulations, 1],dtype=tuple)
    R = np.zeros([n_simulations, 1],dtype=tuple)
    P = np.zeros([n_simulations, 1],dtype=tuple)

    for i in range(n_simulations):
        Time[i][0] = [0.0]
        S[i][0] = [y0[0]]
        I[i][0] = [y0[1]]
        R[i][0] = [y0[2]]
        P[i][0] = [y0[3]]

    for i in range(n_simulations):

        timetime = 0
        while timetime < runtime:

            Rate = propensities(S[i][-1][-1], I[i][-1][-1], P[i][-1][-1], N, beta_d, beta_i, gamma, alpha, mu)

            uT  = np.random.random()
            uSout2  = np.random.random()
            uI  = np.random.random()
            uPin  = np.random.random()
            uPout  = np.random.random()
            
            NT  = poisson.ppf(uT, Rate[0]*tau)
            NSout2  = poisson.ppf(uSout2, Rate[1]*tau)
            NI  = poisson.ppf(uI, Rate[2]*tau)
            NPin  = poisson.ppf(uPin, Rate[3]*tau)
            NPout  = poisson.ppf(uPout, Rate[4]*tau)

            NT  = min(NT + NSout2, S[i][-1][-1])
            NI  = min(NI, I[i][-1][-1])
            NPout  = min(NPout, P[i][-1][-1])

            S[i][0] = np.append(S[i][0], S[i][-1][-1] - NT)
            I[i][0] = np.append(I[i][0], I[i][-1][-1] + NT - NI)
            R[i][0] = np.append(R[i][0], R[i][-1][-1] + NI)
            P[i][0] = np.append(P[i][0], P[i][-1][-1] + NPin - NPout)

            Time[i][0] = np.append(Time[i][0], Time[i][-1][-1] + tau)
            timetime = Time[i][-1][-1]

    Time_ave = Time[0][0]
    S_ave = S.mean(axis=0)[0]
    I_ave = I.mean(axis=0)[0]
    R_ave = R.mean(axis=0)[0]
    P_ave = P.mean(axis=0)[0]

    output = np.column_stack((S_ave, I_ave, R_ave, P_ave))

    return output

def solve_ODE(model, S0s):

    time = np.linspace(0, runtime, runtime)
    tau = time[1] - time[0]
    N = 100

    if model == 'SIR':
        beta = 0.5
        gamma = 0.2
        y0 = [S0s[0], N-S0s[0], 0]

        def SIR_propensities(S, I, N, beta, gamma):
            dSdt = beta * S * I / N
            dIdt = gamma * I
            return [dSdt, dIdt]

        sol = tauleap_SIR(SIR_propensities, y0, tau, runtime, N, beta, gamma)
        return time, sol, ['S', 'I', 'R']

    elif model == 'SIRS':
        beta = 0.5
        gamma = 0.2
        phi = 0.2
        y0 = [S0s[1], N-S0s[1], 0]

        def SIRS_propensities(S, I, R, N, beta, gamma, phi):
            dSdt = beta * S * I / N
            dIdt = gamma * I
            dRdt = phi * R
            return [dSdt, dIdt, dRdt]

        sol = tauleap_SIRS(SIRS_propensities, y0, tau, runtime, N, beta, gamma, phi)
        return time, sol, ['S', 'I', 'R']

    elif model == 'SEIR':
        beta = 0.5
        epsilon = 0.5
        gamma = 0.2
        vital = 0.0
        y0 = [S0s[2], 0, N-S0s[2], 0]  # S, E, I, R

        def SEIR_propensities(S, E, I, N, beta, epsilon, gamma):
            dSdt = beta * S * I / N
            dEdt = epsilon * E
            dIdt = gamma * I
            return [dSdt, dEdt, dIdt]

        sol = tauleap_SEIR(SEIR_propensities, y0, tau, runtime, N, beta, epsilon, gamma)
        return time, sol, ['S', 'E', 'I', 'R']

    elif model == 'SIRP':
        beta_d = 0.25
        beta_i = 0.3 / N
        gamma = 0.05
        alpha = 0.1
        mu = 0.15
        y0 = [S0s[3], N-S0s[3], 0, 10] # S, I, R, P

        def SIRP_propensities(S, I, P, N, beta_d, beta_i, gamma, alpha, mu):
            dSdt = beta_d * S * I / N
            Sout2 = beta_i * S * P
            dIdt = gamma * I
            Pin  = alpha * I
            Pout = mu * P
            return [dSdt, Sout2, dIdt, Pin, Pout]

        sol = tauleap_SIRP(SIRP_propensities, y0, tau, runtime, N, beta_d, beta_i, gamma, alpha, mu)
        return (time), sol, ['S', 'I', 'R', 'P']

    else:
        raise ValueError(f"Unknown model type: {model}")

def load_fixed_simulation_data(base_path, folder_prefix, compartments, env=False):
    column = "States_AfterAnyAddedInds"
    all_data = {c: [] for c in compartments}
    for i in range(100):
        folder_name = f"{folder_prefix}{i}species0"
        folder_path = os.path.join(base_path, folder_name)
        file_path = os.path.join(folder_path, "summary_popAllTime_DiseaseStates.csv")

        if os.path.isfile(file_path):
            df = pd.read_csv(file_path)
            temp = {c: [] for c in compartments}
            for j in range(runtime+1):
                try:
                    if "SEIR" in base_path:
                        values = df[column][j].split('|')[0].split(';')
                        for k, c in enumerate(compartments):
                            temp[c].append(float(values[k]))
                    elif "SIRP" in base_path:
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
            if all(len(temp[c]) == runtime+1 for c in compartments):
                for c in compartments:
                    all_data[c].append(temp[c])

    means = {c: np.mean(np.array(all_data[c]), axis=0) if all_data[c] else np.zeros(runtime+1) for c in compartments}
    return means, all_data

model_configs = {
    "SIR": {
        "base_path": r"Figure_1_summary_data/SIR",
        "prefix": "run0batch0mc",

        "comps": ['S', 'I', 'R']
    },
    "SIRS": {
        "base_path": r"Figure_1_summary_data/SIRS",
        "prefix": "run0batch0mc",
        "comps": ['S', 'I', 'R']
    },
    "SEIR": {
        "base_path": r"Figure_1_summary_data/SEIR",
        "prefix": "run0batch0mc",
        "comps": ['S', 'E', 'I', 'R']
    },
    "SIRP": {
        "base_path": r"Figure_1_summary_data/SIRP",
        "prefix": "run0batch0mc",
        "comps": ['S', 'I', 'R', 'P']
    }
}

fig, axs = plt.subplots(2, 2, figsize=(12, 12), sharex=True, sharey=True)
labels = ['a', 'b', 'c', 'd']
colors = {
    'S': '#56B4E9',     # Sky Blue
    'E': '#E69F00',     # Orange
    'I': '#D55E00',     # Vermillion
    'R': '#009E73',     # Bluish Green
    'P': '#CC79A7'      # Reddish Purple
}

runtime = 50
timesteps = range(runtime+1)

S0s = []
for ax, (model, cfg), label in zip(axs.flat, model_configs.items(), labels):
    means, all_data = load_fixed_simulation_data(cfg["base_path"], cfg["prefix"], cfg["comps"], env=(model == "SIRP"))
    S0s.append(means['S'][0])
    t_ode, sol, comps = solve_ODE(model, S0s)

    file_path = f"{cfg['base_path']}/{model}_results_"

    np.save(file_path+"means.npy", means)
    np.save(file_path+"all_data.npy", all_data)
    np.save(file_path+"t_ode.npy", t_ode)
    np.save(file_path+"sol.npy", sol)
    np.save(file_path+"comps.npy", comps)

    for comp in cfg["comps"]:
        for trace in all_data[comp]:
            ax.plot(timesteps, trace, linestyle='--', color=colors[comp], alpha=0.1, linewidth=0.5)

        idx = comps.index(comp)
        ax.plot(t_ode, sol[:, idx], color=colors[comp], linewidth=2.5)
        ax.plot(timesteps, means[comp], linestyle='--', color=colors[comp], linewidth=2)
        ax.tick_params(axis='both', labelsize=14)  

    ax.text(0.02, 0.9, f"({label})", transform=ax.transAxes, fontsize=18, fontweight='bold')
    ax.grid(True)

color_legend = [
    Line2D([0], [0], color='#56B4E9', lw=2, label='S'),
    Line2D([0], [0], color='#E69F00', lw=2, label='E'),
    Line2D([0], [0], color='#D55E00', lw=2, label='I'),
    Line2D([0], [0], color='#009E73', lw=2, label='R'),
    Line2D([0], [0], color='#CC79A7', lw=2, label='P')
]

style_legend = [
    Line2D([0], [0], color='black', lw=2, linestyle='-', label='ODE'),
    Line2D([0], [0], color='black', lw=2, linestyle='--', label='Simulation')
]

fig.legend(color_legend + style_legend, 
           [h.get_label() for h in color_legend + style_legend],
           loc='upper center', ncol=7, fontsize=16, frameon=False)

fig.supxlabel("Year", fontsize=16, y=0.03)
fig.supylabel("Population size", fontsize=16)
plt.tight_layout(rect=[0, 0.05, 1, 0.93])
plt.savefig("figure_outputs/Figure_1.png")
plt.close()

