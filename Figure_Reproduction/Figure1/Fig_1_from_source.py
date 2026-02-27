import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.lines import Line2D


def solve_ODE(model):
    t_ode = np.linspace(0, 50, 1000)
    N = 100

    if model == 'SIR':
        beta = 0.5
        gamma = 0.2
        vital = 0.0
        y0 = [N - 10, 10, 0]

        def SIR_model(y, t, N, beta, gamma, vital):
            S, I, R = y
            dSdt = vital * N - beta * S * I / N - vital * S
            dIdt = beta * S * I / N - gamma * I - vital * I
            dRdt = gamma * I - vital * R
            return [dSdt, dIdt, dRdt]

        sol = odeint(SIR_model, y0, t_ode, args=(N, beta, gamma, vital))
        return t_ode, sol, ['S', 'I', 'R']

    elif model == 'SIRS':
        beta = 0.5
        gamma = 0.2
        phi = 0.2
        vital = 0.0
        y0 = [N - 10, 10, 0]

        def SIRS_model(y, t, N, beta, gamma, phi, vital):
            S, I, R = y
            dSdt = vital * N - beta * S * I / N + phi * R - vital * S
            dIdt = beta * S * I / N - gamma * I - vital * I
            dRdt = gamma * I - phi * R - vital * R
            return [dSdt, dIdt, dRdt]

        sol = odeint(SIRS_model, y0, t_ode, args=(N, beta, gamma, phi, vital))
        return t_ode, sol, ['S', 'I', 'R']

    elif model == 'SEIR':
        beta = 0.5
        epsilon = 0.5
        gamma = 0.2
        vital = 0.0
        y0 = [90, 0, 10, 0]  # S, E, I, R

        def SEIR_model(y, t, N, beta, epsilon, gamma, vital):
            S, E, I, R = y
            dSdt = -beta * S * I / N
            dEdt = beta * S * I / N - epsilon * E
            dIdt = epsilon * E - gamma * I
            dRdt = gamma * I
            return [dSdt, dEdt, dIdt, dRdt]

        sol = odeint(SEIR_model, y0, t_ode, args=(N, beta, epsilon, gamma, vital))
        return t_ode, sol, ['S', 'E', 'I', 'R']

    elif model == 'SIRP':
        beta_d = 0.25
        beta_i = 0.3 / N
        gamma = 0.05
        alpha = 0.1
        mu = 0.15
        y0 = [90, 10, 0, 10]  # S, I, R, P

        def SIRP_model(y, t, N, beta_d, beta_i, gamma, alpha, mu):
            S, I, R, P = y
            dSdt = -beta_d * S * I / N - beta_i * S * P
            dIdt = beta_d * S * I / N + beta_i * S * P - gamma * I
            dRdt = gamma * I
            dPdt = alpha * I - mu * P
            return [dSdt, dIdt, dRdt, dPdt]

        sol = odeint(SIRP_model, y0, t_ode, args=(N, beta_d, beta_i, gamma, alpha, mu))
        return t_ode, sol, ['S', 'I', 'R', 'P']

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
            for j in range(51):
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
            if all(len(temp[c]) == 51 for c in compartments):
                for c in compartments:
                    all_data[c].append(temp[c])

    means = {c: np.mean(np.array(all_data[c]), axis=0) if all_data[c] else np.zeros(51) for c in compartments}
    return means, all_data

# Plotting
model_configs = {
    "SIR": {
        "base_path": r"SIR",
        "prefix": "run0batch0mc",

        "comps": ['S', 'I', 'R']
    },
    "SIRS": {
        "base_path": r"SIRS",
        "prefix": "run0batch0mc",
        "comps": ['S', 'I', 'R']
    },
    "SEIR": {
        "base_path": r"SEIR",
        "prefix": "run0batch0mc",
        "comps": ['S', 'E', 'I', 'R']
    },
    "SIRP": {
        "base_path": r"SIRP",
        "prefix": "run0batch0mc",
        "comps": ['S', 'I', 'R', 'P']
    }
}



# Setup figure
fig, axs = plt.subplots(2, 2, figsize=(12, 12), sharex=True, sharey=True)

labels = ['a', 'b', 'c', 'd']
colors = {
    'S': '#56B4E9',     # Sky Blue
    'E': '#E69F00',     # Orange
    'I': '#D55E00',     # Vermillion
    'R': '#009E73',     # Bluish Green
    'P': '#CC79A7'      # Reddish Purple
}

timesteps = range(51)

for ax, (model, cfg), label in zip(axs.flat, model_configs.items(), labels):
    t_ode, sol, comps = solve_ODE(model)
    means, all_data = load_fixed_simulation_data(cfg["base_path"], cfg["prefix"], cfg["comps"], env=(model == "SIRP"))

    for comp in cfg["comps"]:
        for trace in all_data[comp]:
            ax.plot(timesteps, trace, linestyle='--', color=colors[comp], alpha=0.1, linewidth=0.5)

        idx = comps.index(comp)
        ax.plot(t_ode, sol[:, idx], color=colors[comp], linewidth=2.5)
        ax.plot(timesteps, means[comp], linestyle='--', color=colors[comp], linewidth=2)
        ax.tick_params(axis='both', labelsize=14)  

    # Only panel label, no model title
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
plt.show()

