import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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

timesteps = range(51)
for ax, (model, cfg), label in zip(axs.flat, model_configs.items(), labels):
    
    file_path = f"{cfg['base_path']}/{model}_results_"
    means = np.load(file_path+'means.npy', allow_pickle=True).item()
    all_data = np.load(file_path+'all_data.npy', allow_pickle=True).item()
    t_ode = np.load(file_path+'t_ode.npy', allow_pickle=True)
    sol = np.load(file_path+'sol.npy', allow_pickle=True)
    comps = list(np.load(file_path+'comps.npy', allow_pickle=True))


    for comp in cfg["comps"]:
        for trace in all_data[comp]:
            ax.plot(timesteps, trace, linestyle='--', color=colors[comp], alpha=0.1, linewidth=0.5)

        idx = comps.index(comp)
        ax.plot(t_ode, sol[:, idx], color=colors[comp], linewidth=2.5)
        if comp == "P": #the time is shifted because CDMetaPOP currently reports P after the first year. (Will be fixed in future update)
            timesteps = np.array(timesteps)[:-1] + 1.0
            means[comp] = means[comp][:-1]
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



