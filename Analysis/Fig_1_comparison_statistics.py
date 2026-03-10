import numpy as np
from permetrics.regression import RegressionMetric

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

    print("\n")
    print(model)
    for comp in cfg["comps"]:
        print(comp)
        idx = comps.index(comp)
        data = sol[:, idx]
        obs = means[comp][:-1]
        if comp == 'P':
            data = sol[:, idx][1:]
            obs = means[comp][:-2]

        evaluator = RegressionMetric(data, obs)
        print("(Pearson’s Correlation Index)**2: ", evaluator.pearson_correlation_coefficient_square())
        print("R2 - Coefficient of Determination: ",evaluator.coefficient_of_determination())
        print("NRMSE: ", evaluator.normalized_root_mean_square_error())
        print("NSE: ", evaluator.nash_sutcliffe_efficiency())
        print("KGE: ", evaluator.kling_gupta_efficiency())

        print("\n")















