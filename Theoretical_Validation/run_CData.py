import os
import pandas as pd
import matplotlib.pyplot as plt

# Path to the base folder
base_path = "C:\GitHub\CDMetaPOP_DiseaseDev\example_files\Defense_SIRD_batlogistic_deaths4_pt8beta_pt1rem_pt1D_1744921684"

plt.figure(figsize=(12, 6))

for folder in sorted(os.listdir(base_path)):
	folder_path = os.path.join(base_path, folder)
	file_path= os.path.join(folder_path, "summary_popAllTime_DiseaseStates.csv")

	if os.path.isfile(file_path):
		df = pd.read_csv(file_path)

		init_vals =str(df.iloc[0, 2]).split('|')
		d_vals =str(df.iloc[0, 5]).split('|')
		initial_S= sum(int(val.split(';')[0]) for val in init_vals if val)
		initial_I= sum(int(val.split(';')[1]) for val in init_vals if val)
		initial_R= sum(int(val.split(';')[2]) for val in init_vals if val)
		initial_D= sum(int(val.split(';')[3]) for val in d_vals if val)

		sus_series = []
		infected_series = []
		r_series = []
		d_series = []
		for val in df.iloc[1:, 6]:	# here skip t=0
			entries = str(val).split('|')
			dentries = str(val).split('|')
			total_S = sum(int(entry.split(';')[0]) for entry in entries if entry)
			total_I = sum(int(entry.split(';')[1]) for entry in entries if entry)
			total_R = sum(int(entry.split(';')[2]) for entry in entries if entry)
			total_D = sum(int(entry.split(';')[3]) for entry in entries if entry)
			sus_series.append(total_S)
			infected_series.append(total_I)
			r_series.append(total_R)
			d_series.append(total_D)
		sus_series = [initial_S] + sus_series
		infected_series = [initial_I] + infected_series
		r_series = [initial_R] + r_series
		d_series = [initial_D] + d_series

		# Plotting here
		plt.plot(range(len(sus_series)), sus_series, 'g--', alpha=0.2, label=folder)
		plt.plot(range(len(infected_series)), infected_series, 'r--', alpha=0.8)
		plt.plot(range(len(r_series)), r_series, 'b--', alpha=0.2)
		plt.plot(range(len(d_series)), d_series, 'k--', alpha=0.2)
		

	else:
		print(f"❌ File notfound: {file_path}")

plt.xlabel("Timestep(t)")
plt.ylabel("Total Population (summed over patches)")
plt.grid(True)
plt.legend(fontsize='small', loc='upper right')
plt.tight_layout()
plt.show()
