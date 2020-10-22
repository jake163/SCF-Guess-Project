import matplotlib.pyplot as plt

water_size = [3,4,5,7,8,9,10]
water = [-1.161377192, -2.87171483, -3.158443689, 31.399617434, 28.013007879, 49.329758406, 70.815126657]
CO2_size = [4,5,6,7,8,9,10]
CO2 = [60.525092363, 142.536129951, 224.71994257, 396.902711391, 625.431005239, 908.067396164, 1319.410159349]
HF_size = [2,4,5,6,8]
HF = [0.093206644, -1.766253471, -2.720202923, -4.283352137, 6.808915615]

fig, ax = plt.subplots()
ax.plot(water_size, water, label='H2O')
ax.plot(CO2_size, CO2, label='CO2')
ax.plot(HF_size, HF, label='HF')
ax.set_xlabel('Cluster Size')
ax.set_ylabel('SCF Convergence Time Reduction (seconds)')
ax.set_title('SCF Convergence Time Reduction vs. Cluster Size')
ax.legend()
plt.savefig('SCF_Convergence_Time_Reduction_vs_Cluster_Size.pdf')