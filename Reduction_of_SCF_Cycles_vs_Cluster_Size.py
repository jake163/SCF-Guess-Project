import matplotlib.pyplot as plt

water_size = [3,4,5,7,8,9,10]
water = [3,2,2,2,1,1,1]
CO2_size = [4,5,6,7,8,9,10]
CO2 = [5,5,4,4,4,4,4]
HF_size = [4,5,6,8]
HF = [2,2,2,1]

fig, ax = plt.subplots()
ax.plot(water_size, water, label='H2O')
ax.plot(CO2_size, CO2, label='CO2')
ax.plot(HF_size, HF, label='HF')
ax.set_xlabel('Cluster Size')
ax.set_ylabel('Reduction of Cycles')
ax.set_title('Reduction of SCF Cycles vs. Cluster Size')
ax.legend()
plt.savefig('Reduction_of_SCF_Cycles_vs_Cluster_Size.pdf')