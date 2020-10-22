import matplotlib.pyplot as plt

water_size = [3,4,5,6,7,8,9,10]
water_di = [1.467879295349121, 2.85368275642395, 4.681039094924927, 7.32571220397949, 10.033047199249268, 
             13.033649206161499, 16.516063451766968, 21.048356771469116]
CO2_size = [4,5,6,7,8,9,10]
CO2_di = [17.913692474365234, 29.259541511535645, 43.61529850959778, 78.68523359298706, 81.20058131217957, 
           101.81348252296448, 128.0662031173706]
HF_size = [2,4,5,6,7,8]
HF_di = [0.2975616455078125, 1.6088576316833496, 2.596971035003662, 3.830887794494629, 5.253161907196045, 
          7.038269519805908]

fig, ax = plt.subplots()
ax.plot(water_size, water_di, label='H2O')
ax.plot(CO2_size, CO2_di, label='CO2')
ax.plot(HF_size, HF_di, label='HF')
ax.set_xlabel('Cluster Size')
ax.set_ylabel('Time to Compute Dimer Densities (seconds)')
ax.set_title('Time to Compute Dimer Densities vs. Cluster Size')
ax.legend()
plt.savefig('Time_to_Compute_Dimer_Densities_vs_Cluster_Size.pdf')