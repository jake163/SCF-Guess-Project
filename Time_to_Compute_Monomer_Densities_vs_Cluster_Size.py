import matplotlib.pyplot as plt

water_size = [3,4,5,6,7,8,9,10]
water_mono = [0.3838937282562256, 0.5223190784454346, 0.6900269985198975, 0.8127787113189697, 0.9290657043457031,
             1.0745518207550049, 1.210066556930542, 1.300506830215454]
CO2_size = [4,5,6,7,8,9,10]
CO2_mono = [1.7547690868377686, 1.635404348373413, 1.9593424797058105, 2.2897789478302, 5.07656455039978, 
           2.9627771377563477, 3.2577767372131348]
HF_size = [2,4,5,6,7,8]
HF_mono = [0.6911334991455078, 0.48136091232299805, 0.5698037147521973, 0.67104172706604, 0.797966480255127, 
          0.9404196739196777]

fig, ax = plt.subplots()
ax.plot(water_size, water_mono, label='H2O')
ax.plot(CO2_size, CO2_mono, label='CO2')
ax.plot(HF_size, HF_mono, label='HF')
ax.set_xlabel('Cluster Size')
ax.set_ylabel('Time to Compute Monomer Densities (seconds)')
ax.set_title('Time to Compute Monomer Densities vs. Cluster Size')
ax.legend()
plt.savefig('Time_to_Compute_Monomer_Densities_vs_Cluster_Size.pdf')