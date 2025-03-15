
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np



df = pd.read_table("GCMC_VS/particle_number.txt")
i = df.iloc[:, 0]
count = df.iloc[:,1]
n_h = df.iloc[:,2]
n_v = df.iloc[:,3]

diff = []
for j in range(len(i)):
	diff.append((n_h[j] - n_v[j]))


fig, ax = plt.subplots()
ax.plot(i,count)
plt.show()

fig, ax = plt.subplots()
ax.plot(i,diff)
plt.show()



