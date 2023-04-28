import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.grid(color='grey', ls='-.', lw=0.2)
yy = range(0, 40, 5)
xx = range(1, 9, 1)

filename = './DATA/rail_info.csv'

with open(filename) as csvfile:
    reader = csv.reader(csvfile)
    rows = [row for row in reader]
    rownum = 16
    colnum = 5
    content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]

for num in range(rownum):
    content[num] = rows[num]

core = [0.0 for i in range(16)]
cub_HA = [0.0 for i in range(16)]
cub_sinple = [0.0 for i in range(16)]
cache_HA = [0.0 for i in range(16)]
cache_sinple = [0.0 for i in range(16)]

for cycle in range(16):
    core[cycle] = content[cycle][0]
    cub_HA[cycle] = (float)(content[cycle][1])
    cache_HA[cycle] = (float)(content[cycle][2])
    cub_sinple[cycle] = (float)(content[cycle][3])
    cache_sinple[cycle] = (float)(content[cycle][4])

#ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)

plt.plot(core, cub_HA, label="cub_HASpGEMM", marker='|', color='red')
plt.plot(core, cache_HA, label="cache_HASpGEMM", marker='*', color='red')
plt.plot(core, cub_sinple, label="cub_normal", marker='|', color='blue')
plt.plot(core, cache_sinple, label="cache_normal", marker='*', color='blue')

plt.subplots_adjust(top=0.85)
plt.legend(fontsize=15, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
plt.xlabel('core Id', fontsize=30)
plt.ylabel('Cache line cost or\n Intermediate products\n(*1e6)', fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
fig.savefig("experience_examplematrix.eps", dpi=600, format='eps', bbox_inches='tight')
plt.show()
plt.close()
