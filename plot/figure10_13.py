#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

yy = range(0 , 12 , 2)
filename = './DATA/output.csv'
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = 22
	colnum = 20
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    

for num in range(rownum):
	content[num] = rows[num]
	
x = [0.0 for i in range(22)]
SPA = [0.0 for i in range(22)]
HASH = [0.0 for i in range(22)]
ESC = [0.0 for i in range(22)]
oneMKL = [0.0 for i in range(22)]
HASpGEMM = [0.0 for i in range(22)]

fig = plt.figure(figsize=(30, 4))
ax1 = fig.add_subplot(1,1,1)
plt.grid(color='grey', ls = '-.', lw = 0.2)

for cycle in range(22): 
	x[cycle] = content[cycle][1]
	SPA[cycle] = (float)(content[cycle][3])
	HASH[cycle] = (float)(content[cycle][8])
	ESC[cycle] = (float)(content[cycle][11])
	oneMKL[cycle] = (float)(content[cycle][15])
	HASpGEMM[cycle] = (float)(content[cycle][19])
		
x_len = np.arange(len(x))
total_width, n = 0.8, 5
width = total_width/n
xticks = x_len - (total_width - width) / 5

plt.bar(xticks, SPA, width=0.9*width, label="SPA", color="blue",edgecolor='black',linewidth = 1.1,  zorder=10)
plt.bar(xticks + 1*width, HASH, width=0.9*width, label="HASH", color="green",edgecolor='black',linewidth = 1.1, zorder=10)
plt.bar(xticks + 2*width, ESC, width=0.9*width, label="ESC", color="red",edgecolor='black',linewidth = 1.1, zorder=10)
plt.bar(xticks + 3*width, oneMKL, width=0.9*width, label="oneMKL", color="orange",edgecolor='black',linewidth = 1.1, zorder=10)
plt.bar(xticks + 4*width, HASpGEMM, width=0.9*width, label="HASpGEMM (our work)", color="gold",edgecolor='black',linewidth = 1.1, zorder=10)

plt.subplots_adjust()
ax1.legend(fontsize = 17,loc='upper right')
plt.xticks(x_len, x, fontsize = 25)
plt.gca().set_xticklabels(x, rotation=60)
plt.yticks(yy,fontsize = 20)
plt.ylabel("Performance\n (Gflops)",fontsize = 34)
fig.savefig("12_example.eps",dpi=600,format='eps',bbox_inches = 'tight')
plt.show()
plt.close()
