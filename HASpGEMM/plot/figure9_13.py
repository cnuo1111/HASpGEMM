#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(16,16))
gs = gridspec.GridSpec(nrows=5, ncols=1, height_ratios=[3, 1, 1, 1, 1], hspace=0.5)

ax = fig.add_subplot(gs[0, :])

gs2 = gridspec.GridSpecFromSubplotSpec(nrows=4, ncols=1, subplot_spec=gs[1:, :], hspace=0.3)
ax1 = fig.add_subplot(gs2[0, 0])

ax2 = fig.add_subplot(gs2[1, 0])

ax3 = fig.add_subplot(gs2[2, 0])

ax4 = fig.add_subplot(gs2[3, 0])

ax.set_position([0.1, 0.7, 0.8, 0.2])
ax1.set_position([0.1, 0.45, 0.8, 0.1])
ax2.set_position([0.1, 0.3, 0.8, 0.1])
ax3.set_position([0.1, 0.15, 0.8, 0.1])
ax4.set_position([0.1, 0.0, 0.8, 0.1])

ax.grid(color='grey', ls = '-.', lw = 0.2)
yy = range(0 , 25 , 5)
xx = range(-8 , 8, 1)

filename = '13-spa.csv'
a = 2585
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = a
	colnum = 7
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

SPA = [0.0 for i in range(a)]
matrix1 = [0.0 for i in range(a)]

for cycle in range(a):
	matrix1[cycle] = math.log(abs((float)(content[cycle][4])),10)
	SPA[cycle] = (float)(content[cycle][5])
	
filename = '13-hash.csv'
b = 2585
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = b
	colnum = 7
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

HASH = [0.0 for i in range(b)]
matrix2 = [0.0 for i in range(b)]

for cycle in range(b):
	matrix2[cycle] = math.log(abs((float)(content[cycle][4])),10)
	HASH[cycle] = (float)(content[cycle][6])
	
filename = '13-esc.csv'
c = 2585
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = c
	colnum = 7
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

ESC = [0.0 for i in range(c)]
matrix3 = [0.0 for i in range(c)]

for cycle in range(c):
	matrix3[cycle] = math.log(abs((float)(content[cycle][4])),10)
	ESC[cycle] = (float)(content[cycle][6])
	
filename = '13-all.csv'
d = 2584
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = d
	colnum = 7
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

HASpGEMM = [0.0 for i in range(d)]
matrix4 = [0.0 for i in range(d)]

for cycle in range(d):
	matrix4[cycle] = math.log(abs((float)(content[cycle][4])),10)
	HASpGEMM[cycle] = (float)(content[cycle][6])
	
filename = '13-mkl.csv'
e = 2627
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = e
	colnum = 5
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

mkl = [0.0 for i in range(e)]
matrix5 = [0.0 for i in range(e)]

for cycle in range(e):
	matrix5[cycle] = math.log(abs((float)(content[cycle][4])),10)
	mkl[cycle] = (float)(content[cycle][6])

ax.scatter(matrix1, SPA, c = '#fd9490', marker = 'o', s = 12,label='SPA')
ax.scatter(matrix2, HASH, c = '#f2cf44', marker = 'o', s = 12,label='HASH')
ax.scatter(matrix3, ESC, c = '#aac28c', marker = 'o', s = 12,label='ESC')
ax.scatter(matrix5, mkl, c = 'green', marker = 'o', s = 12,label='oneMKL')
ax.scatter(matrix4, HASpGEMM, c = 'blue', marker = 'o', s = 12,label='HASpGEMM (our work)')

plt.subplots_adjust()
ax.legend(fontsize = 20,loc='upper left')
ax.set_ylabel('Performance \n(Gflops)', fontsize = 34)


filename = 'result_SPA.csv'
#a = 1864
a = 1708
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = a
	colnum = 6
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

#fig = plt.figure(figsize=(15, 10))
ax1.grid(color='grey', ls = '-.', lw = 0.2)
yy = range(0 , 15, 5)
xx = range(-1 , 8, 1)
SPA = [0.0 for i in range(a)]
matrix = [0.0 for i in range(a)]
ax1.axhline(y=1,c='red',ls="--")

for cycle in range(a):
	matrix[cycle] = math.log(abs((float)(content[cycle][2])),10)
	SPA[cycle] = (float)(content[cycle][5])
ax1.scatter(matrix, SPA, c = 'blue', marker = 'o', s = 12)
ax1.set_ylabel('Speedup\nHASpGEMM\n over SPA', fontsize = 20)

filename = 'result_HASH.csv'
#b = 2482
b = 2394
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = b
	colnum = 6
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

HASH = [0.0 for i in range(b)]
matrix2 = [0.0 for i in range(b)]

for cycle in range(b):
	matrix2[cycle] = math.log(abs((float)(content[cycle][2])),10)
	HASH[cycle] = (float)(content[cycle][5])
ax2.axhline(y=1,c='red',ls="--")
ax2.grid(color='grey', ls = '-.', lw = 0.2)
ax2.scatter(matrix2, HASH, c = 'blue', marker = 'o', s = 12)
ax2.set_ylabel('Speedup\nHASpGEMM\n over HASH', fontsize = 20)

filename = 'result_ESC.csv'
#c = 2467
c = 2156
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = c
	colnum = 6
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

ESC = [0.0 for i in range(c)]
matrix3 = [0.0 for i in range(c)]

for cycle in range(c):
	matrix3[cycle] = math.log(abs((float)(content[cycle][2])),10)
	ESC[cycle] = (float)(content[cycle][5])
ax3.axhline(y=1,c='red',ls="--")
ax3.grid(color='grey', ls = '-.', lw = 0.2)
ax3.scatter(matrix3, ESC, c = 'blue', marker = 'o', s = 12)
ax3.set_ylabel('Speedup\nHASpGEMM\n over ESC', fontsize = 20)

filename = 'result_mkl.csv'
#d = 2084
d = 1881
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = d
	colnum = 6
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    
        
for num in range(rownum):
	content[num] = rows[num]

mkl = [0.0 for i in range(d)]
matrix4 = [0.0 for i in range(d)]
ax4.axhline(y=1,c='red',ls="--")

for cycle in range(d):
	matrix4[cycle] = math.log(abs((float)(content[cycle][2])),10)
	mkl[cycle] = (float)(content[cycle][5])
ax4.axhline(y=1,c='red',ls="--")
ax4.grid(color='grey', ls = '-.', lw = 0.2)
ax4.scatter(matrix4, mkl, c = 'blue', marker = 'o', s = 12)
ax4.set_ylabel('Speedup\nHASpGEMM\n over oneMKL', fontsize = 20)

ax4.set_xlabel('intermediate products (log scale)', fontsize = 30)

ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])
ax.set_xlim([-1, 7])
ax4.set_xlim([-1, 7])
ax.set_ylim([0, 25])
ax1.set_ylim([0, 10])
ax2.set_ylim([0, 10])
ax3.set_ylim([0, 10])
ax4.set_ylim([0, 10])
ax.set_xticklabels([-1,0,1,2,3,4,5,6,7],fontsize = 18)
ax4.set_xticklabels([-1,0,1,2,3,4,5,6,7],fontsize = 18)
ax.set_yticklabels([0,5,10,15,20,25],fontsize = 18)
ax1.set_yticklabels([0,2,4,6,8,10],fontsize = 18)
ax2.set_yticklabels([0,2,4,6,8,10],fontsize = 18)
ax3.set_yticklabels([0,2,4,6,8,10],fontsize = 18)
ax4.set_yticklabels([0,2,4,6,8,10],fontsize = 18)
plt.subplots_adjust()
#ax.legend(fontsize = 10,loc='upper left')
#plt.ylabel("Speedup\nDASP over\ncuSPARSE",fontsize=32,labelpad=20)
#ax.set_xlabel('nonzero of the matrix (log scale)', fontsize = 30)
#ax.set_ylabel('Speedup\nHASpGEMM\n over HASH', fontsize = 30)
fig.savefig("cmpfinal_12.eps",dpi=600,format='eps',bbox_inches = 'tight')
plt.show()
plt.close()

