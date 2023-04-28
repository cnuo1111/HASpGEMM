import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import csv

df = pd.read_csv('12-analyze_time.csv')

df['num1_ratio'] = df['num1'] / (df['num1'] + df['num2'])
df['num2_ratio'] = df['num2'] / (df['num1'] + df['num2'])

fig, ax = plt.subplots()

ax.bar(df['nnz'], df['num2_ratio'], color='#a2bdf6', label='Analysis time')
ax.bar(df['nnz'], df['num1_ratio'], bottom=df['num2_ratio'], color='#edeb77', label='Execution time')

plt.subplots_adjust()
ax.legend(fontsize = 15)
ax.set_ylabel('Ratio of Analysis time\n and Execution time',fontsize = 25)
ax.set_xlabel('matrices',fontsize = 25)
plt.yticks(fontsize = 15)
plt.gca().set_xticklabels('',rotation=60)
fig.savefig("analyse.eps",bbox_inches = 'tight')
plt.show()
