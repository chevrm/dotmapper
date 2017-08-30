import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

xlab = sys.argv[1]
ylab = sys.argv[2]
out = sys.argv[3]
xl = int(sys.argv[4])
yl = int(sys.argv[5])

m = np.loadtxt("./r.csv",delimiter=",",skiprows=1,dtype=np.float) ## data is floats, skip header row
fig,ax = plt.subplots() 
f1 = ax.pcolormesh(m[:,1:], cmap='jet')  ## heatmap generated, : grabs all rows, 1: grabs all but first col 
fig.colorbar(f1) 
plt.xlabel(xlab) 
plt.ylabel(ylab) 
plt.xlim([0,xl]) 
plt.ylim([0,yl])
fig.savefig(out, dpi=600)
