import matplotlib
matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

print(sys.argv)
datbs=pd.read_csv("a.txt");
datbs1=pd.read_csv("traj2.txt");


cmap = plt.cm.get_cmap('RdYlBu');

fig, axs = plt.subplots(2,1);

axs[0].plot(datbs.t/2/np.pi/1e6,datbs.e,label="N-body");
axs[0].plot(datbs1.time/1e6,datbs1.e,label="Secular");
axs[0].legend();

#axs[1].plot(datbs.t/2/np.pi,np.mod(datbs.Omega,2*np.pi),c="orange",label="N-body");
axs[1].plot(datbs.t/2/np.pi/1e6,np.mod(datbs.res,2*np.pi),label="N-body",linewidth=5);
axs[1].plot(datbs1.time/1e6,np.mod(datbs1.sigma,2*np.pi),label="Secular",linewidth=3);

axs[1].legend();


axs[0].set_ylabel(r'$e$')
axs[1].set_ylabel(r'$\sigma$')
axs[1].set_xlabel(r'$time(Myrs)$')


fig.subplots_adjust(wspace=0.4)

sp=[axs[0]]
plt.setp([a.get_xticklabels() for a in sp], visible=False)

filename="comp_int.png";
plt.savefig(filename)

