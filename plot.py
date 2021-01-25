import json
from matplotlib import pyplot as plt

with open("Data.json","r") as infile:
    data = json.load(infile)

Ts =[]
Ms = []
M_error = []
Es = []
E_error = []
xi = []
xi_error = []
C = []
C_error = []
for x in data['measurements'].values():
    Ts.append(x['T'])
    Ms.append(x['M_mean'])
    M_error.append(x['M_error'])
    Es.append(x['E_mean'])
    E_error.append(x['E_error'])
    xi_error.append(x['xi_error'])
    C_error.append(x['C_error'])
    xi.append(x['xi'])
    C.append(x['C'])
    

fig,axs = plt.subplots(2,2,sharex=True)
fig.suptitle("Temperature progression plots")
axs[0,0].errorbar(Ts,Ms,yerr=M_error,ecolor="r",capsize=1,barsabove=True)
axs[0,0].set_title("Average Magnetism")
axs[0,0].set_xlabel("T")
axs[0,0].set_ylabel("M")
axs[0,1].errorbar(Ts,xi,yerr=xi_error,ecolor="r",capsize=1,barsabove=True)
axs[0,1].set_title("Susceptibility")
axs[0,1].set_xlabel("T")
axs[0,1].set_ylabel(r"$\chi$")
axs[1,0].errorbar(Ts,Es,yerr=E_error,ecolor="r",capsize=1,barsabove=True)
axs[1,0].set_title("Average Energy")
axs[1,0].set_xlabel("T")
axs[1,0].set_ylabel("Energy")
axs[1,1].errorbar(Ts,C,yerr=C_error,ecolor="r",capsize=1,barsabove=True)
axs[1,1].set_title("Capacity")
axs[1,1].set_xlabel("T")
axs[1,1].set_ylabel("C")
fig.set_size_inches(8,6)
plt.savefig("TempVariation.png",dpi=100)
plt.show()
