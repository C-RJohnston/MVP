from Ising_Model import Lattice
import copy
import numpy as np
from matplotlib import pyplot as plt
import timeit
import json
#params

def equilibrate(states: Lattice,method: str):
    """
    returns a lattice that has been allowed to reach equilibrium (100 sweeps) using either
    Glauber method or Kawasaki method
    """
    new = copy.copy(states)
    if method.lower() == "g" or method.lower() == "glauber":
        new.sim_glauber(100)
    elif method.lower() == "k" or method.lower() == "kawasaki":
        new.sim_kawasaki(100)
    return new

def susceptibility(av_M,av_M2,N,T):
    return 1/(N*T)*(av_M2-av_M**2)


def capacity(av_E,av_E2,N,T):
    return 1/(N*T**2)*(av_E2-av_E**2)


def main():
    Ts = np.linspace(1,5,20).tolist()
    Ms = []
    xi = []
    Es = []
    C = []
    i=0
    for T in Ts:
        i+=1
        lx,ly = 50,50
        L = Lattice(lx,ly,T,"up")
        L = equilibrate(L,"glauber")
        L.sim_glauber(1000,True)
        E,M = L.get_measurements()
        del L
        E = E[::10]
        M = M[::10]
        av_M = np.mean(M)
        av_M2 = np.mean(np.power(M,2))
        av_E = np.mean(E)
        av_E2 = np.mean(np.power(E,2))
        Ms.append(av_M)
        Es.append(av_E)
        xi.append(susceptibility(av_M,av_M2,lx*ly,T))
        C.append(capacity(av_E,av_E2,lx*ly,T))
        print("completed temp",f"{i}/20")
    Data = {"Temperature":Ts,
            "Magnetisation":Ms,
            "Susceptibility":xi,
            "Energy":Es,
            "Capacity":C}
    with open("Data.json",'w') as outfile:
        json.dump(Data,outfile)
    fig,axs = plt.subplots(2,2,sharex=True)
    fig.suptitle("Temperature progression plots")
    axs[0,0].plot(Ts,Ms)
    axs[0,0].set_title("Average Magnetism")
    axs[0,0].set_xlabel("T")
    axs[0,0].set_ylabel("M")
    axs[0,1].plot(Ts,xi)
    axs[0,1].set_title("Susceptibility")
    axs[0,1].set_xlabel("T")
    axs[0,1].set_ylabel(r"$\chi$")
    axs[1,0].plot(Ts,Es)
    axs[1,0].set_title("Average Energy")
    axs[0,0].set_xlabel("T")
    axs[0,0].set_ylabel("Energy")
    axs[1,1].plot(Ts,C)
    axs[1,1].set_title("Capacity")
    axs[0,0].set_xlabel("T")
    axs[0,0].set_ylabel("C")
    fig.set_size_inches(8,6)
    plt.savefig("TempVariation.png",dpi=100)
    plt.show()

main()



