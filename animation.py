from Ising_Model import Lattice
import sys

def main():
    #retreive user parameters
    lx = int(sys.argv[1])
    ly = int(sys.argv[2])
    T = float(sys.argv[3])
    Dynamic = sys.argv[4]
    sweeps = int(sys.argv[5])
    #instantiate a system at the given temperature and size
    L = Lattice(lx,ly,T)
    if Dynamic.lower()=="g":
        L.sim_glauber(sweeps,True,1)
    elif Dynamic.lower()=="k":
        L.sim_kawasaki(sweeps,True,1)
    else:
        raise ValueError("Please specify a valid model (g or k)")
    L.anim([0,0,1,1],[1,0,0,1],sweeps) #specify the colour for the up state and the down state

main()