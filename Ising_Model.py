import random
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib import animation
import copy
from timeit import timeit
import sys



class Lattice:
    def __init__(self,WIDTH:int,HEIGHT:int,T:float,INITIAL = "random"):
        """
        default constructor for the lattice, initialise with

        WIDTH: int width of the grid

        HEIGHT: height of the grid

        T: Temperature of the system

        INITIAL: how to generate the initial conditions of the lattice (all up, all down, random)
        """
        self._X = WIDTH
        self._Y = HEIGHT        
        #generate a random spin lattice of -1 and 1
        #first generate a random int array of 0s and 1s of size WIDTH,HEIGHT
        self._spins = self.create_lattice(INITIAL)
        #replace 0 with -1
        self._spins = np.where(self._spins==0,-1,self._spins)
        #define the colour map for drawing the array
        self.cache = []
        self.E = []
        self.M = []
        self.T = T

    def create_lattice(self,method):
        if(method.lower() == "up"):
            return np.ones((self._Y,self._X))
        elif(method.lower() == "down"):
            return -np.ones((self._Y,self._X))
        elif(method.lower() == "random"):
            return np.random.randint(2,size=(self._Y,self._X))

    @property
    def size(self):
        return self._X*self._Y
    def __iter__(self):
        return iter(self._spins)
    def __getitem__(self,index):
        return self._spins[index]
    def __repr__(self):
        unique, counts = np.unique(self._spins, return_counts=True)
        u = dict(zip(unique, counts))
        return f"""Lattice of size ({self._X},{self._Y})
        Magnetisation: {self.calc_total_magnetisation()}
        Energy: {self.calc_total_energy()}
        No. Up spins: {u[1]}
        No. Down spins:{u[-1]}"""

    def get_neighbours(self,x,y):
        #gets the on-lattice neighbours for spin(x,y) using periodic boundaries
        neighbours = [self._spins[y-1,x],self._spins[y,x-1]]
        if(x+1==self._X):
            neighbours.append(self._spins[y,0])
        else:
            neighbours.append(self._spins[y,x+1])
        if(y+1==self._Y):
            neighbours.append(self._spins[0,x])
        else:
            neighbours.append(self._spins[y+1,x])
        return neighbours

    def is_neighbour(self,x1,y1,x2,y2):
        """
        Returns true if (x1,y1) and (x2,y2) are neighbours using periodic boundaries
        if the distance between the points == 1 or N and the points or on axis, then they
        are neighbours
        """
        XMAX = self._X-1
        YMAX = self._Y-1
        sy = abs(y2-y1)
        sx = abs(x2-x1)
        #if the distance between x points is 1 (or lx) and the ys are aligned
        onX = (sx==1 or sx==XMAX) and sy==0
        #if the distance between y points is 1 (or ly) and the xs are aligned
        onY = (sy==1 or sy==YMAX) and sx==0
        return onX or onY


    def calc_delta_energy(self,x,y):
        """
        Calculate the change in energy of the state by flipping the given spin
        Returns delta E and flips spin(x,y)

        x: the row of the spin within the state
        y: the column of the spin within the state  
        
        The energy between two neighbouring spins i,j is given by
        E(i,j) = -J*S_i*S_j
        where J is normalised to 1 and S_x is the spin value of x [-1,1]
        The total energy of the system is thus the sum of these energies. The difference between the two states after flipping one spin 
        is given by the change in energy of the neighbourhood of the flipped spin
        """

        spin = self._spins[y,x]
        neighbours = self.get_neighbours(x,y)

        return 2*spin*sum(neighbours)


    def glauber_step(self):
        """obtain a set of states from the Boltzmann distribution using the glauber method"""
        #choose a random spin within the grid
        x = random.randrange(self._X)
        y = random.randrange(self._Y)
        dE = self.calc_delta_energy(x,y)
        #the probability that the spin should flip
        p = min(1,np.exp(-dE/self.T))
        if(random.random()<p):
            #flip the spin based on the determined probability
            self._spins[y,x]*=-1
        return self._spins

    def glauber_sweep(self):
        """completes a whole sweep of the glauber method"""
        for i in range(self._X*self._Y):
            self.glauber_step()
        return self._spins

    def sim_glauber(self,runs,cache=False,interval=10):
        """
        simulates the glauber method and caches the states
        """
        for r in range(runs):
            if cache:
                self.glauber_sweep()
                if(r%interval == 0):
                    self.cache.append(copy.copy(self._spins))
                    self.E.append(self.calc_total_energy())
                    self.M.append(self.calc_total_magnetisation())
                    #print(f"Temperature {self.T}: {r}/{runs}")
            else:
                self.glauber_sweep()
    
    def kawasaki_step(self):
        """obtain a state based off the kawasaki method"""
        #obtain a first random point
        x1 = random.randrange(self._X)
        y1 = random.randrange(self._Y)
        spin1 = self._spins[y1,x1]
        #find a second random point with a different spin
        spin2=spin1
        while(spin2==spin1):
            x2 = random.randrange(self._X)
            y2 = random.randrange(self._Y)
            spin2=self._spins[y2,x2]
        #determine the change in energy
        dE = self.calc_delta_energy(x1,y1)+self.calc_delta_energy(x2,y2)
        #check if the two points are neighbours
        if(self.is_neighbour(x1,y1,x2,y2)):
            #account for the two spins being neighbours
            dE+=4*spin2*spin1
        p = min(1,np.exp(-dE/self.T))
        if(random.random()<p):
            #flip the spin based on the determined probability
            self._spins[y1,x1]*=-1
            self._spins[y2,x2]*=-1
        return self._spins

    def kawasaki_sweep(self):
        """completes a whole sweep of the kawasaki method"""
        for i in range(self._X*self._Y):
            self.kawasaki_step()
        return self._spins

    def sim_kawasaki(self,runs,cache=False,interval=10):
        """
        simulates the kawasaki method and caches the states
        """
        for r in range(runs):
            if cache:
                self.kawasaki_sweep()
                if(r%interval==0):
                    self.cache.append(copy.copy(self._spins))
                    self.E.append(self.calc_total_energy())
                    self.M.append(self.calc_total_magnetisation())
            else:
                self.kawasaki_sweep()

    def calc_total_energy(self):
        """
        Calculates the total energy of the state
        """
        energy = 0
        for x in range(self._X):
            for y in range(self._Y):
                #gets the on-lattice neighbours for spin(x,y) using periodic boundaries
                neighbours = self.get_neighbours(x,y)
                for neighbour in neighbours:
                    energy+=-self._spins[y,x]*neighbour
        return energy

    def calc_total_magnetisation(self):
        """
        calculates the total magnetisation of the current state of the lattice using
        M = sum(spins)
        """
        return sum(sum(self._spins))
    
    def calc_average_magnetisation(self):
        """
        returns the magnetisation per number of spins
        """
        N = self._X*self._Y
        return self.calc_total_magnetisation()/N

    def get_measurements(self):
        return self.E,self.M

    def draw(self,UP_COLOUR:list,DOWN_COLOUR:list):
        """
        Draws the lattice in its current state
        """
        cols = ListedColormap(np.array([DOWN_COLOUR,UP_COLOUR]))
        fig,ax = plt.subplots()
        ax.pcolormesh(self._spins,cmap = cols)
        plt.show()


    def anim(self,UP_COLOUR:list,DOWN_COLOUR:list,steps):
        if(len(self.cache)>0):
            cols = ListedColormap(np.array([DOWN_COLOUR,UP_COLOUR]))
            fig,ax = plt.subplots()
            im = ax.imshow(self.cache[0],cmap=cols)
            def animate(i):
                im.set_array(self.cache[i])
                return im,

            a = animation.FuncAnimation(fig,animate,frames=steps,interval=1)
            plt.show()
            a.save("sim.gif",fps=60)
        else:
            print("Make sure to run a simulation first")



def main():
    if(len(sys.argv[1:])!=4):
        raise TypeError(f"Missing {4-len(sys.argv[1:])} required positional arguments: lx, ly, T, dynamics(G/K)")
    lx = int(sys.argv[1])
    ly = int(sys.argv[2])
    T = float(sys.argv[3])
    Dynamic = sys.argv[4]
    L = Lattice(lx,ly,T)
    if Dynamic=="G":
        L.sim_glauber(360)
    elif Dynamic=="K":
        L.sim_kawasaki(360)
    L.anim([0,0,1,1],[1,0,0,1],360)
    print(L)

if __name__ == "__main__":
    main()