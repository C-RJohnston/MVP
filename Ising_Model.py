import random
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


#params
WIDTH = 50
HEIGHT = 50
RED = [1,0,0,1]
BLUE = [0,0,1,1]
vals = np.array([RED,BLUE])
cmap = ListedColormap(vals)

class Lattice:
    def __init__(self,WIDTH:int,HEIGHT:int,UP_COLOUR:list,DOWN_COLOUR:list):
        """
        default constructor for the lattice, initialise with

        WIDTH: int width of the grid

        HEIGHT: height of the grid

        UP_COLOUR: colour with which the up spin is drawn

        DOWN_COLOUR: colour with which the down spin is drawn
        """
        self._X = WIDTH
        self._Y = HEIGHT        
        #generate a random spin lattice of -1 and 1
        #first generate a random int array of 0s and 1s of size WIDTH,HEIGHT
        self._spins = np.random.randint(2,size=(self._Y,self._X))
        #replace 0 with -1
        self._spins = np.where(self._spins==0,-1,self._spins)
        #define the colour map for drawing the array
        self.cmap = ListedColormap(np.array([DOWN_COLOUR,UP_COLOUR]))

    def __iter__(self):
        return iter(self._spins)
    def __getitem__(self,index):
        return self._spins[index]
    def __str__(self):
        unique, counts = np.unique(self._spins, return_counts=True)
        u = dict(zip(unique, counts))
        return f"Lattice of size ({self._X},{self._Y})\nEnergy (arb units): {self.calc_total_energy()}\nNo. Up spins: {u[1]}\nNo. Down spins:{u[-1]}"

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
        spin = self._spins[y][x]
        neighbours = self.get_neighbours(x,y)

        #determine initial energy pre spin-flip
        E_i = 0
        for neighbour in neighbours:
            E_i+=-spin*neighbour
        #flip the spin
        spin*=-1
        #determine the new energy
        E_f = 0
        for neighbour in neighbours:
            E_f+=-spin*neighbour
        self._spins[y][x]*=-1
        return E_f-E_i

    def metropolis(self):
        """obtain a set of states from the Boltzmann distribution using the Metropolis Algorithm and Markov Chains"""
        #choose a random spin within the grid
        x = random.randrange(self._X)
        y = random.randrange(self._Y)

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

    def draw(self):
        """
        Draws the lattice in its current state
        """
        fig,ax = plt.subplots()
        ax.pcolormesh(self._spins,cmap = self.cmap)
        plt.show()
    


L = Lattice(50,50,[0,0,1,1],[1,0,0,1])
