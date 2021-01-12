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


#generate a random spin lattice of -1 and 1

#first generate a random int array of 0s and 1s of size WIDTH,HEIGHT
spins = np.random.randint(2,size=(HEIGHT,WIDTH))
#replace 0 with -1
spins = np.where(spins==0,-1,spins)



plt.imshow(spins,cmap)
plt.show()