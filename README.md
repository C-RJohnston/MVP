Ising Model
===========
My solution to the Ising Model checkpoint for Modelling and Visualisation in Physics. This solution
uses three separate python files to animate, measure, and plot the model. Premade data and plots are provided in
Plots_data.zip but new data and plots can be generated using the instructions below.

Contents
------------
1. Prerequisites
2. Animating the model
3. Taking measurements
-3.1 File format
4. Plotting results

1.Prerequisites
-------------
This model requires:
Python 3.x
Numpy
Matplotlib

2.Animating the model
-------------------
Visualisation of the model is handled by the file animation.py. animation.py takes 5 command line arguments
lx: grid width
ly: grid height
T: temperature to animate at
Dynamic: model to animate (formatted as g or k)
sweeps: number of sweeps to animate for
e.g.

python animation.py 50 50 1 g 1000

The model will prerender the frames of the animation and will play the number
of sweeps on loop until the window is closed. The file does not save the animation to prevent potential issues
running the file across different operating systems.
For the purposes of animation, the model will cache every sweep - meaning that calculating 10,000 
sweeps for animation may be more resource intensive than calculating 10,000 sweeps during the measurements.

3.Taking measurements
-------------------
Measuring the model for the two different dynamics is handled in the file simulate.py. The file will measure
glauber first and then kawasaki. This file takes 7 command line arguments
lx: grid width
ly: grid height
T0: starting temperature
Tf: end temperature
Ts: number of different temperature measurements
runs: number of full sweeps
tau: autocorrelation time
e.g.

python simulate.py 50 50 1 3 20 10_000 10

On completion, the file will generate two JSON files Glauber_Data.json and Kawasaki_Data.json, the format of
which is explained in 3.1. 
To take measurements, the model first creates a lattice of the given size for each temperature and performs 100
uncached sweeps to bring the system to equilibrium. For the glauber method, the system starts with the initial
condition of all states being up, while Kawasaki starts with all states being random. Though there are initial
conditions which may improve the speed of equilibrating kawasaki, the time to bring the system up to equilibrium
was not a noticable issue.

3.1 File format
---------------
The data is stored in a JSON file with the format
{
parameters:
{
	number of spins,
	autocorrelation time
}
measurements:
{
	variables from the experiment
}
critical temperature:
{
	critical temperature determined from measurements
}
}
n.b. the critical temperature values will not appear until after running plot.py
4. Plotting results
-------------------
the results from running simulate.py can be plotted by running plot.py. This file takes 2 command line arguments:
method: which modelling method to use (g or k)
error_method: which error method to use (b or j)
e.g.

python plot.py g b

the file will present all plots in a single image and save the image as a png on completion. Error bars are marked
on the graph using the user specified and the critical temperature is marked on the graph in green.