from Ising_Model import Lattice
import copy
import numpy as np
from matplotlib import pyplot as plt
import timeit
import json
import time
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

def error(av,av2,N,tau=10):
    """
    returns the uncorrelated error for either the magnetism or the energy
    av = <x>
    av2 = <x^2>
    N = number of measurements
    tau = autocorrelation time (number of sweeps between each measurement), default 10
    """
    return ((av2-av**2)*(2*tau)/(N))**(1/2)

def bootstrap(measurements,value,k,N,T):
    """
    uses the bootstrap method to calculate the error on the capacity or the susceptibility
    measurements: measurements to resample (energy or temperature)
    value: value that is being calculated (capacity or susceptibility)
    k: number of times to resample
    """
    xs = []
    for i in range(k):
        resamples = np.random.choice(measurements,len(measurements))
        mean = np.mean(resamples)
        square_mean = np.mean(np.power(resamples,2))
        xs.append(value(mean,square_mean,N,T))
    x_squares = np.power(xs,2)
    av_x = np.mean(xs)
    av_x2 = np.mean(x_squares)
    return (av_x2-av_x**2)**(0.5)

def jacknife(measurements,value,N,T):
    """
    use the jacknife method to calculate the error on the capacity or the susceptibility
    measurement: measurement to resample (energy or temperature)
    value: value that is being calculated (capacity or susceptibility)
    """
    mean = np.mean(measurements)
    square_mean = np.mean(np.power(measurements,2))
    x = value(mean,square_mean,N,T)
    xs = []
    for i in range(measurements):
        new = np.delete(measurements,i)
        mean = np.mean(new)
        square_mean = np.mean(np.power(new))
        xs.append(value(mean,square_mean,N,T))
    return (sum([(x-xi)**2 for xi in xs]))**0.5
    

def main():
    Ts = np.linspace(1,3,20).tolist()
    lx,ly = 50,50
    experiment = {"params": {"N":lx*ly,"tau":10},"measurements":{}}
    Ms = []
    xi = []
    Es = []
    C = []
    i=0
    for T in Ts:
        L = Lattice(lx,ly,T,"up")
        L = equilibrate(L,"glauber")
        L.sim_glauber(10_000,True)
        E,M = L.get_measurements()
        measurement = {"T":T,"E":np.array(E),"M":np.absolute(np.array(M))}
        experiment["measurements"][i] = measurement
        del L
        i+=1
        print("completed temp",f"{i}/20")
    for index,measurement in experiment["measurements"].items():
        E = measurement["E"]
        E_mean = np.mean(E)
        E_square_mean = np.mean(np.power(E,2))
        M = measurement["M"]
        M_mean = np.mean(M)
        M_square_mean = np.mean(np.power(M,2))
        measurement["E_mean"] = float(E_mean)
        measurement["M_mean"] = float(M_mean)
        measurement["C"] = float(capacity(E_mean,E_square_mean,experiment["params"]["N"],measurement["T"]))
        measurement["xi"] = float(susceptibility(M_mean,M_square_mean,experiment["params"]["N"],measurement["T"]))
        measurement["E_error"] = error(E_mean,E_square_mean,experiment["params"]["N"])
        measurement["M_error"] = error(M_mean,M_square_mean,experiment["params"]["N"])
        measurement["C_berror"]=bootstrap(E,capacity,1000,experiment["params"]["N"],measurement["T"])
        measurement["xi_berror"]=bootstrap(M,susceptibility,1000,experiment["params"]["N"],measurement["T"])
        measurement["C_jerror"]=jacknife(E,capacity,1000,experiment["params"]["N"],measurement["T"])
        measurement["xi_jerror"]=jacknife(M,susceptibility,1000,experiment["params"]["N"],measurement["T"])
        del measurement["E"]
        del measurement["M"]

    with open("Data.json",'w') as outfile:
        json.dump(experiment,outfile)

start = time.time()
main()
print("Time:",time.time()-start)

