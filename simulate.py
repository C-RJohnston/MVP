from Ising_Model import Lattice
import copy
import numpy as np
import json
import sys
import concurrent.futures
import itertools
import time
def equilibrate(states: Lattice,method: str):
    """
    returns a lattice that has been allowed to reach equilibrium (100 sweeps) using either
    Glauber method or Kawasaki method
    """
    new = copy.copy(states)
    if method.lower() == "g" or method.lower() == "glauber":
        new.sim_glauber(100)
    elif method.lower() == "k" or method.lower() == "kawasaki":
        new.sim_kawasaki(200)
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
    for _ in range(k):
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
    for i in range(len(measurements)):
        new = np.delete(measurements,i)
        mean = np.mean(new)
        square_mean = np.mean(np.power(new,2))
        xs.append(value(mean,square_mean,N,T))
    return (sum([(x-chi)**2 for chi in xs]))**0.5
    
def do_glauber(L,runs,tau):
    print(f"Begin Glauber temp {L.T}")
    L = equilibrate(L,"glauber")
    L.sim_glauber(runs,True,tau)
    return L

def mp_glauber(lx,ly,T0,Tf,NT,runs,tau):
    Ts = np.linspace(T0,Tf,NT,False).tolist()
    Ls = [Lattice(lx,ly,T,"up") for T in Ts]
    experiment = {"params": {"N":lx*ly,"tau":tau},"measurements":{}}

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(do_glauber,Ls,itertools.repeat(runs,len(Ls)),itertools.repeat(tau,len(Ls)))
        i=0
        for L in results:
            E,M = L.get_measurements()
            E_mean = float(np.mean(E))
            E_square_mean = float(np.mean(np.power(E,2)))
            M_mean = float(np.mean(np.absolute(M)))
            M_square_mean = float(np.mean(np.power(M,2)))
            C = float(capacity(E_mean,E_square_mean,experiment["params"]["N"],L.T))
            chi = float(susceptibility(M_mean,M_square_mean,experiment["params"]["N"],L.T))
            E_error = error(E_mean,E_square_mean,experiment["params"]["N"])
            M_error = error(M_mean,M_square_mean,experiment["params"]["N"])
            C_berror = bootstrap(E,capacity,1000,experiment["params"]["N"],L.T)
            chi_berror = bootstrap(M,susceptibility,1000,experiment["params"]["N"],L.T)
            C_jerror = jacknife(E,capacity,experiment["params"]["N"],L.T)
            chi_jerror = jacknife(M,susceptibility,experiment["params"]["N"],L.T)
            x = {"T": L.T,
             "E_mean": E_mean,
             "M_mean": M_mean,
             "C":C,
             "chi":chi,
             "E_error":E_error,
             "M_error":M_error,
             "C_berror":C_berror,
             "chi_berror":chi_berror,
             "C_jerror":C_jerror,
             "chi_jerror":chi_jerror}
            experiment['measurements'][i] = x
            i+=1
    with open("Glauber_Data.json",'w') as outfile:
        json.dump(experiment,outfile)

    

def glauber(lx,ly,T0,Tf,NT,runs,tau):
    Ts = np.linspace(T0,Tf,NT,False).tolist()
    experiment = {"params": {"N":lx*ly,"tau":tau},"measurements":{}}
    i=0
    for T in Ts:
        L = Lattice(lx,ly,T,"up")
        L = equilibrate(L,"glauber")
        L.sim_glauber(runs,True,tau)
        E,M = L.get_measurements()
        measurement = {"T":T,"E":np.array(E),"M":np.absolute(np.array(M))}
        experiment["measurements"][i] = measurement
        del L
        i+=1
        print(f"Glauber: {i}/{NT}")
    for measurement in experiment["measurements"].values():
        E = measurement["E"]
        E_mean = np.mean(E)
        E_square_mean = np.mean(np.power(E,2))
        M = measurement["M"]
        M_mean = np.mean(M)
        M_square_mean = np.mean(np.power(M,2))
        measurement["E_mean"] = float(E_mean)
        measurement["M_mean"] = float(M_mean)
        measurement["C"] = float(capacity(E_mean,E_square_mean,experiment["params"]["N"],measurement["T"]))
        measurement["chi"] = float(susceptibility(M_mean,M_square_mean,experiment["params"]["N"],measurement["T"]))
        measurement["E_error"] = error(E_mean,E_square_mean,experiment["params"]["N"])
        measurement["M_error"] = error(M_mean,M_square_mean,experiment["params"]["N"])
        measurement["C_berror"]=bootstrap(E,capacity,1000,experiment["params"]["N"],measurement["T"])
        measurement["chi_berror"]=bootstrap(M,susceptibility,1000,experiment["params"]["N"],measurement["T"])
        measurement["C_jerror"]=jacknife(E,capacity,experiment["params"]["N"],measurement["T"])
        measurement["chi_jerror"]=jacknife(M,susceptibility,experiment["params"]["N"],measurement["T"])
        del measurement["E"]
        del measurement["M"]

    with open("Glauber_Data.json",'w') as outfile:
        json.dump(experiment,outfile)

def do_kawasaki(L,runs,tau):
    print(f"Begin Kawasaki temp {L.T}")
    L = equilibrate(L,"kawasaki")
    L.sim_kawasaki(runs,True,tau)
    return L

def mp_kawasaki(lx,ly,T0,Tf,NT,runs,tau):
    Ts = np.linspace(T0,Tf,NT,False).tolist()
    Ls = [Lattice(lx,ly,T) for T in Ts]
    experiment = {"params": {"N":lx*ly,"tau":tau},"measurements":{}}

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(do_kawasaki,Ls,itertools.repeat(runs,len(Ls)),itertools.repeat(tau,len(Ls)))
        i=0
        for L in results:
            E = L.get_measurements()[0]
            E_mean = float(np.mean(E))
            E_square_mean = float(np.mean(np.power(E,2)))
            C = float(capacity(E_mean,E_square_mean,experiment["params"]["N"],L.T))
            E_error = error(E_mean,E_square_mean,experiment["params"]["N"])
            C_berror = bootstrap(E,capacity,1000,experiment["params"]["N"],L.T)
            C_jerror = jacknife(E,capacity,experiment["params"]["N"],L.T)
            x = {"T": L.T,
             "E_mean": E_mean,
             "C":C,
             "E_error":E_error,
             "C_berror":C_berror,
             "C_jerror":C_jerror,}
            experiment['measurements'][i] = x
            i+=1
    with open("Kawasaki_Data.json",'w') as outfile:
        json.dump(experiment,outfile)

def kawasaki(lx,ly,T0,Tf,NT,runs,tau):
    Ts = np.linspace(T0,Tf,NT,False).tolist()
    experiment = {"params": {"N":lx*ly,"tau":tau},"measurements":{}}
    i=0
    for T in Ts:
        L = Lattice(lx,ly,T,"random")
        L = equilibrate(L,"kawasaki")
        L.sim_kawasaki(runs,True,tau)
        E= L.get_measurements()[0]
        measurement = {"T":T,"E":np.array(E)}
        experiment["measurements"][i] = measurement
        del L
        i+=1
        print(f"Kawasaki: {i}/{NT}")
    for measurement in experiment["measurements"].values():
        E = measurement["E"]
        E_mean = np.mean(E)
        E_square_mean = np.mean(np.power(E,2))
        measurement["E_mean"] = float(E_mean)
        measurement["C"] = float(capacity(E_mean,E_square_mean,experiment["params"]["N"],measurement["T"]))
        measurement["E_error"] = error(E_mean,E_square_mean,experiment["params"]["N"])
        measurement["C_berror"]=bootstrap(E,capacity,1000,experiment["params"]["N"],measurement["T"])
        measurement["C_jerror"]=jacknife(E,capacity,experiment["params"]["N"],measurement["T"])
        del measurement["E"]
    with open("Kawasaki_Data.json",'w') as outfile:
        json.dump(experiment,outfile)
    

def main():
    try: 
        if(sys.argv[8] == "-nomulti"):
            t = time.perf_counter()
            glauber(*[int(x) for x in sys.argv[1:8]])
            kawasaki(*[int(x) for x in sys.argv[1:8]])
            print(f"time to complete: {(time.perf_counter()-t)/60} minutes (which is {(time.perf_counter()-t)/3600} hours)")
    except:
        t = time.perf_counter()
        mp_glauber(*[int(x) for x in sys.argv[1:8]])
        mp_kawasaki(*[int(x) for x in sys.argv[1:8]])
        print(f"time to complete: {(time.perf_counter()-t)/60} minutes (which is {(time.perf_counter()-t)/3600} hours)")

if __name__ == "__main__":
    main()
    