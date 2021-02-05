import json
from matplotlib import pyplot as plt
import sys
import numpy as np

def plot_glauber(error_method = "j"):
    #loading in the data from the json
    with open("Glauber_Data.json","r") as infile:
        data = json.load(infile)
    Ts =[]
    Ms = []
    M_error = []
    Es = []
    E_error = []
    chi = []
    chi_error = []
    C = []
    C_error = []
    #separating the data
    for x in data['measurements'].values():
        Ts.append(x['T'])
        Ms.append(x['M_mean'])
        M_error.append(x['M_error'])
        Es.append(x['E_mean'])
        E_error.append(x['E_error'])
        #extract different errors depending on the method specified
        if(error_method.lower() == "j"):
            filename = "Glauber temperature plot jacknife"
            chi_error.append(x['chi_jerror'])
            C_error.append(x['C_jerror'])
        else:
            filename = "Glauber temperature plot bootstrap"
            chi_error.append(x['chi_berror'])
            C_error.append(x['C_berror'])
        chi.append(x['chi'])
        C.append(x['C'])
    #find the critical temperature by looking at the maximum value for susceptibility and capacity
    chiTc = Ts[np.argmax(chi)]
    CTc = Ts[np.argmax(C)]
    fig,axs = plt.subplots(2,2,sharex=True)
    fig.suptitle("Glauber Dynamics Temperature progression plots")
    axs[0,0].errorbar(Ts,Ms,yerr=M_error,ecolor="r",capsize=1,barsabove=True)
    axs[0,0].set_title("Average Magnetism")
    axs[0,0].set_xlabel("T")
    axs[0,0].set_ylabel("|M|")
    axs[0,1].errorbar(Ts,chi,yerr=chi_error,ecolor="r",capsize=1,barsabove=True)
    axs[0,1].axvline(x=chiTc,color='g',label="Critical Temperature")
    axs[0,1].set_title("Susceptibility")
    axs[0,1].set_xlabel("T")
    axs[0,1].set_ylabel(r"$\chi$")
    axs[0,1].legend()
    axs[1,0].errorbar(Ts,Es,yerr=E_error,ecolor="r",capsize=1,barsabove=True)
    axs[1,0].set_title("Average Energy")
    axs[1,0].set_xlabel("T")
    axs[1,0].set_ylabel("Energy")
    axs[1,1].errorbar(Ts,C,yerr=C_error,ecolor="r",capsize=1,barsabove=True)
    axs[1,1].axvline(x=CTc,color='g',label="Critical Temperature")
    axs[1,1].set_title("Capacity")
    axs[1,1].set_xlabel("T")
    axs[1,1].set_ylabel("C")
    axs[1,1].legend()
    fig.set_size_inches(8,6)
    plt.savefig(filename+".png",dpi=100)
    plt.show()
    print(f"""Critical Temperature measurements:
susceptibility:{round(chiTc,2)} K
Capacity:{round(CTc,2)} K""")
    with open("Glauber_Data.json","w") as outfile:
        data['Critical Temperature'] = {"chi": chiTc,"C":CTc}
        json.dump(data,outfile)

def plot_kawasaki(error_method = "j"):
    with open("Kawasaki_Data.json","r") as infile:
        data = json.load(infile)

    Ts =[]
    Es = []
    E_error = []
    C = []
    C_error = []
    for x in data['measurements'].values():
        Ts.append(x['T'])
        Es.append(x['E_mean'])
        E_error.append(x['E_error'])
        if(error_method.lower() == "j"):
            filename = "Kawasaki temperature plot jacknife"
            C_error.append(x['C_jerror'])
        else:
            filename = "Kawasaki temperature plot bootstrap"
            C_error.append(x['C_berror'])
        C.append(x['C'])
    Tc = Ts[np.argmax(C)]
    fig,axs = plt.subplots(1,2,sharex=True)
    fig.suptitle("Kawasaki Dynamics Temperature progression plots")
    axs[0].errorbar(Ts,Es,yerr=E_error,ecolor="r",capsize=1,barsabove=True)
    axs[0].set_title("Average Energy")
    axs[0].set_xlabel("T")
    axs[0].set_ylabel("Energy")
    axs[1].errorbar(Ts,C,yerr=C_error,ecolor="r",capsize=1,barsabove=True)
    axs[1].axvline(x=Tc,color='g',label="Critical Temperature")
    axs[1].set_title("Capacity")
    axs[1].set_xlabel("T")
    axs[1].set_ylabel("C")
    axs[1].legend()
    fig.set_size_inches(8,6)
    plt.savefig(filename+".png",dpi=100)
    plt.show()
    print(f"Critical Temperature measurement from capacity: {round(Tc,2)} K")
    with open("Kawasaki_Data.json","w") as outfile:
        data['Critical Temperature'] = {"C":Tc}
        json.dump(data,outfile)

def main():
    method = sys.argv[1]
    error_method = sys.argv[2]

    if(method.lower() == "g"):
        plot_glauber(error_method)

    else:
        plot_kawasaki(error_method)

if __name__ == "__main__":
    main()