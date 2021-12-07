import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize

def my_actime(a):
    acf_time = 1.0
    for i in range(1, len(a)):
        if (a[i] > 0):
            acf_time += 2*a[i]
        else:
            break
    return i, acf_time

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = '20'
color_list = np.array(['black', 'crimson', 'dodgerblue'])
# color_list = np.array(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown'])

fig, (ax1, ax2) = plt.subplots(1, 2, sharex = False, figsize = (16, 8))

seed_list = [839301973, 511358380, 857095970, 566521973, 638570023, 924539071, 285032441, 343646014, 956132039, 550502395]

Nmol_ls = [20, 50, 100]
freq_ls = [1, 10, 100]

Nmol_ls = [50]
freq_ls = [10]

freq_ls = np.array(freq_ls)
for nm in range(len(Nmol_ls)):
    Nmol=Nmol_ls[nm]
    tau_nm = []
    for nf in range(len(freq_ls)):
        freq=freq_ls[nf]

        NB=Nmol*(1+1/2)
        icycle = 100000
        iout = 100
        Nframe = int(icycle/iout + 1)

        tau_ls = []
        
        for ns in range(0, len(seed_list)):
            f_log = f"/Users/weiimac/school/ber-multi/Nmol_{Nmol}-freq{freq}/seed_{seed_list[ns]}_bondlist.txt"             
            f = np.loadtxt(f_log, dtype=(int))
            mcstep = np.linspace(0, icycle, num=Nframe)
            
            nb_inf = 1.0/float(NB)
            bls = []
            binitial = f[0:Nmol, :]
            for i in range(0, Nframe):
                barray = f[i*Nmol:i*Nmol+Nmol, :]
                btb0 = np.multiply(binitial, barray)
                nb_t = np.sum(btb0)/np.sum(binitial)
                acf = (nb_t - nb_inf)/(1-nb_inf)
                bls.append(acf)
            bls = np.array(bls)
            ax1.plot(mcstep, bls)
        
            i, tau = my_actime(bls)
            print(f"Bond relaxation time = {tau*iout} MC steps")
            tau_ls.append(iout*tau)
           
        print(f"Avg bond relaxation time = {np.mean(np.array(tau_ls))}")
        
        tau_nm.append(np.mean(np.array(tau_ls)))
    
    ax2.scatter(1/freq_ls, tau_nm, marker='o', s=80, facecolors=color_list[nm], edgecolors=color_list[nm], alpha=1, label=f"Nmol={Nmol}") 

    
ax1.set_xlabel('MC steps',fontsize=24)
ax1.set_ylabel("Bond correlation function",fontsize=24)
ax2.set_xlabel("BER rate",fontsize=24)
ax2.set_ylabel("log($\u03C4_{}$)",fontsize=24)
ax1.legend(fontsize=20)
ax1.set_xscale('log')
fig.tight_layout()
