#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D



cmap = mpl.colors.LinearSegmentedColormap.from_list("tasetful rbg", colors=[u'#0072bb',  u'#d62728', u'#2ca02c'])


mpl.rcParams["font.size"] = 16
def plot(x, y, spins, ns, name, color=True, bar=True, sz=False):
    fig, ax = plt.subplots()
    for i in range(ns):
        points = np.concatenate([x.reshape(-1,1), y[:,i].reshape(-1,1)], axis=1).reshape(-1,1,2)
        segs = np.concatenate([points[:-1], points[1:]],axis=1)

        if(color):
            if(sz):
                if(i <= ns/3):
                    lcol = u'#8c564b'
                elif(i < 2*ns/3 ):
                    lcol = u'#17becf'
                else:
                    lcol = u'#ff7f0e'
                lc = LineCollection(segs, color = lcol)
            else:
                lc = LineCollection(segs, array=spins[:-1,i], cmap=cmap)
        else:
            lc = LineCollection(segs, color=u'#9467bd')
        lc.set_capstyle('round')
        lc.set_linewidth(3)
        ax.add_collection(lc)
        lc.set_clim(minSpin, maxSpin)

    ax.set_xlabel("r ($\AA$)")
    ax.set_ylabel("E (a.u.)")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)
    if(color and bar):
        fig.colorbar(lc).set_label(r'$\langle S^2 \rangle$')
    plt.tight_layout()
    plt.savefig("plots/" + name, format="eps", dpi=1000)
    plt.clf()
    plt.close()








# ----------------- HYDROGEN STO-3G --------------------
ymin = -1.2
ymax = -0.2
xmin = 0.5
xmax = 3.5
nStates = 3
CIS_Edata = np.genfromtxt(open("H2_min/H2_CISEnergyPlot.csv", "r"), delimiter=",")
CIS_Sdata = np.genfromtxt(open("H2_min/H2_s2Plot.csv", "r"), delimiter=",")
spinsCIS = CIS_Sdata[:,1:]
maxSpin = np.max(spinsCIS)
minSpin = np.min(spinsCIS)

# Plot Normal CIS
x = CIS_Edata[:,0]
EnergiesCIS = CIS_Edata[:,1:]
plot(x, EnergiesCIS, spinsCIS, nStates, "H2_min.eps")





# plot TDHF
TDHF_Edata = np.genfromtxt(open("H2_min/H2_TDHFEnergyPlot.csv", "r"), delimiter=",")
maxSpin = np.max(spinsCIS)
minSpin = np.min(spinsCIS)

# Plot Normal CIS
x = TDHF_Edata[:,0]
EnergiesTDHF = TDHF_Edata[:,1:]
plot(x, EnergiesTDHF, spinsCIS, nStates, "H2_min_TDHF.eps")







# Plot sf-t-CIS
sf_t_CIS_Edata = np.genfromtxt(open("H2_min/H2_sf_t_CISEnergyPlot.csv", "r"), delimiter=",")
sf_t_Sdata = np.genfromtxt(open("H2_min/H2_sf_t_s2Plot.csv", "r"), delimiter=",")

x = sf_t_CIS_Edata[:,0]
Energies_sftCIS = sf_t_CIS_Edata[:,2:]
spins = sf_t_Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_sftCIS, spins, nStates, "H2_min_sft.eps")





# Plot sf-s-CIS
nStates = 5
sf_s_CIS_Edata = np.genfromtxt(open("H2_min/H2_sf_CISEnergyPlot.csv", "r"), delimiter=",")
sf_s_Sdata = np.genfromtxt(open("H2_min/H2_sf_s2Plot.csv", "r"), delimiter=",")

x = sf_s_CIS_Edata[:,0]
Energies_sfsCIS = np.hstack((sf_s_CIS_Edata[:,1:], EnergiesCIS))
spins = np.hstack((sf_s_Sdata[:,1:], spinsCIS))

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_sfsCIS, spins, nStates, "H2_min_sfs.eps")




# Plot exact
nStates = 4
exact_CIS_Edata = np.genfromtxt(open("H2_min/H2_exact_CISEnergyPlot.csv", "r"), delimiter=",")
exact_Sdata = np.genfromtxt(open("H2_min/H2_exact_s2Plot.csv", "r"), delimiter=",")

x = exact_CIS_Edata[:,0]
Energies_exactCIS = exact_CIS_Edata[:,1:]
spins = exact_Sdata[:,1:]


maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_exactCIS, spins, nStates, "H2_min_exact.eps", bar=False)





















# ----------------- HYDROGEN TZ --------------------
# Triple Zeta
nSinglets = 6
nTriplets = 5

ymax=-0.5
xmax = 5

rCF = 100
cutoff = 0.95


CIS_Edata = np.genfromtxt(open("H2/augtz/H2_CISEnergyPlot.csv", "r"), delimiter=",")
CIS_Sdata = np.genfromtxt(open("H2/augtz/H2_s2Plot.csv", "r"), delimiter=",")


x = CIS_Edata[:,0]
Energies_CIS = CIS_Edata[:,3:]
spins = CIS_Sdata[:,3:]
plot(x, Energies_CIS, spins, 10, "test.eps")


# separate into singlets and triplets up to CF point
spins_1 = spins[-rCF:,:]
energies_1 = Energies_CIS[-rCF:,:]

energies_s_1 = np.zeros((energies_1.shape[0],nSinglets))
energies_t_1 = np.zeros((energies_1.shape[0],nTriplets))
spins_s_1 = np.zeros((energies_1.shape[0],nSinglets))
spins_t_1 = np.zeros((energies_1.shape[0],nTriplets))

for i in range(energies_1.shape[0]):
    energies_s_1[i,:] = ((energies_1[i,:])[spins_1[i,:] <= cutoff])[:nSinglets]
    energies_t_1[i,:] = ((energies_1[i,:])[spins_1[i,:] >= 2 - cutoff])[:nTriplets]
    spins_s_1[i,:] = ((spins_1[i,:])[spins_1[i,:] <= cutoff])[:nSinglets]
    spins_t_1[i,:] = ((spins_1[i,:])[spins_1[i,:] >= 2-cutoff])[:nTriplets]


# follow the states past CF
s_inds = spins_1[0,:] <= cutoff
t_inds = spins_1[0,:] >= 2 - cutoff
s_inds[np.cumsum(s_inds) > nSinglets] = False
t_inds[np.cumsum(t_inds) > nTriplets] = False
energies_s_2 = Energies_CIS[:-rCF, s_inds]
spins_s_2 = spins[:-rCF,s_inds]
energies_t_2 = Energies_CIS[:-rCF, t_inds]
spins_t_2 = spins[:-rCF, t_inds]

energies_s = np.concatenate((energies_s_2, energies_s_1), axis=0)
spins_s = np.concatenate((spins_s_2, spins_s_1), axis=0)
energies_t = np.concatenate((energies_t_2, energies_t_1), axis=0)
spins_t = np.concatenate((spins_t_2, spins_t_1), axis=0)

# add in lowest triplet
energies_t = np.concatenate((CIS_Edata[:,2].reshape((-1,1)), energies_t), axis=1)
spins_t = np.concatenate((CIS_Sdata[:,2].reshape((-1,1)), spins_t), axis=1)

maxSpin = np.max(spins)
minSpin = np.min(spins)

# plot singlets
ymax = -0.45
ymin = -0.8
plot(x, energies_s, spins_s, nSinglets, "H2_tz_s.eps")

# plot triplets
ymin = -1.05
plot(x, energies_t, spins_t, nTriplets, "H2_tz_t.eps")




# TDHF
ymin = -1.15
ymax = -0.55
nStates = 6
TDHF_Edata = np.genfromtxt(open("H2/augtz/H2_TDHFEnergyPlot.csv", "r"), delimiter=",")

maxSpin = np.max(spins)
minSpin = np.min(spins)

x = CIS_Edata[:,0]
EnergiesTDHF = TDHF_Edata[:,1:]
spins = CIS_Sdata[:,1:]
plot(x, EnergiesTDHF, spins, nStates, "H2_tz_tdhf.eps")







# Triple Zeta SF CIS
nStates = 8

ymax = -0.45
ymin = -1.18

sf_t_CIS_Edata = np.genfromtxt(open("H2/augtz/H2_sf_t_CISEnergyPlot.csv", "r"), delimiter=",")
sf_t_Sdata = np.genfromtxt(open("H2/augtz/H2_sf_t_s2Plot.csv", "r"), delimiter=",")


x = sf_t_CIS_Edata[:,0]
Energies_sftCIS = sf_t_CIS_Edata[:,2:]
spins = sf_t_Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_sftCIS, spins, nStates, "H2_tz_sf.eps")





# Triple Zeta SFS CIS
nStates = 5

ymax = -0.45
ymin = -1.05

maxSpin = 2 #np.max(spins)
minSpin = 0 #np.min(spins)

sf_s_CIS_Edata = np.genfromtxt(open("H2/augtz/H2_sf_CISEnergyPlot.csv", "r"), delimiter=",")
sf_s_Sdata = np.genfromtxt(open("H2/augtz/H2_sf_s2Plot.csv", "r"), delimiter=",")


x = sf_s_CIS_Edata[:,0]
Energies_sfsCIS = sf_s_CIS_Edata[:,2:]
spins = sf_s_Sdata[:,2:]


plot(x, Energies_sfsCIS, spins, nStates, "H2_tz_sfs.eps", bar=False)







# Triple Zeta XCIS
nStates = 8
x_CIS_Edata = np.genfromtxt(open("H2/augtz/H2_x_t_CISEnergyPlot.csv", "r"), delimiter=",")
x_Sdata = np.genfromtxt(open("H2/augtz/H2_x_t_s2Plot.csv", "r"), delimiter=",")


x = x_CIS_Edata[:,0]
Energies_xCIS = x_CIS_Edata[:,2:]
spins = x_Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

fig, ax = plt.subplots()
for i in range(nStates):
    points = np.concatenate([x.reshape(-1,1), Energies_xCIS[:,i].reshape(-1,1)], axis=1).reshape(-1,1,2)
    segs = np.concatenate([points[:-1], points[1:]],axis=1)
    colors = [(1-s,0,s) for s in spins[:-1,i]]

    lc = LineCollection(segs, array=spins[:-1,i], cmap=cmap)
    lc.set_linewidth(3)
    ax.add_collection(lc)
    lc.set_clim(minSpin, maxSpin)


ax.set_xlabel("r ($\AA$)")
ax.set_ylabel("E (a.u.)")
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
fig.colorbar(lc).set_label(r'$\langle S^2 \rangle$')

plt.tight_layout()
plt.savefig("plots/H2_tz_x.eps", format="eps", dpi=1000)
plt.clf()






# Triple Zeta RCIS
nStates = 6

maxSpin = 2 #np.max(spins)
minSpin = 0 #np.min(spins)

RCIS_Edata = np.genfromtxt(open("H2/augtz/rH2_CISEnergyPlot.csv", "r"))
r_Sdata = np.genfromtxt(open("H2/augtz/rH2_s2Plot.csv", "r"), delimiter=",")

x = RCIS_Edata[:,0]
Energies_RCIS = RCIS_Edata[:,2:]
spins = r_Sdata[:,2:]

t_inds = spins[200, :] == 2.0
s_inds = spins[200,:] == 0


energies_t = []
energies_s = []
spins_t = []
spins_s = []
for dst in range(Energies_RCIS.shape[0]):
    energies_t.append(Energies_RCIS[dst, spins[dst,:] == 2.0][:nStates])
    energies_s.append(Energies_RCIS[dst, spins[dst,:] == 0][:nStates])
    spins_t.append(spins[dst, spins[dst,:] == 2.0][:nStates])
    spins_s.append(spins[dst, spins[dst,:] == 0][:nStates])

Energies_RCIS_t = np.array(energies_t)
Energies_RCIS_s = np.array(energies_s)
spins_t = np.array(spins_t)
spins_s = np.array(spins_s)


ymin = -1.05
ymax = -0.45
plot(x, Energies_RCIS_t, spins_t, nStates, "H2_tz_rcis_t.eps", bar=False)

ymin = -0.8
ymax = -0.45
plot(x, Energies_RCIS_s, spins_s, nStates, "H2_tz_rcis_s.eps", bar=False)











# Triple Zeta Exact
nStates = 5
nSinglets = 6
nTriplets = 6
cutoff_state = 20
edata = np.genfromtxt(open("H2_exact/aug-cc-pVTZ_larger_energyPlot.csv", "r"), delimiter=",")
sdata = np.genfromtxt(open("H2_exact/aug-cc-pVTZ_larger_s2Plot.csv", "r"), delimiter=",")


x_ex = edata[:,0]
energies = edata[:,2:cutoff_state]
spins = sdata[:,2:cutoff_state]

energies_s_ex = np.zeros((energies.shape[0],nSinglets))
energies_t_ex = np.zeros((energies.shape[0],nTriplets))
spins_s_ex = np.zeros((energies.shape[0],nSinglets))
spins_t_ex = np.zeros((energies.shape[0],nTriplets))

for i in range(energies.shape[0]):
    energies_s_ex[i,:] = ((energies[i,:])[spins[i,:] <= cutoff])[:nSinglets]
    energies_t_ex[i,:] = ((energies[i,:])[spins[i,:] >= 2 - cutoff])[:nTriplets]
    spins_s_ex[i,:] = ((spins[i,:])[spins[i,:] <= cutoff])[:nSinglets]
    spins_t_ex[i,:] = ((spins[i,:])[spins[i,:] >= 2-cutoff])[:nTriplets]

maxSpin = np.max(spins)
minSpin = np.min(spins)



# plot singlets
ymax = -0.45
ymin = -0.8
plot(x_ex, energies_s_ex, spins_s_ex, nSinglets, "H2_tz_exact_s.eps", bar=False)



# plot triplets
ymin = -1.05
plot(x_ex, energies_t_ex, spins_t_ex, nTriplets, "H2_tz_exact_t.eps", bar=False)




























# ----------------- H2 TDDFT --------------------
# PBE
# augtz
nStates = 4
ymin = -0.99
ymax = -0.6

TDA_Edata = np.genfromtxt(open("H2/augtz/H2_pbe_CISEnergyPlot.csv", "r"), delimiter=",")
TDDFT_Edata = np.genfromtxt(open("H2/augtz/H2_pbe_TDHFEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("H2/augtz/H2_pbe_s2Plot.csv", "r"), delimiter=",")


x = x_CIS_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
Energies_TDDFT = TDDFT_Edata[:,2:]
spins = Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_TDA, spins, nStates, "H2_tz_pbe.eps")
plot(x, Energies_TDDFT, spins, nStates, "H2_tz_pbe_rpa.eps")



# sto-3g
nStates = 3
ymin = -0.99
ymax = -0.6

TDA_Edata = np.genfromtxt(open("H2_min/H2_pbe_CISEnergyPlot.csv", "r"), delimiter=",")
TDDFT_Edata = np.genfromtxt(open("H2_min/H2_pbe_TDHFEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("H2_min/H2_pbe_s2Plot.csv", "r"), delimiter=",")


x = x_CIS_Edata[:,0]
Energies_TDA = TDA_Edata[:,1:]
Energies_TDDFT = TDDFT_Edata[:,1:]
spins = Sdata[:,1:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_TDA, spins, nStates, "H2_min_pbe.eps")
plot(x, Energies_TDDFT, spins, nStates, "H2_min_pbe_rpa.eps")


















# sfs
nStates = 4
ymin = -0.99
ymax = -0.6

TDA_Edata = np.genfromtxt(open("H2/augtz/H2_pbe_sfs_CISEnergyPlot.csv", "r"))
Sdata = np.genfromtxt(open("H2/augtz/H2_pbe_sfs_s2Plot.csv", "r"), delimiter=",")

x = TDA_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
spins = Sdata[:,2:]

plot(x, Energies_TDA, spins, nStates, "H2_tz_pbe_sfs.eps")










# PBE0
# TDA
TDA_Edata = np.genfromtxt(open("H2/augtz/H2_pbe0_CISEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("H2/augtz/H2_pbe0_s2Plot.csv", "r"), delimiter=",")

x = x_CIS_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
spins = Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_TDA, spins, nStates, "H2_tz_pbe0.eps")







# sfs
nStates = 4
ymin = -0.99
ymax = -0.6

TDA_Edata = np.genfromtxt(open("H2/augtz/H2_pbe0_sfs_CISEnergyPlot.csv", "r"))
Sdata = np.genfromtxt(open("H2/augtz/H2_pbe0_sfs_s2Plot.csv", "r"), delimiter=",")

x = TDA_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
spins = Sdata[:,2:]

plot(x, Energies_TDA, spins, nStates, "H2_tz_pbe0_sfs.eps")











# LRC-wPBEh
TDA_Edata = np.genfromtxt(open("H2/augtz/H2_wpbe_CISEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("H2/augtz/H2_wpbe_s2Plot.csv", "r"), delimiter=",")


x = x_CIS_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
spins = Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_TDA, spins, nStates, "H2_tz_wpbe.eps")



# sfs
nStates = 4
ymin = -0.99
ymax = -0.6

TDA_Edata = np.genfromtxt(open("H2/augtz/H2_wpbe_sfs_CISEnergyPlot.csv", "r"))
Sdata = np.genfromtxt(open("H2/augtz/H2_wpbe_sfs_s2Plot.csv", "r"), delimiter=",")

x = TDA_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
spins = Sdata[:,2:]

plot(x, Energies_TDA, spins, nStates, "H2_tz_wpbe_sfs.eps")






















# ----------------- ETHANE --------------------


ymin = -79.3
ymax = -78.6
xmin = 1
xmax = 5




# Plot Normal CIS
nStates = 7
CIS_Edata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_CISEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_s2Plot.csv", "r"), delimiter=",")

x = CIS_Edata[:,0]
EnergiesCIS = CIS_Edata[:,1:]
spins = Sdata[:,1:]

plot(x, EnergiesCIS, spins, nStates, "ethane_augdz.eps")




# Plot t CIS
ymin = -79.2
ymax = -78.7
xmin = 1
xmax = 5

nStates = 7
CIS_Edata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_t_CISEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_t_s2Plot.csv", "r"), delimiter=",")

x = CIS_Edata[:,0]
EnergiesCIS = CIS_Edata[:,1:]
spins = Sdata[:,1:]

plot(x, EnergiesCIS, spins, nStates, "ethane_augdz_t.eps", bar=False)






# Plot sf CIS
ymin = -79.3
ymax = -78.6
xmin = 1
xmax = 5

sf_Edata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_sf_t_CISEnergyPlot.csv", "r"), delimiter=",")
sf_Sdata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_sf_t_s2Plot.csv", "r"), delimiter=",")

x = sf_Edata[:,0]
sf_Energies = sf_Edata[:,2:]
spins = sf_Sdata[:,2:]

plot(x, sf_Energies, spins, nStates, "ethane_augdz_sf.eps")




# Plot sfs CIS
sf_Edata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_sf_CISEnergyPlot.csv", "r"), delimiter=",")
sf_Sdata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_sf_s2Plot.csv", "r"), delimiter=",")

x = sf_Edata[:,0]
sf_Energies = sf_Edata[:,2:]
spins = sf_Sdata[:,2:]

plot(x, sf_Energies, spins, nStates, "ethane_augdz_sfs.eps")








# Plot sf CIS with MOM
mom_Edata = np.genfromtxt(open("jobs/CH3CH3_mom_lesstight_CISEnergyPlot.csv", "r"), delimiter=",")
mom_Sdata = np.genfromtxt(open("jobs/CH3CH3_mom_lesstight_s2Plot.csv", "r"), delimiter=",")

x = mom_Edata[:,0]
mom_Energies = mom_Edata[:,2:]
spins = mom_Sdata[:,2:]

plot(x, mom_Energies, spins, nStates, "ethane_augdz_mom.eps")






# Plot x CIS
nSinglets = 5
nTriplets = 5
x_Edata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_x_t_CISEnergyPlot.csv", "r"))
x_Sdata = np.genfromtxt(open("CH3CH3/augdz/CH3CH3_x_t_s2Plot.csv", "r"), delimiter=',')

x = x_Edata[:,0]
x_Energies = x_Edata[:,2:]
spins = x_Sdata[:,2:]

x_Energies_s = np.zeros((x_Energies.shape[0], nSinglets))
x_Energies_t = np.zeros((x_Energies.shape[0], nTriplets))
x_spins_s = np.zeros((x_Energies.shape[0], nSinglets))
x_spins_t = np.zeros((x_Energies.shape[0], nTriplets))

# get lowest singlet
for d in range(x_Energies.shape[0]):
    x_Energies_s[d,:] = x_Energies[d, spins[d,:] <= cutoff][:nSinglets]
    x_Energies_t[d,:] = x_Energies[d, spins[d,:] >= 2 - cutoff][:nTriplets]
    x_spins_s[d,:] = spins[d, spins[d,:] <= cutoff][:nSinglets]
    x_spins_t[d,:] = spins[d, spins[d,:] >= 2 - cutoff][:nTriplets]

x_Energies = np.concatenate((x_Energies_s, x_Energies_t), axis=1)
spins = np.concatenate((x_spins_s, x_spins_t), axis=1)

plot(x, x_Energies, spins, nStates, "ethane_augdz_x.eps", bar=False)
















# ----------------------- SILANE ------------------------------
xmin = 1.7
xmax = 5.0
ymin = -581.4
ymax = -580.8

xmin_in = 2.3
xmax_in = 2.6
ymin_in = -581.337
ymax_in = -581.33
ymin_in2 = -581.0825
ymax_in2 = -581.0725
in_loc = 0.75
in_size = 0.22

nStates = 3

sf_t_CIS_Edata = np.genfromtxt(open("silane/daugdz/silane_sf_t_CISEnergyPlot.csv", "r"), delimiter=",")
sf_t_Sdata = np.genfromtxt(open("silane/daugdz/silane_sf_t_s2Plot.csv", "r"), delimiter=",")

# Plot SF-CIS
x = sf_t_CIS_Edata[:,0]






spins = sf_t_Sdata[:,2:]
Energies_sf_CIS = sf_t_CIS_Edata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

# axes and inset axes
fig, ax = plt.subplots()
ax_in = ax.inset_axes([in_loc, 1 - in_loc - in_size, in_size, in_size])
ax_in2 = ax.inset_axes([in_loc, in_loc, in_size, in_size])

for i in range(nStates):
    points = np.concatenate([x.reshape(-1,1), Energies_sf_CIS[:,i].reshape(-1,1)], axis=1).reshape(-1,1,2)
    segs = np.concatenate([points[:-1], points[1:]],axis=1)

    lc = LineCollection(segs, array=spins[:-1,i], cmap=cmap)
    lc.set_linewidth(3)
    lc.set_clim(minSpin, maxSpin)
    lc.set_capstyle('round')
    ax.add_collection(lc)

    #lc2 = LineCollection(segs, array=spins[:-1,i], cmap=cmap)
    #lc2.set_linewidth(3)
    #lc2.set_clim(minSpin, maxSpin)
    #lc2.set_capstyle('round')
    #ax_in.add_collection(lc2)
    ax_in.plot(points[:,0,0], points[:,0,1], marker='o', markersize=2, color=cmap(spins[0,i]/2.0))
    ax_in2.plot(points[:,0,0], points[:,0,1], marker='o', markersize=2,  color=cmap(spins[0,i]/2.0))


ax.set_xlabel("r ($\AA$)")
ax.set_ylabel("E (a.u.)")
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
plt.colorbar(lc).set_label(r'$\langle S^2 \rangle$')

ax_in.set_xlim(xmin_in, xmax_in)
ax_in.set_ylim(ymin_in, ymax_in)
ax_in.set_xticklabels('')
ax_in.set_yticklabels('')
ax.indicate_inset_zoom(ax_in)

ax_in2.set_xlim(xmin_in, xmax_in)
ax_in2.set_ylim(ymin_in2, ymax_in2)
ax_in2.set_xticklabels('')
ax_in2.set_yticklabels('')
ax.indicate_inset_zoom(ax_in2)

plt.tight_layout()
plt.savefig("plots/silane_sf.eps", format="eps", dpi=1000)
plt.clf()






# plot MOM
xmin_in = 2.3
xmax_in = 2.7
ymin_in = -581.333
ymax_in = -581.31
ymin_in2 = -581.0825
ymax_in2 = -581.0725
in_loc = 0.75
in_size = 0.22

mom_Edata = np.genfromtxt(open("silane/daugdz_mom/sf-t-MOM-CIS_CISEnergyPlot.csv", "r"), delimiter=",")
mom_Sdata = np.genfromtxt(open("silane/daugdz_mom/sf-t-MOM-CIS_s2Plot.csv", "r"), delimiter=",")

x2 = mom_Edata[:,0]
spins = mom_Sdata[:,2:]
mom_CIS = mom_Edata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

# axes and inset axes
fig, ax = plt.subplots()
ax_in = ax.inset_axes([in_loc, 1 - in_loc - in_size, in_size, in_size])
ax_in2 = ax.inset_axes([in_loc, in_loc, in_size, in_size])


for i in range(nStates):
    points = np.concatenate([x2.reshape(-1,1), mom_CIS[:,i].reshape(-1,1)], axis=1).reshape(-1,1,2)
    segs = np.concatenate([points[:-1], points[1:]],axis=1)
    colors = [(1-s,0,s) for s in spins[:-1,i]]

    lc = LineCollection(segs, array=spins[:-1,i], cmap=cmap)
    lc.set_linewidth(3)
    ax.add_collection(lc)
    lc.set_capstyle('round')
    lc.set_clim(minSpin, maxSpin)

    ax_in.plot(points[:,0,0], points[:,0,1], marker='o', markersize=2, color=cmap(spins[0,i]/2.0))
    ax_in2.plot(points[:,0,0], points[:,0,1], marker='o', markersize=2,  color=cmap(spins[0,i]/2.0))


ax.set_xlabel("r ($\AA$)")
ax.set_ylabel("E (a.u.)")
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
plt.colorbar(lc).set_label(r'$\langle S^2 \rangle$')

ax_in.set_xlim(xmin_in, xmax_in)
ax_in.set_ylim(ymin_in, ymax_in)
ax_in.set_xticklabels('')
ax_in.set_yticklabels('')
ax.indicate_inset_zoom(ax_in)

ax_in2.set_xlim(xmin_in, xmax_in)
ax_in2.set_ylim(ymin_in2, ymax_in2)
ax_in2.set_xticklabels('')
ax_in2.set_yticklabels('')
ax.indicate_inset_zoom(ax_in2)

plt.tight_layout()
plt.savefig("plots/silane_sf_mom.eps", format="eps", dpi=1000)
plt.clf()















# ----------------------- LiH ------------------------------
xmin = 1.2
xmax = 5
ymin = -8
ymax = -7.75

# plot CIS
nStates = 8

lih_edata = np.genfromtxt(open("LiH/augtz/LiH_CISEnergyPlot.csv", "r"), delimiter=",")
lih_sdata = np.genfromtxt(open("LiH/augtz/LiH_s2Plot.csv", "r"), delimiter=",")


x = lih_edata[:,0]
spins = lih_sdata[:,1:]
lih_CIS = lih_edata[:,1:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, lih_CIS, spins, nStates, "LiH_tz.eps")






# plot sfs-CIS
xmin = 1.2
xmax = 5
ymin = -7.95
ymax = -7.45

nStates = 8

lih_spin1_edata = np.genfromtxt(open("LiH/spin/LiH_spin1_CISEnergyPlot.csv", "r"))
lih_spin2_edata = np.genfromtxt(open("LiH/spin/LiH_spin2_CISEnergyPlot.csv", "r"))


x = lih_edata[:,0]
lih_CIS = np.concatenate((lih_spin1_edata[:,2:nStates+2], lih_spin2_edata[:,2:nStates+2]), axis=1)

ns = lih_CIS.shape[1]
fig, ax = plt.subplots()
for i in range(ns):
    if(i < ns/2 ):
        lcol = u'#17becf'
    else:
        lcol = u'#ff7f0e'
    ax.plot(x, lih_CIS[:,i], color=lcol, linewidth=3)


ax.set_xlabel("r ($\AA$)")
ax.set_ylabel("E (a.u.)")
ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)

custom_lines = [Line2D([0], [0], color=u'#17becf', lw=3),
                Line2D([0], [0], color=u'#ff7f0e', lw=3)]
plt.legend(custom_lines, ['$M_s=-1$','$M_s=+1$'])
plt.tight_layout()
plt.savefig("plots/LiH_tz_sfs.eps", format="eps", dpi=1000)
plt.clf()
plt.close()














# ----------------------- NH3 ------------------------------
xmin = 0.5
xmax = 5
ymin = -56.25
ymax = -55.7

# plot CIS
nStates = 5

nh3_edata = np.genfromtxt(open("NH3/NH3_CISEnergyPlot.csv", "r"), delimiter=",")
nh3_sdata = np.genfromtxt(open("NH3/NH3_s2Plot.csv", "r"), delimiter=",")


x = nh3_edata[:,0]
spins = nh3_sdata[:,1:nStates+1]
nh3_CIS = nh3_edata[:,1:nStates+1]

maxSpin = 2 #np.max(spins)
minSpin = 0 #np.min(spins)

plot(x, nh3_CIS, spins, nStates, "nh3_tz.eps")








# plot SFS-CIS
ymax = -55.9
nStates = 2

nh3_s1_edata = np.genfromtxt(open("NH3/spin1_CISEnergyPlot.csv", "r"))
nh3_s1_sdata = np.genfromtxt(open("NH3/spin1_s2Plot.csv", "r"), delimiter=",")
nh3_s2_edata = np.genfromtxt(open("NH3/spin2_CISEnergyPlot.csv", "r"))
nh3_s2_sdata = np.genfromtxt(open("NH3/spin2_s2Plot.csv", "r"), delimiter=",")


x = nh3_s1_edata[:,0]
spins = np.concatenate((nh3_sdata[:,1:nStates+1], nh3_s1_sdata[:,2:nStates+1], nh3_s2_sdata[:,2:nStates+1]), axis=1)
nh3_CIS = np.concatenate((nh3_edata[:,1:nStates+1], nh3_s1_edata[:,2:nStates+1], nh3_s2_edata[:,2:nStates+1]), axis=1)

maxSpin = 2 #np.max(spins)
minSpin = 0 #np.min(spins)

ns = nh3_CIS.shape[1]
fig, ax = plt.subplots()
for i in range(ns):
    if(i <= ns/3):
        lcol = u'#8c564b'
    elif(i < 2*ns/3 ):
        lcol = u'#17becf'
    else:
        lcol = u'#ff7f0e'
    ax.plot(x, nh3_CIS[:,i], color=lcol, linewidth=3)


ax.set_xlabel("r ($\AA$)")
ax.set_ylabel("E (a.u.)")
ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)

custom_lines = [Line2D([0], [0], color=u'#8c564b', lw=3),
                Line2D([0], [0], color=u'#17becf', lw=3),
                Line2D([0], [0], color=u'#ff7f0e', lw=3)]
plt.legend(custom_lines, ['$M_s=0$', '$M_s=-1$','$M_s=+1$'])
plt.tight_layout()
plt.savefig("plots/nh3_tz_sfs.eps", format="eps", dpi=1000)
plt.clf()
plt.close()












# plot XCIS
nStates = 5

nh3_edata = np.genfromtxt(open("NH3/NH3_XCIS_CISEnergyPlot.csv", "r"), delimiter=",")
nh3_sdata = np.genfromtxt(open("NH3/NH3_XCIS_s2Plot.csv", "r"), delimiter=",")


x = nh3_edata[:,0]
spins = nh3_sdata[:,1:nStates+1]
nh3_CIS = nh3_edata[:,1:nStates+1]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, nh3_CIS, spins, nStates, "nh3_tz_x.eps")
