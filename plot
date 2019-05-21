#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection




cmap = mpl.colors.LinearSegmentedColormap.from_list("tasetful rbg", colors=[u'#0072bb',  u'#d62728', u'#2ca02c'])


mpl.rcParams["font.size"] = 16
def plot(x, y, spins, ns, name, color=True, bar=True):
    fig, ax = plt.subplots()
    for i in range(ns):
        points = np.concatenate([x.reshape(-1,1), y[:,i].reshape(-1,1)], axis=1).reshape(-1,1,2)
        segs = np.concatenate([points[:-1], points[1:]],axis=1)

        if(color):
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
spins = CIS_Sdata[:,1:]
maxSpin = np.max(spins)
minSpin = np.min(spins)

# Plot Normal CIS
x = CIS_Edata[:,0]
EnergiesCIS = CIS_Edata[:,1:]
plot(x, EnergiesCIS, spins, nStates, "H2_min.eps")







# Plot sf-t-CIS
sf_t_CIS_Edata = np.genfromtxt(open("H2_min/H2_sf_t_CISEnergyPlot.csv", "r"), delimiter=",")
sf_t_Sdata = np.genfromtxt(open("H2_min/H2_sf_t_s2Plot.csv", "r"), delimiter=",")

x = sf_t_CIS_Edata[:,0]
Energies_sftCIS = sf_t_CIS_Edata[:,2:]
spins = sf_t_Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_sftCIS, spins, nStates, "H2_min_sft.eps")




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
nStates = 5
TDHF_Edata = np.genfromtxt(open("H2/augtz/H2_TDHFEnergyPlot.csv", "r"), delimiter=",")
maxSpin = np.max(spins)
minSpin = np.min(spins)

# Plot Normal CIS
x = CIS_Edata[:,0]
EnergiesTDHF = TDHF_Edata[:,2:]
plot(x, EnergiesTDHF, spins, nStates, "H2_tz_tdhf.eps", color=False)












# Triple Zeta SF
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








# Triple Zeta Exact
nStates = 5
nSinglets = 6
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







# PBE0
TDA_Edata = np.genfromtxt(open("H2/augtz/H2_pbe0_CISEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("H2/augtz/H2_pbe0_s2Plot.csv", "r"), delimiter=",")


x = x_CIS_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
spins = Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_TDA, spins, nStates, "H2_tz_pbe0.eps")









# LRC-wPBEh
TDA_Edata = np.genfromtxt(open("H2/augtz/H2_wpbe_CISEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("H2/augtz/H2_wpbe_s2Plot.csv", "r"), delimiter=",")


x = x_CIS_Edata[:,0]
Energies_TDA = TDA_Edata[:,2:]
spins = Sdata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, Energies_TDA, spins, nStates, "H2_tz_wpbe.eps")

























# ----------------- ETHANE --------------------


ymin = -79.3
ymax = -78.8
xmin = 1
xmax = 5




# Plot Normal CIS
nStates = 3
CIS_Edata = np.genfromtxt(open("CH3CH3/dz/CH3CH3_CISEnergyPlot.csv", "r"), delimiter=",")
Sdata = np.genfromtxt(open("CH3CH3/dz/CH3CH3_s2Plot.csv", "r"), delimiter=",")

x = CIS_Edata[:,0]
EnergiesCIS = CIS_Edata[:,1:]
spins = Sdata[:,1:]

plot(x, EnergiesCIS, spins, nStates, "ethane_dz.eps")




# Plot sf CIS
sf_Edata = np.genfromtxt(open("CH3CH3/dz/CH3CH3_sf_t_CISEnergyPlot.csv", "r"), delimiter=",")
sf_Sdata = np.genfromtxt(open("CH3CH3/dz/CH3CH3_sf_t_s2Plot.csv", "r"), delimiter=",")

x = sf_Edata[:,0]
sf_Energies = sf_Edata[:,2:]
spins = sf_Sdata[:,2:]

plot(x, sf_Energies, spins, nStates, "ethane_dz_sf.eps")






# Plot x CIS
nSinglets = 1
nTriplets = 2
x_Edata = np.genfromtxt(open("CH3CH3/dz/CH3CH3_x_t_CISEnergyPlot.csv", "r"), delimiter=",")
x_Sdata = np.genfromtxt(open("CH3CH3/dz/CH3CH3_x_t_s2Plot.csv", "r"), delimiter=",")

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

plot(x, x_Energies, spins, nStates, "ethane_dz_x.eps", bar=False)
















# ----------------------- SILANE ------------------------------
xmin = 1.7
xmax = 5.0
ymin = -581.4
ymax = -580.8
nStates = 3

sf_t_CIS_Edata = np.genfromtxt(open("silane/daugdz/silane_sf_t_CISEnergyPlot.csv", "r"), delimiter=",")
sf_t_Sdata = np.genfromtxt(open("silane/daugdz/silane_sf_t_s2Plot.csv", "r"), delimiter=",")

# Plot SF-CIS
x = sf_t_CIS_Edata[:,0]
spins = sf_t_Sdata[:,2:]
Energies_sf_CIS = sf_t_CIS_Edata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

fig, ax = plt.subplots()

for i in range(nStates):
    points = np.concatenate([x.reshape(-1,1), Energies_sf_CIS[:,i].reshape(-1,1)], axis=1).reshape(-1,1,2)
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
plt.savefig("plots/silane_sf.eps", format="eps", dpi=1000)
plt.clf()






# plot MOM
mom_Edata = np.genfromtxt(open("silane/daugdz_mom/sf-t-MOM-CIS_CISEnergyPlot.csv", "r"), delimiter=",")
mom_Sdata = np.genfromtxt(open("silane/daugdz_mom/sf-t-MOM-CIS_s2Plot.csv", "r"), delimiter=",")

# Plot SF-CIS
x2 = mom_Edata[:,0]
spins = mom_Sdata[:,2:]
mom_CIS = mom_Edata[:,2:]

maxSpin = np.max(spins)
minSpin = np.min(spins)

fig, ax = plt.subplots()

for i in range(nStates):
    points = np.concatenate([x2.reshape(-1,1), mom_CIS[:,i].reshape(-1,1)], axis=1).reshape(-1,1,2)
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

















# ----------------------- NH3 ------------------------------
xmin = 0.5
xmax = 5
ymin = -56.3
ymax = -55.7

# plot CIS
nStates = 5

nh3_edata = np.genfromtxt(open("NH3/NH3_CISEnergyPlot.csv", "r"), delimiter=",")
nh3_sdata = np.genfromtxt(open("NH3/NH3_s2Plot.csv", "r"), delimiter=",")


x = nh3_edata[:,0]
spins = nh3_sdata[:,1:nStates+1]
nh3_CIS = nh3_edata[:,1:nStates+1]

maxSpin = np.max(spins)
minSpin = np.min(spins)

plot(x, nh3_CIS, spins, nStates, "nh3_tz.eps")










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
