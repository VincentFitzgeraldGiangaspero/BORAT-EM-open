import matplotlib.pyplot as plt
import math
import numpy as np
import SolutionPLK
import POPlots
import os
from matplotlib import rc
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

#########################################################################
# region

c = 299792458.0  # Speed of light in m/s
Frequency = 400e6
Rho = 1e-2

outputFolder= 'ThesisSolutions'

if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
    
#endregion
#########################################################################

importFile = 'ThesisSolutions/ExoMarsDebug_SPHERE_NONionized_80_coarse20k_phiSBR.plk'

solution = SolutionPLK.importSolutionList(importFile)


thetadeg = solution[0]
phideg = solution[1]
Ethscat = solution[2]

EthTOT = np.sqrt(abs(Ethscat[0, :, :])**2 + abs(Ethscat[1, :, :])**2)

#########################################################################
""" 
#PLots
"""
# region

phi4plot = 0

POPlots.SimplePlot(
    phideg[phi4plot], Ethscat[:, :, phi4plot], outputFolder)

POPlots.PlotApertureExoMars(
    phideg[phi4plot], Ethscat[:, :, phi4plot], outputFolder)


# endregion
#########################################################################


fig, axs = plt.subplots(1, 3, figsize=(10, 5))

axs[0].pcolormesh(phideg, thetadeg, EthTOT, cmap=plt.cm.jet)

axs[1].pcolormesh(phideg, thetadeg, abs(Ethscat[0, :, :]), cmap=plt.cm.jet)

axs[2].pcolormesh(phideg, thetadeg, abs(Ethscat[1, :, :]), cmap=plt.cm.jet)
fig.colorbar()

plt.figure()
plt.pcolormesh(phideg, thetadeg, EthTOT, cmap=plt.cm.jet)
plt.gca().invert_yaxis() #theta=0  on top
plt.colorbar()
plt.title('Directivity, linear scale')


phi4plot  = 0


EthTOT = np.sqrt(Ethscat[0, :, :]**2 + Ethscat[1, :, :]**2)

anglerad = np.deg2rad(thetadeg[phi4plot])
angle = thetadeg[phi4plot]

gEtot = np.where(abs(EthTOT[0, :]) == 0,
                np.finfo(float).eps, abs(EthTOT[0, :]))
gdBEtot = 10*np.log10(gEtot)


fig = plt.figure(figsize=(10, 5))
fig.suptitle('Farfield Radiation Amplitude Φ component')
gs = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])

ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(angle, gdBEtot,label='E')
ax1.grid(True)
ax1.set_xlabel('θ [deg]')
ax1.set_ylabel('dBV/m')
ax1.set_xlim([min(angle), max(angle)])
ax1.legend(loc='upper right')

ax2 = fig.add_subplot(gs[0, 1], projection='polar')
ax2.set_theta_zero_location('N')
ax2.set_theta_direction(-1)
lineObjects = ax2.plot(anglerad, gdBEtot, label='E')

ax2.grid(True)
ax2.set_rlabel_position(-30)  # Move radial labels away from plotted line
angleLegend = np.deg2rad(67.5)
ax2.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))

fig.tight_layout()

plt.savefig(outputFolder+"/Etot.png")





