import matplotlib.pyplot as plt
import numpy as np

#check this
from matplotlib import rc
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

#mpl.rcParams['figure.dpi'] = 300
mpl.rc('text', usetex=False)

###################################################################

def sph2cart(r, theta, phi):
  '''
  Function transform points stored in spherical coordinate system to points on Catersian coordinate system  
  r is the radial distance 
  theta is the vector of polar angles
  phi is the vector of azimuth angles
  '''
  x = r * np.cos(phi)*np.sin(theta)
  y = r * np.sin(phi)*np.sin(theta)
  z = r * np.cos(theta)
  return x, y, z

###################################################################
def Plot_ThetaRCS_Efield(angle,Ethscat,outputFolder):

  anglerad=np.deg2rad(angle)

  fig = plt.figure(figsize=(10, 5))
  ax1 = fig.add_subplot(122, projection='polar')


  ax1.set_theta_zero_location("N")
  ax1.set_theta_direction(-1)
  lineObjects = ax1.plot(anglerad, (Ethscat[0, :]), label='Eθ')
  lineObjects = ax1.plot(anglerad, (Ethscat[1, :]), label='EΦ')
  ax1.set_title('E')
  angleLegend = np.deg2rad(67.5)
  ax1.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))
  # ax1.set_rticks([0, 0.5, 1])  # Less radial ticks
  ax1.set_rlabel_position(-30)  # Move radial labels away from plotted line
  
  ax2 = fig.add_subplot(121)
  ax2.plot(angle, (Ethscat[0,:]),label='Eθ')
  ax2.plot(angle, (Ethscat[1,:]), label='EΦ')
  ax2.grid(True)
  ax2.set_xlabel('θ [deg]')
  ax2.set_ylabel('E')
  ax2.set_xlim([min(angle), max(angle)])
  ax2.legend(loc='upper right')

  plt.savefig(outputFolder+"/ThetaRCS_Efield.png")

  return

###################################################################
def Plot_ThetaRCS(angle,Ethscat,outputFolder):

  anglerad=np.deg2rad(angle)

  fig = plt.figure(figsize=(10, 5))
  ax1 = fig.add_subplot(122, projection='polar')


  ax1.set_theta_zero_location("N")
  ax1.set_theta_direction(-1)
  lineObjects = ax1.plot(anglerad, (Ethscat[0, :]), label='RCSθ')
  lineObjects = ax1.plot(anglerad, (Ethscat[1, :]), label='RCSΦ')
  ax1.set_title('RCS [dBsm]')
  angleLegend = np.deg2rad(67.5)
  ax1.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))
  ax1.set_rmax(30)
  ax1.set_rmin(-50)
  # ax1.set_rticks([0, 0.5, 1])  # Less radial ticks
  ax1.set_rlabel_position(-30)  # Move radial labels away from plotted line
  
  ax2 = fig.add_subplot(121)
  ax2.plot(angle, (Ethscat[0,:]),label='RCSθ')
  ax2.plot(angle, (Ethscat[1,:]), label='RCSΦ')
  ax2.grid(True)
  ax2.set_xlabel('θ [deg]')
  ax2.set_ylabel('RCS[dBsm]')
  ax2.set_xlim([min(angle), max(angle)])
  ax2.set_ylim([-50,30])
  ax2.legend(loc='upper right')

  plt.savefig(outputFolder+"/ThetaRCS.png")

  return
###################################################################
def PolarPlot_ThetaRCS(angle,Ethscat,outputFolder):
    
  angle=np.deg2rad(angle)
  # Ethscat=np.flip(Ethscat,axis=1)

  fig, (ax1,ax2) = plt.subplots(1,2, subplot_kw={'projection': 'polar'}, figsize=(10, 5))
  ax1.set_theta_zero_location("N")
  ax1.set_theta_direction(-1)
  # ax1.set_theta_offset(np.pi/2.0)
  lineObjects = ax1.plot(angle, (Ethscat[0, :]), label='Aθ')
  lineObjects = ax1.plot(angle, (Ethscat[1, :]), label='AΦ')
  ax1.set_title('Farfield Radiation Amplitude A(θ,Φ)')
  angleLegend = np.deg2rad(67.5)
  ax1.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))  
  # ax1.set_rmax(1.2)
  # ax1.set_rmin(0)
  # ax1.set_rticks([0, 0.5, 1])  # Less radial ticks
  ax1.set_rlabel_position(-30)  # Move radial labels away from plotted line

  # ax1.set_thetamin(min(angle))
  # ax1.set_thetamax(max(angle))

  # To avoid divide by zeros error in taking log10 replace 0 values by machine eps
  gtheta = np.where(abs(Ethscat[0, :]) == 0,
                    np.finfo(float).eps, abs(Ethscat[0, :]))
  gdBtheta = 10*np.log10(gtheta)
  # gain clipped to minimum -50 dB level for plot
  gdBtheta = np.where(gdBtheta <= -50, -50, gdBtheta)

  gphi = np.where(abs(Ethscat[1, :]) == 0,
                    np.finfo(float).eps, abs(Ethscat[1, :]))
  gdBphi = 10*np.log10(gphi)
  # gain clipped to minimum -50 dB level for plot
  gdBphi = np.where(gdBphi <= -50, -50, gdBphi)


  ax2.set_theta_zero_location("N")
  lineObjects = ax2.plot(angle,gdBtheta, label='Gain θ')
  lineObjects = ax2.plot(angle,gdBphi, label='Gain Φ')
  ax2.set_theta_direction(-1)
  # ax2.set_theta_offset(np.pi/2.0)
  ax2.set_title('Dipole gain (dB units)')
  ax2.set_rmax(1)
  ax2.set_rmin(-50)
  # ax2.set_rticks([-50, -40, -20, 0])  # Less radial ticks
  ax2.set_rlabel_position(-30)  # Move radial labels away from plotted line

  angleLegend = np.deg2rad(67.5)
  ax2.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))  
  # ax2.set_thetamin(min(angle))
  # ax2.set_thetamax(max(angle))

  plt.savefig(outputFolder+"/PolarPlot_ThetaRCS.png")

  return

###################################################################

def CartesianPlot_ThetaRCS(angle, Ethscat, outputFolder):
    
  
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    fig.suptitle('Farfield Radiation Amplitude A(θ,Φ)')

    axs[0, 0].plot(angle, np.real(Ethscat[0,:]),
                'o-', label='Real')
    axs[0, 0].plot(angle, np.imag(Ethscat[0,:]),
                'ro-', label='Imag')
    axs[0, 0].grid(True)
    axs[0, 0].set_xlabel('θ [deg]')
    axs[0, 0].set_ylabel('Aθ')
    axs[0, 0].set_xlim([min(angle), max(angle)])
    axs[0, 0].legend(loc='upper right')
    axs[0, 0].set_ylim([-50,30])

    axs[0, 1].plot(angle, (Ethscat[0, :]), 'o-')
    axs[0, 1].grid(True)
    axs[0, 1].set_xlabel('θ [deg]')
    axs[0, 1].set_ylabel('|Aθ|')
    axs[0, 1].set_xlim([min(angle), max(angle)])
    axs[0, 1].set_ylim([-50,30])
    # axs[0, 1].legend(loc='upper right')

    axs[1, 0].plot(angle, np.real(Ethscat[1,:]),
                'o-', label='Real')
    axs[1, 0].plot(angle, np.imag(Ethscat[1,:]),
                'ro-', label='Imag')
    axs[1, 0].grid(True)
    axs[1, 0].set_xlabel('θ [deg]')
    axs[1, 0].set_ylabel('AΦ')
    axs[1, 0].set_xlim([min(angle), max(angle)])
    # axs[1, 0].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])    
    axs[1, 0].legend(loc='upper right')

    axs[1, 1].plot(angle, (Ethscat[1, :]), 'o-')
    axs[1, 1].grid(True)
    axs[1, 1].set_xlabel('θ [deg]')
    axs[1, 1].set_ylabel('|AΦ|')
    axs[1, 1].set_xlim([min(angle), max(angle)])
    # axs[1, 1].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])
    # axs[1, 1].legend(loc='upper right')

    fig.tight_layout()
        
    plt.savefig(outputFolder+"/CartesianPlot_ThetaRCS.png")
    
    
    return
  
###################################################################

def PolarPlot_Theta(angle,Ethscat,outputFolder):
    
  angle=np.deg2rad(angle)
  # Ethscat=np.flip(Ethscat,axis=1)

  fig, (ax1,ax2) = plt.subplots(1,2, subplot_kw={'projection': 'polar'}, figsize=(10, 5))
  ax1.set_theta_zero_location("N")
  ax1.set_theta_direction(-1)
  # ax1.set_theta_offset(np.pi/2.0)
  lineObjects = ax1.plot(angle, abs(Ethscat[0, :]), label='Aθ')
  lineObjects = ax1.plot(angle, abs(Ethscat[1, :]), label='AΦ')
  ax1.set_title('Farfield Radiation Amplitude A(θ,Φ)')
  angleLegend = np.deg2rad(67.5)
  ax1.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))
  
  
  # ax1.set_rmax(1.2)
  # ax1.set_rmin(0)
  # ax1.set_rticks([0, 0.5, 1])  # Less radial ticks
  ax1.set_rlabel_position(-30)  # Move radial labels away from plotted line

  # ax1.set_thetamin(min(angle))
  # ax1.set_thetamax(max(angle))

  # To avoid divide by zeros error in taking log10 replace 0 values by machine eps
  gtheta = np.where(abs(Ethscat[0, :]) == 0,
                    np.finfo(float).eps, abs(Ethscat[0, :]))
  gdBtheta = 10*np.log10(gtheta)
  # gain clipped to minimum -50 dB level for plot
  gdBtheta = np.where(gdBtheta <= -50, -50, gdBtheta)

  gphi = np.where(abs(Ethscat[1, :]) == 0,
                    np.finfo(float).eps, abs(Ethscat[1, :]))
  gdBphi = 10*np.log10(gphi)
  # gain clipped to minimum -50 dB level for plot
  gdBphi = np.where(gdBphi <= -50, -50, gdBphi)


  ax2.set_theta_zero_location("N")
  lineObjects = ax2.plot(angle,gdBtheta, label='Gain θ')
  lineObjects = ax2.plot(angle,gdBphi, label='Gain Φ')
  ax2.set_theta_direction(-1)
  # ax2.set_theta_offset(np.pi/2.0)
  ax2.set_title('Dipole gain (dB units)')
  ax2.set_rmax(1)
  ax2.set_rmin(-50)
  # ax2.set_rticks([-50, -40, -20, 0])  # Less radial ticks
  ax2.set_rlabel_position(-30)  # Move radial labels away from plotted line
  angleLegend = np.deg2rad(67.5)
  ax2.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))
  # ax2.set_thetamin(min(angle))
  # ax2.set_thetamax(max(angle))

  plt.savefig(outputFolder+"/PolarPlot_Theta.png")

  return

###################################################################

def CartesianPlot_Theta(angle, Ethscat, outputFolder):
    
      
    ExactSolution = np.zeros_like(np.array(angle))
    ExactSolution = np.sin(np.array(angle)*np.pi/180)
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    fig.suptitle('Farfield Radiation Amplitude A(θ,Φ)')

    axs[0, 0].plot(angle, np.real(Ethscat[0,:]),
                'o-', label='Real')
    axs[0, 0].plot(angle, np.imag(Ethscat[0,:]),
                'ro-', label='Imag')
    axs[0, 0].grid(True)
    axs[0, 0].set_xlabel('θ [deg]')
    axs[0, 0].set_ylabel('Aθ')
    axs[0, 0].set_xlim([min(angle), max(angle)])
    # axs[0, 0].set_ylim([-1, 1])
    axs[0, 0].legend(loc='upper right')

    axs[0, 1].plot(angle, abs(Ethscat[0, :]), 'o-')
    # axs[0, 1].plot(angle, abs(ExactSolution), 'r', label='Exact')

    axs[0, 1].grid(True)
    axs[0, 1].set_xlabel('θ [deg]')
    axs[0, 1].set_ylabel('|Aθ|')
    axs[0, 1].set_xlim([min(angle), max(angle)])
    # axs[0, 1].set_ylim([0, 1])
    # axs[0, 1].set_ylim([0, 1.1*max((Ethscat[0,:]))])
    # axs[0, 1].legend(loc='upper right')

    axs[1, 0].plot(angle, np.real(Ethscat[1,:]),
                'o-', label='Real')
    axs[1, 0].plot(angle, np.imag(Ethscat[1,:]),
                'ro-', label='Imag')
    axs[1, 0].grid(True)
    axs[1, 0].set_xlabel('θ [deg]')
    axs[1, 0].set_ylabel('AΦ')
    axs[1, 0].set_xlim([min(angle), max(angle)])
    # axs[1, 0].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])    
    axs[1, 0].legend(loc='upper right')

    axs[1, 1].plot(angle, abs(Ethscat[1, :]), 'o-')
    axs[1, 1].grid(True)
    axs[1, 1].set_xlabel('θ [deg]')
    axs[1, 1].set_ylabel('|AΦ|')
    axs[1, 1].set_xlim([min(angle), max(angle)])
    # axs[1, 1].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])
    # axs[1, 1].legend(loc='upper right')

    fig.tight_layout()
        
    plt.savefig(outputFolder+"/CartesianPlot_Theta.png")
    
    
    return

###################################################################  

def CartesianPlot_Phi(angle,Ethscat,outputFolder):
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    fig.suptitle('Farfield Radiation Amplitude A(θ,Φ)')

    axs[0, 0].plot(angle, np.real(Ethscat[0,:]),
                'o-', label='Real')
    axs[0, 0].plot(angle, np.imag(Ethscat[0,:]),
                'ro-', label='Imag')
    axs[0, 0].grid(True)
    axs[0, 0].set_xlabel('Φ [deg]')
    axs[0, 0].set_ylabel('Aθ')
    axs[0, 0].set_xlim([min(angle), max(angle)])
    axs[0, 0].set_ylim([min(Ethscat[0,:]), 1.1*max((Ethscat[0,:]))])
    axs[0, 0].legend(loc='upper right')

    axs[0, 1].plot(angle, abs(Ethscat[0, :]), 'o-')
    axs[0, 1].grid(True)
    axs[0, 1].set_xlabel('Φ [deg]')
    axs[0, 1].set_ylabel('|Aθ|')
    axs[0, 1].set_xlim([min(angle), max(angle)])
    axs[0, 1].set_ylim([0, 1.1*max((Ethscat[0,:]))])
    axs[0, 1].legend(loc='upper right')

    axs[1, 0].plot(angle, np.real(Ethscat[1,:]),
                'o-', label='Real')
    axs[1, 0].plot(angle, np.imag(Ethscat[1,:]),
                'ro-', label='Imag')
    axs[1, 0].grid(True)
    axs[1, 0].set_xlabel('Φ [deg]')
    axs[1, 0].set_ylabel('AΦ')
    axs[1, 0].set_xlim([min(angle), max(angle)])
    axs[1, 0].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])    
    axs[1, 0].legend(loc='upper right')

    axs[1, 1].plot(angle, abs(Ethscat[1, :]), 'o-')
    axs[1, 1].grid(True)
    axs[1, 1].set_xlabel('Φ [deg]')
    axs[1, 1].set_ylabel('|AΦ|')
    axs[1, 1].set_xlim([min(angle), max(angle)])
    axs[1, 1].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])
    axs[1, 1].legend(loc='upper right')

    fig.tight_layout()
        
    plt.savefig(outputFolder+"/CartesianPlot_Phi.png")
    
    
    return

###################################################################

"""
Radiation Pattern

"""


def RadiationPattern3D(outputFolder):
    
    numPts=100
    theta = np.linspace(start = 0, stop = np.pi, num = numPts) # N point vector for polar angle theta ranging from o to pi
    phi = np.linspace(start = 0, stop = 2*np.pi, num = numPts) # N point vectors for azimuth angle ranging from o to 2pi
    [P,T] = np.meshgrid(phi,theta) # NxN coordinate meshgrid matrix
    G = dipole(T,L) # recompute the normalized power gain with these new vectors

    [X,Y,Z] = sph2cart(G, T, P) # transform points stored in spherical coordinate system to points on Catersian coordinate system

    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(12, 10))
    ax.plot_surface(X,Y,Z, rstride=1, cstride=1,cmap='winter')
    ax.set_title('Normalized power gain of Dipole (L = ' + str(L) + ') antenna')
    ax.set_xlabel('X axis');ax.set_ylabel('Y axis');ax.set_zlabel('Z axis')

    
    
    plt.savefig(outputFolder+"/RadiationPattern3D.png")

    return


def dipole(theta, L=0):
  '''
  Calculate the normalized power gain of dipole antenna of length L specified in terms of wavelength lambda and theta is an numpy array of
  equally spaced polar angles (angle with respect to polar axis) in the range 0 to 2pi. Note: L=0 implies infinitesimally small dipole which is the Hertzian dipole.
  '''
  if L == 0:
    g = np.sin(theta)**2  # power gain for Hertzian dipole
  else:
    with np.errstate(divide='ignore', invalid='ignore'):
      g = ((np.cos(np.pi*L*np.cos(theta)) - np.cos(np.pi*L)) /
           np.sin(theta))**2  # power gain for ordinary dipole

  g[theta == 0] = 0  # handle nan points
  g[theta == np.pi] = 0  # handle nan points

  g = g/np.max(g)  # normalize the power gain
  return g


###################################################################


def PlotApertureExoMars(angle,Ethscat,outputFolder):

  gtheta = np.where(abs(Ethscat[0, :]) == 0,
                    np.finfo(float).eps, abs(Ethscat[0, :]))
  gdBtheta = 10*np.log10(gtheta)


  gphi = np.where(abs(Ethscat[1, :]) == 0,
                  np.finfo(float).eps, abs(Ethscat[1, :]))
  gdBphi = 10*np.log10(gphi)

  anglerad = np.deg2rad(angle)

  fig = plt.figure(figsize=(10, 5))
  fig.suptitle('Farfield Radiation Amplitude A(θ,Φ)')
  gs = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])

  ax1 = fig.add_subplot(gs[0, 0])
  ax1.plot(angle, gdBtheta,label='Aθ')
  ax1.plot(angle, gdBphi,label='AΦ')
  ax1.grid(True)
  ax1.set_xlabel('θ [deg]')
  ax1.set_ylabel('dBV/m')
  ax1.set_xlim([min(angle), max(angle)])
  ax1.legend(loc='upper right')
  
  ax2 = fig.add_subplot(gs[0, 1], projection='polar')
  ax2.set_theta_zero_location('N')
  ax2.set_theta_direction(-1)
  lineObjects = ax2.plot(anglerad, gdBtheta, label='Gain θ')
  lineObjects = ax2.plot(anglerad, gdBphi, label='Gain Φ')
  
  ax2.grid(True)
  ax2.set_rlabel_position(-30)  # Move radial labels away from plotted line
  angleLegend = np.deg2rad(67.5)
  ax2.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))

  fig.tight_layout()
  
  plt.savefig(outputFolder+"/Plotall.png")

  fig = plt.figure(figsize=(10, 5))
  fig.suptitle('Farfield Radiation Theta Component')
  gs = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])

  ax1 = fig.add_subplot(gs[0, 0])
  ax1.plot(angle, gdBtheta, label='Aθ')
  ax1.grid(True)
  ax1.set_xlabel('θ [deg]')
  ax1.set_ylabel('dBV/m')
  ax1.set_xlim([min(angle), max(angle)])
  ax1.legend(loc='upper right')

  ax2 = fig.add_subplot(gs[0, 1], projection='polar')
  ax2.set_theta_zero_location('N')
  ax2.set_theta_direction(-1)
  lineObjects = ax2.plot(anglerad, gdBtheta, label='Gain θ')

  ax2.grid(True)
  ax2.set_rlabel_position(-30)  # Move radial labels away from plotted line
  angleLegend = np.deg2rad(67.5)
  ax2.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))

  fig.tight_layout()
  

  plt.savefig(outputFolder+"/PlotTheta.png")
  
  


  fig = plt.figure(figsize=(10, 5))
  fig.suptitle('Farfield Radiation Amplitude Φ component')
  gs = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])

  ax1 = fig.add_subplot(gs[0, 0])
  ax1.plot(angle, gdBphi,label='AΦ')
  ax1.grid(True)
  ax1.set_xlabel('θ [deg]')
  ax1.set_ylabel('dBV/m')
  ax1.set_xlim([min(angle), max(angle)])
  ax1.legend(loc='upper right')
  
  ax2 = fig.add_subplot(gs[0, 1], projection='polar')
  ax2.set_theta_zero_location('N')
  ax2.set_theta_direction(-1)
  lineObjects = ax2.plot(anglerad, gdBphi, label='Gain Φ')
  
  ax2.grid(True)
  ax2.set_rlabel_position(-30)  # Move radial labels away from plotted line
  angleLegend = np.deg2rad(67.5)
  ax2.legend(loc="lower left",
             bbox_to_anchor=(.5 + np.cos(angleLegend)/2, .5 + np.sin(angleLegend)/2))

  fig.tight_layout()
  

  plt.savefig(outputFolder+"/PlotPhi.png")
  
  
def SimplePlot(angle, Ethscat, outputFolder):
  
  gtheta = np.where(abs(Ethscat[0, :]) == 0,
                    np.finfo(float).eps, abs(Ethscat[0, :]))
  gdBtheta = 10*np.log10(gtheta)


  gphi = np.where(abs(Ethscat[1, :]) == 0,
                  np.finfo(float).eps, abs(Ethscat[1, :]))
  gdBphi = 10*np.log10(gphi)

  fig, axs = plt.subplots(2, 2, figsize=(10, 10))
  fig.suptitle('Farfield Radiation Amplitude A(θ,Φ)')

  axs[0, 0].plot(angle, np.real(Ethscat[0, :]),
                label='Real')
  axs[0, 0].plot(angle, np.imag(Ethscat[0, :]),
                label='Imag')
  axs[0, 0].grid(True)
  axs[0, 0].set_xlabel('θ [deg]')
  axs[0, 0].set_ylabel('dBV/m')
  axs[0, 0].set_xlim([min(angle), max(angle)])
  axs[0, 0].legend(loc='upper right')

  axs[0, 1].plot(angle, gdBtheta)
  axs[0, 1].grid(True)
  axs[0, 1].set_xlabel('θ [deg]')
  axs[0, 1].set_ylabel('dBV/m')
  axs[0, 1].set_xlim([min(angle), max(angle)])
  # axs[0, 1].set_ylim([0, 1.1*max((Ethscat[0,:]))])
  # axs[0, 1].legend(loc='upper right')

  axs[1, 0].plot(angle, np.real(Ethscat[1, :]),
                label='Real')
  axs[1, 0].plot(angle, np.imag(Ethscat[1, :]),
                label='Imag')
  axs[1, 0].grid(True)
  axs[1, 0].set_xlabel('θ [deg]')
  axs[1, 0].set_ylabel('AΦ')
  axs[1, 0].set_xlim([min(angle), max(angle)])
  # axs[1, 0].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])
  axs[1, 0].legend(loc='upper right')

  axs[1, 1].plot(angle, gdBphi)
  axs[1, 1].grid(True)
  axs[1, 1].set_xlabel('θ [deg]')
  axs[1, 1].set_ylabel('dBV/m')
  axs[1, 1].set_xlim([min(angle), max(angle)])
  # axs[1, 1].set_ylim([min(Ethscat[0, :]), 1.1*max((Ethscat[0, :]))])
  # axs[1, 1].legend(loc='upper right')

  fig.tight_layout()

  plt.savefig(outputFolder+"/CartesianPlotSimple.png")