import numpy as np


def Cart2sphere(x, y, z):
    '''
    Function transform points stored in Catersian coordinate system to points on spherical coordinate system  
    r is the radial distance 
    theta is the vector of polar angles
    phi is the vector of azimuth angles
    '''
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    return r, theta, phi


def cartesian2spherical_Vector(V_x,V_y,V_z):

    
    V=np.array([V_x,V_y,V_z])
    
    r = np.sqrt(V_x**2 + V_y**2 + V_z**2)
    theta = np.arccos(V_z/r)
    
  
    phi = np.arctan2(V_y, V_x)
    
    st = np.sin(theta) 	
    ct = np.cos(theta)
    cp = np.cos(phi)		
    sp = np.sin(phi)
    
    
    r_versor = np.array([st*cp, st*sp, ct])
    theta_versor = np.array([ct*cp,ct*sp,-st])


    phi_versor = np.array([-sp,cp,0])

    T=[[st*cp, ct*cp, -st],
       [st*sp, ct*sp, cp],
       [ct,-st,0]]
    
    T2=np.transpose(T)
    
    Vnew=np.dot(T2,V)
    
    
    V_r = np.dot(r_versor,V)
    
    V_theta=np.dot(theta_versor,V)
    
    V_phi = np.dot(phi_versor,V)
    
    
    
    return V_r,V_theta,V_phi