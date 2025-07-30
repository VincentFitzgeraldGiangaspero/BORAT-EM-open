import numpy as np


def fact(m):
    #factorial of a  number m
    f = 1
    if m >= 1:
        for n in range(1, m+1):
            f = f*n
    return f
  
 

def G(n,w):
    #recursive function G.
    jw = 1j*w
    g = (np.exp(jw) - 1)/jw
    if n > 0:
        for m in range(1, n+1):
            go = g
            g = (np.exp(jw) - n*go)/jw
            
    return g

###################################################################

def LudwigShapeFactor(bk,x,y,z,vind,Ri,Robs,Area,Co):
    
    Lt = 1e-5
    #in python range is +1
    Nt=5+1
    
    u=Robs[0]
    v=Robs[1]
    w=Robs[2]
        
    ui=Ri[0]
    vi=Ri[1]    
    wi=Ri[2]    
    
    # s= bk(u+ui...)
    
    Dp = bk*((x[vind[0]] - x[vind[2]])*(u+ui) + (y[vind[0]] - y[vind[2]])*(v+vi) + (z[vind[0]] - z[vind[2]])*(w+wi))
    Dq = bk*((x[vind[1]] - x[vind[2]])*(u+ui) + (y[vind[1]] - y[vind[2]])*(v+vi) + (z[vind[1]] - z[vind[2]])*(w+wi))
    Do = bk * (x[vind[2]]*(u+ui) + y[vind[2]]*(v+vi) + z[vind[2]]*(w+wi))



    # Area integral for general case
    DD = Dq - Dp
    expDo = np.exp(1j*Do)
    expDp = np.exp(1j*Dp)
    expDq = np.exp(1j*Dq)



    # Special case 1
    if abs(Dp) < Lt and abs(Dq) >= Lt:
        sic=0
        for n in range(Nt):
            sic = sic + (((1j*Dp)**n)/ fact(n))*(-Co/(n+1)+expDq*(Co*G(n,-Dq)))
            
        Ic=sic*2*Area*expDo/(1j*Dq)
        
    # Special case 2
    elif np.abs(Dp) < Lt and np.abs(Dq) < Lt:
        sic = 0.
        for n in range(Nt):
            for nn in range(Nt):
                sic = sic+Co*((1j*Dp)**n)*((1j*Dq)**nn)/fact(nn+n+2)
            
        Ic = sic*2*Area*expDo
        
    # Special case 3
    elif np.abs(Dp) >= Lt and np.abs(Dq) < Lt:
        sic = 0.
        for n in range(Nt):
            sic = sic+(((1j*Dq)**n)/fact(n))*(Co/(n+1))*G(n+1,-Dp)
            
        Ic = sic*2*Area*expDo*expDp
    # Special case 4
    elif np.abs(Dp) >= Lt and np.abs(Dq) >= Lt and np.abs(DD) < Lt:
        sic = 0.
        for n in range(Nt):
            sic = sic+(((1j*DD)**n)/fact(n))*(-Co*G(n,Dq)+expDq*Co/(n+1))
        
        Ic = sic*2*Area*expDo/(1j*Dq)
            
    else:
            Ic = 2*Area*expDo*(expDp*Co/(Dp*DD)-expDq*Co/(Dq*DD)-Co/(Dp*Dq))

    return Ic

###################################################################


def GisbonShapeFactor(bk, x, y, z, vind, Ri, Robs, Area):
    
    u = Robs[0]
    v = Robs[1]
    w = Robs[2]

    ui = Ri[0]
    vi = Ri[1]
    wi = Ri[2]
    
    s=bk*(Robs+Ri)
    
    v1 = np.array([x[vind[0]],y[vind[0]],z[vind[0]]])
    v2 = np.array([x[vind[1]],y[vind[1]],z[vind[1]]])
    v3 = np.array([x[vind[2]],y[vind[2]],z[vind[2]]])

    a1=np.dot(s,v1)
    a2=np.dot(s,v2)
    a3=np.dot(s,v3)
    
    d_A2A1=a2-a1
    d_A3A1=a3-a1
    d_A3A2=a3-a2
    
    tolerance=1e-5
    
    if abs(d_A2A1) <  tolerance  and abs(d_A3A1) > tolerance and abs(d_A3A2) > tolerance :
        
        Ic = 2*Area/(a3-a2)*(1j*np.exp(1j*a1) - (np.exp(1j*a1)-np.exp(1j*a3))/(a1-a3 ))
    
    elif abs(d_A3A1) < tolerance and abs(d_A2A1) > tolerance and abs(d_A3A2) > tolerance:

        Ic = 2*Area/(a3-a2)*(-1j*np.exp(1j*a1) + (np.exp(1j*a1)-np.exp(1j*a2))/(a1-a2 ))

    elif abs(d_A3A2) < tolerance and abs(d_A2A1) > tolerance and abs(d_A3A1) > tolerance:

        Ic = 2*Area/(a1-a2)*(1j*np.exp(1j*a3) - (np.exp(1j*a1)-np.exp(1j*a3))/(a1-a3 ))

    elif abs(d_A3A2) < tolerance and abs(d_A2A1) < tolerance and abs(d_A3A1) < tolerance:

        Ic = Area*np.exp(1j*a1)
        
    else:
        
        Ic  =  2*Area/(a3-a2)*( (np.exp(1j*a1)-np.exp(1j*a2))/(a1-a2)  - (np.exp(1j*a1)-np.exp(1j*a3))/(a1-a3))
    
    return Ic

###################################################################


def GisbonShapeFactorModified(bk, x, y, z, vind, Ri, Robs, Area,Rc):

    u = Robs[0]
    v = Robs[1]
    w = Robs[2]

    ui = Ri[0]
    vi = Ri[1]
    wi = Ri[2]
    
    Xc=Rc[0]
    Yc=Rc[1]
    Zc=Rc[2]

    s = bk*(Robs+Ri)

    v1 = np.array([x[vind[0]]-Xc, y[vind[0]]-Yc, z[vind[0]]-Zc])
    v2 = np.array([x[vind[1]]-Xc, y[vind[1]]-Yc, z[vind[1]]-Zc])
    v3 = np.array([x[vind[2]]-Xc, y[vind[2]]-Yc, z[vind[2]]-Zc])

    a1 = np.dot(s, v1)
    a2 = np.dot(s, v2)
    a3 = np.dot(s, v3)

    d_A2A1 = a2-a1
    d_A3A1 = a3-a1
    d_A3A2 = a3-a2

    tolerance = 1e-5

    if abs(d_A2A1) < tolerance and abs(d_A3A1) > tolerance and abs(d_A3A2) > tolerance:

        Ic = 2*Area/(a3-a2)*(1j*np.exp(1j*a1) -
                             (np.exp(1j*a1)-np.exp(1j*a3))/(a1-a3))

    elif abs(d_A3A1) < tolerance and abs(d_A2A1) > tolerance and abs(d_A3A2) > tolerance:

        Ic = 2*Area/(a3-a2)*(-1j*np.exp(1j*a1) +
                             (np.exp(1j*a1)-np.exp(1j*a2))/(a1-a2))

    elif abs(d_A3A2) < tolerance and abs(d_A2A1) > tolerance and abs(d_A3A1) > tolerance:

        Ic = 2*Area/(a1-a2)*(1j*np.exp(1j*a3) -
                             (np.exp(1j*a1)-np.exp(1j*a3))/(a1-a3))

    elif abs(d_A3A2) < tolerance and abs(d_A2A1) < tolerance and abs(d_A3A1) < tolerance:

        Ic = Area*np.exp(1j*a1)

    else:

        Ic = 2*Area/(a3-a2)*((np.exp(1j*a1)-np.exp(1j*a2)) /
                             (a1-a2) - (np.exp(1j*a1)-np.exp(1j*a3))/(a1-a3))

    return Ic
