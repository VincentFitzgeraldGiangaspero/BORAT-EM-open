import numpy as np
import pyvista as pv

import Functions

###################################################################

def FinalRay(RaySolutions,Ray0,rho):
    
    totalRays = Ray0.n_points

    pointsRay = np.column_stack((RaySolutions[0, -1, :],
                                RaySolutions[1, -1, :],
                                RaySolutions[2, -1, :]))

    RayPolyData = pv.PolyData(pointsRay)

    # lines = np.column_stack(
    #     (np.full(nSteps-1, 2), np.arange(nSteps - 1), np.arange(nSteps - 1) + 1))

    # RayPolyData = pv.PolyData(pointsRay, lines=lines)


    RayPolyData['phase'] = RaySolutions[9, -1, :]
    RayPolyData['ds'] = RaySolutions[10, -1, :]

    RayPolyData['n'] = RaySolutions[11, -1, :]

    # direction versor
    RayPolyData['kx'] = RaySolutions[3, -1, :]/RayPolyData['n']
    RayPolyData['ky'] = RaySolutions[4, -1, :]/RayPolyData['n']
    RayPolyData['kz'] = RaySolutions[5, -1, :]/RayPolyData['n']

    RayPolyData['ex'] = RaySolutions[6, -1, :]
    RayPolyData['ey'] = RaySolutions[7, -1, :]
    RayPolyData['ez'] = RaySolutions[8, -1, :]

    SphericalE = np.zeros(shape=(3, totalRays))

    RayPolyData['K'] = np.column_stack(
        [RayPolyData['kx'], RayPolyData['ky'], RayPolyData['kz']])

    RayPolyData['Edir'] = np.column_stack([RayPolyData['ex'], RayPolyData['ey'], RayPolyData['ez']])

    RayPolyData['DF'] = rho/(rho+RayPolyData['ds'])

    RayPolyData['Efield'] = np.column_stack([Ray0['E0'][:,0]*RayPolyData['DF'],
                                            Ray0['E0'][:,1]*RayPolyData['DF'],
                                            Ray0['E0'][:,2]*RayPolyData['DF']])    

    for iRay in range(totalRays):
        
        # how do we define a spherical system at the center tube? 
        # for now we use the Xc,Yc,Zc
        
        _, theta, phi = Functions.Cart2sphere(
            RaySolutions[3, -1, iRay], RaySolutions[4, -1, iRay], RaySolutions[5, -1,iRay])
        
        st = np.sin(theta) 	
        ct = np.cos(theta)
        cp = np.cos(phi)		
        sp = np.sin(phi)
        
        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp,ct*sp,-st])
        phi_versor = np.array([-sp,cp,0])
        
        SphericalE[0,iRay] = np.dot(RayPolyData['Efield'][iRay],r_versor)
        SphericalE[1,iRay] =np.dot(RayPolyData['Efield'][iRay],theta_versor)
        SphericalE[2, iRay] = np.dot(RayPolyData['Efield'][iRay], phi_versor)
        

    RayPolyData['Er'] = SphericalE[0,:]
    RayPolyData['Et'] = SphericalE[1,:]
    RayPolyData['Ep'] = SphericalE[2,:]


    RayPolyData['|E|']=np.linalg.norm(RayPolyData['Efield'],axis=1)
        
        
    #if you use modules you can get NAN if Efield is zero
    # RayPolyData['Area'] = (np.linalg.norm(Ray0['E0'],axis=1)**2) /(np.linalg.norm(RayPolyData['Efield'],axis=1)**2) * A0_per_ray
    # RayPolyData['Area']=A0_per_ray/(RayPolyData['DF']**2)
        
        
    return RayPolyData

###################################################################

def CompleteRay(RaySolutions, Ray0, rho):
    
    RayBlock = pv.MultiBlock()

    totalRays = Ray0.n_points
    
    nSteps=np.shape(RaySolutions)[1]

    SphericalE = np.zeros(shape=(3, nSteps))

    for iRay in range (totalRays):

       
        pointsRay = np.column_stack((RaySolutions[0, :, iRay], 
                                     RaySolutions[1, :, iRay], 
                                     RaySolutions[2, :, iRay]))
        
        
        lines = np.column_stack(
            (np.full(nSteps-1, 2), np.arange(nSteps - 1), np.arange(nSteps - 1) + 1))

        RayPolyData = pv.PolyData(pointsRay, lines=lines)


        RayPolyData['phase'] = RaySolutions[9, :, iRay]
        RayPolyData['ds'] = RaySolutions[10, :, iRay]

        RayPolyData['n'] = RaySolutions[11, :, iRay]

        # direction versor
        RayPolyData['kx'] = RaySolutions[3, :, iRay]/RayPolyData['n']
        RayPolyData['ky'] = RaySolutions[4, :, iRay]/RayPolyData['n']
        RayPolyData['kz'] = RaySolutions[5, :, iRay]/RayPolyData['n']

        RayPolyData['ex'] = RaySolutions[6, :, iRay]
        RayPolyData['ey'] = RaySolutions[7, :, iRay]
        RayPolyData['ez'] = RaySolutions[8, :, iRay]

        RayPolyData['K'] = np.column_stack(
            [RayPolyData['kx'], RayPolyData['ky'], RayPolyData['kz']])

        RayPolyData['Edir'] = np.column_stack([RayPolyData['ex'], RayPolyData['ey'], RayPolyData['ez']])

        RayPolyData['DF'] = rho/(rho+RayPolyData['ds'])

        RayPolyData['Efield'] = np.column_stack([Ray0['E0'][iRay,0]*RayPolyData['DF'],
                                                 Ray0['E0'][iRay,1]*RayPolyData['DF'],
                                                 Ray0['E0'][iRay,2]*RayPolyData['DF']])    

        for s in range (nSteps):

            _, theta, phi = Functions.Cart2sphere(
                RaySolutions[3, s, iRay], RaySolutions[4, s, iRay], RaySolutions[5, s,iRay])
            
            st = np.sin(theta) 	
            ct = np.cos(theta)
            cp = np.cos(phi)		
            sp = np.sin(phi)
            
            r_versor = np.array([st*cp, st*sp, ct])
            theta_versor = np.array([ct*cp,ct*sp,-st])
            phi_versor = np.array([-sp,cp,0])
            
            SphericalE[0,s] = np.dot(RayPolyData['Efield'][s,:],r_versor)
            SphericalE[1,s] =np.dot(RayPolyData['Efield'][s,:],theta_versor)
            SphericalE[2, s] = np.dot(RayPolyData['Efield'][s, :], phi_versor)
        

        RayPolyData['Er'] = SphericalE[0,:]
        RayPolyData['Et'] = SphericalE[1,:]
        RayPolyData['Ep'] = SphericalE[2,:]


        RayPolyData['|E|']=np.linalg.norm(RayPolyData['Efield'],axis=1)
        
        RayBlock.append(RayPolyData, name=str(iRay))
        
    return RayBlock

###################################################################

def FinalTubes(RaySolutions,RayTube0,rho):
    
    totalTubes = RayTube0.n_faces
    
    # quantities to be calculated for the tubes
    TubesFinalArea = np.zeros(totalTubes)
    TubesFinalEvector = np.zeros(shape=(3, totalTubes))
    TubesFinalKvector = np.zeros(shape=(3, totalTubes))
    TubesFinalPhase = np.zeros(totalTubes)
    TubesFinalds = np.zeros(totalTubes)
    TubesFinalN = np.zeros(totalTubes)
    TubesDF = np.zeros(totalTubes)
    TubesDF2 = np.zeros(totalTubes)


    TubesInitialArea = RayTube0['Area']


    for nTube in range(totalTubes):

        ray1 = RayTube0.faces[nTube*4+1]
        ray2 = RayTube0.faces[nTube*4+2]
        ray3 = RayTube0.faces[nTube*4+3]

        TubesFinalN[nTube] = (RaySolutions[11, -1, ray1] +
                            RaySolutions[11, -1, ray2] +
                            RaySolutions[11, -1, ray3])/3

        TubesFinalPhase[nTube] = (RaySolutions[9, -1, ray1] +
                                RaySolutions[9, -1, ray2] +
                                RaySolutions[9, -1, ray3])/3

        TubesFinalds[nTube] = (RaySolutions[10, -1, ray1] +
                            RaySolutions[10, -1, ray2] +
                            RaySolutions[10, -1, ray3])/3

        TubesFinalEvector[0, nTube] = (RaySolutions[6, -1, ray1] +
                                    RaySolutions[6, -1, ray2] +
                                    RaySolutions[6, -1, ray3])/3

        TubesFinalEvector[1, nTube] = (RaySolutions[7, -1, ray1] +
                                    RaySolutions[7, -1, ray2] +
                                    RaySolutions[7, -1, ray3])/3

        TubesFinalEvector[2, nTube] = (RaySolutions[8, -1, ray1] +
                                    RaySolutions[8, -1, ray2] +
                                    RaySolutions[8, -1, ray3])/3

        TubesFinalKvector[0, nTube] = (RaySolutions[3, -1, ray1] +
                                    RaySolutions[3, -1, ray2] +
                                    RaySolutions[3, -1, ray3])/(3*TubesFinalN[nTube])

        TubesFinalKvector[1, nTube] = (RaySolutions[4, -1, ray1] +
                                    RaySolutions[4, -1, ray2] +
                                    RaySolutions[4, -1, ray3])/(3*TubesFinalN[nTube])

        TubesFinalKvector[2, nTube] = (RaySolutions[5, -1, ray1] +
                                    RaySolutions[5, -1, ray2] +
                                    RaySolutions[5, -1, ray3])/(3*TubesFinalN[nTube])

        # this is important, check it all the time otherwise the triangle is not the right one

        triangle = pv.Triangle([RaySolutions[:3, -1, ray1],
                                RaySolutions[:3, -1, ray2],
                                RaySolutions[:3, -1, ray3]])

        TubesFinalArea[nTube] = triangle.area

        if TubesFinalArea[nTube] == 0:
            
            TubesDF[nTube] = 0
        else:
            TubesDF[nTube] = np.sqrt(
                TubesInitialArea[nTube]/TubesFinalArea[nTube])

        TubesDF2[nTube] = rho/(rho+TubesFinalds[nTube])

        # print('A0 = ', TubesInitialArea[nTube], '\n')
        # print('Af = ', TubesFinalArea[nTube],'\n' )
        # print('DF = ', TubesDF[nTube], '\n')
        # print('rho/(rho+s) = ', rho/(rho+TubesFinalds[nTube]),'\n' )


    FinalRayTube = pv.PolyData(np.transpose(
        RaySolutions[0:3, -1, :]),  RayTube0.faces)

    FinalRayTube= FinalRayTube.compute_cell_sizes(area=True)
    
    FinalRayTube.cell_data['DF'] = TubesDF
    FinalRayTube.cell_data['DF2'] = TubesDF2

    FinalRayTube.cell_data['n'] = TubesFinalN
    FinalRayTube.cell_data['phase'] = TubesFinalPhase
    FinalRayTube.cell_data['ds'] = TubesFinalds

    FinalRayTube.cell_data['ex'] = TubesFinalEvector[0, :]
    FinalRayTube.cell_data['ey'] = TubesFinalEvector[1, :]
    FinalRayTube.cell_data['ez'] = TubesFinalEvector[2, :]

    FinalRayTube.cell_data['kx'] = TubesFinalKvector[0, :]
    FinalRayTube.cell_data['ky'] = TubesFinalKvector[1, :]
    FinalRayTube.cell_data['kz'] = TubesFinalKvector[2, :]

    FinalRayTube.cell_data['K'] = np.column_stack(
        [FinalRayTube['kx'], FinalRayTube['ky'], FinalRayTube['kz']])

    FinalRayTube.cell_data['Edir'] = np.column_stack(
        [FinalRayTube['ex'], FinalRayTube['ey'], FinalRayTube['ez']])

    FinalRayTube['DF3'] = np.sqrt(
        TubesInitialArea[:]/FinalRayTube['Area'])
    
    """
    fix this
    """
    # FinalRayTube.cell_data['|E|'] = FinalRayTube['DF']*RayTube0['|E0|']
    # FinalRayTube.cell_data['|E|'] = FinalRayTube['DF2']*RayTube0['|E0|']
    FinalRayTube.cell_data['|E|'] = FinalRayTube['DF'] * \
        RayTube0['|E0|']*np.sqrt(FinalRayTube.cell_data['n'])

    FinalRayTube.cell_data['Efield'] = np.column_stack([FinalRayTube['|E|']*FinalRayTube['ex'],
                                                        FinalRayTube['|E|'] *
                                                        FinalRayTube['ey'],
                                                        FinalRayTube['|E|']*FinalRayTube['ez']])

    FinalRayTube.cell_data['Area'] = TubesFinalArea


    xiyizi = FinalRayTube.cell_centers()

    SphericalE = np.zeros(shape=(3, totalTubes))

    for iTube in range(totalTubes):

        # how do we define a spherical system at the center tube?
        # for now we use the Xc,Yc,Zc

        _, theta, phi = Functions.Cart2sphere(
            xiyizi.points[iTube, 0], xiyizi.points[iTube, 1], xiyizi.points[iTube, 2])

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        SphericalE[0, iTube] = np.dot(
            FinalRayTube.cell_data['Efield'][iTube], r_versor)
        SphericalE[1, iTube] = np.dot(
            FinalRayTube.cell_data['Efield'][iTube], theta_versor)
        SphericalE[2, iTube] = np.dot(
            FinalRayTube.cell_data['Efield'][iTube], phi_versor)


    FinalRayTube.cell_data['Er'] = SphericalE[0, :]
    FinalRayTube.cell_data['Etheta'] = SphericalE[1, :]
    FinalRayTube.cell_data['Ephi'] = SphericalE[2, :]

    FinalRayTube.compute_normals()
    FinalRayTube.cell_data['Normal'] = FinalRayTube.face_normals
    
    return FinalRayTube

###################################################################
