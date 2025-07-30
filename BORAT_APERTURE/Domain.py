import numpy as np
import pyvista as pv
import pytecio


def create(Solver):

    iVarRefractive = 0
    iVarGradient = 0

    data = pytecio.read(Solver.CFDsolution)

    print('CFD Solution : ', Solver.CFDsolution)
    print('Number of Zones : ' ,data.numZones)
    print('Name of Zones : ',data.nameZones)
    print('Number of Variables : ', data.numVars)
    print('Name of Variables : ', data.nameVars)
    print('Input Dimension : ', Solver.Dimension)
    print("\n")
    
    dimension = Solver.Dimension
    
    # refractive index and gradient in input
    for i in range(len(data.nameVars)):
        

        if data.nameVars[i] == 'Ri' or  data.nameVars[i] == 'ri':
            iVarRefractive = i
        
        if data.nameVars[i] == 'X Density Gradient':
            iVarGradient = i
            
        # if data.nameVars[i] == 'gradientX':
        #     iVarGradient = i


    # generate domain mesh
    pvMesh = pv.UnstructuredGrid()
    
    pvMesh.add_field_data([dimension],'dimension')
    
   
    for z in range(data.numZones):
        
        print("*\n")
        print('Name of Zone : ',data.zone_info[z]['name'])

        nNode=data.zone_info[z]['IJK'][0]
        nElem=data.zone_info[z]['IJK'][1]
        eleType=data.zone_info[z]['IJK'][2]
        
        print('N of Points : ',nNode)
        print('N of Elements : ',nElem)

        celltypes = np.empty(nElem, dtype=np.uint8)

        if eleType == 4 and dimension == '2D':
            celltypes[:] = pv.CellType.QUAD
            print('Element Type : ',eleType, 'QUAD')
        elif eleType == 4 and dimension == '3D':
            celltypes[:] = pv.CellType.TETRA
            print('Element Type : ',eleType, 'TETRA')
        elif eleType == 3:
            celltypes[:] = pv.CellType.TRIANGLE
            print('Element Type : ',eleType, 'TRIANGLE')
        elif eleType == 8:
            celltypes[:] = pv.CellType.HEXAHEDRON
            print('Element Type : ',eleType, 'HEXAHEDRON')
        else:
            print('Element Type : ',eleType, 'UNKNOWN')

        # read connectivity
        arrayNodeMap = np.zeros(shape=(nElem,eleType+1), dtype='<i4')  
                                        
        nodeMap=data._retrieve_zone_node_map(z+1)
        arrayNodeMap[:,0]=eleType
        arrayNodeMap[:,1:eleType+1]=nodeMap[:]-1
        
        if dimension == '3D':
            pvMesh_temp=pv.UnstructuredGrid(arrayNodeMap.ravel(),celltypes,np.transpose([data._read_zone_var(z+1,1),data._read_zone_var(z+1,2),data._read_zone_var(z+1,3)]))
        elif dimension == '2D':
            pvMesh_temp=pv.UnstructuredGrid(arrayNodeMap.ravel(),celltypes,np.transpose([data._read_zone_var(z+1,1),data._read_zone_var(z+1,2),np.zeros(nNode)]))
        else:
            print('Wrong Dimension : ',dimension)

        if iVarRefractive!=0:
            # arraySolution[:,0]=data._read_zone_var(z+1,iVarRefractive+1) #refractive index
            pvMesh_temp['RefractiveIndex']=data._read_zone_var(z+1,iVarRefractive+1)

        else:
            print('Refractive index read problem: ',iVarRefractive)

        if iVarGradient !=0:
            
            if Solver.Dimension == '3D':
            
                pvMesh_temp['Gradient X Import']=data._read_zone_var(z+1,iVarGradient+1) 
                pvMesh_temp['Gradient Y Import']=data._read_zone_var(z+1,iVarGradient+2) 
                pvMesh_temp['Gradient Z Import']=data._read_zone_var(z+1,iVarGradient+3) 
            else:
                pvMesh_temp['Gradient X Import']=data._read_zone_var(z+1,iVarGradient+1) 
                pvMesh_temp['Gradient Y Import']=data._read_zone_var(z+1,iVarGradient+2)
                pvMesh_temp['Gradient Z Import']=np.zeros(nNode)
        else:
            
            print("No gradients in the solution ")
     
        print("\n")

        pvMesh=pvMesh.merge(grid=pvMesh_temp, merge_points=True, inplace=True, main_has_priority=False)
    
    
    if Solver.precomputedGrad:
        pvMesh["gradientX"]=pvMesh['Gradient X Import']
        pvMesh["gradientY"]=pvMesh['Gradient Y Import']
        pvMesh["gradientZ"]=pvMesh['Gradient Z Import']
        
        print('Using pre-computed Gradients   \n') 

    elif not Solver.precomputedGrad:

        pvMesh=pvMesh.compute_derivative(scalars='RefractiveIndex')
        pvMesh["gradientX"]=pvMesh["gradient"].T[0,:]
        pvMesh["gradientY"]=pvMesh["gradient"].T[1,:]
        pvMesh["gradientZ"]=pvMesh["gradient"].T[2,:]
        
        print('Computing Gradients   \n') 


    else:
        print('What about gradients?  \n')
        
    print("********************************************\n")


    return pvMesh


def create_boundary(Solver):


    # generate domain mesh
    pvMesh = pv.UnstructuredGrid()

    dimension = Solver.Dimension
    
    data = pytecio.read(Solver.CFDsolutionSurface)

    for z in range(data.numZones):

        nNode = data.zone_info[z]['IJK'][0]
        nElem = data.zone_info[z]['IJK'][1]
        eleType = data.zone_info[z]['IJK'][2]

        celltypes = np.empty(nElem, dtype=np.uint8)

        if eleType == 1 and dimension == '2D':
            celltypes[:] = pv.CellType.POLY_LINE
        if eleType == 4 and dimension == '2D':
            celltypes[:] = pv.CellType.QUAD
            # print('Element Type : ', eleType, 'QUAD')
        elif eleType == 4 and dimension == '3D':
            celltypes[:] = pv.CellType.TETRA
            # print('Element Type : ', eleType, 'TETRA')
        elif eleType == 3:
            celltypes[:] = pv.CellType.TRIANGLE
            # print('Element Type : ', eleType, 'TRIANGLE')
        elif eleType == 8:
            celltypes[:] = pv.CellType.HEXAHEDRON
            # print('Element Type : ', eleType, 'HEXAHEDRON')
        else:
            print('Element Type : ', eleType, 'UNKNOWN')

        # read connectivity
        arrayNodeMap = np.zeros(shape=(nElem, eleType+1), dtype='<i4')

        nodeMap = data._retrieve_zone_node_map(z+1)
        arrayNodeMap[:, 0] = eleType
        arrayNodeMap[:, 1:eleType+1] = nodeMap[:]-1

        if dimension == '3D':
            pvMesh_temp = pv.UnstructuredGrid(arrayNodeMap.ravel(), celltypes, np.transpose(
                [data._read_zone_var(z+1, 1), data._read_zone_var(z+1, 2), data._read_zone_var(z+1, 3)]))
        elif dimension == '2D':
            pvMesh_temp = pv.UnstructuredGrid(arrayNodeMap.ravel(), celltypes, np.transpose(
                [data._read_zone_var(z+1, 1), data._read_zone_var(z+1, 2), np.zeros(nNode)]))
        else:
            print('Wrong Dimension : ', dimension)

        pvMesh = pvMesh.merge(
            grid=pvMesh_temp, merge_points=True, inplace=True, main_has_priority=False)

    #     else:
    #         nodes = np.zeros(shape=(nNode, 3), dtype=float)
    #         nodes[:,0] = data._read_zone_var(z+1, 1)[0][0][:]
    #         nodes[:,1] = data._read_zone_var(z+1, 2)[0][0][:]
    #         nodes[:,2] = 0
            
    # if eleType != 1:
    if dimension == '2D':
        pvMesh = pvMesh.extract_feature_edges()
        pvMesh = pvMesh.extrude([0, 0, 1], capping=False)

    else:
        pvMesh=pvMesh.extract_surface()

    # if not pvMesh.is_manifold:
    #     print('Non manifold boundary surface : ', Solver.CFDsolutionSurface)

    return pvMesh
