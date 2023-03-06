import numpy as np

def diffuser(Envir,Biomass,MetDiffCoeff,BCType,BCs):
    
    # global BCType,BCs
    
    nLoop = np.shape(Envir)[2]
    
    [row,col] = np.where(Biomass[:,:,-1]>0.001)[0],np.where(Biomass[:,:,-1]>0.001)[1]
    SpaceDiffInd=np.zeros((np.shape(Biomass)[0],np.shape(Biomass)[1]))
    
    # [SpaceDiffInd[row[i],col[i]] for i in range(len(row))] = [1]*len(row)
    
    for i in range(len(row)):
        SpaceDiffInd[row[i],col[i]]=1
    
    for i in range(nLoop):
        SpaceDiffCoffs = SpaceDiffInd * MetDiffCoeff[i] * 3600 * 1e-4 #from cm^2/s to m^2/h
        Envir[:,:,i,-1] = diffuser2(Envir[:,:,i,-1],SpaceDiffCoffs,BCType[i,:],BCs[i,:])
        # Envir(:,:,i,end) = MakeDiffuse1M2(Envir(:,:,i,end), 
        # SpaceDiffCoffs,BCType(i,:),BCs(i,:));

    return(Envir)


def diffuser2(Envir, SpaceDiffCoffs,BCType,BCs):
    # space = um
    # time = hour
    # DiffCoff = cm**2/s
    # 0 for Neumann
    # 1 for Dirichlet
    
    SpaceUnit = 1e+4 
    nx = np.shape(Envir)[0]
    ny = np.shape(Envir)[1]
    dx = dy = SpaceUnit

    dt=1;
    nt=1;
    t=0;
    Tn=Envir
    
    DiffCoffs = SpaceDiffCoffs
    
    for n in range(nt): 
        Tc = Tn
        t = t+dt
        for i in range(1,ny-1):
            for j in range(1,nx-1):
                Tn[i,j] = Tc[i,j]+DiffCoffs[i,j]*dt*(Tc[i+1,j]-2*Tc[i,j]+Tc[i-1,j])/dx**2+(Tc[i,j+1]-2*Tc[i,j]+Tc[i,j-1])/dy**2
    
    
        # BCs
        # Top
        
        if BCType[0] == 0:
            Tn[0,:]=Tn[1+1,:]
        elif BCType[0] == 1:
            Tn[0,:]=BCs[1]
        else:
            print('Wrong type of BC is given.')
        
        # Right
        
        if BCType[1] == 0:
            Tn[:,-1]=Tn[:,-2]
        elif BCType[1] == 1:
            Tn[:,-1]=BCs[2]
        else:
            print('Wrong type of BC is given.')
    
        # Left
        
        if BCType[2] == 0:
            Tn[-1,:]=Tn[-2,:]
        elif BCType[2] == 1:
            Tn[-1,:]=BCs[2]
        else:
            print('Wrong type of BC is given.')
            
        # Bottom
        
        if BCType[3] == 0:
            Tn[:,0]=Tn[:,1+1]
        elif BCType[3] == 1:
            Tn[:,0]=BCs[3]
        else:
            print('Wrong type of BC is given.')

    return(Tn)
