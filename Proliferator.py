import numpy as np

def MakeProliferate(BiomassIn, mu):
    fluxes = mu # cm^2/s to m^2/hour
    
    dt = 1 # Diffusion at the end of an hour
    
    dim_x,dim_y = np.shape(BiomassIn)[0],np.shape(BiomassIn)[1]
    
    dx,dy = 1/dim_x,1/dim_y #0.01m scar section
    
    
    u = np.zeros((dim_x,dim_y))
    
    
    # Biomass_to_Ones = np.where(Biomass[:,:,-1]> 0.0001,1,0)
    
    # max_iter = int(sum(Biomass_to_Ones[:,Biomass_to_Ones.sum(axis=0).argmax()]))
    
    
    u = BiomassIn
    max_iter = 3
    for iter in range(max_iter): #60 minutes
        u0 = u
        for i in range(1, dim_x-1):
            for j in range(1, dim_y-1):
                uxx = (u0[i+1,j] - u0[i,j]) / dx**2
                uyy = (u0[i,j+1] - 2*u0[i,j]+ u0[i,j-1]) / dy**2
                u[i,j] = round(u0[i,j] + dt/max_iter * fluxes * (uxx + uyy),4)
    BiomassIn = u

    return BiomassIn
