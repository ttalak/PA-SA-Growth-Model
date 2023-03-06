import numpy as np


def diffuser(max_iter,Envir,fluxes_no_QS):
    
    fluxes = fluxes_no_QS * 3600 * 1e-4 # cm^2/s to m^2/hour
    
    dt = 1 # Diffusion at the end of an hour
    
    dim_x,dim_y = np.shape(Envir)[0],np.shape(Envir)[1]
    
    dx,dy = 1/dim_x,1/dim_y #0.01m scar section
    
    
    u = np.zeros((dim_x,dim_y))
    

    
    for metabolite in range(len(fluxes)): # One iteration per metabolite for its flux
        u = Envir[:,:,metabolite,-1]
        if Envir[:,:,metabolite,-1].any() != 0:
            for iter in range(max_iter): #60 minutes
                u0 = u
                for i in range(1, dim_x-1):
                    for j in range(1, dim_y-1):
                        uxx = (u0[i+1,j] - u0[i,j]) / dx**2
                        uyy = (u0[i,j+1] - 2*u0[i,j] + u0[i,j-1]) / dy**2
                        u[i,j] = round(u0[i,j] + dt/max_iter * fluxes[metabolite] * (uxx + uyy),4)
            Envir[:,:,metabolite,-1] = u
        else:
            continue
        
    return Envir