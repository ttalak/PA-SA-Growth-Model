import numpy as np

def TSFBA(Biomass,Envir,model,obj_value,All_Fluxes,Biomass2):
    
    row,col = np.where(Biomass[:,:,-1]>0.001)[0],np.where(Biomass[:,:,-1]>0.001)[1]
    
    Biomass_new = np.zeros((np.shape(Biomass)[0],np.shape(Biomass)[1],np.shape(Biomass)[2]+1))
    Envir_new = np.zeros((np.shape(Envir)[0],np.shape(Envir)[1],np.shape(Envir)[2],np.shape(Envir)[3]+1))
    
    Biomass_new[:,:,:-1] = Biomass
    Envir_new[:,:,:,:-1] = Envir
    
    Envir_new[:,:,:,-1] = Envir[:,:,:,0] #to set nutritious surface values
    
    
    All_Fluxes_new = np.zeros((np.shape(All_Fluxes)[0],np.shape(All_Fluxes)[1],np.shape(All_Fluxes)[2],np.shape(All_Fluxes)[3]+1))
    All_Fluxes_new[:,:,:,:-1] = All_Fluxes
    
    obj_value_new = np.zeros((np.shape(Biomass)[0],np.shape(Biomass)[1],np.shape(Biomass)[2]+1))
    obj_value_new[:,:,:-1] = obj_value
    
    for i in range(len(row)):
        
        [Envir_new[row[i],col[i],:,-1],Biomass_new[row[i],col[i],-1],obj_value_new[row[i],col[i],-1],
         All_Fluxes_new[row[i],col[i],:,-1]] = FBA(Biomass[row[i],col[i],-1],Envir[row[i],col[i],:,-1],
                                                   model,obj_value[row[i],col[i],-1],All_Fluxes[row[i],col[i],:,-1],
                                                   Biomass2[row[i],col[i],-1])
    
    return Biomass_new, Envir_new, obj_value_new, All_Fluxes_new

def FBA(Biomass,Envir,model,obj_value,All_Fluxes,Biomass2):
    if Biomass + Biomass2 <200:
        dt = 1 #1 for growth per hour
        
        DeathRate = 0.8 #to be set later, no death for now #it used to be 1
    
        # Setting exchanges according to what is available in the environment
        if Biomass > 20:
            for i in range(len(model.exchanges)):
                model.exchanges[i].lower_bound = Envir[i]*(-1/4)*Biomass/(Biomass+Biomass2) #lower diffusion due high-QS
        else:
            for i in range(len(model.exchanges)):
                model.exchanges[i].lower_bound = Envir[i]*(-1)*Biomass/(Biomass+Biomass2)
                
        Opt_Res = model.optimize()
        
        summary = model.summary()
        
        obj_value = Opt_Res.objective_value
        
        BiomassOut = Biomass * np.exp(Opt_Res.objective_value * DeathRate * dt)  #Population growth formula
        
        #Metabolite Consumption
       
        Consumed_Mets_index = [model.exchanges.index(summary.uptake_flux.reaction[i]) for i in range(len(summary.uptake_flux))]
        
        # Reduction of consumed metabolites according to summary of FBA
        
        All_Fluxes = Opt_Res.fluxes
        
        for i in range(len(Consumed_Mets_index)):
            if summary.uptake_flux.flux[i] >= Envir[Consumed_Mets_index[i]]:
                Envir[Consumed_Mets_index] = 0
            else:
                Envir[Consumed_Mets_index] = Envir[Consumed_Mets_index] - summary.uptake_flux.flux[i]
    else:
        Envir = Envir
        BiomassOut = Biomass
        obj_value = 0
        
    return Envir , BiomassOut, obj_value, All_Fluxes


def TSFBA_SA(Biomass,Envir,model,obj_value,All_Fluxes,Biomass2):
    
    row,col = np.where(Biomass[:,:,-1]>0.001)[0],np.where(Biomass[:,:,-1]>0.001)[1]
    
    Biomass_new = np.zeros((np.shape(Biomass)[0],np.shape(Biomass)[1],np.shape(Biomass)[2]+1))
    Envir_new = np.zeros((np.shape(Envir)[0],np.shape(Envir)[1],np.shape(Envir)[2],np.shape(Envir)[3]+1))
    
    Biomass_new[:,:,:-1] = Biomass
    Envir_new[:,:,:,:-1] = Envir
    
    Envir_new[:,:,:,-1] = Envir[:,:,:,0] #to set nutritious surface values
    
    
    All_Fluxes_new = np.zeros((np.shape(All_Fluxes)[0],np.shape(All_Fluxes)[1],np.shape(All_Fluxes)[2],np.shape(All_Fluxes)[3]+1))
    All_Fluxes_new[:,:,:,:-1] = All_Fluxes
    
    obj_value_new = np.zeros((np.shape(Biomass)[0],np.shape(Biomass)[1],np.shape(Biomass)[2]+1))
    obj_value_new[:,:,:-1] = obj_value
    
    for i in range(len(row)):
        
        [Envir_new[row[i],col[i],:,-1],Biomass_new[row[i],col[i],-1],obj_value_new[row[i],col[i],-1],
        All_Fluxes_new[row[i],col[i],:,-1]] = FBA_SA(Biomass[row[i],col[i],-1],Envir[row[i],col[i],:,-1],
                                                   model,obj_value[row[i],col[i],-1],All_Fluxes[row[i],col[i],:,-1],
                                                   Biomass2[row[i],col[i],-1])
    
    return Biomass_new, Envir_new, obj_value_new, All_Fluxes_new

def FBA_SA(Biomass,Envir,model,obj_value,All_Fluxes,Biomass2):
    if Biomass + Biomass2 < 200:        
        dt = 1 #1 for growth per hour
        
        DeathRate = 0.8-(0.2*0.4*(Biomass2/(Biomass+Biomass2))) #to be set later, no death for now #it used to be 1
    
        # Setting exchanges according to what is available in the environment
        if Biomass > 20:
            for i in range(len(model.exchanges)):
                model.exchanges[i].upper_bound = Envir[i]*(1/4)*Biomass/(Biomass+Biomass2) #lower diffusion due high-QS
        else:
            for i in range(len(model.exchanges)):
                model.exchanges[i].upper_bound = Envir[i]*(1)*Biomass/(Biomass+Biomass2)
        Opt_Res = model.optimize()
        
        summary = model.summary()
        
        obj_value = Opt_Res.objective_value
        
        BiomassOut = Biomass * np.exp(Opt_Res.objective_value * DeathRate * dt)  #Population growth formula
        
        #Metabolite Consumption
       
        Consumed_Mets_index = [model.exchanges.index(summary.uptake_flux.reaction[i]) for i in range(len(summary.uptake_flux))]
        
        # Reduction of consumed metabolites according to summary of FBA
        
        All_Fluxes = Opt_Res.fluxes
        
        for i in range(len(Consumed_Mets_index)):
            if summary.uptake_flux.flux[i] >= Envir[Consumed_Mets_index[i]]:
                Envir[Consumed_Mets_index] = 0
            else:
                Envir[Consumed_Mets_index] = Envir[Consumed_Mets_index] - summary.uptake_flux.flux[i]
    else:
        Envir = Envir
        BiomassOut = Biomass
        obj_value = 0
    
    return Envir , BiomassOut, obj_value, All_Fluxes
