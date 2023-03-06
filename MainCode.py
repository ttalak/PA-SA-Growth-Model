"""
The temporospatial model is a three-phased model, each following one another to create a dynamic, 
living environment where the bacterium thrives. The first step updates the environment, based on 
the available constraints, diffused from the nutrient-available boundary of the lattice, taking 
diffusion constants into consideration, up to where biomass is present. This new environment 
profile is than used to determine upper and lower boundaries of the model for each individual 
compartment of the lattice. As diffusion into a biofilm, under QS effects, would be lower in 
numbers, this setting is conveyed accordingly considering the impact of biofilm formation. As 
the second step, growth rates of biomasses in compartments as well as their production/consumption 
rates (fluxes) are calculated via Flux Balance Analysis, followed by the deduction of consumed 
metabolites through exchange reactions from the environment. The third step is the proliferation 
of biomasses after growth. The new biomasses migrate to adjacent according to the specie specific 
diffusion coefficient. After the loop is complete, a new iteration starts to repeat the cycle, 
creating a growth model of bacterium.
"""

import numpy as np
from cobra.io import load_model
from TSFBA import TSFBA,TSFBA_SA
from Diffusion import diffuser
from Proliferator import MakeProliferate


model_SA = load_model("MODEL1507180072")
model_PA = load_model('MODEL1507180020')
# About the model:
# unit: mmol / g*h(3600s)
# area: m^2
# length: m
# time: s

model_SA.solver = 'gurobi'
model_PA.solver = 'gurobi'

dim_x = 10 #length of the nutrient surface
dim_y = 20 #height of the nutrient surface

dt = 1 #hour

#Growth Environment (Bloody scar surface at the bottom)
Envir_PA = np.zeros((dim_x,dim_y,len(model_PA.exchanges),dt))
Envir_SA = np.zeros((dim_x,dim_y,len(model_SA.exchanges),dt))



#Initial biomass on the scar surface
Biomass_SA = np.zeros((dim_x,dim_y,dt))

Biomass_SA[8,9] = 3
Biomass_SA[8,10] = 3
Biomass_SA[8,11] = 3
Biomass_SA[7,10] = 3

Biomass_PA = np.zeros((dim_x,dim_y,dt))

Biomass_PA[8,9] = 3
Biomass_PA[8,10] = 3
Biomass_PA[8,11] = 3
Biomass_PA[7,10] = 3



obj_value_PA = np.zeros(np.shape(Biomass_PA))
obj_value_SA = np.zeros(np.shape(Biomass_SA))


#From surface to environment fluxes

fluxes_no_QS_PA = np.ones(len(model_PA.exchanges))*6.7e-6
fluxes_no_QS_SA = np.ones(len(model_SA.exchanges))*6.7e-6
#Glucose
# fluxes_no_QS[46] = 6.4e-6 #Glucose from blood to skin rate

# Blood Concentrations of Exchanges
# https://hmdb.ca/
# The Human Metabolome Database
Blood_Conc_PA = (0,0.5,2.447,0,0,0,41.9,1,0.64,0,
              8.9,3,444.61,427.2,0,2.1,129.5,60.8,32.9,0,
              0,0,0.019,0,14.5,0,(21000+31000)/2,(0.008+0.00067)/2,43,6.4,
              (12+22)/2,212,109,0,0,40,9766,16.7,0,0,
              0,48,1.5,0,0,0,11100,3.295,657.9,(0.2+0.014)/2, #5400 glc is normal, 11100 is Diabetes mellitus 
              65.19,0.8,329.9,82,53,30,0,0,55000000,1064.643,
              0,(120+40)/2,11.02,77.5,4200,1510,160,198,0,12,
              15,64,0,33.4,833,0.17292,34,0,0,0,
              0,25,0,(14500+13500)/2,127,0.64,32,41.8,6960,467.524, #oxygen is 6960 (88)
              0,93.8,5,0,0,57,820.4,((1.12+1.45)/2)*1000,0,198.3,
              10,0.12,63,2.3,2.28,159.8,490,10.3,23.5,162,
              0.096,127.7,0,78.4,72,2.1,4530,233,0.7,(1797.5+3088.7)/2,
              13.5)
# Dictionary of Concentrations
# exch_list = [model.exchanges[i].id for i in range(len(model.exchanges))]
# exc = [exch_list[i].strip('EX_') for i in range(len(exch_list))]
# exc2 = [exc[i].strip('_e') for i in range(len(exch_list))]

Blood_Conc_SA = [0, 0, 1510, 0, 12, 14000.0, 8.9, 41.9, 0, 0, 0, 0, 4200, 
              43, 14.5, 833, 26000.0, 0.004335, 444.61, 10.2, 0, 0, 0, 
              29.1, 0.10700000000000001, 0, 2.3, 2.28, 0, 40, 16.7, 53, 
              329.9, 7.9, 30.0, 55000000, 0.00518, 9766, 427.2, 129.5, 
              93.8, 65.19, 657.9, 77.5, 160, 198, 198.3, 233, 0, 0, 
              0.17292, 0, 35.0, 0.64, 0, 1000000, 41.8, 32, 1285.0000000000002, 6960, 0, 0, 48, #Oxygen is [59] , Nicotinate [55]
              11100, 0, 34, 64, 0, 0, 0, 0, 0.12, 0, 4.8, 
              0, 10.3, 1.8, 0, 0, 0, 0, 0.096, 1000000, 2.1, 60.8, 0.7, 13.5] #Thiosulfate is available [82]


model_SA.exchanges.EX_Thiosulfate_e.upper_bound = 10
# model_SA.exchanges.EX_Nicotinate_e.upper_bound = 10 #id 56
model_SA.exchanges.EX_H2O_e.lower_bound = - 1000 


LBs_SA = [round(Blood_Conc_SA[i]/1000,4) for i in range(len(Blood_Conc_SA))] # umol to mmol

LBs_PA = [round(Blood_Conc_PA[i]/1000,4) for i in range(len(Blood_Conc_PA))]

for i in range(len(Blood_Conc_SA)):
    Envir_SA[(-1,-2),:,i] = LBs_SA[i] #bottom two lines make a blood stream
    

for i in range(len(Blood_Conc_PA)):
    Envir_PA[(-1,-2),:,i] = LBs_PA[i]*(-1) #bottom two lines make a blood stream
    

All_Fluxes_PA = np.zeros((dim_x,dim_y,int(len(model_PA.reactions)),dt))


All_Fluxes_SA = np.zeros((dim_x,dim_y,int(len(model_SA.reactions)),dt))


for i in range(120): #30 hours of growth and diffusion
    
    Biomass = Biomass_SA[:,:,i]+Biomass_PA[:,:,i]
    Biomass_to_Ones = np.where(Biomass[:,:]> 0.0001,1,0)
    max_iter = int(sum(Biomass_to_Ones[:,Biomass_to_Ones.sum(axis=0).argmax()]))
    
    
    # for k in range(dim_x):
    #     for l in range(dim_y):
    #         if Biomass_to_Ones[k,l] != 1:
    Envir_PA[:,:,88,i] = 1000
    Envir_SA[:,:,59,i] = 1000
    
    for j in range(len(Blood_Conc_SA)):
        Envir_SA[(-1,-2),:,j,-1] = LBs_SA[j]
    
    
    for j in range(len(Blood_Conc_PA)):
        Envir_PA[(-1,-2),:,j,-1] = LBs_PA[j]

    
    Envir_SA = diffuser(max_iter,Envir_SA,fluxes_no_QS_SA)
    Envir_PA = diffuser(max_iter,Envir_PA,fluxes_no_QS_PA)
    
        
    [Biomass_SA,Envir_SA,obj_value_SA,All_Fluxes_SA] = TSFBA_SA(Biomass_SA,Envir_SA,model_SA,obj_value_SA,All_Fluxes_SA,Biomass_PA)
    [Biomass_PA,Envir_PA,obj_value_PA,All_Fluxes_PA] = TSFBA(Biomass_PA,Envir_PA,model_PA,obj_value_PA,All_Fluxes_PA,Biomass_SA)
    
    Biomass_SA[:,:,-1] = MakeProliferate(Biomass_SA[:,:,-1],0.00003)
    Biomass_PA[:,:,-1] = MakeProliferate(Biomass_PA[:,:,-1],0.00003)

