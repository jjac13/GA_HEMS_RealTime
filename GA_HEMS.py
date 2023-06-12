# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################
import time
from numpy import append, array, cos, pi, arange

###############################################################################
###########################   Functions definition  ###########################

###############################################################################
##########################   House model definition  ##########################

def heat_loss(U, Area, T_in, T_out):
    return U * Area * (T_out - T_in)
    
def new_house_Temperature(T_0, Qdot_roof, Qdot_windows, Qdot_wall, mc_T, Qdot_HP = 0, Qdot_TESS = 0, Qdot_SC = 0, dt=1):
    T_new = T_0 + dt*(Qdot_roof + Qdot_windows + Qdot_wall + Qdot_HP + Qdot_TESS + Qdot_SC)/(mc_T)
    
    return T_new

def update_TESS(active, T_0, T_soil, Qdot_SC = 0, Qdot_HP = 0, Tmax = 95 + 273, Tmin = 50 + 273, mdot = 0.1, T_network = 40 + 273, Qdot_SD = 100, efficiency = 0.8, m = 4000, c = 4200, dt=1*3600):

    if active and T_0 >= Tmin:
        Qdot_TESS = mdot*c*(T_0*0.97 - T_network)
    else:
        Qdot_TESS = 0

    
    if T_0 <= Tmin:         # TESS discharged, only charge
        if T_0 <= T_soil:   # Code for heat dissipation/absoprtion through the soil needed
            T_new = T_0
            
        else:
            T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS)*dt/(m*c)
            
    elif T_0 <= Tmax:       # TESS available for charge and discharge
        
        T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS)*dt/(m*c)
        
    else:                   # TESS fully charged, only discharge
        
        T_new = T_0 + (- Qdot_SD - Qdot_TESS)*dt/(m*c)        
        
    return [T_new, Qdot_TESS*efficiency]          

def Qdot_SolarCollector(active, G, A = 6, SC_eff = 0.45, dt = 1*3600):
    
    if active:
        return A*SC_eff*G/dt
    else:
        return 0

# Enphase IQ3 https://enphase.com/download/iq-battery-3-data-sheet

def update_BESS(SoC_0, P_BESS, P_Load, P_PV = 0, SoCmax = 0.9, SoCmin = 0.2, P_BESS_max = 1.28, P_Grid_max = 0, Capacity_BESS = 3.36, charge_efficiency = 0.943, discharge_efficiency = 0.943, P_SD = 0, dt = 0.25):

    
    E_BESS_0 = Capacity_BESS*SoC_0
    
    if P_BESS > 0:
        E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
    elif P_BESS <= 0:
        E_BESS = E_BESS_0 - P_BESS*dt/charge_efficiency
    
    E_BESS = E_BESS*(1-P_SD)
    SoC_BESS = E_BESS /  Capacity_BESS
    
    return SoC_BESS
    

def BESS_perm_min(SoC, Capacity_BESS = 3.36, SoCmax = 0.9, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmax)/dt, -P_BESS_max, P_BESS_max)
    
    
def BESS_perm_max(SoC, Capacity_BESS = 3.36, SoCmin = 0.2, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmin)/dt, -P_BESS_max, P_BESS_max)
    
        
def HP_Power(active, P_in = 2.7*1000, COP = 4.1):
    if active:
        return [P_in/1000, P_in*COP]
    else:
        return [0,0]
    
def Thermal_Electrical_model(T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G, SC_active, TESS_active, HP_active, PV_curtail, P_BESS, T_0_TESS, SoC_0_BESS, dt = 0.25):

    #############################   Paremeters ################################
    
    # Convective heat transfer coefficients [W/m^2K]
    h_air_wall   = 0.9#24       # Indoor air -> walls, scaled to R-value of a C-label house
    h_wall_atm   = 0.9#34       # Walls -> atmosphere, scaled to R-value of a C-label house
    h_air_window = 25            # Indoor air -> windows
    h_window_atm = 32            # Windows -> atmosphere
    h_air_roof   = 12            # Indoor air -> roof
    h_roof_atm   = 38            # Roof -> atmosphere
    
    ## House
    
    # Air
    c_air        = 1005.4        # Specific heat of air at 273 K [J/kgK]
    airDensity   = 1.025         # Densiity of air at 293 K [kg/m^3]
    kAir         = 0.0257        # Thermal conductivity of air at 293 K [W/mK]
    
    
    # Windows (glass)
    n1_window     = 3            # Number of windows in room 1
    n2_window     = 2            # Number of windows in room 2
    n3_window     = 2            # Number of windows in room 3
    n4_window     = 1            # Number of windows in room 4
    htWindows     = 1            # Height of windows [m]
    widWindows    = 1            # Width of windows [m]
    windows_area  = (n1_window + n2_window + n3_window + n4_window) * htWindows * widWindows
    LWindow       = 0.004        # Thickness of a single window pane [m]
    LCavity       = 0.014        # Thickness of the cavity between the double glass window [m]  
    windowDensity = 2500         # Density of glass [kg/m^3]
    c_window      = 840          # Specific heat of glass [J/kgK]
    kWindow       = 0.8          # Thermal conductivity of glass [W/mK]
    U_windows = ((1/h_air_window) + (LWindow/kWindow) + (LCavity/kAir) + (LWindow/kWindow) + (1/h_window_atm))**-1
    m_windows = windowDensity * windows_area * LWindow
    
    # Walls (concrete)
    lenHouse    = 15             # House length [m]
    widHouse    = 8              # House width [m]
    htHouse     = 2.6            # House height [m]
    LWall       = 0.25           # Wall thickness [m]
    wallDensity = 2400           # Density [kg/m^3]
    c_wall      = 750            # Specific heat [J/kgK]
    kWall       = 0.14           # Thermal conductivity [W/mK]
    walls_area = 2*(lenHouse + widHouse) * htHouse - windows_area
    U_wall = ((1/h_air_wall) + (LWall /kWall) + (1/h_wall_atm))**-1
    m_walls = wallDensity * walls_area * LWall
    
    # Roof (glass fiber)
    pitRoof     = 40/180/pi      # Roof pitch (40 deg)
    LRoof       = 0.2            # Roof thickness [m]
    roofDensity = 2440           # Density of glass fiber [kg/m^3]
    c_roof      = 835            # Specific heat of glass fiber [J/kgK]
    kRoof       = 0.04           # Thermal conductivity of glass fiber [W/mK]
    roof_Area = 2 * (widHouse/(2*cos(pitRoof))*lenHouse)
    U_roof = ((1/h_air_roof) + (LRoof/kRoof) + (1/h_roof_atm))**-1
    m_roof = roofDensity * roof_Area * LRoof
    
    
    m_air = airDensity * lenHouse * widHouse * htHouse
    
    mc_T = m_air*c_air + m_roof*c_roof + m_windows*c_window + m_walls*c_wall
     
    
    ############################   Calculations ###############################    
    ################### Thermal carrier ######################
    
    ## Thermal losses
    # Roof
    Qdot_roof = heat_loss(U_roof, roof_Area, T_0_in, T_amb)
    
    # Windows
    Qdot_windows = heat_loss(U_windows, windows_area, T_0_in, T_amb)
    
    # Walls
    Qdot_wall = heat_loss(U_wall, walls_area, T_0_in, T_amb)
    
#    Qdot_losses = Qdot_roof + Qdot_windows + Qdot_wall


    # Solar Collector
    Qdot_SC = Qdot_SolarCollector(SC_active, G)
    
    # TESS     
    Qdot_SC_TESS = Qdot_SolarCollector(not SC_active, G)
    TESS_state = update_TESS(TESS_active, T_0_TESS, T_soil = T_soil, Qdot_SC = Qdot_SC_TESS, dt = dt)
    T_TESS = TESS_state[0]
    Qdot_TESS =  TESS_state[1]

    
    # HP        
    HP_state = HP_Power(HP_active)
    P_HP = HP_state[0]
    Qdot_HP = HP_state[1]        
    
    # Thermal demand
    T_in = new_house_Temperature(T_0_in, Qdot_roof, Qdot_windows, Qdot_wall, mc_T, Qdot_HP = Qdot_HP, Qdot_TESS = Qdot_TESS, Qdot_SC = Qdot_SC, dt = dt)  
    Qdot_Losses = -(Qdot_roof + Qdot_windows + Qdot_wall)


    ################### Electric carrier ######################

    # PV
    P_PV = P_PV_av*PV_curtail

    # BESS
#   [SoC_BESS_state, P_Grid_state, P_BESS_state] = update_BESS(SoC_BESS, P_Load + P_HP/1000, P_PV)
    SoC_BESS = update_BESS(SoC_0_BESS, P_BESS, P_Load + P_HP/1000)
    
    
    # Grid balance
    
    P_Grid = P_Load + P_HP - P_PV - P_BESS
    
    return [Qdot_SC, Qdot_SC_TESS, Qdot_TESS, Qdot_HP, Qdot_Losses, T_TESS, T_in, P_HP, P_PV, P_BESS, P_Grid, SoC_BESS]

###############################################################################
#############################   GA EMS definition  ############################

## The function ---------------------------------------------------------------

def chromosome(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS):
# Thermal_Electrical_model(T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G, SC_active, TESS_active, HP_active, PV_curtail, P_BESS, T_0_TESS, SoC_0_BESS, dt = 0.25):    
    from numpy import append, random

# For discrete PV curtailment.    
#    chromosome = array([random.randint(0,2), random.randint(0,2), random.randint(0,2), random.randint(0,2), round(random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1])),1)])
    
# For continous PV curtailment.
    
    chromosome = array([random.randint(0,2), random.randint(0,2), random.randint(0,2), round(random.uniform(0, 1),1), round(random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1])),1)])
    
    chromosome = append(chromosome, cost_function(chromosome, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS))
    
    return chromosome

## The function ---------------------------------------------------------------
    
def initial_population(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, individuals = 200):
    
    population = array([chromosome(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS) for i in range(individuals)])
    
    return population

## The function ---------------------------------------------------------------

## Considerations:
##  - Consider the size of the population in case the elite_rate and 
##    non_elite_rate variables are modified, as their sum should result into an
##    even number when multiplied by the size of the population.


def create_new_population(population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, elite_rate = 0.05, non_elite_rate = 0.15):    
    from numpy import append, random
    
    # Elite selection
    population = population[population[:,4].argsort()]
    elite = population[int((1-elite_rate)*len(population)): len(population), :]
    
    # Non-elite selection
    non_elite = population[0:int((1-elite_rate)*len(population)), :]
    non_elite = non_elite[random.choice(non_elite.shape[0], int(len(population)*non_elite_rate), replace=False), :]
    
    # Mating
    mating_population = append(elite, non_elite, axis=0)
    mating_population = mating_population[random.choice(list(range(0,len(mating_population))),len(mating_population), replace=False)]    
    new_generation = create_new_generation(mating_population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS)
    
    # Inclusion of new generation
    
    new_population = append(population, new_generation, axis=0)
    
    # Selecting survivors
    
    new_population = new_population[random.choice(new_population.shape[0], len(population), replace=False), :]
    
    
    return new_population

## The function ---------------------------------------------------------------
    
## Considerations:
##  - The number of genes to be crossed has to be at least 1, with a maximum of
##    1 less than the number of genes per chromosome.

def create_new_generation(mating_population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, rand_gene_number = True):
    from numpy import append, ceil, random, reshape
    
    new_generation = array([])   
    number_of_genes = len(mating_population[0,:])
    
    if rand_gene_number:        
        cross_genes = random.choice(list(range(0,number_of_genes)),random.randint(1,number_of_genes), replace=False)
    else:
        cross_genes = random.choice(list(range(0,number_of_genes)),int(ceil(number_of_genes/2)), replace=False)

    for i in range(int(len(mating_population)/2)):         
        
        chromosome_1 = mating_population[i]
        chromosome_2 = mating_population[int(len(mating_population)/2)+i]
        
        new_chromosome_1 = array([i for i in chromosome_1])
        new_chromosome_2 = array([i for i in chromosome_2])
        
        for i in cross_genes:            
            new_chromosome_1[i] = chromosome_2[i]
            new_chromosome_2[i] = chromosome_1[i]
            
        
        new_chromosome_1[-1] = cost_function(new_chromosome_1, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS)
        new_chromosome_2[-1] = cost_function(new_chromosome_2, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS)
        
        new_chromosomes = append(new_chromosome_1, new_chromosome_2, axis=0)
        new_generation = append(new_generation, new_chromosomes, axis=0)
    
    
    new_generation = reshape(new_generation, (-1,number_of_genes))
    return new_generation
    

def best_individual(population):    
    from numpy import argmin
#    from statistics import mean
    
#    print('Min costs', min([i[-1] for i in population]))
#    a = population[argmin(population[:,-1])]
#    
#    if min([i[-1] for i in population]) == a[-1]:
#        print('True')

    return population[argmin(population[:,-1])]
    
    
## The function ---------------------------------------------------------------

def cost_function(chromosome, energy_costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, beta = 6, theta_E = 1/26, theta_T = 1/1, theta_CO2 = 1/4):
#    from numpy import random
    CO2_costs = [0.071, 0.00155, 0.1, 0.041, 0.1916, 0.325]
    
    model_state = Thermal_Electrical_model(T_amb = T_amb, T_0_in = T_0_in, T_set = T_set, T_soil = T_amb, P_PV_av = P_PV_av, P_Load = P_Load, G = G, SC_active = chromosome[0], TESS_active = chromosome[1], HP_active = chromosome[2], PV_curtail = chromosome[3], P_BESS = chromosome[4], T_0_TESS = T_0_TESS, SoC_0_BESS = SoC_0_BESS, dt = dt)
    
#    Qdot_SC = append(Qdot_SC, model_state[0])
#    Qdot_SC_TESS = append(Qdot_SC_TESS, model_state[1])
#    Qdot_TESS = append(Qdot_TESS, model_state[2])
#    Qdot_HP = append(Qdot_HP, model_state[3])
#    Qdot_Losses = append(Qdot_Losses, model_state[4])
#    T_TESS = append(T_TESS, model_state[5])
#    T_in = append(T_in, model_state[6])
#    P_HP = append(P_HP, model_state[7])
#    P_PV = append(P_PV, model_state[8])
#    P_BESS = append(P_BESS, model_state[9])
#    P_Grid = append(P_Grid, model_state[10])
#    SoC_BESS = append(SoC_BESS, model_state[11])
    
    if model_state[9] < 0:
    
        Variable_powers = [(model_state[0] + model_state[1])/1000, model_state[2]/1000, model_state[3]/1000, model_state[8], model_state[9], model_state[10]]
    
    else:
        Variable_powers = [model_state[0]/1000, model_state[2]/1000, model_state[3]/1000, model_state[8], 0, model_state[10]]
        
#    
#    electric_cost = sum([Variable_powers[i]*costs[i] for i in range(len(costs))])
#    thermal_cost = beta*abs(T_set - T_0_in)
#    
    cost = theta_E*sum([Variable_powers[i]*energy_costs[i] for i in range(len(energy_costs))]) + theta_T*beta*abs(T_set - model_state[6]) + theta_CO2*sum([Variable_powers[i]*CO2_costs[i] for i in range(len(CO2_costs))])
#    print('Electric cost: ', electric_cost, ', Thermal cost: ', thermal_cost, ', Total cost: ', cost)
    
    
    return cost 
    
## The function ---------------------------------------------------------------

def GA_Optimization(T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, costs = [0.15, 1.4, 0.11, 0.18, 0.13, 0.483], consecutive_generations = 5):

    population = initial_population(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS)
    best_candidate = best_individual(population)
    
    generation = 0
    number_of_candidates=1
    iteration = 0
    
    while iteration < consecutive_generations:
        new_population = create_new_population(population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS)
        new_candidate = best_individual(new_population)
        
        if new_candidate[-1] > best_candidate[-1]:
            best_candidate = [i for i in new_candidate]
            iteration = 0
            number_of_candidates+=1
#            print('-----New best candidate found-----')
#        
#        print('Current generation: ', generation)
#        print('Current iteration: ', iteration)
        generation+=1  
        iteration+=1
    
#    print('The best candidate from the initial population is: ', best_individual(population))
#    print('The best candidate from the last generation is: ', best_candidate)   
#    print('The number of generations was: ', generation)
    
    return best_candidate    


###############################################################################
################################   Simulations ################################

import matplotlib.pyplot as plt
import csvreader
from numpy import savetxt


start = time.time()    # The timer is initializad.

start_day = 0
end_day = 365
# [SC, TESS, HP, PV, BESS, Grid]
energy_costs = [0.15, 1.4, 0.11, 0.18, 0.13, 0.483]
CO2_costs = [0.071, 0.00155, 0.1, 0.041, 0.1916, 0.325]


# PV
n_modules = 10
module_power_ref = 0.315
module_power = 0.400   

CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
CSVDataP_Load = csvreader.read_data(csv='Load_Profile_15min.csv', address='', delim=',')
CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
CSVDataPV.data2array()
CSVDataTamb.data2array()
CSVDataP_Load.data2array()
CSVDataRad.data2array()
P_PV_av = [i[0]*n_modules*module_power/module_power_ref/1000 for i in CSVDataPV.ar]
T_amb = [i[0]+273 for i in CSVDataTamb.ar[start_day*24*4:end_day*24*4]]
P_Load = [i for i in CSVDataP_Load.ar[0]]
a = arange(0,len(CSVDataRad.ar),15)
G = array([CSVDataRad.ar[i][0] for i in a])


# Initial conditions
t = array([0])                      # In s
dt = 60*15                          # In s
t_final = int((end_day - start_day)*24*3600/dt)         # In s

# TESS
T_TESS = array([75 + 273])          # In K
Qdot_TESS = array([0, 0])              # In W
TESS_active = False

# HP
HP_active = False
HP_status = HP_Power(HP_active)
P_HP = array([HP_status[0],HP_status[0],HP_status[0]])       # In kW
Qdot_HP = array([HP_status[1],HP_status[1],HP_status[1]])    # In kW


# Thermal demand
T_0 = 20 + 273
T_in = array([T_0, T_0, T_0])          # In K
T_set_day = [17+273]*int((6-0)*4) + [20+273]*int((22-6)*4)+ [17+273]*int((24-22)*4)
T_set = array(T_set_day*(end_day-start_day))
Qdot_Losses = array([0, 0, 0])
Qdot_SC_TESS = array([0])

# Solar Collectors
SC_active = False
Qdot_SC = array([0, 0, 0])                # In W

# PV

P_PV = array([0]) 

# BESS
SoC_BESS = array([.50])                     # In %
P_BESS = array([0])                         # In kW
SoCmax = 0.9
SoCmin = 0.2
charge_efficiency = 0.943
discharge_efficiency = 0.943


# Grid
P_Grid = array([0])                         # In kW


i = 0
day = 0

for step in range(t_final-1):
    
#    from numpy import random
#    SC_active = random.randint(0,2)
#    TESS_active = random.randint(0,2)
#    HP_active = random.randint(0,2)
#    PV_curtail = random.uniform(0,1)
#    P_BESS_active = random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1]))
    
    
    best_candidate = GA_Optimization(T_amb = T_amb[step], T_0_in = T_in[-1], T_set = T_set[step], T_soil = T_amb[step], P_PV_av = P_PV_av[step], P_Load = P_Load[step], G = G[step],  T_0_TESS = T_TESS[-1], SoC_0_BESS = SoC_BESS[-1], costs = energy_costs, consecutive_generations = 50)

    SC_active = best_candidate[0]
    TESS_active = best_candidate[1]
    HP_active = best_candidate[2]
    PV_curtail = best_candidate[3]
    P_BESS_active = best_candidate[4]

    model_state = Thermal_Electrical_model(T_amb = T_amb[step], T_0_in = T_in[-1], T_set = T_set[step], T_soil = T_amb[step], P_PV_av = P_PV_av[step], P_Load = P_Load[step], G = G[step], SC_active = SC_active, TESS_active = TESS_active, HP_active = HP_active, PV_curtail = PV_curtail, P_BESS = P_BESS_active, T_0_TESS = T_TESS[-1], SoC_0_BESS = SoC_BESS[-1], dt = dt)


#    P_BESS_active = random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1]))
#    SoC_BESS_state = update_BESS(SoC_BESS[-1], P_BESS[-1], P_Load[-1] + P_HP[-1]/1000)
#    SoC_BESS = append(SoC_BESS[-1], SoC_BESS_state)    
    
#    return [Qdot_SC, Qdot_TESS, Qdot_HP, Qdot_losses, T_TESS, T_in, P_HP, P_PV, P_BESS, P_Grid, SoC_BESS]
    
    Qdot_SC = append(Qdot_SC, model_state[0])
    Qdot_SC_TESS = append(Qdot_SC_TESS, model_state[1])
    Qdot_TESS = append(Qdot_TESS, model_state[2])
    Qdot_HP = append(Qdot_HP, model_state[3])
    Qdot_Losses = append(Qdot_Losses, model_state[4])
    T_TESS = append(T_TESS, model_state[5])
    T_in = append(T_in, model_state[6])
    P_HP = append(P_HP, model_state[7])
    P_PV = append(P_PV, model_state[8])
    P_BESS = append(P_BESS, model_state[9])
    P_Grid = append(P_Grid, model_state[10])
    SoC_BESS = append(SoC_BESS, model_state[11])

    t = append(t, dt*step/3600)
    
#    SC_active = random.randint(0,2)
#    TESS_active = random.randint(0,2)
#    HP_active = random.randint(0,2)
#    PV_curtail = random.uniform(0,1)
#    P_BESS_active = random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1]))
    
    i += 1
    if i == 96:
        i = 0
        day += 1
        print('Current day: ', day, '. Current temperature: ', T_in[-1]-273)
    
    


end = time.time()    # The timer is initializad.
totalelapsed = end - start  # The total time is calculated.
###################################   CSVs  ###################################

#savetxt('T_in_15min.csv', T_in, delimiter =", ", fmt ='% s')
#savetxt('Thermal_load_15min.csv', Qdot_Losses, delimiter =", ", fmt ='% s')


##################################   Plots  ###################################

plt.rcParams.update({
#    "text.usetex": True,
    "font.family": "Times New Roman"
})


plt.figure(1)
plt.plot(t, [i-273 for i in T_TESS])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Temperature in the TESS, $T_{TESS}$, [°C]')
#plt.title('Temperature in the TESS')
plt.show()

plt.figure(2)
plt.plot(t, [i/1000 for i in Qdot_TESS[1:]])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal power of the TESS, $\dot{Q}_{TESS}$, [kW]')
#plt.title('Thermal power provided by the TESS')
plt.show()

plt.figure(3)
plt.plot(t, [i-273 for i in T_amb], 'r', label='Ambient temperature, $T_{amb}$')
plt.plot(t, [i-273 for i in T_in[2:]], 'b', label='Temperature inside the house, $T_{in}$')
plt.plot(t, [i-273 for i in T_set], 'g', label='Setpoint temperature inside the house, $T_{set}$')
plt.legend(loc='lower center')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Temperature [°C]')
#plt.title('Temperature inside the house')
plt.show()    
#
#plt.figure(4)
#plt.plot(t, [i/1000 for i in Qdot_Losses[2:]])
#plt.grid()
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
##plt.xlabel('Time [h]')
#plt.ylim([0, 1.8])
#plt.ylabel('Thermal losses, $\dot{Q}_{L}$, [kW]')
##plt.title('Thermal losses in the house')
#plt.show()    
#
plt.figure(5)
plt.plot(t, [i/1000 for i in Qdot_HP[2:]])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal power of the HP, $\dot{Q}_{HP}$, [kW]')
#plt.title('Heat provided by the HP')
plt.show()

plt.figure(6)
plt.plot(t, [i/1000 for i in Qdot_SC[2:]], 'b', label='To the thermal load')
plt.plot(t, [i/1000 for i in Qdot_SC_TESS], 'r', label='To the TESS')
plt.legend(loc='lower center')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Thermal power of the SC, $\dot{Q}_{SC}$, [kW]')
#plt.title('Thermal power provided by the SC')
plt.show()
#
plt.figure(7)
plt.plot(t, P_Load[0:end_day*4*24])
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Electric load, $P_{L}$, [kW]')
#plt.title('Electric load')
plt.show()  
#
#plt.figure(8)
#plt.plot(t, [i/1000 for i in P_HP[2:]])
#plt.grid()
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
##plt.xlabel('Time [h]')
#plt.ylabel('Power of the HP, $P_{HP}$, [kW]')
##plt.title('Power consumed by the HP')
#plt.show()

plt.figure(9)
plt.plot(t, P_PV, 'r', label='PV power curtailed')
plt.plot(t, P_PV_av[1:end_day*4*24+1], 'b', label='PV power available')
plt.legend(loc='lower center')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('PV power load, $P_{PV}$, [kW]')
#plt.title('Electric load')
plt.show()  

plt.figure(10)
plt.plot(t, SoC_BESS)
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('State-of-charge of the BESS, $SoC_{BESS}$ [%]')
#plt.title('SoC of the BESS')
plt.show()

plt.figure(11)
plt.plot(t, P_BESS, 'b', label='Power delivered by the BESS')
#plt.legend(loc='center right')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Power of the BESS, $P_{BESS}$, [kW]')
#plt.title('Power of the BESS')
plt.show()

plt.figure(12)
plt.plot(t, P_Grid, 'b', label='Power consumed from the grid')
#plt.legend(loc='center right')
plt.grid()
plt.xlim([0, end_day*24])
plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
#plt.xlabel('Time [h]')
plt.ylabel('Power from the grid, $P_{G}$, [kW]')
#plt.title('Power from the grid')
plt.show()