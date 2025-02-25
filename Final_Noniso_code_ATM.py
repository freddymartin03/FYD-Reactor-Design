#%% ==================== Importing necessary libraries ==================== 

import numpy as np
import matplotlib.pyplot as plt

#%% ==================== Initial Parameters ====================

# Kinetic and model parameters
k0_best = 11.9 #( min^-1) could also be 156
Ea_best = 40.4 # ( kJ/mol ) could also be 63.67
delta2_best = 0.1
delta3_best = 0.00000005
lambda_decay_best = 1e10

# Fixed Constants
R_val = 0.008314 # ( Kj/mol*K )

T_initial = 25 + 273.15 # ( K )
T_reaction = 900 + 273.15 # ( K )
T_target = 1500 + 273.15  #Cap for temperature (1500 °C)

# Species data: molar mass (kg/kmol) and initial mass (kg)

species = {
       "MnO2":  {"molar_mass": 86.94,   "initial_mass": 2380},
       "Mn2O3": {"molar_mass": 157.88,  "initial_mass": 2740},
       "Mn3O4": {"molar_mass": 228.82,  "initial_mass": 0},
       "MnO":   {"molar_mass": 70.94,   "initial_mass": 2380},
       "Fe2O3": {"molar_mass": 159.70,  "initial_mass": 980},
       "Fe3O4": {"molar_mass": 231.55,  "initial_mass": 0},
       "FeO":   {"molar_mass": 71.85,   "initial_mass": 0},
       "Fe":    {"molar_mass": 55.85,   "initial_mass": 0},
       "SiO2":  {"molar_mass": 60.0843, "initial_mass": 790},
       "Al2O3": {"molar_mass": 101.961, "initial_mass": 30},
       "CaO":   {"molar_mass": 56.0774, "initial_mass": 1830},
       "MgO":   {"molar_mass": 40.3044, "initial_mass": 480},
       "C":     {"molar_mass": 12.011,  "initial_mass": 10}
    }


# Define initial and final weight fractions (for the ore)
w0 = 1.0    # initial weight fraction
wh = 0.913  # final weight fraction (at complete conversion)

# Experimental Data for 900°C
experimental_time = np.array([0, 6, 10, 15, 20, 25, 30, 35, 40, 65, 92.5, 120])
experimental_alpha = np.array([0, 0.75, 0.78, 0.85, 0.87, 0.9, 0.914, 0.925, 0.94, 0.96, 1, 1])



# Set a time frame in minutes
Heating_time = 1 # for heating stage
reduction_time = np.linspace(0, 120, 121) # For reduction

# Reactor Geometry parameters
rho_ore = 4462.17 # ( Kg/m^3 )
porosity = 0.253 # Also known as void fraction
design_factor = 0.75


# Heat transfer coefficients
h_conv = 25
A_solid = 34.8

# Reaction Enthalpies ( Kj / mol )
Delta_Hr = np.array([-163.7, -135.1, -16.6, -50.0, 16.6, -13.8])

# --- Gas stream parameters based on H₂ ---
H2_total_needed_full = 62.65   # kmol (total H₂ needed at full conversion)
H2_inlet_rate = 3.16  / 3       # kmol/min (assumed constant inlet rate)
MW_H2 = 2.016    # kg/kmol
Cp_H2 = 15400   # J/(kg·K)
T_in_H2 = T_reaction
P = 5
#%% ==================== Kinetic Modelling ====================

# Kinetic rate constant using the Arrhenius equation
def rate_constant(k0, Ea, R, T):
    return k0 * np.exp(-Ea / (R * T))

# Dynamic Avrami exponent as a function of conversion (alpha)
def n_alpha(alpha, delta2, delta3, lambda_decay):
    """
    Returns a dynamic Avrami exponent based on the current conversion alpha.
    Conversion thresholds are defined as:
      - Region 1 (reaction-controlled): α < 0.54
      - Region 2 (transition): 0.54 <= α < 0.74
      - Region 3 (diffusion-controlled): 0.74 <= α < 0.99
      - Region 4 (final phase): α >= 0.99
    """
    # Define conversion thresholds
    alpha1 = 0.54
    alpha2 = 0.74
    alpha3 = 0.99
   
    # Define regions based on conversion
    cond1 = alpha < alpha1
    cond2 = (alpha >= alpha1) & (alpha < alpha2)
    cond3 = (alpha >= alpha2) & (alpha < alpha3)
    cond4 = alpha >= alpha3
   
    # Region 1: Reaction-controlled
    n1 = lambda a: 1
    # Region 2: Transition from reaction to diffusion control
   
    n2 = lambda a: 1 - 0.5 * (1 - np.exp(-(a - alpha1) / delta2))
   
    # Region 3: Diffusion-controlled
    n3 = lambda a: 0.5 + (4.0 - 0.5) * (1 - (1 / (1 + delta3 * (a - alpha2))))
   
    # Region 4: Final phase
    n4 = lambda a: 0.5 + (4.0 - 0.5) * np.exp(-lambda_decay * (a - alpha3))
   
    return np.piecewise(alpha, [cond1, cond2, cond3, cond4], [n1, n2, n3, n4])

# Avrami Model with a dynamic exponent.
def update_kinetics(T_current, alpha_prev, cum_time, dt_min, params):
    """
    Update the conversion α using a simple Avrami-like model.
   
    Inputs:
      T_current: current temperature (K)
      alpha_prev: conversion from previous time step
      cum_time: cumulative time in minutes (for kinetics)
      dt_min: time step in minutes
      params: dictionary containing kinetic parameters:
          - k0, Ea, delta2, delta3, lambda_decay
         
    Outputs:
      alpha_new: updated conversion (clipped to <=1)
      d_alpha_dt: rate of conversion in 1/s (using dt in seconds)
      k_current: current rate constant based on T_current
      cum_time_new: updated cumulative time (minutes)
    """
    # Update cumulative time
    cum_time_new = cum_time + dt_min

    # Compute current rate constant using Arrhenius
    k_current = rate_constant(params['k0'], params['Ea'], R_val, T_current)
   
    # Determine the dynamic exponent using previous conversion
    n_current = n_alpha(alpha_prev, params['delta2'], params['delta3'], params['lambda_decay'])
   
    # Compute new conversion using a simplified Avrami model:
    #   α = 1 - exp[-(k_current * cum_time_new)^(n_current)]
    alpha_new = 1 - np.exp( - (k_current * cum_time_new)**(n_current) )
    if alpha_new > 1.0:
        alpha_new = 1.0

    # Compute dα/dt in 1/s (dt_min converted to seconds)
    dt_sec = dt_min * 60.0
    d_alpha_dt = (alpha_new - alpha_prev) / dt_sec

    return alpha_new, d_alpha_dt, k_current, cum_time_new


#%% ==================== Mass Transfer Modelling ====================

def compute_total_feed_mass(species_dict):
    """
    Given a dictionary of species with 'initial_mass' for each,
    return the total feed mass (in kg).
    """
    total_mass = 0.0
    for sp_data in species_dict.values():
        total_mass += sp_data["initial_mass"]
    return total_mass

def compute_xi_for_alpha(a, species, w0, wh):
    """
    Partition the overall conversion 'a' into six reaction stages with thresholds:
      Stage 1: 0      to 0.219 (MnO2 → Mn2O3)
      Stage 2: 0.219  to 0.382 (Mn2O3 → Mn3O4)
      Stage 3: 0.382  to 0.709 (Mn3O4 → MnO)
      Stage 4: 0.709  to 0.741 (Fe2O3 → Fe3O4)
      Stage 5: 0.741  to 0.806 (Fe3O4 → FeO)
      Stage 6: 0.806  to 1     (FeO → Fe)
   
    Each stage’s extent is computed as:
        ξ = (Δα) * (w0 – wh) * im / ΔM.
    """
   
    total_feed_mass = compute_total_feed_mass(species)
    xi = np.zeros(6)
    conv_factor = (w0 - wh)
    a1, a2, a3, a4, a5 = 0.219, 0.382, 0.709, 0.741, 0.806

    # Stage 1: MnO2 → Mn2O3
    delta1 = min(a, a1)
    xi[0] = delta1 * conv_factor * total_feed_mass / (2 * species["MnO2"]["molar_mass"] - species["Mn2O3"]["molar_mass"])
   
    # Stage 2: Mn2O3 → Mn3O4
    if a > a1:
        delta2_val = min(a, a2) - a1
        xi[1] = delta2_val * conv_factor * total_feed_mass / (3 * species["Mn2O3"]["molar_mass"] - 2 * species["Mn3O4"]["molar_mass"])
   
    # Stage 3: Mn3O4 → MnO
    if a > a2:
        delta3_val = min(a, a3) - a2
        xi[2] = delta3_val * conv_factor * total_feed_mass / (species["Mn3O4"]["molar_mass"] - 3 * species["MnO"]["molar_mass"])
   
    # Stage 4: Fe2O3 → Fe3O4
    if a > a3:
        delta4 = min(a, a4) - a3
        xi[3] = delta4 * conv_factor * total_feed_mass / (3 * species["Fe2O3"]["molar_mass"] - 2 * species["Fe3O4"]["molar_mass"])
   
    # Stage 5: Fe3O4 → FeO
    if a > a4:
        delta5 = min(a, a5) - a4
        xi[4] = delta5 * conv_factor * total_feed_mass / (species["Fe3O4"]["molar_mass"] - 3 * species["FeO"]["molar_mass"])
   
    # Stage 6: FeO → Fe
    if a > a5:
        delta6 = a - a5
        xi[5] = delta6 * conv_factor * total_feed_mass / (species["FeO"]["molar_mass"] - species["Fe"]["molar_mass"])
   
    return xi

def compute_all_xi(alpha_values, species, w0, wh):
    xi_all = np.zeros((len(alpha_values), 6))
    for i, a in enumerate(alpha_values):
        xi_all[i, :] = compute_xi_for_alpha(a, species, w0, wh)
    return xi_all

def compute_mass_balances(xi, species):
    m_MnO2 = species["MnO2"]["initial_mass"] - 2 * xi[:, 0] * species["MnO2"]["molar_mass"]
    m_Mn2O3 = species["Mn2O3"]["initial_mass"] + (xi[:, 0] - 3 * xi[:, 1]) * species["Mn2O3"]["molar_mass"]
    m_Mn3O4 = species["Mn3O4"]["initial_mass"] + (2 * xi[:, 1] - xi[:, 2]) * species["Mn3O4"]["molar_mass"]
    m_MnO   = species["MnO"]["initial_mass"]   + 3 * xi[:, 2] * species["MnO"]["molar_mass"]
    m_Fe2O3 = species["Fe2O3"]["initial_mass"] - 3 * xi[:, 3] * species["Fe2O3"]["molar_mass"]
    m_Fe3O4 = species["Fe3O4"]["initial_mass"] + (2 * xi[:, 3] - xi[:, 4]) * species["Fe3O4"]["molar_mass"]
    m_FeO   = species["FeO"]["initial_mass"]   + (3 * xi[:, 4] - xi[:, 5]) * species["FeO"]["molar_mass"]
    m_Fe    = species["Fe"]["initial_mass"]    + xi[:, 5] * species["Fe"]["molar_mass"]
    m_SiO2 = species["SiO2"]["initial_mass"]
    m_Al2O3 = species["Al2O3"]["initial_mass"]
    m_CaO = species["CaO"]["initial_mass"]
    m_MgO = species["MgO"]["initial_mass"]
    m_C = species["C"]["initial_mass"]
    m_non_reactive = m_SiO2 + m_Al2O3 + m_CaO + m_MgO + m_C
   
    m_total = m_SiO2 + m_Al2O3 + m_CaO + m_MgO + m_C + m_MnO2 + m_Mn2O3 + m_Mn3O4 + m_MnO + m_Fe2O3 + m_Fe3O4 + m_FeO + m_Fe
    return {"MnO2": m_MnO2, "Mn2O3": m_Mn2O3, "Mn3O4": m_Mn3O4, "MnO": m_MnO,
            "Fe2O3": m_Fe2O3, "Fe3O4": m_Fe3O4, "FeO": m_FeO, "Fe": m_Fe, "SiO2": m_SiO2, "Al2O3": m_Al2O3, "CaO": m_CaO, "MgO": m_MgO, "total": m_total, "non_reactive": m_non_reactive}

def compute_mass_fractions(xi, species):
    """
    Compute the mass fractions for each species based on the reaction extents xi.
    """
    mass_data = compute_mass_balances(xi, species)
    total_mass = mass_data["total"]
    mass_fractions = {key: mass_data[key] / total_mass for key in mass_data if key != "total"}
    mass_fractions["total"] = np.ones_like(total_mass)  # total mass fraction is always 1
    return mass_fractions

def gas_phase_calculations(mass_total, d_alpha_dt, H2_inlet_rate):
    """
    Compute gas-phase outlet properties based on the solid-phase reaction progress.
    """
    o2_mass = mass_total[0] - mass_total[-1]  # kg of oxygen removed
   
   
    H2_total_needed = o2_mass / 16.0  # kmol of H2 needed = kmol of O2 removed
   
    H2_consumption_rate = H2_total_needed * d_alpha_dt  # kmol/min
   
    H2_inlet_mols = H2_total_needed * 2 # kmol
   
    H2_mass_in = H2_inlet_rate * (2 / 60) # kg/s
   
    H2_outlet_rate = H2_inlet_rate - H2_consumption_rate
   
    H2O_rate = H2_consumption_rate  # kmol/min
   
    total_gas_rate = H2_outlet_rate + H2O_rate # kmol/min
   
    H2_mass_out = H2_outlet_rate * 2    # kg/min
    H2O_mass_out = H2O_rate * 18        # kg/min
   
    total_outlet_mass_flow = H2_mass_out + H2O_mass_out # kg/min
    H2_mass_fraction = H2_mass_out / total_outlet_mass_flow
    H2O_mass_fraction = H2O_mass_out / total_outlet_mass_flow
   
    results = {
            "H2_total_needed": H2_total_needed, # kmol
            "H2_consumption_rate": H2_consumption_rate, # kmol / min
            "H2_inlet_rate": H2_inlet_rate, # kmol / min
            "H2_outlet_rate": H2_outlet_rate, # kmol / min
            "H2O_rate": H2O_rate, # kmol / min
            "total_gas_rate": total_gas_rate, # kmol / min
            "H2_mass_in": H2_mass_in,
            "H2_mass_out": H2_mass_out, # kg / min
            "H2O_mass_out": H2O_mass_out, # kg / min
            "H2_mass_fraction": H2_mass_fraction,
            "H2O_mass_fraction": H2O_mass_fraction,
        }
    return results
 
#%% ==================== Energy Balance Modelling ====================

def cp_Mn3O4(T):

    # Validate temperature range
    if T < 298:
        raise ValueError("Below valid range for Mn3O4.")
    if T > 2000:
        raise ValueError("Above valid range for Mn3O4.")


    if 298.15 <= T <= 2000:
       
        A, B, C, D, E = (120.48, 0.078, 1.16e-5, -1.527e-8, 0.003)

    # Shomate equation for Cp:
    # Cp(T) = A + B*T + C*T^2 + D*T^3 + E/T^2
    cp = A + B*T + C*(T**2) + D*(T**3) + E/(T**2) # J / (K * mol)
    cp_Kg = cp / 0.22882 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)


def cp_Mn2O3(T):

    # Validate temperature range
    if T < 298:
        raise ValueError("Below valid range for Mn2O3.")
    if T > 2000:
        raise ValueError("Above valid range for Mn2O3.")


    if 298.15 <= T <= 2000:
       
        A, B, C, D, E = (104.482, -0.002, 5.98e-5, -2.48e-8, 0.003)

    # Shomate equation for Cp:
    # Cp(T) = A + B*T + C*T^2 + D*T^3 + E/T^2
    cp = A + B*T + C*(T**2) + D*(T**3) + E/(T**2) # J / (K * mol)
    cp_Kg = cp / 0.15788 # J / (K * Kg)
    return cp_Kg

def cp_MnO2(T):

    # Validate temperature range
    if T < 298:
        raise ValueError("Below valid range for MnO2.")
    if T > 2000:
        raise ValueError("Above valid range for MnO2.")

    if 298.15 <= T <= 2000:
       
        A, B, C, D, E = (3.923, 0.252, -0.000312, 1.32e-7, 0.000179)

    # Shomate equation for Cp:
    # Cp(T) = A + B*T + C*T^2 + D*T^3 + E/T^2
    cp = A + B*T + C*(T**2) + D*(T**3) + E/(T**2) # J / (K * mol)
    cp_Kg = cp / 0.08694 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)

def cp_MnO(T):

    # Ensure T is within the valid range for alpha-Fe (298 K to 1183 K).
    if T < 298:
        raise ValueError("Below valid range for MnO.")
    if T > 2000:
        raise ValueError("Above valid range for MnO.")
   
    # Piecewise selection of coefficients
    if 298 <= T <= 623.21:
        A, B, C= ( 0.611448e1, 0.280163e-2, -0.598956e5)
    elif 623.21 < T <= 980.21:
        A, B, C = (0.573001e1, 0.335735e-2, -0.450511e5)
    else:  
        A, B, C = (0.818319e1, 0.904175e-3, -0.820377e5)
   
    # NASA - JANAF Aerotherm
    cp = A + B*T + C/(T**2) # cal / (K * mol)
    cp_J = cp * 4.184 # J / (K * mol)
    cp_Kg = cp_J / 0.07094 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)

def cp_Fe(T):


    if T < 298:
        raise ValueError("Below valid range for Fe.")
    if T > 2000:
        raise ValueError("Above valid range for Fe.")
   
    # Piecewise selection of coefficients
    if 298 <= T <= 385.21:
        A, B, C= (0.44113e1, 0.53881e-2, -0.175800e4 )
       
    elif 385.21 < T <= 500.21:
        A, B, C = (0.281069e1, 0.787171e-2, 0.937615e5)
       
    elif 500.21 < T <= 900.21:
        A, B, C = (-0.195896e1, 0.126295e-1, 0.692575e6)
       
    elif 900.21 < T <= 1042.21:
        A, B, C = ( -0.101379e4, -0.72813e0, 0.298666e9)
       
    else:  
        A, B, C = (-0.110678e4, 0.647210e0, 0.491183e9)
   
    # NASA - JANAF Aerotherm
    cp = A + B*T + C/(T**2) # cal / (K * mol)
    cp_J = cp * 4.184 # J / (K * mol)
    cp_Kg = cp_J / 0.05585 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)

def cp_Fe2O3(T):


    if T < 298:
        raise ValueError("Below valid range for Fe2O3.")
    if T > 1490.23:
        raise ValueError("Above valid range for Fe2O3.")
   
    # Piecewise selection of coefficients
    if 298 <= T <= 729.21:
        A, B, C= (0.245068e2, 0.173975e-1, -0.435099e6 )
       
    elif 729.21 < T <= 950.21:
        A, B, C = (0.238212e2, 0.183501e-1, -0.439763e6)
       
    elif 950.21 < T <= 1050.22:
        A, B, C = (0.36e2, 0, -0.353343e-3)
       
    else:  
        A, B, C = (0.317225e2, 0.175296e-2, -0.564358e4)
   
    # NASA - JANAF Aerotherm
    cp = A + B*T + C/(T**2) # cal / (K * mol)
    cp_J = cp * 4.184 # J / (K * mol)
    cp_Kg = cp_J / 0.15970 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)

def cp_SiO2(T):

    # Validate temperature range
    if T < 298:
        raise ValueError("Below valid range for SiO2.")
    if T > 1699.21:
        raise ValueError("Above valid range for SiO2.")


    if 298.15 <= T <= 1699.21:
       
        A, B, C = (0.173941e2, 0.307699e-3, -0.989658e6)

    # NASA - JANAF Aerotherm
    cp = A + B*T + C/(T**2) # cal / (K * mol)
    cp_J = cp * 4.184 # J / (K * mol)
    cp_Kg = cp_J / 0.0600843 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)


def cp_Al2O3(T):

    # Validate temperature range
    if T < 298:
        raise ValueError("Below valid range for Al2O3.")
    if T > 2000:
        raise ValueError("Above valid range for Al2O3.")


    if 298.15 <= T <= 489.21:
       
        A, B, C = (0.226264e2, 0.103223e-1, -0.606477e6)
       
    else:
       
        A, B, C = (0.283346e2, 0.254438e-2, -0.106193e7)  
       
    # NASA - JANAF Aerotherm
    cp = A + B*T + C/(T**2) # cal / (K * mol)
    cp_J = cp * 4.184 # J / (K * mol)
    cp_Kg = cp_J / 0.101961 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)

def cp_CaO(T):


    if T < 298:
        raise ValueError("Below valid range for CaO.")
    if T > 3200.21:
        raise ValueError("Above valid range for CaO.")
   
    # Piecewise selection of coefficients
    if 298 <= T <= 400.21:
        A, B, C= (0.121433e2, 0.657842e-3, -0.202015e6)
       
    elif 400.21 < T <= 603.21:
        A, B, C = (0.119380e2, 0.113622e-2, -0.199772e6)
       
    elif 603.21 < T <= 1040.21:
        A, B, C = (0.120737e2, 0.984742e-3, -0.215918e6)
       
    else:  
        A, B, C = (0.121138e2, 0.959254e-3, -0.230653e6)
   
    # NASA - JANAF Aerotherm
    cp = A + B*T + C/(T**2) # cal / (K * mol)
    cp_J = cp * 4.184 # J / (K * mol)
    cp_Kg = cp_J / 0.0560774 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)

def cp_MgO(T):


    if T < 298:
        raise ValueError("Below valid range for MgO.")
    if T > 2000:
        raise ValueError("Above valid range for MgO.")
   
    # Piecewise selection of coefficients
    if 298 <= T <= 480.21:
        A, B, C= (0.105616e2, 0.235351e-2, -0.212914e6)
       
    elif 480.21 < T <= 816.21:
        A, B, C = (0.115232e2, 0.101658e-2, -0.286606e6)        
       
    else:  
        A, B, C = (0.117167e2, 0.843867e-3, -0.321617e6)
       
    # NASA - JANAF Aerotherm
    cp = A + B*T + C/(T**2) # cal / (K * mol)
    cp_J = cp * 4.184 # J / (K * mol)
    cp_Kg = cp_J /  0.0403044 # J / (K * Kg)
    return cp_Kg # J / (K * Kg)


def cp_species(species, T):
    """
    Switch or dictionary that calls the correct function
    for each species at temperature T.
    Returns Cp in J/(kg*K).
    """
    if species == "MnO":
        return cp_MnO(T)
    elif species == "MnO2":
        return cp_MnO2(T)
    elif species == "Mn2O3":
        return cp_Mn2O3(T)
    elif species == "Mn3O4":
        return cp_Mn3O4(T)
    elif species == "Fe":
        return cp_Fe(T)
    elif species == "Fe2O3":
        return cp_Fe2O3(T)
    elif species == "Al2O3":
        return cp_Al2O3(T)
    elif species == "SiO2":
        return cp_SiO2(T)
    elif species == "CaO":
        return cp_CaO(T)
    elif species == "MgO":
        return cp_MgO(T)
    else:
        return 0.0

# Function to compute total heating time from T_initial to T_final:
def heating_time_fixed_power(species, T_initial, T_reaction, Heating_time, n_steps=10000):



    T_array = np.linspace(T_initial, T_reaction, n_steps+1)
    total_heat_required = 0.0  # J


    for i in range(n_steps):
        T1 = T_array[i]
        T2 = T_array[i+1]
        dT = (T2 - T1)

        for sp_name, sp_data in species.items():
            mass_kg = sp_data["initial_mass"]  # total mass in kg
           
            # average Cp for species from T1->T2
            cp1 = cp_species(sp_name, T1)
            cp2 = cp_species(sp_name, T2)
            cp_avg = 0.5*(cp1 + cp2)  # J/(kg*K)

            Q_part = mass_kg * cp_avg * dT
            total_heat_required += Q_part
   
    # Power = Q / time
    power = total_heat_required / (Heating_time * 3600)
    return power

def reaction_cumulative_at_alpha(a, species, w0, wh, Delta_Hr_J):
    """Compute the cumulative reaction enthalpy (J) at conversion a."""
    xi = compute_xi_for_alpha(a, species, w0, wh)  # 6-element array
    xi_mol = xi * 1000  # convert kmol to mol
    return np.sum(xi_mol * Delta_Hr_J)

def effective_Cp_at_alpha(a, species, w0, wh, T):

    xi = compute_xi_for_alpha(a, species, w0, wh)
    mass_data = compute_mass_balances(np.array([xi]), species)
    effective_cp = np.sum( mass_data["MnO2"] * cp_species("MnO2", T) +
                     mass_data["Mn2O3"] * cp_species("Mn2O3", T) +
                     mass_data["Mn3O4"] * cp_species("Mn3O4", T) +
                     mass_data["MnO"]   * cp_species("MnO", T) +
                     mass_data["Fe2O3"] * cp_species("Fe2O3", T) +
                     mass_data["Fe3O4"] * cp_species("Fe3O4", T) +
                     mass_data["FeO"]   * cp_species("FeO", T) +
                     mass_data["Fe"]    * cp_species("Fe", T) +
                     mass_data["SiO2"] * cp_species("SiO2", T) +
                     mass_data["Al2O3"] * cp_species("Al203", T) +
                     mass_data["CaO"] * cp_species("MgO", T)
                      )
    return effective_cp.item()  # return as scalar

def update_temperature(T_prev, a_prev, a_new, d_alpha_dt, dt_sec, species, w0, wh, Delta_Hr_J,
                       H2_total_needed_full, H2_inlet_rate, MW_H2, Cp_H2, T_in_H2, T_target, h_conv, A_solid, wt):
   
   
    Q_total = reaction_cumulative_at_alpha(1, species, w0, wh, Delta_Hr_J) *  (w0 - wh)
    # Compute the total initial mass (kg)
    dH = Q_total * d_alpha_dt
    # Compute the effective heat capacity (J/K) of the solids at the current conversion and temperature.
    Cp_eff = effective_Cp_at_alpha(a_prev, species, w0, wh, T_prev)
    if Cp_eff == 0:
        Cp_eff = 1  # safeguard
    # Gas-phase heat loss:
    # H2_consumption_rate in kmol/min based on the reaction rate.
    H2_consumption_rate = (H2_total_needed_full / 60) * d_alpha_dt
    # H2_out_rate (kmol/min):
    H2_out_rate = max(H2_inlet_rate - H2_consumption_rate, 0)
    # Convert to mass flow rate in kg/min:
    m_dot_H2 = H2_out_rate * MW_H2
    # Q_loss_H2 (J) over dt_min minutes:
    Q_loss_H2 = m_dot_H2 * Cp_H2 * (T_prev - T_in_H2)
    # Solids-to-gas convective heat loss (J)
    Q_loss_solids = 0#h_conv * A_solid * (T_prev - T_in_H2)   # W * s = J

    # Energy balance: the net energy available to heat the solids is (-dH) (because dH is negative for exothermic)
    # minus the energy lost to the gas.
    dT = (((-dH) - Q_loss_H2 - Q_loss_solids) / Cp_eff) * dt_sec
    T_new = T_prev + dT
    if T_new > T_target:
        T_new = T_target
    return T_new

#%% ==================== Process Operation ====================

# Calculate power required to heat up the ore from 25°C to 900°C
heat_power = heating_time_fixed_power(species, T_initial, T_reaction, Heating_time)
print("Estimated Power:", heat_power)
print("In MW: ", heat_power / 1e6)

# Set up time grid: use minutes for kinetics.
n_steps = len(reduction_time)
dt_min = reduction_time[1] - reduction_time[0]   # time step in minutes
dt_sec = dt_min * 60                     # time step in seconds

# Initialise arrays
alpha_coupled = np.zeros(n_steps)
T_energy = np.zeros(n_steps)
k_values = np.zeros(n_steps)
wt_dynamic = np.zeros(n_steps)
d_alpha_dt = np.zeros(n_steps)

# Initial boundaries
alpha_coupled[0] = 0.0
T_energy[0] = T_reaction
k_values[0] = rate_constant(k0_best, Ea_best, R_val, T_reaction)
cum_time = 0.0  # cumulative kinetics time in minutes
wt_dynamic[0] = w0

# Reaction enthalpy (kJ/mol) for each stage; convert to J/mol:
Delta_Hr_J = Delta_Hr * 1000

# Create a dictionary for kinetic parameters.
kinetic_params = {
    'k0': k0_best,
    'Ea': Ea_best,
    'delta2': delta2_best,
    'delta3': delta3_best,
    'lambda_decay': lambda_decay_best
}

# Main coupled loop:
for i in range(1, n_steps):
    # Update kinetics: get new conversion and dα/dt using current T
    alpha_new, d_alpha_dt_new, k_current, cum_time = update_kinetics(T_energy[i-1],
                                                                  alpha_coupled[i-1],
                                                                  cum_time,
                                                                  dt_min,
                                                                  kinetic_params)
    d_alpha_dt[i] = d_alpha_dt_new
    alpha_coupled[i] = alpha_new
    k_values[i] = k_current
    # Calculate the weight fraction over time
    wt_dynamic_new = w0 - alpha_coupled[i - 1] * (w0 - wh)
    wt_dynamic[i] = wt_dynamic_new
   
    # Update temperature using energy balance.
    # We use the lumped average Cp (Cp_avg) computed earlier.
    T_new = update_temperature(T_energy[i-1], alpha_coupled[i-1], alpha_coupled[i], d_alpha_dt[i], dt_sec, species, w0, wh, Delta_Hr_J, H2_total_needed_full, H2_inlet_rate, MW_H2, Cp_H2, T_in_H2, T_target, h_conv, A_solid, wt_dynamic[i-1])
   
    # Optionally, if T exceeds T_target, cap it (and record excess cooling if desired)
    if T_new <= T_target:
        T_energy[i] = T_new
    else:
        T_energy[i] = T_target


# --- Reactor Geometry Calculations --- #

# from your alpha vs time:
desired_alpha = 0.99
idx = np.where(alpha_coupled >= desired_alpha)[0]
if len(idx) == 0:
    t_res = reduction_time[-1]
    print("Warning: desired alpha not reached within time range!")
else:
    t_res = reduction_time[idx[0]]

M_solids   = compute_total_feed_mass(species) # kg per batch
rho_bulk   = rho_ore *(1 - porosity)  # kg/m^3, bulk density = solid denisty( 1 - porosity)
V_bulk   = M_solids / rho_bulk
V_solids_only = M_solids / rho_ore
V_void_only = V_bulk - V_solids_only
V_reactor = V_bulk / 0.75
V_gas = V_reactor - V_solids_only
print("Volume of Reactor in m^3", V_reactor)


# --- Mass Transfer Calculations --- #

# Calculate the reaction extents xi and mass data based on conversion
xi = compute_all_xi(alpha_coupled, species, w0, wh)
mass_data_new = compute_mass_balances(xi, species)
mass_fractions = compute_mass_fractions(xi, species)

# Gas stream
gas_results = gas_phase_calculations(mass_data_new["total"], d_alpha_dt, H2_inlet_rate)

H2_mass_fraction = gas_results["H2_mass_fraction"]
H2O_mass_fraction = gas_results["H2O_mass_fraction"]

H2_inlet_rate = gas_results["H2_inlet_rate"]
total_gas_rate = gas_results["total_gas_rate"]

#%% ==================== Plotting Results ====================

# Alpha vs Time
plt.figure(figsize=(10, 6))
plt.plot(reduction_time, alpha_coupled, label="900 °C (Dynamic Model)", linewidth=2)
plt.scatter(experimental_time, experimental_alpha, color='red', zorder=5, label='Experimental Data')
plt.xlabel('Time (min)')
plt.ylabel('Alpha (α)')
plt.yticks(np.arange(0, 1.05, 0.05))
plt.title('Model Fit at 900°C')
plt.ylim(-0.05, 1.05)
plt.grid(True)
plt.legend()
plt.show()

# Plot temperature evolution
plt.figure(figsize=(10,6))
plt.plot(reduction_time, T_energy - 273.15, label="Furnace Temperature (K)", linewidth=2)
plt.xlabel("Time (min)")
plt.ylabel("Temperature (K)")
plt.title("Temperature Evolution from Reaction Heat")
plt.grid(True)
plt.legend()
plt.show()


# Plot evolution of rate constant k(T)
plt.figure(figsize=(10,6))
plt.plot(reduction_time, k_values, label="k(T) from Arrhenius", linewidth=2)
plt.xlabel("Time (min)")
plt.ylabel("Rate Constant k")
plt.title("Evolution of k(T) with Temperature")
plt.grid(True)
plt.legend()
plt.show()

# Mass fraction vs Alpha for each species (and total mass)
plt.figure(figsize=(12, 7))
plt.plot(alpha_coupled, mass_fractions["MnO2"], label="MnO2", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["Mn2O3"], label="Mn2O3", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["Mn3O4"], label="Mn3O4", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["MnO"], label="MnO", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["Fe2O3"], label="Fe2O3", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["Fe3O4"], label="Fe3O4", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["FeO"], label="FeO", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["Fe"], label="Fe", linewidth=2)
plt.plot(alpha_coupled, mass_fractions["total"], label="Total Mass", linewidth=2)
plt.xlabel("Alpha (α)", fontsize=14)
plt.ylabel("Mass Fraction", fontsize=14)
plt.title("Mass Fraction vs Alpha", fontsize=16)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.show()


# Mass fraction vs Time
plt.figure(figsize=(12, 7))
plt.plot(reduction_time, mass_fractions["MnO2"], label="MnO2", linewidth=2)
plt.plot(reduction_time, mass_fractions["Mn2O3"], label="Mn2O3", linewidth=2)
plt.plot(reduction_time, mass_fractions["Mn3O4"], label="Mn3O4", linewidth=2)
plt.plot(reduction_time, mass_fractions["MnO"], label="MnO", linewidth=2)
plt.plot(reduction_time, mass_fractions["Fe2O3"], label="Fe2O3", linewidth=2)
plt.plot(reduction_time, mass_fractions["Fe3O4"], label="Fe3O4", linewidth=2)
plt.plot(reduction_time, mass_fractions["FeO"], label="FeO", linewidth=2)
plt.plot(reduction_time, mass_fractions["Fe"], label="Fe", linewidth=2)
plt.plot(reduction_time, mass_fractions["total"], label="Total Mass", linewidth=2)
plt.xlabel("Time (min)", fontsize=14)
plt.ylabel("Mass Fraction", fontsize=14)
plt.title("Mass Fraction vs Time", fontsize=16)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.show()

# Outlet gas mass fraction vs. time:
plt.figure(figsize=(10, 6))
plt.plot(reduction_time, H2_mass_fraction, label='Outlet H₂ Mass Fraction', linewidth=2)
plt.plot(reduction_time, H2O_mass_fraction, label='Outlet H₂O Mass Fraction', linewidth=2)
plt.xlabel('Time (min)', fontsize=14)
plt.ylabel('H₂ Mass Fraction', fontsize=14)
plt.yticks(np.arange(0, 1.05, 0.1))
plt.title('Outlet H₂ Mass Fraction vs Time', fontsize=16)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.show()
