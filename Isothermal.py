#%% ==================== Importing necessary libraries ====================
import numpy as np
import matplotlib.pyplot as plt

#==================== Initial Parameters ====================
# Kinetic and model parameters

k0_best = 11.9 # could also be 156
Ea_best = 40.4 # could also be 63.67
delta2_best = 0.12
delta3_best = 0.00000005
lambda_decay_best = 1e10

# Fixed Constants

R_val = 0.008314
T = 900 + 273.15
# Species data: molar mass (kg/kmol) and initial mass (kg)

species = {
    "MnO2": {"molar_mass": 86.94, "initial_mass": 7140},
    "Mn2O3": {"molar_mass": 157.88, "initial_mass": 8230},
    "Mn3O4": {"molar_mass": 228.82, "initial_mass": 0},
    "MnO":   {"molar_mass": 70.94,  "initial_mass": 7150},
    "Fe2O3": {"molar_mass": 159.7,  "initial_mass": 2950},
    "Fe3O4": {"molar_mass": 231.55, "initial_mass": 0},
    "FeO":   {"molar_mass": 71.85,  "initial_mass": 0},
    "Fe":    {"molar_mass": 55.85,  "initial_mass": 0},
    "SiO2": {"molar_mass": 10 , "initial_mass": 2360},
    "Al2O3": {"molar_mass": 10, "initial_mass": 80},
    "CaO": {"molar_mass": 10, "initial_mass": 5490},
    "MgO": {"molar_mass": 10, "initial_mass": 1430},
    "C": {"molar_mass": 10, "initial_mass": 40}        
}


# Define initial and final weight fractions (for the ore)
w0 = 1.0    # initial weight fraction
wh = 0.913  # final weight fraction (at complete conversion)

# Experimental Data for 900°C
experimental_time = np.array([0, 6, 10, 15, 20, 25, 30, 35, 40, 65, 92.5, 120])
experimental_alpha = np.array([0, 0.75, 0.78, 0.85, 0.87, 0.9, 0.914, 0.925, 0.94, 0.96, 1, 1])


# ==================== Kinetic Modelling ====================

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

def alpha_dynamic(t, k, delta2, delta3, lambda_decay):
    alpha_vals = np.zeros_like(t)
    for i, time in enumerate(t):
        if i == 0:
            alpha_vals[i] = 0
        else:
            # Determine the exponent based on the previous conversion value.
            exponent = n_alpha(alpha_vals[i-1], delta2, delta3, lambda_decay)
            # Compute the conversion using the Avrami equation.
            alpha_vals[i] = 1 - np.exp(- (k * time)**exponent)
    return alpha_vals

#==================== Mass Transfer Modelling ====================

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

        ξ = (Δα) * (w0 – wh) * im / ΔM,

    where ΔM is taken as ~16 kg/kmol.
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
    m_total = m_SiO2 + m_Al2O3 + m_CaO + m_MgO + m_C + m_MnO2 + m_Mn2O3 + m_Mn3O4 + m_MnO + m_Fe2O3 + m_Fe3O4 + m_FeO + m_Fe

    return {"MnO2": m_MnO2, "Mn2O3": m_Mn2O3, "Mn3O4": m_Mn3O4, "MnO": m_MnO,
            "Fe2O3": m_Fe2O3, "Fe3O4": m_Fe3O4, "FeO": m_FeO, "Fe": m_Fe, "total": m_total}

 

def compute_mass_fractions(xi, species):
    """
    Compute the mass fractions for each species based on the reaction extents xi.
    """
    mass_data = compute_mass_balances(xi, species)
    total_mass = mass_data["total"]
    mass_fractions = {key: mass_data[key] / total_mass for key in mass_data if key != "total"}
    mass_fractions["total"] = np.ones_like(total_mass)  # total mass fraction is always 1
    return mass_fractions

 

def gas_phase_calculations(mass_total, d_alpha_dt):
    """
    Compute gas-phase outlet properties based on the solid-phase reaction progress.
    """

    o2_mass = mass_total[0] - mass_total[-1]  # kg of oxygen removed
    H2_total_needed = o2_mass / 16.0  # kmol of H2 needed = kmol of O2 removed
    H2_consumption_rate = H2_total_needed * d_alpha_dt  # kmol/min
    H2_inlet_mols = H2_total_needed * 2 # kmol
    H2_inlet_rate = H2_inlet_mols / 120 # kmol/min
    H2_mass_in = H2_inlet_rate * (2 / 60) # kg/s
    H2_outlet_rate = np.maximum(H2_inlet_rate - H2_consumption_rate, 0) # kmol/min
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

 

def pressure_calc(molar_flowrate, Vr, R, T):
    moles = molar_flowrate * (1000 / 60) # mol / s
    P_Pa = (moles * R * T) / Vr # Pa
    P_kPa = P_Pa / 1000 # kPa
    P_bar = P_kPa * 0.01 # bar

    results = {
            "P_Pa": P_Pa,
            "P_kPa": P_kPa,
            "P_bar": P_bar
        }
    return results
#==================== Energy Balance Modelling ====================


#==================== Process Operation ====================

# --- Kinetic Calculations --- #

# Set a time frame (in minutes)

dense_time = np.linspace(0, 120, 121)

# Calculate the kinetic rate constant k

k_val = rate_constant(k0_best, Ea_best, R_val, T)

# Compute the conversion (alpha) dynamically using the conversion-based dynamic exponent

alpha_model = alpha_dynamic(dense_time, k_val, delta2_best, delta3_best, lambda_decay_best)

gradient = np.gradient(alpha_model, dense_time)


# --- Reactor Geometry Calculations --- #


# from your alpha vs time:

desired_alpha = 0.99
idx = np.where(alpha_model >= desired_alpha)[0]
if len(idx) == 0:
    t_res = dense_time[-1]
    print("Warning: desired alpha not reached within time range!")
else:
    t_res = dense_time[idx[0]]

 

M_solids   = compute_total_feed_mass(species) # kg per batch
rho_solid = 4462.17
porosity = 0.253
rho_bulk   = rho_solid *(1 - porosity)  # kg/m^3, bulk density = solid denisty( 1 - porosity)
V_bulk   = M_solids / rho_bulk
V_solids_only = M_solids / rho_solid
V_void_only = V_bulk - V_solids_only

design_factor = 0.75
V_reactor = V_bulk / 0.75

 

# --- Mass Transfer Calculations --- #

 

# Calculate the weight fraction over time
wt_dynamic = w0 - alpha_model * (w0 - wh)

 

# Calculate the reaction extents xi and mass data based on conversion

xi = compute_all_xi(alpha_model, species, w0, wh)
mass_data = compute_mass_balances(xi, species)
mass_fractions = compute_mass_fractions(xi, species)

# Gas stream

gas_results = gas_phase_calculations(mass_data["total"], gradient)
H2_mass_fraction = gas_results["H2_mass_fraction"]
H2O_mass_fraction = gas_results["H2O_mass_fraction"]
H2_inlet_rate = gas_results["H2_inlet_rate"]
total_gas_rate = gas_results["total_gas_rate"]

#Pressure of inlet

pH2 = pressure_calc(H2_inlet_rate, V_reactor, 8.314, T)
P = pressure_calc(total_gas_rate, V_reactor, 8.314, T)
P_bar = P["P_bar"]

# Alpha vs Time

plt.figure(figsize=(10, 6))
plt.plot(dense_time, alpha_model, label="900 °C (Dynamic Model)", linewidth=2)
plt.scatter(experimental_time, experimental_alpha, color='red', zorder=5, label='Experimental Data')
plt.xlabel('Time (min)')
plt.ylabel('Alpha (α)')
plt.yticks(np.arange(0, 1.05, 0.05))
plt.title('Model Fit at 900°C')
plt.ylim(-0.05, 1.05)
plt.grid(True)
plt.legend()
plt.show()

 

# Weight Fraction vs Time

plt.figure(figsize=(10, 6))
plt.plot(dense_time, wt_dynamic, label="Weight Fraction (wt)", linewidth=2, color='purple')
plt.xlabel("Time (min)", fontsize=14)
plt.ylabel("Weight Fraction (wt)", fontsize=14)
plt.title("Continuous Weight Fraction vs Time", fontsize=16)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.show()

 

# Mass fraction vs Alpha for each species (and total mass)

plt.figure(figsize=(12, 7))
plt.plot(alpha_model, mass_fractions["MnO2"], label="MnO2", linewidth=2)
plt.plot(alpha_model, mass_fractions["Mn2O3"], label="Mn2O3", linewidth=2)
plt.plot(alpha_model, mass_fractions["Mn3O4"], label="Mn3O4", linewidth=2)
plt.plot(alpha_model, mass_fractions["MnO"], label="MnO", linewidth=2)
plt.plot(alpha_model, mass_fractions["Fe2O3"], label="Fe2O3", linewidth=2)
plt.plot(alpha_model, mass_fractions["Fe3O4"], label="Fe3O4", linewidth=2)
plt.plot(alpha_model, mass_fractions["FeO"], label="FeO", linewidth=2)
plt.plot(alpha_model, mass_fractions["Fe"], label="Fe", linewidth=2)
plt.plot(alpha_model, mass_fractions["total"], label="Total Mass", linewidth=2)
plt.xlabel("Alpha (α)", fontsize=14)
plt.ylabel("Mass Fraction", fontsize=14)
plt.title("Mass Fraction vs Alpha", fontsize=16)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.show()

 

 

# Mass fraction vs Time

plt.figure(figsize=(12, 7))
plt.plot(dense_time, mass_fractions["MnO2"], label="MnO2", linewidth=2)
plt.plot(dense_time, mass_fractions["Mn2O3"], label="Mn2O3", linewidth=2)
plt.plot(dense_time, mass_fractions["Mn3O4"], label="Mn3O4", linewidth=2)
plt.plot(dense_time, mass_fractions["MnO"], label="MnO", linewidth=2)
plt.plot(dense_time, mass_fractions["Fe2O3"], label="Fe2O3", linewidth=2)
plt.plot(dense_time, mass_fractions["Fe3O4"], label="Fe3O4", linewidth=2)
plt.plot(dense_time, mass_fractions["FeO"], label="FeO", linewidth=2)
plt.plot(dense_time, mass_fractions["Fe"], label="Fe", linewidth=2)
plt.plot(dense_time, mass_fractions["total"], label="Total Mass", linewidth=2)
plt.xlabel("Time (min)", fontsize=14)
plt.ylabel("Mass Fraction", fontsize=14)
plt.title("Mass Fraction vs Time", fontsize=16)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.show()

# Outlet gas mass fraction vs. time:

plt.figure(figsize=(10, 6))
plt.plot(dense_time, H2_mass_fraction, label='Outlet H₂ Mass Fraction', linewidth=2)
plt.plot(dense_time, H2O_mass_fraction, label='Outlet H₂O Mass Fraction', linewidth=2)
plt.xlabel('Time (min)', fontsize=14)
plt.ylabel('H₂ Mass Fraction', fontsize=14)
plt.yticks(np.arange(0, 1.05, 0.1))
plt.title('Outlet H₂ Mass Fraction vs Time', fontsize=16)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.show()

 

# Pressure vs. time:
plt.figure(figsize=(10, 6))
plt.plot(dense_time, P_bar, label='Pressure vs Time', linewidth=2)
plt.xlabel('Time (min)', fontsize=14)
plt.ylabel('Pressure (bar)', fontsize=14)
plt.title('Pressure vs Time', fontsize=16)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.show()