from Ores import Ore

class Calcinator:
    """
    Example: calcination in air, including oxidation of carbon to CO2.
    Also transforms Mn-oxides, removes H2O, etc.
    """

    def __init__(self, 
                 heating_rate: float = 10.0,
                 final_temp: float = 900.0,
                 hold_time: float = 60.0,
                 gas_flow_rate: float = 1.0,
                 phase_densities: dict = None):
        """
        :param heating_rate: ramp rate (°C/min)
        :param final_temp: final calcination temperature (°C)
        :param hold_time: soak time (minutes)
        :param gas_flow_rate: air flow rate, if relevant to advanced modeling
        :param phase_densities: for recomputing bulk density after transformations
        """
        self.heating_rate = heating_rate
        self.final_temp = final_temp
        self.hold_time = hold_time
        self.gas_flow_rate = gas_flow_rate
        self.phase_densities = phase_densities if phase_densities else {}

        # Track total CO2 released (kg)
        self.total_co2_released = 0.0

    def process(self, ore: Ore) -> Ore:
        """
        1) Heat from ore.temperature to final_temp
        2) Perform transformations
        3) Recompute density if desired
        """
        if ore.temperature < self.final_temp:
            delta_t = self.final_temp - ore.temperature
            time_to_heat = delta_t / self.heating_rate
            ore.temperature = self.final_temp
        else:
            time_to_heat = 0.0

        total_calcination_time = time_to_heat + self.hold_time

        # Perform transformations
        self.perform_phase_transformations(ore, total_calcination_time)

        # Recompute density (optional)
        self.recompute_bulk_density(ore)

        return ore

    def perform_phase_transformations(self, ore: Ore, calcination_time: float) -> None:
        """
        Typical steps in air:
          1) Remove H2O
          2) Burn some fraction of elemental carbon -> CO2
          3) e.g., MnO2 -> Mn2O3, Mn2O3 -> Mn3O4 partial transformations
        """
        # (1) Remove H2O if present
        if "H2O" in ore.phases:
            water_mass = ore.get_phase_mass("H2O")
            ore.remove_phase_mass("H2O", water_mass)

        # (2) Oxidize carbon -> CO2 (assume 100% or partial)
        self.transform_carbon_to_co2(ore, fraction=1)  # burn all carbon

        # (3) If you have Mn-oxide transformations, do them
        # e.g. partial MnO2->Mn2O3, Mn2O3->Mn3O4, etc.
        # self.transform_mno2_to_mn2o3(ore, fraction=0.60)
        # self.transform_mn2o3_to_mn3o4(ore, fraction=0.20)

    def transform_carbon_to_co2(self, ore: Ore, fraction: float):
        """
        Oxidize 'C' in the ore to CO2:
        C + O2 -> CO2
        M(C)=12 g/mol, M(CO2)=44 g/mol => 1 kg C forms 3.67 kg CO2
        We remove the carbon from the ore mass, and track how much CO2 is released.
        """
        if "C" not in ore.phases:
            return

        old_carbon_mass = ore.get_phase_mass("C")
        burn_mass = old_carbon_mass * fraction  # how much C we actually burn
        if burn_mass <= 0:
            return

        # Amount of CO2 produced from burned carbon
        produced_co2 = burn_mass * (44.0 / 12.0)  # 3.67 factor

        # Remove that burned carbon from the ore
        ore.remove_phase_mass("C", burn_mass)

        # Subtract carbon from total ore mass (oxygen is from air, not from the ore)
        ore.mass -= burn_mass
        ore.re_normalize_phases()

        # Track how much CO2 was released
        self.total_co2_released += produced_co2

    # ------------------------------------------------
    # Optional transformations for Mn-oxides, if needed
    # ------------------------------------------------
    def transform_mno2_to_mn2o3(self, ore: Ore, fraction: float):
        """
        Example:
          2 MnO2 -> Mn2O3 + 1/2 O2
        """
        pass  # omitted here for brevity

    def transform_mn2o3_to_mn3o4(self, ore: Ore, fraction: float):
        """
        Example:
          3 Mn2O3 -> 2 Mn3O4 + 1/2 O2
        """
        pass  # omitted here for brevity

    # ------------------------------------------------
    # Bulk Density Recalculation (optional)
    # ------------------------------------------------
    def recompute_bulk_density(self, ore: Ore):
        if not self.phase_densities:
            return
        denom = 0.0
        for phase_name, frac in ore.phases.items():
            rho_phase = self.phase_densities.get(phase_name, 3000.0)
            denom += (frac / rho_phase)
        if denom > 0:
            new_density = 1.0 / denom
            ore.update_bulk_density(new_density)