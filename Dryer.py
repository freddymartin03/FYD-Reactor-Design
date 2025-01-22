from Ores import Ore

class Dryer:
    """
    Drying step to remove free H2O.
    """
    def __init__(self, drying_efficiency: float = 1.0):
        self.drying_efficiency = drying_efficiency
        self.total_water_removed = 0.0

    def process(self, ore: Ore) -> Ore:
        if "H2O" not in ore.phases:
            return ore

        old_water = ore.get_phase_mass("H2O")
        if old_water <= 0:
            return ore

        remove_amt = old_water * self.drying_efficiency
        ore.remove_phase_mass("H2O", remove_amt)
        # also reduce ore.mass
        ore.mass -= remove_amt
        ore.re_normalize_phases()

        self.total_water_removed += remove_amt
        return ore
