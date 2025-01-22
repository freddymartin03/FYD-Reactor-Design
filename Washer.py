from Ores import Ore

class Washer:
    """
    Washing to remove clay/fines/impurities.
    Example: remove some fraction of 'Impurities' from the ore.
    """
    def __init__(self, removal_efficiency: float = 1):
        self.removal_efficiency = removal_efficiency
        self.total_removed = 0.0

    def process(self, ore: Ore) -> Ore:
        if "Impurities" not in ore.phases:
            return ore

        old_imp_mass = ore.get_phase_mass("Impurities")
        if old_imp_mass <= 0:
            return ore

        removed = old_imp_mass * self.removal_efficiency
        ore.remove_phase_mass("Impurities", removed)

        # remove that from ore.mass as well
        ore.mass -= removed
        ore.re_normalize_phases()

        self.total_removed += removed
        return ore

