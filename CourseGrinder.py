from Ores import Ore

class CoarseCrusher:
    """
    Coarse crushing:
      Moves mass from >50 mm -> 9.5-50 mm
      leaving 2.5-9.5 mm and <2.5 mm bins alone.
    """
    def __init__(self,
                 xl_bin: str = ">50 mm",
                 large_bin: str = "9.5-50 mm",
                 medium_bin: str = "2.5-9.5 mm",
                 fines_bin: str = "<2.5 mm",
                 dust_loss_fraction: float = 0.0):
        self.xl_bin = xl_bin
        self.large_bin = large_bin
        self.medium_bin = medium_bin
        self.fines_bin = fines_bin
        self.dust_loss_fraction = dust_loss_fraction

    def process(self, ore: Ore) -> Ore:
        # Optional dust removal
        if self.dust_loss_fraction > 0:
            self.remove_dust(ore, self.dust_loss_fraction)

        # Hammer lumps from >50 mm -> 9.5-50 mm
        lumps_mass = ore.get_size_bin_mass(self.xl_bin)
        if lumps_mass > 0:
            med_mass = ore.get_size_bin_mass(self.large_bin)
            new_med = med_mass + lumps_mass
            ore.set_size_bin_mass(self.xl_bin, 0.0)
            ore.set_size_bin_mass(self.large_bin, new_med)

        ore.re_normalize_size_dist()
        return ore

    def remove_dust(self, ore: Ore, fraction: float):
        old_mass = ore.mass
        dust_loss = old_mass * fraction
        new_mass = old_mass - dust_loss

        if dust_loss <= 0:
            return

        # Remove proportionally from all phases
        for phase_name, frac in ore.phases.items():
            if frac > 0:
                old_phase_mass = frac * old_mass
                new_phase_mass = old_phase_mass * (1.0 - fraction)
                new_frac = new_phase_mass / new_mass
                ore.phases[phase_name] = new_frac

        ore.mass = new_mass
        ore.re_normalize_phases()