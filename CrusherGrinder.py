from Ores import Ore

class CrusherGrinder:
    """
    Final crush:
      Moves 9.5-50 mm -> 2.5-9.5 mm
      Removes <2.5 mm
      (No mention of 10-50 mm anymore!)
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
        # 1) optional dust removal
        if self.dust_loss_fraction > 0:
            self.remove_dust(ore, self.dust_loss_fraction)

        # 2) lumps from 9.5-50 -> 2.5-9.5
        lumps_mass = ore.get_size_bin_mass(self.large_bin)
        if lumps_mass > 0:
            med_mass = ore.get_size_bin_mass(self.medium_bin)
            new_med = med_mass + lumps_mass
            ore.set_size_bin_mass(self.large_bin, 0.0)
            ore.set_size_bin_mass(self.medium_bin, new_med)

        # 3) remove <2.5 mm
        fines_mass = ore.get_size_bin_mass(self.fines_bin)
        if fines_mass > 0:
            old_mass = ore.mass
            new_mass = old_mass - fines_mass
            frac_removed = fines_mass / old_mass

            # proportionally remove from phases
            for phase_name, frac in ore.phases.items():
                if frac <= 0:
                    continue
                old_phase_mass = frac * old_mass
                new_phase_mass = old_phase_mass * (1.0 - frac_removed)
                new_frac = new_phase_mass / new_mass if new_mass>0 else 0
                ore.phases[phase_name] = new_frac

            ore.mass = new_mass
            ore.re_normalize_phases()
            ore.set_size_bin_mass(self.fines_bin, 0.0)

        ore.re_normalize_size_dist()
        return ore

    def remove_dust(self, ore: Ore, fraction: float):
        old_mass = ore.mass
        dust_loss = old_mass * fraction
        new_mass = old_mass - dust_loss

        if dust_loss <= 0:
            return

        for phase_name, frac in ore.phases.items():
            if frac>0:
                old_phase_mass = frac * old_mass
                new_phase_mass = old_phase_mass * (1.0 - fraction)
                new_frac = new_mass and new_phase_mass / new_mass
                ore.phases[phase_name] = new_frac
        ore.mass = new_mass
        ore.re_normalize_phases()