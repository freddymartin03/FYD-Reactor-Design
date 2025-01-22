from dataclasses import dataclass, field
from typing import Dict

@dataclass
class Ore:
    """
    Represents a batch of ore in terms of:
      - total mass
      - phases (each as a fraction of total mass)
      - temperature
      - bulk density
      - particle-size distribution
    """
    mass: float
    temperature: float
    density: float
    phases: Dict[str, float] = field(default_factory=dict)
    size_distribution: Dict[str, float] = field(default_factory=dict)

    def re_normalize_phases(self) -> None:
        total_fraction = sum(self.phases.values())
        if total_fraction > 0:
            for phase in self.phases:
                self.phases[phase] /= total_fraction
        else:
            for phase in self.phases:
                self.phases[phase] = 0.0

    def re_normalize_size_dist(self) -> None:
        total_sd = sum(self.size_distribution.values())
        if total_sd > 0:
            for bin_name in self.size_distribution:
                self.size_distribution[bin_name] /= total_sd
        else:
            for bin_name in self.size_distribution:
                self.size_distribution[bin_name] = 0.0

    def get_phase_mass(self, phase_name: str) -> float:
        frac = self.phases.get(phase_name, 0.0)
        return frac * self.mass

    def set_phase_mass(self, phase_name: str, new_mass: float) -> None:
        new_fraction = new_mass / self.mass
        self.phases[phase_name] = new_fraction
        self.re_normalize_phases()

    def remove_phase_mass(self, phase_name: str, remove_mass: float) -> None:
        current_mass = self.get_phase_mass(phase_name)
        if remove_mass > current_mass:
            remove_mass = current_mass

        new_mass = current_mass - remove_mass
        self.set_phase_mass(phase_name, new_mass)

    def update_bulk_density(self, new_density: float) -> None:
        self.density = new_density

    def get_size_bin_mass(self, size_bin: str) -> float:
        fraction = self.size_distribution.get(size_bin, 0.0)
        return fraction * self.mass

    def set_size_bin_mass(self, size_bin: str, new_bin_mass: float) -> None:
        new_fraction = new_bin_mass / self.mass
        self.size_distribution[size_bin] = new_fraction
        self.re_normalize_size_dist()

