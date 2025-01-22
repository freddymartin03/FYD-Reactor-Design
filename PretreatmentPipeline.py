from Ores import Ore
from CourseGrinder import CoarseCrusher
from Washer import Washer
from CrusherGrinder import CrusherGrinder
from Dryer import Dryer
from Calcinator import Calcinator

class PretreatmentPipeline:
    def __init__(self,
                 throughput: float,
                 coarse_crusher: CoarseCrusher = None,
                 washer: Washer = None,
                 dryer: Dryer = None,
                 final_crusher: CrusherGrinder = None,
                 calcinator: Calcinator = None):
        """
        :param throughput: overall feed rate (kg/h) for each step
        """
        self.throughput = throughput
        self.coarse_crusher = coarse_crusher
        self.washer = washer
        self.dryer = dryer
        self.final_crusher = final_crusher
        self.calcinator = calcinator

        # We'll store the time spent in each step
        self.time_coarse_crush = 0.0
        self.time_washing = 0.0
        self.time_drying = 0.0
        self.time_final_crush = 0.0
        self.time_calcination = 0.0

    def run_pipeline(self, ore: Ore) -> Ore:
        # 1) Coarse Crush
        if self.coarse_crusher:
            self.time_coarse_crush = self.compute_time(ore.mass)
            ore = self.coarse_crusher.process(ore)

        # 2) Washing
        if self.washer:
            self.time_washing = self.compute_time(ore.mass)
            ore = self.washer.process(ore)

        # 4) Drying
        if self.dryer:
            self.time_drying = self.compute_time(ore.mass)
            ore = self.dryer.process(ore)


        # 3) Final Crusher
        if self.final_crusher:
            self.time_final_crush = self.compute_time(ore.mass)
            ore = self.final_crusher.process(ore)

        # 5) Calcination
        if self.calcinator:
            self.time_calcination = self.compute_time(ore.mass)
            ore = self.calcinator.process(ore)

        return ore

    def compute_time(self, mass_kg: float) -> float:
        if self.throughput <= 0:
            return 0.0
        return mass_kg / self.throughput

    def total_time(self) -> float:
        return (self.time_coarse_crush + self.time_washing +
                self.time_final_crush + self.time_drying +
                self.time_calcination)
