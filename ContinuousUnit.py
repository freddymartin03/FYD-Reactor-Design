class ContinuousUnit:
    """
    A mixin to define throughput (kg/h) and compute processing time for 'mass' of ore.
    """
    def __init__(self, throughput: float = 1000.0):
        """
        :param throughput: kg/h
        """
        self.throughput = throughput

    def compute_time(self, ore_mass: float) -> float:
        """
        Time (hours) to process 'ore_mass' at self.throughput (kg/h).
        """
        if self.throughput <= 0:
            return 0.0
        return ore_mass / self.throughput
