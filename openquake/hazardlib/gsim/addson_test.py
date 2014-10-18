
import unittest

from openquake.hazardlib.gsim.addson import MeanSACalculator
from openquake.hazardlib.gsim.boore_atkinson_2008 import BooreAtkinson2008


class MeanSACalculatorTest(unittest.TestCase):
    """
    """

    def setUp(self):
        gsim = BooreAtkinson2008()
        calculator = MeanSACalculator(gsim)

