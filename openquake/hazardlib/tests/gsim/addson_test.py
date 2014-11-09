
import numpy
import unittest

from openquake.hazardlib import const
from openquake.hazardlib.imt import SA
from openquake.hazardlib.gsim.addson import AvgSaCalc
from openquake.hazardlib.gsim.boore_atkinson_2008 import BooreAtkinson2008
from openquake.hazardlib.gsim.base import (RuptureContext, DistancesContext,
                                           SitesContext)


class MeanSACalculatorTest(unittest.TestCase):
    """
    """

    def setUp(self):

        self.gsim = BooreAtkinson2008()

        # We must define the correlation model

        self.rctx = RuptureContext()
        self.rctx.mag = 6.0
        self.rctx.rake = 0.0

        self.dctx = DistancesContext()
        self.dctx.rjb = numpy.array([10.])

        self.sctx = SitesContext()
        self.sctx.vs30 = numpy.array([800.])
        self.sctx.vs30measured = numpy.array([True])
        self.sctx.z1pt0 = numpy.array([100.])
        self.sctx.z2pt5 = numpy.array([1.])

        self.stddev_types = [const.StdDev.TOTAL]

    def get_mean_test1(self):
        """
        AAA
        """
        gsim = AvgSaCalc(self.gsim)

        print gsim.gsim
        print self.stddev_types

        lgm, std = gsim.get_mean_and_stddevs(self.sctx,
                                             self.rctx,
                                             self.dctx,
                                             SA(damping=5, period=0.1),
                                             self.stddev_types)
