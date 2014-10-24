"""
Module :mod:`openquake.hazardlib.gsim.addson` provide additional
functionalities to the class :class:`ground shaking intensity models
<GroundShakingIntensityModel>`
"""

import numpy
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib.gsim.base import GMPE


class AvgSaCalc(GMPE):
    """
    :param gsim:
    :param t_low:
    :param t_upp:
    :param t_step:
    """

    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        PGV,
        SA
    ])
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL,
        const.StdDev.INTER_EVENT,
        const.StdDev.INTRA_EVENT
    ])
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL
    REQUIRES_SITES_PARAMETERS = set(('vs30', ))
    REQUIRES_RUPTURE_PARAMETERS = set(('rake', 'mag'))
    REQUIRES_DISTANCES = set(('rjb', ))

    def __init__(self, gsim, t_low=None, t_upp=None, t_step=1.0):

        super(AvgSaCalc, self).__init__()
        self.gsim = gsim

        self.t_low = t_low if t_low is not None else \
            self.gsim.COEFFS.sa_coeffs.keys()[0].period
        self.t_upp = t_upp if t_upp is not None else \
            self.gsim.COEFFS.sa_coeffs.keys()[-1].period
        self.t_step = t_step

        self.imt = []
        for per in numpy.arange(self.t_low, self.t_upp+self.t_step*0.1,
                                self.t_step):
            print per
            self.imt.append(SA(period=per, damping=5))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        This methos has a similar behaviour to the
        """

        if len(imt) < 1:
            imt = self.imt
