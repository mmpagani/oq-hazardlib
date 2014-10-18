"""
Module :mod:`openquake.hazardlib.gsim.addson` provide additional
functionalities to the class :class:`ground shaking intensity models
<GroundShakingIntensityModel>`
"""

import numpy
from openquake.hazardlib.gsim.base import (SitesContext, RuptureContext,
                                           DistancesContext)


class MeanSpectralIntensityCalculator(object):
    """
    :param gsim:
    :param t_low:
    :param t_upp:
    :param t_step:
    """

    def __init__(self, gsim, t_low=None, t_upp=None, t_step=0.1):
        self.gsim = gsim
        self.t_low = t_low if t_low is not None else 0.0
        self.t_upp = t_upp if t_upp is not None else \
            self.gsim.COEFFS.sa_coeffs.keys()[-1].period
        self.t_step = t_step

    def get_mean_and_stddevs(self, sites, rup, dists, stddev_types):
        """
        This methos has a similar behaviour to the
        """
        print self.gsim.COEFFS.sa_coeffs.keys()

        assert isinstance(sites, SitesContext)
        assert isinstance(rup, RuptureContext)
        assert isinstance(dists, DistancesContext)

        for imt in numpy.arange(self.t_low, self.t_upp, self.t_step):
            print imt
