"""
Module :mod:`openquake.hazardlib.gsim.addson` provide additional
functionalities to the class :class:`ground shaking intensity models
<GroundShakingIntensityModel>`
"""

import numpy
from scipy.constants import pi
from openquake.hazardlib import const
from openquake.hazardlib.imt import AvgSA, SA
from openquake.hazardlib.gsim.base import GMPE


class AvgSaCalc(GMPE):
    """
    :param gsim:
    :param t_low:
    :param t_upp:
    :param t_step:
    """

    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([AvgSA])
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([const.StdDev.TOTAL])
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL
    REQUIRES_SITES_PARAMETERS = set(('vs30', ))
    REQUIRES_RUPTURE_PARAMETERS = set(('rake', 'mag'))
    REQUIRES_DISTANCES = set(('rjb', ))

    def __init__(self, gsim, periods):

        assert len(periods) > 0

        super(AvgSaCalc, self).__init__()
        self.gsim = gsim
        self.periods = periods

        # Create the vector of IMTs

        self.imts = []
        for per in periods:
            self.imts.append(SA(period=per, damping=5))
            
    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        This methos has a similar behaviour to the
        """

        # Compute median and standard deviation for each IMT
        stddev_types = [const.StdDev.TOTAL]
        means = numpy.zeros(len(self.imts))
        stds = numpy.zeros(len(self.imts))
        for i, imt in enumerate(self.imts):
            mean, std = self.gsim.get_mean_and_stddevs(sites, rup, dists,
                                                       imt, stddev_types)
            means[i] = mean
            stds[i] = std[0]
            
        moc = self._get_correlation_matrix()
        avgmean = 1./(len(self.periods))*sum(means)
        avgstd = 0.
        for i, per1 in enumerate(self.periods):
            for j, per2 in enumerate(self.periods):
                avgstd = (avgstd+(1./len(self.periods))**2 *
                          (moc[i, j]*stds[i]*stds[j]))

        print avgmean, avgstd
        return numpy.array(avgmean), [avgstd]

    def _get_correlation_matrix(self):
        moc = numpy.zeros((len(self.periods),len(self.periods)))
        for i, per1 in enumerate(self.periods):
            for j, per2 in enumerate(self.periods):
                rho = jb_correlation(per1, per2)
                moc[i,j] = rho
        return moc
                
#
#    
def jb_correlation(t1, t2):
    """
    parameter 
    """
    T_min = min(t1, t2)
    T_max = max(t1, t2)

    C1 = (1-numpy.cos(pi/2 - numpy.log(T_max/max(T_min, 0.109)) * 0.366 ))
    
    if T_max < 0.2:
        C2 = (1 - 0.105*(1 - 1./(1+numpy.exp(100*T_max-5)))*
              (T_max-T_min)/(T_max-0.0099))
    
    if T_max < 0.109:
        C3 = C2
    else:
        C3 = C1
    
    C4 = (C1 + 0.5 * (numpy.sqrt(C3) - C3) * 
            (1 + numpy.cos(pi*(T_min)/(0.109))))

    if T_max <= 0.109:
        rho = C2
    elif T_min > 0.109:
        rho = C1
    elif T_max < 0.2:
        rho = min(C2, C4)
    else:
        rho = C4
    return rho
