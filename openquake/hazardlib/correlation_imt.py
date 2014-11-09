# The Hazard Library
# Copyright (C) 2012-2014, GEM Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Module :mod:`openquake.hazardlib.imt_correlation` defines correlation models
between different intensity measure types
"""

import abc
import numpy

from scipy.constants import pi

from openquake.hazardlib.imt import SA, PGA


class BaseIMTCorrelationModel(object):
    """
    Base class for correlation models for different intensity measure types.
    """
    __metaclass__ = abc.ABCMeta

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = abc.abstractproperty()

    @abc.abstractmethod
    def get_correlation_coefficient(self, imt_1, imt_2):
        """
        Calculate the correlation coefficient between the epsilon values of two 
        intensity measure types.

        :parameter imt_1:
            The first intensity measure type. Can be an instance of
            the ones defined in :mod:`openquake.hazardlib.imt`
        :parameter imt_2:
            The first intensity measure type. Can be an instance of
            the ones defined in :mod:`openquake.hazardlib.imt`
        :returns:
            The correlation coefficient
        """


class BJ2008cmEpsilonIMT(BaseIMTCorrelationModel):
    """
    This implements the correlation model proposed by Baker and Jayaram in
    2008 (Reference: Baker, J.W. and Jayaram, N. (2008), "Correlation of
    spectral acceleration values from NGA ground motion models," Earthquake
    Spectra, 24 (1), 299-317.
    """

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGA, SA])

    def get_correlation_coefficient(self, imt_1, imt_2):
        """
        """

        # Checks if the two IMTS are supported by the current correlation
        # model
        assert 

        # Compute the first coefficient 
        c1 = (1 - numpy.cos(pi/2 - log(T_max/max(T_min, 0.109)) * 0.366 ))

        if T_max < 0.2
            c2 = 1 - 0.105*(1 - 1./(1+exp(100*T_max-5)))*(T_max-T_min)/(T_max-0.0099)
        end
        if T_max < 0.109
            c3 = c2;
        else
            c3 = c1;
        end
        c4 = c1 + 0.5 * (sqrt(c3) - c3) * (1 + cos(pi*(T_min)/(0.109)));

        if T_max <= 0.109
            rho = c2;
        elseif T_min > 0.109
            rho = c1;
        elseif T_max < 0.2
            rho = min(c2, c4);
        else
            rho = c4;
        end

        

