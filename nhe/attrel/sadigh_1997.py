"""
Module exports :class:`SadighEtAl1997`.
"""
from __future__ import division
from math import log, exp

from nhe.attrel.base import GMPE, CoeffsTable, SACoeffsTable
from nhe import const
from nhe.imt import PGA, SA


class SadighEtAl1997(GMPE):
    """
    Implements GMPE developed by Sadigh, K., C. -Y. Chang, J. A. Egan,
    F. Makdisi, and R. R. Youngs (1997) and published as "Attenuation
    relationships for shallow crustal earthquakes based on California
    strong motion data", Seismological Research Letters, 68(1), 180-189.
    """
    #: Supported tectonic region type is only active shallow crust,
    #: since data consists of California earthquakes mainly.
    DEFINED_FOR_TECTONIC_REGION_TYPES = set([
        const.TRT.ACTIVE_SHALLOW_CRUST
    ])

    #: Supported intensity measure types are spectral acceleration,
    #: and peak ground acceleration, see page 180.
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        SA
    ])

    #: Supported intensity measure component is the geometric mean
    #: of two horizontal components :attr:`~nhe.const.IMC.AVERAGE_HORIZONTAL`,
    #: see page 180.
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENTS = set([
        const.IMC.AVERAGE_HORIZONTAL
    ])

    #: Supported standard deviation type is only total, see table 3.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: Required site parameter is only Vs30 (used to distinguish rock
    #: and deep soil).
    REQUIRES_SITE_PARAMETERS = set(('vs30', ))

    #: Required rupture parameters are magnitude and rake (eq. 1).
    REQUIRES_RUPTURE_PARAMETERS = set(('rake', 'mag'))

    #: Required distance measure is RRup (eq. 1).
    REQUIRES_DISTANCES = set(('rrup', ))

    #: If site vs30 is more than 750 m/s -- treat the soil as rock.
    #: See page 180.
    ROCK_VS30 = 750

    #: Magnitude value to separate coefficients table because of near field
    #: saturation effect is 6.5. See page 184.
    NEAR_FIELD_SATURATION_MAG = 6.5

    def get_mean_and_stddevs(self, ctx, imt, stddev_types, component_type):
        """
        See :meth:`superclass method
        <nhe.attrel.base.AttenuationRelationship.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        assert component_type in self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENTS
        assert all(stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
                   for stddev_type in stddev_types)

        # GMPE differentiates strike-slip, reverse and normal ruptures,
        # but combines normal and strike-slip into one category. See page 180.
        is_reverse = (45 <= ctx.rup_rake <= 135)

        if ctx.site_vs30 > self.ROCK_VS30:
            # site is rock
            fmean = self._get_mean_rock
            fstddev = self._get_stddev_rock
        else:
            # site is deep soil
            fmean = self._get_mean_deep_soil
            fstddev = self._get_stddev_deep_soil

        mean = fmean(ctx.rup_mag, ctx.rup_rake, ctx.dist_rrup, is_reverse, imt)
        stddevs = [fstddev(ctx.rup_mag, imt) for stddev_type in stddev_types]
        return mean, stddevs

    def _get_mean_deep_soil(self, mag, rake, rrup, is_reverse, imt):
        """
        Calculate and return the mean intensity for deep soil sites.

        Implements an equation from table 4.
        """
        if mag <= self.NEAR_FIELD_SATURATION_MAG:
            c4 = self.COEFFS_SOIL_IMT_INDEPENDENT['c4lowmag']
            c5 = self.COEFFS_SOIL_IMT_INDEPENDENT['c5lowmag']
        else:
            c4 = self.COEFFS_SOIL_IMT_INDEPENDENT['c4himag']
            c5 = self.COEFFS_SOIL_IMT_INDEPENDENT['c5himag']
        c2 = self.COEFFS_SOIL_IMT_INDEPENDENT['c2']
        c3 = self.COEFFS_SOIL_IMT_INDEPENDENT['c3']
        C = self.COEFFS_SOIL[imt]
        if is_reverse:
            c1 = self.COEFFS_SOIL_IMT_INDEPENDENT['c1r']
            c6 = C['c6r']
        else:
            c1 = self.COEFFS_SOIL_IMT_INDEPENDENT['c1ss']
            c6 = C['c6ss']
        return (c1 + c2 * mag - c3 * log(rrup + c4 * exp(c5 * mag))
                + c6 + C['c7'] * ((8.5 - mag) ** 2.5))

    def _get_mean_rock(self, mag, rake, rrup, is_reverse, imt):
        """
        Calculate and return the mean intensity for rock sites.

        Implements an equation from table 2.
        """
        if mag <= self.NEAR_FIELD_SATURATION_MAG:
            C = self.COEFFS_ROCK_LOWMAG[imt]
        else:
            C = self.COEFFS_ROCK_HIMAG[imt]
        mean = (
            C['c1'] + C['c2'] * mag + C['c3'] * ((8.5 - mag) ** 2.5)
            + C['c4'] * log(rrup + exp(C['c5'] + C['c6'] * mag))
            + C['c7'] * log(rrup + 2)
        )
        if is_reverse:
            # footnote in table 2 says that for reverse ruptures
            # the mean amplitude value should be multiplied by 1.2
            mean += 0.1823215567939546  # == log(1.2)
        return mean

    def _get_stddev_rock(self, mag, imt):
        """
        Calculate and return total standard deviation for rock sites.

        Implements formulae from table 3.
        """
        C = self.COEFFS_ROCK_STDDERR[imt]
        if mag > C['maxmag']:
            return C['maxsigma']
        else:
            return C['sigma0'] + C['magfactor'] * mag

    def _get_stddev_deep_soil(self, mag, imt):
        """
        Calculate and return total standard deviation for deep soil sites.

        Implements formulae from the last column of table 4.
        """
        # footnote from table 4 says that stderr for magnitudes over 7
        # is equal to one of magnitude 7.
        if mag > 7:
            mag = 7
        C = self.COEFFS_SOIL[imt]
        return C['sigma0'] + C['magfactor'] * mag

    #: Coefficients tables for rock sites (table 2), for magnitude
    #: values of :attr:`NEAR_FIELD_SATURATION_MAG` and below. Damping
    #: for spectral acceleration here and in other SA-tables is 5%,
    #: see "introduction" section.
    #:
    #: Here we use lowest SA period of 0.075 s instead of 0.07 which
    #: is listed in table 2 for consistency with table 4.
    COEFFS_ROCK_LOWMAG = CoeffsTable("""\
        IMT        c1     c2    c3      c4      c5       c6      c7
        PGA       -0.624  1.0   0.000  -2.100   1.29649  0.250   0.0
    """).chain(SACoeffsTable("""\
    period damping c1     c2    c3      c4      c5       c6      c7
     0.075    5    0.110  1.0   0.006  -2.128   1.29649  0.250  -0.082
     0.10     5    0.275  1.0   0.006  -2.148   1.29649  0.250  -0.041
     0.20     5    0.153  1.0  -0.004  -2.080   1.29649  0.250   0.0
     0.30     5   -0.057  1.0  -0.017  -2.028   1.29649  0.250   0.0
     0.40     5   -0.298  1.0  -0.028  -1.990   1.29649  0.250   0.0
     0.50     5   -0.588  1.0  -0.040  -1.945   1.29649  0.250   0.0
     0.75     5   -1.208  1.0  -0.050  -1.865   1.29649  0.250   0.0
     1.0      5   -1.705  1.0  -0.055  -1.800   1.29649  0.250   0.0
     1.5      5   -2.407  1.0  -0.065  -1.725   1.29649  0.250   0.0
     2.0      5   -2.945  1.0  -0.070  -1.670   1.29649  0.250   0.0
     3.0      5   -3.700  1.0  -0.080  -1.610   1.29649  0.250   0.0
     4.0      5   -4.230  1.0  -0.100  -1.570   1.29649  0.250   0.0
    """))

    #: Coefficients tables for rock sites (table 2), for magnitude
    #: values above :attr:`NEAR_FIELD_SATURATION_MAG`.
    COEFFS_ROCK_HIMAG = CoeffsTable("""\
        IMT        c1     c2    c3      c4      c5       c6      c7
        PGA       -1.274  1.1   0.000  -2.100  -0.48451  0.524   0.0
    """).chain(SACoeffsTable("""\
    period damping c1     c2    c3      c4      c5       c6      c7
     0.075    5   -0.540  1.1   0.006  -2.128  -0.48451  0.524  -0.082
     0.10     5   -0.375  1.1   0.006  -2.148  -0.48451  0.524  -0.041
     0.20     5   -0.497  1.1  -0.004  -2.080  -0.48451  0.524   0.0
     0.30     5   -0.707  1.1  -0.017  -2.028  -0.48451  0.524   0.0
     0.40     5   -0.948  1.1  -0.028  -1.990  -0.48451  0.524   0.0
     0.50     5   -1.238  1.1  -0.040  -1.945  -0.48451  0.524   0.0
     0.75     5   -1.858  1.1  -0.050  -1.865  -0.48451  0.524   0.0
     1.0      5   -2.355  1.1  -0.055  -1.800  -0.48451  0.524   0.0
     1.5      5   -3.057  1.1  -0.065  -1.725  -0.48451  0.524   0.0
     2.0      5   -3.595  1.1  -0.070  -1.670  -0.48451  0.524   0.0
     3.0      5   -4.350  1.1  -0.080  -1.610  -0.48451  0.524   0.0
     4.0      5   -4.880  1.1  -0.100  -1.570  -0.48451  0.524   0.0
    """))

    #: Coefficient tables for standard error on rock sites (table 3).
    COEFFS_ROCK_STDDERR = CoeffsTable("""\
        IMT        sigma0  magfactor maxsigma maxmag
        PGA        1.39   -0.14      0.38     7.21
    """).chain(SACoeffsTable("""\
    period damping sigma0  magfactor maxsigma maxmag
     0.075    5    1.40   -0.14      0.39     7.21
     0.10     5    1.41   -0.14      0.40     7.21
     0.20     5    1.43   -0.14      0.42     7.21
     0.30     5    1.45   -0.14      0.44     7.21
     0.40     5    1.48   -0.14      0.47     7.21
     0.50     5    1.50   -0.14      0.49     7.21
     0.75     5    1.52   -0.14      0.51     7.21
     1.0      5    1.53   -0.14      0.52     7.21
     4.0      5    1.53   -0.14      0.52     7.21
    """))

    #: Coefficient tables for deep soil sites (table 4).
    COEFFS_SOIL = CoeffsTable("""\
        IMT        c6ss     c6r      c7     sigma0  magfactor maxmag
        PGA        0.0000   0.0000   0.0    1.52   -0.16      7
    """).chain(SACoeffsTable("""\
    period damping c6ss     c6r      c7     sigma0  magfactor maxmag
     0.075    5    0.4572   0.4572   0.005  1.54   -0.16      7
     0.1      5    0.6395   0.6395   0.005  1.54   -0.16      7
     0.2      5    0.9187   0.9187  -0.004  1.565  -0.16      7
     0.3      5    0.9547   0.9547  -0.014  1.58   -0.16      7
     0.4      5    0.9251   0.9005  -0.024  1.595  -0.16      7
     0.5      5    0.8494   0.8285  -0.033  1.61   -0.16      7
     0.75     5    0.7010   0.6802  -0.051  1.635  -0.16      7
     1.0      5    0.5665   0.5075  -0.065  1.66   -0.16      7
     1.5      5    0.3235   0.2215  -0.090  1.69   -0.16      7
     2.0      5    0.1001  -0.0526  -0.108  1.70   -0.16      7
     3.0      5   -0.2801  -0.4905  -0.139  1.71   -0.16      7
     4.0      5   -0.6274  -0.8907  -0.160  1.71   -0.16      7
     """))

    #: IMT-independent coefficients for deep soil sites (table 4).
    COEFFS_SOIL_IMT_INDEPENDENT = {
        'c1ss': -2.17,
        'c1r': -1.92,
        'c2': 1.0,
        'c3': 1.7,
        'c4lowmag': 2.1863,
        'c5lowmag': 0.32,
        'c4himag': 0.3825,
        'c5himag': 0.5882
    }
