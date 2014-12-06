# coding: utf-8
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
:mod:`openquake.hazardlib.calc.hazard_curve` implements
:func:`hazard_curves`.
"""
import numpy

from openquake.hazardlib.source.base import SeismicSourceCluster
from openquake.hazardlib.calc import filters
from openquake.hazardlib.imt import from_string
from openquake.hazardlib.gsim.base import deprecated


@deprecated('Use calc_hazard_curves instead')
def hazard_curves(
        sources, sites, imts, gsims, truncation_level,
        source_site_filter=filters.source_site_noop_filter,
        rupture_site_filter=filters.rupture_site_noop_filter):
    """
    Deprecated. It does the same job of
    :func:`openquake.hazardlib.calc.hazard_curve.calc_hazard_curves`,
    with the only difference that the intensity measure types in input
    and output are hazardlib objects instead of simple strings.
    """
    imtls = {str(imt): imls for imt, imls in imts.iteritems()}
    curves_by_imt = calc_hazard_curves(
        sources, sites, imtls, gsims, truncation_level,
        source_site_filter=filters.source_site_noop_filter,
        rupture_site_filter=filters.rupture_site_noop_filter)
    return {from_string(imt): curves
            for imt, curves in curves_by_imt.iteritems()}


def calc_hazard_curves(
        sources, sites, imtls, gsims, truncation_level,
        source_site_filter=filters.source_site_noop_filter,
        rupture_site_filter=filters.rupture_site_noop_filter):
    """
    Compute hazard curves on a list of sites, given a set of seismic sources
    and a set of ground shaking intensity models (one per tectonic region type
    considered in the seismic sources).


    Probability of ground motion exceedance is computed using the following
    formula ::

        P(X≥x|T) = 1 - ∏ ∏ Prup_ij(X<x|T)

    where ``P(X≥x|T)`` is the probability that the ground motion parameter
    ``X`` is exceeding level ``x`` one or more times in a time span ``T``, and
    ``Prup_ij(X<x|T)`` is the probability that the j-th rupture of the i-th
    source is not producing any ground motion exceedance in time span ``T``.
    The first product ``∏`` is done over sources, while the second one is done
    over ruptures in a source.

    The above formula computes the probability of having at least one ground
    motion exceedance in a time span as 1 minus the probability that none of
    the ruptures in none of the sources is causing a ground motion exceedance
    in the same time span. The basic assumption is that seismic sources are
    independent, and ruptures in a seismic source are also independent.

    :param sources:
        An iterator of seismic sources objects (instances of subclasses
        of :class:`~openquake.hazardlib.source.base.BaseSeismicSource`).
    :param sites:
        Instance of :class:`~openquake.hazardlib.site.SiteCollection` object,
        representing sites of interest.
    :param imtls:
        Dictionary mapping intensity measure type strings
        to lists of intensity measure levels.
    :param gsims:
        Dictionary mapping tectonic region types (members
        of :class:`openquake.hazardlib.const.TRT`) to
        :class:`~openquake.hazardlib.gsim.base.GMPE` or
        :class:`~openquake.hazardlib.gsim.base.IPE` objects.
    :param trunctation_level:
        Float, number of standard deviations for truncation of the intensity
        distribution.
    :param source_site_filter:
        Optional source-site filter function. See
        :mod:`openquake.hazardlib.calc.filters`.
    :param rupture_site_filter:
        Optional rupture-site filter function. See
        :mod:`openquake.hazardlib.calc.filters`.

    :returns:
        Dictionary mapping intensity measure type strings (same keys
        as in parameter ``imtls``) to 2d numpy arrays of float, where
        first dimension differentiates sites (the order and length
        are the same as in ``sites`` parameter) and the second one
        differentiates IMLs (the order and length are the same as
        corresponding value in ``imts`` dict).
    """
    imts = {from_string(imt): imls for imt, imls in imtls.iteritems()}

    # This is the object where we cumulate the contributions
    curves = dict((imt, numpy.ones([len(sites), len(imtls[imt])]))
                  for imt in imtls)

    # This is for back compatibility
    if not isinstance(sources[0], SeismicSourceCluster):
        sources = [SeismicSourceCluster(sources, numpy.ones(len(sources)),
                                        'indep', 'indep')]





    return curves
