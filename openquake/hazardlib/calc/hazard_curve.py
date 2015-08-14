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
import sys
import numpy

from openquake.hazardlib.source.base import SeismicSourceCluster, \
    SeismicSourceClusterCollection
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
        src_clusts, sites, imtls, gsims, truncation_level,
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

    :param src_clusts:
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
    :param truncation_level:
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

    # This is for backward compatibility
    if not isinstance(src_clusts, SeismicSourceClusterCollection):
        cluster = SeismicSourceCluster(src_clusts, 1, 'indep', 'indep')
        src_clusts = SeismicSourceClusterCollection([cluster])

    # Compute hazard curves
    for cluster in src_clusts:
        tot_wei = 0.0

        # Source-site tuple
        sources_sites = ((source, sites) for source in cluster.clusters)

        # Temporary collection of curves used within the cluster of
        # sources
        if cluster.src_indep is 'indep':
            tc_clu = dict((imt, numpy.ones([len(sites), len(imtls[imt])]))
                          for imt in imtls)
        else:
            tc_clu = dict((imt, numpy.zeros([len(sites), len(imtls[imt])]))
                          for imt in imtls)

        # Processing sources in the cluster
        for source, s_sites in source_site_filter(sources_sites):

            try:

                # if cluster.rup_indep is 'indep':
                tc_src = dict((imt, numpy.ones([len(sites),
                              len(imtls[imt])])) for imt in imtls)
                # else:
                #    tc_src = dict((imt, numpy.zeros([len(sites),
                #                  len(imtls[imt])])) for imt in imtls)

                ruptures_sites = ((rupture, s_sites)
                                  for rupture in source.iter_ruptures())

                for rupture, r_sites in rupture_site_filter(ruptures_sites):
                    gsim = gsims[rupture.tectonic_region_type]
                    sctx, rctx, dctx = gsim.make_contexts(r_sites, rupture)
                    for imt in imts:
                        poes = gsim.get_poes(sctx, rctx, dctx, imt, imts[imt],
                                             truncation_level)
                        pno = rupture.get_probability_no_exceedance(poes)

                        # Update the object we use to cumulate contributions
                        # from the different ruptures in a single source
                        if cluster.rup_indep is 'indep':
                            tc_src[str(imt)] *= r_sites.expand(pno,
                                                               placeholder=1)
                        else:
                            tc_src[str(imt)] += r_sites.expand(pno,
                                                               placeholder=1)

                if cluster.src_indep is 'indep':
                    for imt in imtls:
                        tc_clu[imt] *= tc_src[imt]
                else:
                    for imt in imtls:
                        # Here we update the curves for a specific intensity
                        # measure type.
                        print 'Source ID:', source.source_id
                        tc_clu[imt] += (tc_src[imt] *
                                        cluster.weights[source.source_id])
                        
                        tot_wei += cluster.weights[source.source_id]
                del tc_src

            except Exception, err:
                etype, err, tb = sys.exc_info()
                msg = 'An error occurred with source id=%s. Error: %s'
                msg %= (source.source_id, err.message)
                raise etype, msg, tb

        # Update the probability of non-exceedance assuming that
        # SeismicSourceClusters are independent. Note that in the case of
        # mutually exclusive events with associated weights we must cope with
        # the case when some of the sources are filtered out.
        for imt in imtls:
            if tot_wei < 1. and cluster.src_indep is 'mutex':
                curves[imt] *= (tc_clu[imt] + 1. - tot_wei)
            else:
                curves[imt] *= tc_clu[imt]
        del tc_clu

    for imt in imtls:
        curves[imt] = 1. - curves[imt]

    return curves
