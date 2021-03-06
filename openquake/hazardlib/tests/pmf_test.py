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
from decimal import Decimal

import unittest

from openquake.hazardlib.pmf import PMF


class PMFTestCase(unittest.TestCase):
    def test_creation(self):
        pmf = PMF((Decimal('0.1'), i) for i in xrange(10))
        self.assertEqual(pmf.data, [(0.1, i) for i in xrange(10)])

    def test_wrong_sum(self):
        data = [(0.1, i) for i in xrange(10)]
        self.assertRaises(ValueError, PMF, data, 1E-16)

    def test_empty_data(self):
        self.assertRaises(ValueError, PMF, [])

    def test_negative_or_zero_prob(self):
        # negative probs are refused
        data = [(-1, 0)] + [(Decimal('1.0'), 1), (Decimal('1.0'), 2)]
        self.assertRaises(ValueError, PMF, data)

        # 0 probs are accepted
        data = [(0, 0)] + [(Decimal('0.5'), 1), (Decimal('0.5'), 2)]
        PMF(data)
