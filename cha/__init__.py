##############################################################################
# cha: Common Hits Approach
#
# Copyright 2010-2017 University Vienna, Department of Pharmaceutical Chemistry
#
# Authors: Marcus Wieder, Arthur Garon
# #
# cha is free software: you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# MIT License for more details.
#
# You should have received a copy of the MIT License along with pymbar.
##############################################################################

"""The common hits aproach.
"""

__author__ = "Marcus Wieder, Arthur Garon"
__license__ = "MIT"
__maintainer__ = "Marcus Wieder"
__email__ = "marcus.wieder@univie.ac.at"

from cha import generate_pdbs
from pymbar.mbar import MBAR
from pymbar.bar import BAR, BARzero
from pymbar.exp import EXP, EXPGauss
import pymbar.old_mbar

try:
    from pymbar import version
except:
    # Fill in information manually.
    # TODO: See if we can at least get the git revision info in here.
    version = 'dev'
    full_version = 'dev'
    git_revision = 'dev'
    isrelease = False
