"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:21 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2016-05-23 10:02:47 +0100 (Mon, 23 May 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

import math
from numba import jit

@jit(nopython = True)
def qScore(value1: float, value2: float):
    if value1 + value2 == 0: return 0 # otherwise is a ZeroDivision error.
    return math.sqrt(((value1 - value2) ** 2) / ((value1 + value2) ** 2))



def averageQScore(valueLists):
    if len(valueLists[0]) == 0: return 0 # otherwise is a ZeroDivision error.
    score = sum([qScore(valueList[0], valueList[1]) for valueList in valueLists]) / len(valueLists[0])
    return score



def euclidean(valueList):
    score = sum([(scoringValue[0] - scoringValue[1]) ** 2 for scoringValue in valueList])
    return math.sqrt(score)


functionDict = {
    'averageQScore': averageQScore,
    'euclidean': euclidean,
    }


def getNmrResidueMatches(queryShifts, matchNmrResiduesDict, scoringMethod, isotopeCode='13C'):
    scoringMatrix = {}
    for res, mShifts in matchNmrResiduesDict.items():

        if res.isDeleted or res._flaggedForDelete:
            continue

        scoringValues = []
        mShifts2 = [shift for shift in mShifts if shift and shift.nmrAtom.isotopeCode == isotopeCode]
        for mShift in mShifts2:
            qShifts2 = [shift for shift in queryShifts if shift and shift.nmrAtom.isotopeCode == isotopeCode]
            for qShift in queryShifts:
                if qShift and mShift and qShift != mShift:
                    if mShift.nmrAtom.name == qShift.nmrAtom.name:
                        scoringValues.append((mShift.value, qShift.value))
        if scoringValues and len(scoringValues) == len(qShifts2):
            score = functionDict[scoringMethod](scoringValues)
            scoringMatrix[score] = res

    return scoringMatrix
