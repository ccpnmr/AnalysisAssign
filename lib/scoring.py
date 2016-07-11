"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - : 2014-06-04 18:13:10 +0100 (Wed, 04 Jun 2014) $"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon P Skinner, Geerten W Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/license"
               "or ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/license"
                 " or ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author: skinnersp $"
__date__ = "$Date: 2016-05-23 10:02:47 +0100 (Mon, 23 May 2016) $"
__version__ = "$Revision: 9395 $"
#=========================================================================================
# Start of code
#=========================================================================================

import math

def qScore(value1:float, value2:float):
  return math.sqrt(((value1-value2)**2)/((value1+value2)**2))

def averageQScore(valueList):
  score = sum([qScore(scoringValue[0], scoringValue[1]) for scoringValue in valueList])/len(valueList)
  return score

def euclidean(valueList):
  score = sum([(scoringValue[0]-scoringValue[1])**2 for scoringValue in valueList])
  return math.sqrt(score)


def getNmrResidueMatches(queryShifts, matchNmrResiduesDict):
  scoringMatrix = {}
  scores = []
  isotopeCode = '13C'
  for res, mShifts in matchNmrResiduesDict.items():
    scoringValues = []
    mShifts2 = [shift for shift in mShifts if shift and shift.nmrAtom.isotopeCode == isotopeCode]
    for mShift in mShifts2:
      qShifts2 = [shift for shift in queryShifts if shift and shift.nmrAtom.isotopeCode == isotopeCode]
      for qShift in queryShifts:
          if qShift and mShift:
              if mShift.nmrAtom.name == qShift.nmrAtom.name:
                  scoringValues.append((mShift.value, qShift.value))
    if scoringValues and len(scoringValues) == len(qShifts2):
      score = euclidean(scoringValues)
      scoringMatrix[score] = res
      scores.append(score)

  return scoringMatrix, scores



