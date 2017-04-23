"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2017"
__credits__ = ("Wayne Boucher, Ed Brooksbank, Rasmus H Fogh, Luca Mureddu, Timothy J Ragan"
               "Simon P Skinner & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for licence text")
__reference__ = ("For publications, please use reference from http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2017-04-07 11:40:21 +0100 (Fri, April 07, 2017) $"
__version__ = "$Revision: 3.0.b1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"

__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import typing
from collections import OrderedDict

from ccpn.AnalysisAssign.lib.scoring import getNmrResidueMatches
from ccpn.core.ChemicalShift import ChemicalShift
from ccpn.core.NmrResidue import NmrResidue
from ccpn.ui.gui.lib.SpectrumDisplay import makeStripPlot

from ccpn.ui.gui.lib.Strip import navigateToNmrAtomsInStrip, matchAxesAndNmrAtoms
from ccpn.ui.gui.lib.Strip import navigateToNmrResidueInDisplay

from ccpn.ui.gui.lib.Window import markPositions
from ccpn.ui.gui.modules.GuiStrip import GuiStrip
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTable, NmrResidueTableModule
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.CompoundWidgets import ListCompoundWidget
from ccpn.ui.gui.widgets.PulldownListsForObjects import ChemicalShiftListPulldown
from ccpn.ui.gui.widgets.MessageDialog import showWarning


from ccpn.core.lib.Notifiers import Notifier
from ccpn.ui.gui.lib.GuiNotifier import GuiNotifier
from ccpn.ui.gui.widgets.DropBase import DropBase

from ccpn.util.Logging import getLogger
logger = getLogger()

ALL = '<all>'


class BackboneAssignmentModule(NmrResidueTableModule):

  includeSettingsWidget = True
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsOnTop = True
  settingsMinimumSizes = (500, 200)

  className = 'BackboneAssignmentModule'

  def __init__(self, parent, application):

    super(BackboneAssignmentModule, self).__init__(parent=parent, application=application,
                                                   name='Backbone Assignment'
                                                  )
    # derive project, current, and mainWindow from application
    self.application = application

    self.nmrChains = self.application.project.nmrChains
    self.matchCheckBoxWidget = CheckBox(self.nmrResidueTable._widget,
                                        grid=(0,2), checked=True, text='Find matches')

    ### Settings ###

    # change some of the defaults setting inherited from NmrResidueTableModule
    self.sequentialStripsWidget.checkBox.setChecked(True)
    self.displaysWidget.addPulldownItem(0)

    minWidth = 150  # for labels of the compound widgets
    self.numberOfMatches = 5

    row = -1
    col = 1
    # Match module selection
    row += 1
    # cannot set a notifier for displays, as these are not (yet?) implemented
    self.matchWidget = ListCompoundWidget(self.settingsWidget, grid=(row,col), vAlign='top',
                                          minimumWidths=(minWidth, 0, 0),
                                          orientation='left',
                                          labelText="Match module(s):",
                                          texts=[display.pid for display in self.application.ui.mainWindow.spectrumDisplays]
                                         )

    # Chemical shift list selection
    row += 1
    self.shiftListWidget = ChemicalShiftListPulldown(self.settingsWidget, self.application.project,
                                                     grid=(row,col), vAlign='top', minimumWidths=(minWidth,0),
                                                     callback=self._setupShiftDicts, default=0
                                                    )
    self._setupShiftDicts()

    # for compatibility with previous implementation
    self.moduleList = self.matchWidget.listWidget

    self._stripNotifiers = []  # list to store GuiNotifiers for strips

  def _getDisplays(self):
    "return list of displays to navigate"
    displays = []
    dGids = self.displaysWidget.getTexts() # gid's of displays
    if len(dGids) == 0: return displays
    mGids = self.matchWidget.getTexts() # gid of the match displays
    if ALL in dGids:
        displays = [dp for dp in self.application.ui.mainWindow.spectrumDisplays if dp.pid not in mGids]
    else:
        displays = [self.application.getByGid(gid) for gid in dGids if (gid != ALL and gid not in mGids)]
    return displays

  def navigateToNmrResidue(self, nmrResidue, row=None, col=None):
    """
    Navigate in selected displays to nmrResidue; skip if no displays defined defined
    If matchCheckbox is checked, also call findAndDisplayMatches
    """
    displays = self._getDisplays()
    if len(displays) == 0:
      logger.warn('Undefined display module(s); select in settings first')
      showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
      return

    if self.matchCheckBoxWidget.isChecked() and len(self.matchWidget.getTexts()) == 0:
      logger.warn('Undefined match module; select in settings first or unselect "Find matches"')
      showWarning('startAssignment', 'Undefined match module;\nselect in settings first or unselect "Find matches"')
      return


    self.application._startCommandBlock(
        'BackboneAssignmentModule.navigateToResidue(project.getByPid(%r))' % nmrResidue.pid)
    try:
      # optionally clear the marks
      if self.autoClearMarksWidget.checkBox.isChecked():
        self.application.ui.mainWindow.clearMarks()

      # clear any notifiers of previous strips
      for notifier in self._stripNotifiers:
        notifier.unRegister()
        del(notifier)
      self._stripNotifiers = []

      nr = nmrResidue.mainNmrResidue
      # navigate the displays
      for display in displays:
        if len(display.strips) > 0:
          strips = navigateToNmrResidueInDisplay(nr, display, stripIndex=0,
                                      widths=['full']*len(display.strips[0].axisCodes),
                                      showSequentialResidues=(len(display.axisCodes) > 2) and
                                                             self.sequentialStripsWidget.checkBox.isChecked(),
                                      markPositions=self.markPositionsWidget.checkBox.isChecked()
                                      )
          # activate a callback notifiers; allow dropping onto the NmrResidueLabel
          for strip in strips:
            notifier = GuiNotifier(strip.planeToolbar.spinSystemLabel,
                                   [GuiNotifier.DROPEVENT], [DropBase.IDS],
                                   self._processDroppedNmrResidue, nmrResidue=nr)
            self._stripNotifiers.append(notifier)

      if self.matchCheckBoxWidget.isChecked():
        self.findAndDisplayMatches(nmrResidue)

      # update current (should trigger SequenceGraph)
      self.application.current.nmrResidue = nmrResidue
      self.application.current.nmrChain = nmrResidue.nmrChain

    finally:
      self.application._endCommandBlock()

  def findAndDisplayMatches(self, nmrResidue):
    "Find and displays the matches to nmrResidue"

    # If NmrResidue is a -1 offset NmrResidue, set queryShifts as value from self.interShifts dictionary
    # Set matchShifts as self.intraShifts
    if nmrResidue.sequenceCode.endswith('-1'):
      direction = '-1'
      iNmrResidue = nmrResidue.mainNmrResidue
      queryShifts = [shift for shift in self.interShifts[nmrResidue] if shift.nmrAtom.isotopeCode == '13C']
      matchShifts = self.intraShifts

    # If NmrResidue is not an offset NmrResidue, set queryShifts as value from self.intraShifts dictionary
    # Set matchShifts as self.interShifts
    else:
      direction = '+1'
      iNmrResidue = nmrResidue
      queryShifts = [shift for shift in self.intraShifts[nmrResidue] if shift.nmrAtom.isotopeCode == '13C']
      matchShifts = self.interShifts

    assignMatrix = getNmrResidueMatches(queryShifts, matchShifts, 'averageQScore')
    if not assignMatrix.values():
      logger.info('No matches found for NmrResidue: %s' % nmrResidue.pid)
      return
    self._createMatchStrips(assignMatrix)
    return

  def _processDroppedNmrResidue(self, data, nmrResidue):
    "Process the dropped NmrResidue id"

    droppedNmrResidue = None
    if DropBase.IDS in data and len(data[DropBase.IDS]) > 0:
      droppedNmrResidue = self.application.project.getByPid('NR:'+data[DropBase.IDS][0])
    if droppedNmrResidue is None:
      logger.info('Backbone assigned: invalid "id" of dropped item')
      raise Warning('Backbone assigned: invalid "id" of dropped item')

    logger.debug('nmrResidue:%s, droppedNmrResidue:%s', nmrResidue, droppedNmrResidue)
    if droppedNmrResidue == nmrResidue:
      logger.warning('Cannot connect residue to itself')
      return

    # silence the updat of the nmrResidueTable as we will to an explicit update later
    self.nmrResidueTable.setUpdateSilence(True)
    if data['shiftLeftMouse']:
      # leftShift drag; connect to previous
      nmrResidue.connectPrevious(droppedNmrResidue)
    else:
      nmrResidue.connectNext(droppedNmrResidue)
    self.nmrResidueTable.setUpdateSilence(False)

    self.nmrResidueTable.displayTableForNmrChain(droppedNmrResidue.nmrChain)
    self.navigateToNmrResidue(droppedNmrResidue)

  def _centreStripForNmrResidue(self, nmrResidue, strip):
    """
    Centre y-axis of strip based on chemical shifts of from NmrResidue.nmrAtoms
    """
    if not nmrResidue:
      logger.warn('No NmrResidue specified')
      return

    if not strip:
      logger.warn('No Strip specified')
      return

    yShifts = matchAxesAndNmrAtoms(strip, nmrResidue.nmrAtoms)[strip.axisOrder[1]]
    yShiftValues = [x.value for x in yShifts]
    yPosition = (max(yShiftValues) + min(yShiftValues))/2
    yWidth = max(yShiftValues)-min(yShiftValues)+10
    strip.orderedAxes[1].position = yPosition
    strip.orderedAxes[1].width = yWidth

  def _setupShiftDicts(self):
    """
    Creates two ordered dictionaries for the inter residue and intra residue CA and CB shifts for
    all NmrResidues in the project.
    """
    self.intraShifts = OrderedDict()
    self.interShifts = OrderedDict()
    chemicalShiftList = self.application.project.getByPid(self.shiftListWidget.pulldownList.currentText())

    for nmrResidue in self.application.project.nmrResidues:
      nmrAtoms = [nmrAtom for nmrAtom in nmrResidue.nmrAtoms]
      shifts = [chemicalShiftList.getChemicalShift(atom.id) for atom in nmrAtoms]
      if nmrResidue.sequenceCode.endswith('-1'):
        self.interShifts[nmrResidue] = shifts
      else:
        self.intraShifts[nmrResidue] = shifts

  def _createMatchStrips(self, assignMatrix:typing.Tuple[typing.Dict[NmrResidue, typing.List[ChemicalShift]], typing.List[float]]):
    """
    Creates strips in match module corresponding to the best assignment possibilities
    in the assignMatrix.
    """
    if not assignMatrix:
      logger.warn('No assignment matrix specified')
      return

    # Assignment score has format {score: nmrResidue} where score is a float
    # assignMatrix[0] is a dict {score: nmrResidue} assignMatrix[1] is a concurrent list of scores
    assignmentScores = sorted(list(assignMatrix.keys()))[0:self.numberOfMatches]
    nmrAtomPairs = []
    for assignmentScore in assignmentScores:
      matchResidue = assignMatrix[assignmentScore]
      if matchResidue.sequenceCode.endswith('-1'):
        iNmrResidue = matchResidue.mainNmrResidue
      else:
        iNmrResidue = matchResidue
      nmrAtomPairs.append((iNmrResidue.fetchNmrAtom(name='N'), iNmrResidue.fetchNmrAtom(name='H')))

    for modulePid in self.matchWidget.getTexts():
      module = self.application.project.getByPid(modulePid)
      makeStripPlot(module, nmrAtomPairs)

      for ii, strip in enumerate(module.strips):
        nmrResidueId = nmrAtomPairs[ii][0].nmrResidue._id
        strip.planeToolbar.spinSystemLabel.setText(nmrResidueId)

      self._centreStripForNmrResidue(assignMatrix[assignmentScores[0]], module.strips[0])

  def _closeModule(self):
    """
    Re-implementation of the closeModule method of the CcpnModule class required 
    """
    # TODO: use proper subclassing
    for notifier in self._stripNotifiers:
      notifier.unRegister()
    self._stripNotifiers = []
    self.close()



#=====  Just some code to 'save' =====
def hasNmrResidue(nmrChain, residueCode):
  "Simple function to check if sequenCode is found within the nmrResidues of nmrChain"
  resCodes = [res.sequenceCode for res in nmrChain.nmrResidues]
  return (residueCode in resCodes)


def endOfchain(nmrResidue):
  # changes to end of connected chain; not a good idea
  if nmrResidue.nmrChain.isConnected:
    if nmrResidue.sequenceCode.endswith('-1'):
      nmrResidue = nmrResidue.nmrChain.mainNmrResidues[0].getOffsetNmrResidue(-1)
    else:
      nmrResidue = nmrResidue.nmrChain.mainNmrResidues[-1]
  return nmrResidue


def getPids(fromObject, attributeName):
  "Get a list of pids fromObject.attributeName or None on error"
  if not hasattr(fromObject, attributeName): return None
  return [obj.pid for obj in getattr(fromObject, attributeName)]
#===== end code save =====