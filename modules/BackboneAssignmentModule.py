"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date$"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Geerten W Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/license"
              "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/license"
                " or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author: Geerten Vuister $"
__date__ = "$Date: 2017-04-18 15:19:26 +0100 (Tue, April 18, 2017) $"

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
from ccpn.ui.gui.DropBase import DropBase

from ccpn.util.Logging import getLogger
logger = getLogger()

ALL = '<all>'


class BackboneAssignmentModule(NmrResidueTableModule):

  includeSettingsWidget = True
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsOnTop = True

  className = 'BackboneAssignmentModule'

  def __init__(self, parent=None):

    super(BackboneAssignmentModule, self).__init__(parent=parent, name='Backbone Assignment',
                                                   callback=self.navigateToNmrResidue)
    # project, current, application and mainWindow are inherited from
    # BackBoneAssignmentModule/CcpnModule

    self.nmrChains = self.project.nmrChains
    self.matchCheckBoxWidget = CheckBox(self.nmrResidueTable, grid=(0,2), checked=False, text='Find matches')

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
                                          texts=[display.pid for display in self.mainWindow.spectrumDisplays]
                                         )

    # Chemical shift list selection
    row += 1
    self.shiftListWidget = ChemicalShiftListPulldown(self.settingsWidget, self.project,
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
        displays = [dp for dp in self.mainWindow.spectrumDisplays if dp.pid not in mGids]
    else:
        displays = [self.application.getByGid(gid) for gid in dGids if (gid != ALL and gid not in mGids)]
    return displays

  def navigateToNmrResidue(self, nmrResidue, row=None, col=None):
    """
    Navigate in selected displays to nmrResidue; skip if no displays defined defined
    If matchCheckbox is checked, also call findAndDisplayMatches
    """
    displays = self._getDisplays()
    if len(displays) == 0: return

    self.application._startCommandBlock(
        'BackboneAssignmentModule.navigateToResidue(project.getByPid(%r)' % nmrResidue.pid)
    try:
      # optionally clear the marks
      if self.autoClearMarksWidget.checkBox.isChecked():
        self.mainWindow.clearMarks()

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
      self.current.nmrResidue = nmrResidue
      self.current.nmrChain = nmrResidue.nmrChain

    finally:
      self.application._endCommandBlock()

  def findAndDisplayMatches(self, nmrResidue):
    "Find and displays the matches to nmrResidue"

    if len(self.matchWidget.getTexts()) == 0:
      logger.warn('Undefined match module; select in settings first')
      showWarning('startAssignment', 'Undefined match module; select in settings first')
      return

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

    #print('_processDroppedNmrResidue>>', data)
    droppedNmrResidue = self.project.getByPid('NR:'+data['ids'][0])
    print('>_processDroppedNmrResidue>', nmrResidue, droppedNmrResidue)
    if droppedNmrResidue == nmrResidue:
      logger.warning('Cannot connect residue to itself')
      return

    if data['shiftLeftMouse']:
      # leftShift drag; connect to previous
      nmrResidue.connectPrevious(droppedNmrResidue)
    else:
      nmrResidue.connectNext(droppedNmrResidue)

    # this should trigger the update of the Sequence Graph
    self.current.nmrResidue = droppedNmrResidue
    self.current.nmrChain = droppedNmrResidue.nmrChain
    self._updateNmrResidueTable(droppedNmrResidue.nmrChain)

  # def _updateNmrChainList(self, nmrChain):
  #   """
  #   Convenience function for notifiers to update the NmrResidueTable when notifier is called in
  #   response to creation, deletion and changes to NmrChain objects.
  #   """
  #   if not nmrChain:
  #     logger.warn('No NmrChain specified')
  #     return
  #   self.nmrResidueTable.nmrResidueTable.objectLists.append(nmrChain)

  # def _updateNmrChainPulldown(self, nmrChain):
  #   """
  #   Convenience function for notifiers to update the NmrResidueTable when notifier is called in
  #   response to creation, deletion and changes to NmrChain objects.
  #   """
  #   if not nmrChain:
  #     logger.warn('No NmrChain specified')
  #     return
  #   self.nmrResidueTable.nmrResidueTable.objectLists = self.project.nmrChains
  #   self.nmrResidueTable.nmrResidueTable._updateSelectorContents()

  def _updateNmrResidueTable(self, nmrChain):
    """
    Convenience function to update the nmrResidueTable.
    """
    if not nmrChain:
      logger.warn('No NmrChain specified')
      return
#    self.nmrResidueTable.selectNmrChain(nmrChain)
    self.nmrResidueTable.ncWidget.pulldownList.select(nmrChain)
    self.nmrResidueTable.nmrResidueTable.objectLists = self.project.nmrChains
    self.nmrResidueTable.nmrResidueTable._updateSelectorContents()
#    self.nmrResidueTable.updateTable()

  # def _matchModules(self):
  #
  #   return [self.moduleList.item(i).text() for i in range(self.moduleList.count())]

  # def _selectMatchModule(self, text):
  #   """
  #   Call back to assign modules as match modules in response to a signal from the modulePulldown
  #   above. Adds the item to a list containing match modules and adds the Pid of the module to the
  #   moduleList ListWidget object.
  #   """
  #   if not text:  # blank row
  #     return
  #   if text not in self.matchWidget.getTexts():
  #     self.moduleList.addItem(text)
  #   self.modulePulldown.setIndex(0)

  # def _startAssignment(self, nmrResidue:NmrResidue, row:int=None, col:int=None):
  #   """
  #   Initiates assignment procedure when triggered by selection of an NmrResidue from the nmrResidueTable
  #   inside the module.
  #   """
  #   if not nmrResidue:
  #     logger.warn('No NmrResidue specified')
  #     return
  #   if len(self.matchWidget.getTexts()) == 0:
  #     logger.warn('Undefined match module; select in settings first')
  #     showWarning('startAssignment', 'Undefined match module; select in settings first')
  #     return
  #
  #   self.project._startFunctionCommandBlock('_startAssignment', nmrResidue)
  #   try:
  #     self._setupShiftDicts()
  #
  #     # if hasattr(self, 'sequenceGraph'):
  #     #   self.sequenceGraph.clearAllItems()
  #     #   self.sequenceGraph.nmrChainPulldown.select(self.current.nmrChain.pid)
  #
  #     if nmrResidue.nmrChain.isConnected:
  #       if nmrResidue.sequenceCode.endswith('-1'):
  #         nmrResidue = nmrResidue.nmrChain.mainNmrResidues[0].getOffsetNmrResidue(-1)
  #       else:
  #         nmrResidue = nmrResidue.nmrChain.mainNmrResidues[-1]
  #
  #     self._navigateTo(nmrResidue, row, col)
  #     # update current (should trigger SequenceGraph)
  #     self.current.nmrResidue = nmrResidue
  #     self.current.nmrChain = nmrResidue.nmrChain
  #
  #   finally:
  #     self.project._appBase._endCommandBlock()
  #
  #
  # def _navigateTo(self, nmrResidue:NmrResidue, row:int=None, col:int=None, strip:GuiStrip=None):
  #   """
  #   Takes an NmrResidue and an optional GuiStrip and changes z position(s) of all available displays
  #   to chemical shift value NmrAtoms in the NmrResidue. Takes corresponding value from inter-residual
  #   or intra-residual chemical shift dictionaries, using the NmrResidue pid as the key.
  #   Determines which nmrResidue(s) match the query NmrResidue and creates up to five strips, one for
  #   each of the match NmrResidue, and marks the carbon positions.
  #   """
  #
  #   if not nmrResidue:
  #     logger.warn('No NmrResidue specified')
  #     return
  #
  #   self.project._startFunctionCommandBlock('_navigateTo', nmrResidue, strip)
  #   try:
  #     mainWindow = self.mainWindow
  #     mainWindow.clearMarks()
  #     self.nmrResidueTable.nmrResidueTable.updateTable()
  #     selectedDisplays = [display for display in self.project.spectrumDisplays
  #                         if display.pid not in self.matchWidget.getTexts()]
  #
  #
  #     # If NmrResidue is a -1 offset NmrResidue, set queryShifts as value from self.interShifts dictionary
  #     # Set matchShifts as self.intraShifts
  #     if nmrResidue.sequenceCode.endswith('-1'):
  #       direction = '-1'
  #       iNmrResidue = nmrResidue.mainNmrResidue
  #       queryShifts = [shift for shift in self.interShifts[nmrResidue] if shift.nmrAtom.isotopeCode == '13C']
  #       matchShifts = self.intraShifts
  #
  #     # If NmrResidue is not an offset NmrResidue, set queryShifts as value from self.intraShifts dictionary
  #     # Set matchShifts as self.interShifts
  #     else:
  #       direction = '+1'
  #       iNmrResidue = nmrResidue
  #       queryShifts = [shift for shift in self.intraShifts[nmrResidue] if shift.nmrAtom.isotopeCode == '13C']
  #       matchShifts = self.interShifts
  #
  #     self.current.nmrResidue = iNmrResidue
  #     # If a strip is not specified, use the first strip in the each of the spectrumDisplays in selectedDisplays.
  #     if not strip:
  #       strips = [display.strips[0] for display in selectedDisplays]
  #     else:
  #       strips = [strip]
  #     for strip in strips:
  #       self._displayNmrResidueInStrip(iNmrResidue, strip)
  #       if len(strip.axisCodes) > 2:
  #         self._centreStripForNmrResidue(nmrResidue, strip)
  #       amidePair = [iNmrResidue.fetchNmrAtom(name='N'), iNmrResidue.fetchNmrAtom(name='H')]
  #       carbonAtoms = [x for x in nmrResidue.nmrAtoms if x.isotopeCode == '13C']
  #       axisCodePositionDict = matchAxesAndNmrAtoms(strip, nmrAtoms=set((amidePair+carbonAtoms)))
  #       markPositions(self.project, list(axisCodePositionDict.keys()), list(axisCodePositionDict.values()))
  #
  #     assignMatrix = getNmrResidueMatches(queryShifts, matchShifts, 'averageQScore')
  #     if not assignMatrix.values():
  #       logger.info('No matches found for NmrResidue: %s' % nmrResidue.pid)
  #       return
  #     self._createMatchStrips(assignMatrix)
  #
  #     # if hasattr(self, 'sequenceGraph'):
  #     #   if self.sequenceGraph.nmrChainPulldown.currentText() != nmrResidue.nmrChain.pid:
  #     #     self.sequenceGraph.nmrChainPulldown.select(nmrResidue.nmrChain.pid)
  #     #   elif not nmrResidue.nmrChain.isConnected:
  #     #     self.sequenceGraph.addResidue(iNmrResidue, direction)
  #     #   else:
  #     #     self.sequenceGraph.setNmrChainDisplay(nmrResidue.nmrChain.pid)
  #   finally:
  #     self.project._appBase._endCommandBlock()

  # def _displayNmrResidueInStrip(self, nmrResidue, strip):
  #   """
  #   navigate strip position to position specified by nmrResidue and set spinSystemLabel to nmrResidue id
  #   """
  #   if not nmrResidue:
  #     logger.warn('No NmrResidue specified')
  #     return
  #
  #   if not strip:
  #     logger.warn('No Strip specified')
  #     return
  #
  #   navigateToNmrAtomsInStrip(strip=strip, nmrAtoms=nmrResidue.nmrAtoms, widths=['default']*len(strip.axisCodes))
  #   strip.planeToolbar.spinSystemLabel.setText(nmrResidue._id)

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
    chemicalShiftList = self.project.getByPid(self.shiftListWidget.pulldownList.currentText())

    for nmrResidue in self.project.nmrResidues:
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
      module = self.project.getByPid(modulePid)
      makeStripPlot(module, nmrAtomPairs)

      for ii, strip in enumerate(module.strips):
        nmrResidueId = nmrAtomPairs[ii][0].nmrResidue._id
        strip.planeToolbar.spinSystemLabel.setText(nmrResidueId)

      self._centreStripForNmrResidue(assignMatrix[assignmentScores[0]], module.strips[0])

  # def _connectSequenceGraph(self, sequenceGraph:CcpnModule):
  #   """
  #   # CCPN INTERNAL - called in showSequenceGraph method of GuiMainWindow.
  #   Connects Sequence Graph to this module.
  #   """
  #   self.sequenceGraph = sequenceGraph
  #   self.project._appBase.current.assigner = sequenceGraph
  #   self.sequenceGraph.nmrResidueTable = self.nmrResidueTable
  #   self.sequenceGraph.setMode('fragment')

  def _closeModule(self):
    """
    Re-implementation of the closeModule method of the CcpnModule class required 
    """
    # TODO: use proper subclassing
    for notifier in self._notifiers:
      notifier.unRegister()
    self._notifiers = []
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