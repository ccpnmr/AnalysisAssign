"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date$"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon P Skinner, Geerten W Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/license"
              "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/license"
                " or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Revision$"

#=========================================================================================
# Start of code
#=========================================================================================

from collections import OrderedDict
from functools import partial

import typing

from ccpn.AnalysisAssign.lib.scoring import getNmrResidueMatches

from ccpn.core.ChemicalShift import ChemicalShift
from ccpn.core.NmrResidue import NmrResidue

from ccpn.ui.gui.lib.SpectrumDisplay import makeStripPlot
from ccpn.ui.gui.lib.Strip import navigateToNmrAtomsInStrip, matchAxesAndNmrAtoms
from ccpn.ui.gui.lib.Window import markPositions



from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTable
from ccpn.ui.gui.modules.GuiStrip import GuiStrip

from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.Module import CcpnModule
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget



class BackboneAssignmentModule(CcpnModule):

  includeSettingsWidget = True

  def __init__(self, parent, project):

    super(BackboneAssignmentModule, self).__init__(parent=None, name='Backbone Assignment', logger=project._logger)

    self.project = project
    self.current = project._appBase.current
    self.numberOfMatches = 5
    self.nmrChains = project.nmrChains
    self.matchModules = []
    self.nmrResidueTable = NmrResidueTable(self.mainWidget, project, callback=self._startAssignment)
    self.parent = parent

    self.placeSettingsButton(buttonParent=self.nmrResidueTable, buttonGrid=(0, 5))

    self.layout.addWidget(self.nmrResidueTable, 0, 0, 1, 3)

    modules = [display.pid for display in project.spectrumDisplays]
    modules.insert(0, '  ')
    modulesLabel = Label(self.settingsWidget, text="Selected Modules", grid=(1, 0))
    self.modulePulldown = PulldownList(self.settingsWidget, grid=(2, 0), callback=self._selectMatchModule)
    self.modulePulldown.setData(modules)
    self.moduleList = ListWidget(self.settingsWidget, grid=(3, 0), )
    chemicalShiftListLabel = Label(self.settingsWidget, text='Chemical Shift List', grid=(1, 1))
    self.chemicalShiftListPulldown = PulldownList(self.settingsWidget, grid=(2, 1), callback=self._setupShiftDicts)
    self.chemicalShiftListPulldown.setData([shiftlist.pid for shiftlist in project.chemicalShiftLists])
    self.project.registerNotifier('NmrResidue', 'rename', self._updateNmrResidueTable)
    self.project.registerNotifier('NmrChain', 'delete', self._updateNmrChainPulldown)
    self.project.registerNotifier('NmrChain', 'create', self._updateNmrChainList)


  def _updateNmrChainList(self, nmrChain):
    """
    Convenience function for notifiers to update the NmrResidueTable when notifier is called in
    response to creation, deletion and changes to NmrChain objects.
    """
    if not nmrChain:
      self.project._logger.warn('No NmrChain specified')
      return

    self.nmrResidueTable.nmrResidueTable.objectLists.append(nmrChain)

  def _updateNmrChainPulldown(self, nmrChain):
    """
    Convenience function for notifiers to update the NmrResidueTable when notifier is called in
    response to creation, deletion and changes to NmrChain objects.
    """
    if not nmrChain:
      self.project._logger.warn('No NmrChain specified')
      return

    self.nmrResidueTable.nmrResidueTable.objectLists = self.project.nmrChains
    self.nmrResidueTable.nmrResidueTable._updateSelectorContents()

  def _updateNmrResidueTable(self, nmrResidue, oldPid=None):
    """
    Convenience function for notifiers to update the NmrResidueTable when notifier is called in
    response to creation, deletion and changes to the current.nmrResidue object.
    """

    if not nmrResidue:
      self.project._logger.warn('No NmrResidue specified')
      return

    if nmrResidue == self.current.nmrResidue:
      self.nmrResidueTable.nmrResidueTable._updateSelectorContents()
      self.nmrResidueTable.nmrResidueTable.selector.select(nmrResidue.nmrChain)
      self.nmrResidueTable.updateTable()


  def _selectMatchModule(self, item):
    """
    Call back to assign modules as match modules in reponse to a signal from the modulePulldown
    above. Adds the item to a list containing match modules and adds the Pid of the module to the
    moduleList ListWidget object.
    """
    self.moduleList.addItem(item)
    self.matchModules.append(item)


  def _startAssignment(self, nmrResidue:NmrResidue, row:int=None, col:int=None):
    """
    Initiates assignment procedure when triggered by selection of an NmrResidue from the nmrResidueTable
    inside the module.
    """
    if not nmrResidue:
      self.project._logger.warn('No NmrResidue specified')
      return

    self.project._startFunctionCommandBlock('_startAssignment', nmrResidue)
    try:
      self._setupShiftDicts()

      self.current.nmrChain = nmrResidue.nmrChain
      if hasattr(self, 'sequenceGraph'):
        self.sequenceGraph.clearAllItems()
        self.sequenceGraph.nmrChainPulldown.select(self.current.nmrChain.pid)
      if self.current.nmrChain.isConnected:
        if nmrResidue.sequenceCode.endswith('-1'):
          nmrResidue = self.current.nmrChain.mainNmrResidues[0].getOffsetNmrResidue(-1)
        else:
          nmrResidue = self.current.nmrChain.mainNmrResidues[-1]

      self._navigateTo(nmrResidue, row, col)

    finally:
      self.project._appBase._endCommandBlock()


  def _navigateTo(self, nmrResidue:NmrResidue, row:int=None, col:int=None, strip:GuiStrip=None):
    """
    Takes an NmrResidue and an optional GuiStrip and changes z position(s) of all available displays
    to chemical shift value NmrAtoms in the NmrResidue. Takes corresponding value from inter-residual
    or intra-residual chemical shift dictionaries, using the NmrResidue pid as the key.
    Determines which nmrResidue(s) match the query NmrResidue and creates up to five strips, one for
    each of the match NmrResidue, and marks the carbon positions.
    """

    if not nmrResidue:
      self.project._logger.warn('No NmrResidue specified')
      return

    self.project._startFunctionCommandBlock('_navigateTo', nmrResidue, strip)
    try:
      if self.project._appBase.ui.mainWindow is not None:
        mainWindow = self.project._appBase.ui.mainWindow
      else:
        mainWindow = self.project._appBase._mainWindow
      mainWindow.clearMarks()
      self.nmrResidueTable.nmrResidueTable.updateTable()
      selectedDisplays = [display for display in self.project.spectrumDisplays
                          if display.pid not in self.matchModules]


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

      self.current.nmrResidue = iNmrResidue
      # If a strip is not specified, use the first strip in the each of the spectrumDisplays in selectedDisplays.
      if not strip:
        strips = [display.strips[0] for display in selectedDisplays]
      else:
        strips = [strip]
      for strip in strips:
        self._displayNmrResidueInStrip(iNmrResidue, strip)
        if len(strip.axisCodes) > 2:
          self._centreStripForNmrResidue(nmrResidue, strip)
        amidePair = [iNmrResidue.fetchNmrAtom(name='N'), iNmrResidue.fetchNmrAtom(name='H')]
        carbonAtoms = [x for x in nmrResidue.nmrAtoms if x.isotopeCode == '13C']
        axisCodePositionDict = matchAxesAndNmrAtoms(strip, nmrAtoms=set((amidePair+carbonAtoms)))
        markPositions(self.project, list(axisCodePositionDict.keys()), list(axisCodePositionDict.values()))

      assignMatrix = getNmrResidueMatches(queryShifts, matchShifts, 'averageQScore')
      if not assignMatrix.values():
        self.project._logger.info('No matches found for NmrResidue: %s' % nmrResidue.pid)
        return
      self._createMatchStrips(assignMatrix)

      if hasattr(self, 'sequenceGraph'):
        if self.sequenceGraph.nmrChainPulldown.currentText() != nmrResidue.nmrChain.pid:
          self.sequenceGraph.nmrChainPulldown.select(nmrResidue.nmrChain.pid)
        elif not nmrResidue.nmrChain.isConnected:
          self.sequenceGraph.addResidue(iNmrResidue, direction)
        else:
          self.sequenceGraph.setNmrChainDisplay(nmrResidue.nmrChain.pid)
    finally:
      self.project._appBase._endCommandBlock()



  def _displayNmrResidueInStrip(self, nmrResidue, strip):
    """
    navigate strip position to position specified by nmrResidue and set spinSystemLabel to nmrResidue id
    """
    if not nmrResidue:
      self.project._logger.warn('No NmrResidue specified')
      return

    if not strip:
      self.project._logger.warn('No Strip specified')
      return

    navigateToNmrAtomsInStrip(strip=strip, nmrAtoms=nmrResidue.nmrAtoms, widths=['default']*len(strip.axisCodes))
    strip.planeToolbar.spinSystemLabel.setText(nmrResidue._id)

  def _centreStripForNmrResidue(self, nmrResidue, strip):
    """
    Centre y-axis of strip based on chemical shifts of from NmrResidue.nmrAtoms
    """
    if not nmrResidue:
      self.project._logger.warn('No NmrResidue specified')
      return

    if not strip:
      self.project._logger.warn('No Strip specified')
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
    chemicalShiftList = self.project.getByPid(self.chemicalShiftListPulldown.currentText())

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
      self.project._logger.warn('No assignment matrix specified')
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

    for modulePid in self.matchModules:
      module = self.project.getByPid(modulePid)
      makeStripPlot(module, nmrAtomPairs)

      for ii, strip in enumerate(module.strips):
        nmrResidueId = nmrAtomPairs[ii][0].nmrResidue._id
        strip.planeToolbar.spinSystemLabel.setText(nmrResidueId)

      self._centreStripForNmrResidue(assignMatrix[assignmentScores[0]], module.strips[0])

  def _connectSequenceGraph(self, sequenceGraph:CcpnModule):
    """
    # CCPN INTERNAL - called in showSequenceGraph method of GuiMainWindow.
    Connects Sequence Graph to this module.
    """
    self.sequenceGraph = sequenceGraph
    self.project._appBase.current.assigner = sequenceGraph
    self.sequenceGraph.nmrResidueTable = self.nmrResidueTable
    self.sequenceGraph.setMode('fragment')

  def _closeModule(self):
    """
    Re-implementation of the closeModule method of the CcpnModule class required to remove backboneModule
    attribute from mainWindow when the module closes.
    """
    print(self.parent)
    delattr(self.parent, 'backboneModule')
    self.close()



