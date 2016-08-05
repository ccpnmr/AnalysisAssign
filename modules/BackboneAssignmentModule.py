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
import typing

from ccpn.AnalysisAssign.lib.scoring import qScore

from ccpn.core.ChemicalShift import ChemicalShift
from ccpn.core.NmrResidue import NmrResidue

from ccpn.ui.gui.lib.Window import navigateToNmrAtomsInStrip, matchAxesAndNmrAtoms, makeStripPlot #markPositionsInStrip

from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTable
from ccpn.ui.gui.modules.GuiStrip import GuiStrip

from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.Module import CcpnModule
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget

from ccpnmodel.ccpncore.lib.spectrum import Spectrum as spectrumLib


class BackboneAssignmentModule(CcpnModule):

  def __init__(self, parent, project):

    super(BackboneAssignmentModule, self).__init__(parent=None, name='Backbone Assignment')

    self.project = project
    self.current = project._appBase.current
    self.numberOfMatches = 5
    self.nmrChains = project.nmrChains
    self.matchModules = []
    self.nmrResidueTable = NmrResidueTable(self.widget1, project, callback=self._startAssignment)
    self.parent = parent
    self.settingsButton = Button(self.nmrResidueTable, icon='icons/applications-system',
                                grid=(0, 5), hPolicy='fixed', toggle=True)

    self.layout.addWidget(self.nmrResidueTable, 0, 0, 1, 3)
    self.settingsButton.toggled.connect(self._toggleWidget2)
    modules = [display.pid for display in project.spectrumDisplays]
    modules.insert(0, '  ')
    modulesLabel = Label(self.widget2, text="Selected Modules", grid=(1, 0))
    self.modulePulldown = PulldownList(self.widget2, grid=(2, 0), callback=self._selectMatchModule)
    self.modulePulldown.setData(modules)
    self.moduleList = ListWidget(self.widget2, grid=(3, 0), )
    chemicalShiftListLabel = Label(self.widget2, text='Chemical Shift List', grid=(1, 1))
    self.chemicalShiftListPulldown = PulldownList(self.widget2, grid=(2, 1), callback=self._setupShiftDicts)
    self.chemicalShiftListPulldown.setData([shiftlist.pid for shiftlist in project.chemicalShiftLists])
    self.settingsButton.setChecked(False)

    self.project.registerNotifier('NmrResidue', 'rename', self._updateNmrResidueTable)
    self.project.registerNotifier('NmrChain', 'delete', self._updateNmrChainPulldown)
    self.project.registerNotifier('NmrChain', 'create', self._updateNmrChainList)

  def _updateNmrChainList(self, nmrChain):
    self.nmrResidueTable.nmrResidueTable.objectLists.append(nmrChain)

  def _updateNmrChainPulldown(self, nmrChain):
    self.nmrResidueTable.nmrResidueTable.objectLists = self.project.nmrChains
    self.nmrResidueTable.nmrResidueTable._updateSelectorContents()

  def _updateNmrResidueTable(self, nmrResidue, oldPid=None):
    if nmrResidue == self.current.nmrResidue:
      self.nmrResidueTable.nmrResidueTable._updateSelectorContents()
      self.nmrResidueTable.nmrResidueTable.selector.select(nmrResidue.nmrChain)
      self.nmrResidueTable.updateTable()


  def _toggleWidget2(self):
    if self.settingsButton.isChecked():
      self.widget2.show()
    else:
      self.widget2.hide()

  def _selectMatchModule(self, item):
    self.moduleList.addItem(item)
    self.matchModules.append(item)


  def _startAssignment(self, nmrResidue:NmrResidue, row:int=None, col:int=None):
    """
    Initiates assignment procedure when triggered by selection of an NmrResidue from the nmrResidueTable
    inside the module.
    """
    self.project._startFunctionCommandBlock('_startAssignment', nmrResidue)
    try:
      self._setupShiftDicts()

      self.current.nmrChain = nmrResidue.nmrChain
      if hasattr(self, 'assigner'):
        self.assigner.clearAllItems()
        self.assigner.nmrChainPulldown.select(self.current.nmrChain.pid)
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
    to chemical shift value NmrAtoms in the NmrResidue. Creates assignMatrix for strip matching and
    add strips to matchModule(s) corresponding to assignment matches.
    """
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

      if nmrResidue.sequenceCode.endswith('-1'):
        direction = '-1'
        iNmrResidue = nmrResidue.mainNmrResidue
        queryShifts = self.interShifts[nmrResidue]
        matchShifts = self.intraShifts

      else:
        direction = '+1'
        iNmrResidue = nmrResidue
        queryShifts = self.intraShifts[nmrResidue]
        matchShifts = self.interShifts

      self.current.nmrResidue = iNmrResidue

      for display in selectedDisplays:
        if not strip:
          strip = display.strips[0]
        self._displayNmrResidueInStrip(iNmrResidue, strip)
        self._centreStripForNmrResidue(nmrResidue, strip)

      assignMatrix = self._buildAssignmentMatrix(queryShifts, matchShifts)
      if not assignMatrix[1]:
        self.project._logger.info('No matches found for NmrResidue: %s' % nmrResidue.pid)
        return
      self._createMatchStrips(assignMatrix)
      if hasattr(self, 'assigner'):
        if self.assigner.nmrChainPulldown.currentText() != nmrResidue.nmrChain.pid:
          self.assigner.nmrChainPulldown.select(nmrResidue.nmrChain.pid)
        elif not nmrResidue.nmrChain.isConnected:
          self.assigner.addResidue(iNmrResidue, direction)
        else:
          self.assigner.setNmrChainDisplay(nmrResidue.nmrChain.pid)
    finally:
      self.project._appBase._endCommandBlock()



  def _displayNmrResidueInStrip(self, nmrResidue, strip):
    navigateToNmrAtomsInStrip(nmrResidue.nmrAtoms, strip=strip, widths=['default']*len(strip.axisCodes))
    strip.planeToolbar.spinSystemLabel.setText(nmrResidue._id)
    # self._centreStripForNmrResidue(nmrResidue, strip)

  def _centreStripForNmrResidue(self, nmrResidue, strip):
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
    # Assignment score has format {score: nmrResidue} where score is a float
    # assignMatrix[0] is a dict {score: nmrResidue} assignMatrix[1] is a concurrent list of scores
    assignmentScores = sorted(assignMatrix[1])[0:self.numberOfMatches]
    nmrAtomPairs = []
    for assignmentScore in assignmentScores:
      matchResidue = assignMatrix[0][assignmentScore]
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
      self._centreStripForNmrResidue(nmrAtomPairs[0][0].nmrResidue, module.strips[0])

  def _connectSequenceGraph(self, assigner:CcpnModule):
    """
    # CCPN INTERNAL - called in showSequenceGraph method of GuiMainWindow.
    Connects Sequence Graph to this module.
    """
    self.assigner = assigner
    self.project._appBase.current.assigner = assigner
    self.assigner.nmrResidueTable = self.nmrResidueTable
    self.assigner.setMode('fragment')


  def _buildAssignmentMatrix(self, queryShifts:typing.List[ChemicalShift],
                             matchShifts:typing.Dict[NmrResidue, ChemicalShift]) -> typing.Tuple[typing.Dict[NmrResidue, typing.List[ChemicalShift]], typing.List[float]]:
    """
    Creates a dictionary of NmrResidues and qScores between queryShifts and matching shifts.
    Returns dictionary and a list of the qScores.
    """
    scores = []
    matrix = OrderedDict()
    for res, shift in matchShifts.items():

      if len(queryShifts) > 1 and len(shift) > 1:

        if qScore(queryShifts[0].value, shift[0].value) is not None and qScore(queryShifts[1].value, shift[1].value) is not None:

          score = (qScore(queryShifts[0].value, shift[0].value)+qScore(queryShifts[1].value, shift[1].value))/2
          scores.append(score)
          matrix[score] = res
      elif len(queryShifts) == 1 and shift:
        if qScore(queryShifts[0].value, shift[0].value) is not None:

          score = qScore(queryShifts[0].value, shift[0].value)
          scores.append(score)
          matrix[score] = res

    return matrix, scores

  def closeModule(self):
    delattr(self.parent, 'backboneModule')
    self.close()



