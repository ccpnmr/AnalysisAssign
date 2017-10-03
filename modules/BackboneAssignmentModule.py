"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2017"
__credits__ = ("Wayne Boucher, Ed Brooksbank, Rasmus H Fogh, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license",
               "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for licence text")
__reference__ = ("For publications, please use reference from http://www.ccpn.ac.uk/v3-software/downloads/license",
               "or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:21 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b2 $"
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
from PyQt5 import QtGui, QtWidgets

from ccpn.AnalysisAssign.lib.scoring import getNmrResidueMatches
from ccpn.core.ChemicalShift import ChemicalShift
from ccpn.core.NmrResidue import NmrResidue
from ccpn.ui.gui.lib.SpectrumDisplay import makeStripPlot

from ccpn.ui.gui.lib.Strip import matchAxesAndNmrAtoms
from ccpn.ui.gui.lib.Strip import navigateToNmrResidueInDisplay

from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTableModule

from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.CompoundWidgets import ListCompoundWidget, PulldownListCompoundWidget
from ccpn.ui.gui.widgets.MessageDialog import showWarning, progressManager
from ccpn.ui.gui.widgets.PulldownListsForObjects import ChemicalShiftListPulldown
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.lib.GuiNotifier import GuiNotifier
from ccpn.ui.gui.widgets.DropBase import DropBase

from ccpn.util.Logging import getLogger
from ccpn.core.NmrAtom import NmrAtom

ALL = '<all>'


class BackboneAssignmentModule(NmrResidueTableModule):

  className = 'BackboneAssignmentModule'

  includeSettingsWidget = True
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsPosition = 'left'
  settingsMinimumSizes = (500, 200)

  def __init__(self, mainWindow, name='Backbone Assignment'):

    super(BackboneAssignmentModule, self).__init__(mainWindow=mainWindow, name=name)

    # Derive application, project, and current from mainWindow
    self.mainWindow = mainWindow
    self.application = mainWindow.application
    self.project = mainWindow.application.project
    self.current = mainWindow.application.current

    self.nmrChains = self.application.project.nmrChains
    self.matchCheckBoxWidget = CheckBox(self.nmrResidueTable._widget,
                                        grid=(0,2), checked=True, text='Find matches')

    ### Settings ###

    # change some of the defaults setting inherited from NmrResidueTableModule
    self.sequentialStripsWidget.checkBox.setChecked(True)
    self.displaysWidget.addPulldownItem(0)

    colWidth = 200  # for labels of the compound widgets
    colWidth2 = 120  # for the numberOfMatchesWidget

    row = 4 ## Number of widgets of NmrResidueTable
    col = 0

    # Number of matches to show
    row += 1
    self.numberOfMatchesWidget = PulldownListCompoundWidget(self.settingsWidget,
                                          grid=(row,col), vAlign='top', hAlign='left',
                                          fixedWidths=(colWidth2, colWidth2, colWidth2),
                                          orientation='left',
                                          labelText="Matches to show:",
                                          texts=[str(tt) for tt in range(3,7)]
                                          )

    # Match module selection
    row += 1
    # cannot set a notifier for displays, as these are not (yet?) implemented
    self.matchWidget = ListCompoundWidget(self.settingsWidget,
                                          grid=(row,col), vAlign='top', hAlign='left',
                                          fixedWidths=(colWidth, colWidth, colWidth),
                                          orientation='left',
                                          labelText="Match module(s):",
                                          texts=[display.pid for display in self.mainWindow.spectrumDisplays]
                                          )
    self.matchWidget.setFixedHeigths((None, None, 40))

    # Chemical shift list selection
    row += 1
    self.shiftListWidget = ChemicalShiftListPulldown(self.settingsWidget, self.application.project,
                                                     grid=(row,col), vAlign='top', hAlign='left',
                                                     fixedWidths=(colWidth, colWidth),
                                                     callback=self._setupShiftDicts, default=0
                                                     )
    self._setupShiftDicts()
    self._spacer = Spacer(self.settingsWidget, 5, 5
                         , QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
                         , grid=(row+1,10), gridSpan=(1,1))

    # for compatibility with previous implementation
    #self.moduleList = self.matchWidget.listWidget

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
    Navigate in selected displays to nmrResidue; skip if no displays defined
    If matchCheckbox is checked, also call findAndDisplayMatches
    """
    displays = self._getDisplays()
    if len(displays) == 0:
      getLogger().warning('Undefined display module(s); select in settings first')
      showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
      return

    if self.matchCheckBoxWidget.isChecked() and len(self.matchWidget.getTexts()) == 0:
      getLogger().warning('Undefined match module; select in settings first or unselect "Find matches"')
      showWarning('startAssignment', 'Undefined match module;\nselect in settings first or unselect "Find matches"')
      return


    self.application._startCommandBlock(
        'BackboneAssignmentModule.navigateToNmrResidue(project.getByPid(%r))' % nmrResidue.pid)
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
                                      markPositions=False    #self.markPositionsWidget.checkBox.isChecked()
                                      )
          # activate a callback notifiers; allow dropping onto the NmrResidueLabel
          for st, strip in enumerate(strips):
            if strip is not None:
              # NB connections are made as connectPrevious / connectNext to passed-in NmrResidue
              # It follows that it IS the mainNMr Residue that should be passed in here
              # Note, though, that you get teh same connections WHICHEVER strip you drop on
              notifier = GuiNotifier(strip.getStripLabel(),
                                     [GuiNotifier.DROPEVENT], [DropBase.TEXT],
                                     self._processDroppedNmrResidue, nmrResidue=nr)
              self._stripNotifiers.append(notifier)

            strip.spectrumDisplay.setColumnStretches(True)


          # layout.setColumnStretch(col, colStr
          # strips[0].spectrumDisplay.stripFrame.setStretch(1,1)

      # ejb
      # if 'i-1' residue, take CA CB, and take H, N from the 'i' residue (.mainNmrResidue)
      # check if contains '-1' in pid, is this robust? no :)

      if nmrResidue.relativeOffset == -1:
        # -1 residue so need to split the CA, CB from the N, H
        nmrAtomsMinus = nmrAtomsFromResidue(nmrResidue)
        nmrAtomsCentre = nmrAtomsFromResidue(nmrResidue.mainNmrResidue)

        nmrAtoms=[]
        # this should check the experiment type and choose the correct atoms from there
        for naMinus in nmrAtomsMinus:
          if naMinus.name.startswith('CA') or naMinus.name.startswith('CB'):
            nmrAtoms.append(naMinus)
        for naCentre in nmrAtomsCentre:
          if naCentre.name.startswith('N') or naCentre.name.startswith('H'):
            nmrAtoms.append(naCentre)

        markNmrAtoms(mainWindow=self.mainWindow, nmrAtoms=nmrAtoms)
      else:
        nmrAtoms = nmrAtomsFromResidue(nmrResidue.mainNmrResidue)
        markNmrAtoms(mainWindow=self.mainWindow, nmrAtoms=nmrAtoms)

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
    # if nmrResidue.sequenceCode.endswith('-1'):
    if nmrResidue.relativeOffset == -1:
      # direction = '-1'
      # iNmrResidue = nmrResidue.mainNmrResidue
      queryShifts = [shift for shift in self.interShifts[nmrResidue] if shift.nmrAtom.isotopeCode == '13C']
      matchShifts = self.intraShifts

    elif nmrResidue.relativeOffset:
      getLogger().warning(
        "Assignment matching not supported for NmrResidue offset %s. Matching display skipped"
        % nmrResidue.relativeOffset
      )

    # If NmrResidue is not an offset NmrResidue, set queryShifts as value from self.intraShifts dictionary
    # Set matchShifts as self.interShifts
    else:
      # relative offset is None or 0
      # direction = '+1'
      # iNmrResidue = nmrResidue
      queryShifts = [shift for shift in self.intraShifts[nmrResidue] if shift.nmrAtom.isotopeCode == '13C']
      matchShifts = self.interShifts

    assignMatrix = getNmrResidueMatches(queryShifts, matchShifts, 'averageQScore')
    if not assignMatrix.values():
      getLogger().info('No matches found for NmrResidue: %s' % nmrResidue.pid)
      return
    self._createMatchStrips(assignMatrix)
    return

  def _processDroppedNmrResidue(self, data, nmrResidue):
    """Process the dropped NmrResidue id"""

    droppedNmrResidue = None
    if DropBase.TEXT in data and len(data[DropBase.TEXT]) > 0:
      droppedNmrResidue = self.application.project.getByPid(data[DropBase.TEXT])
    if droppedNmrResidue is None:
      getLogger().info('Backbone assignment: invalid "pid" of dropped item')

    getLogger().debug('nmrResidue:%s, droppedNmrResidue:%s', nmrResidue, droppedNmrResidue)
    if droppedNmrResidue == nmrResidue:
      getLogger().warning('Cannot connect residue to itself')
      return

    # silence the update of the nmrResidueTable as we will to an explicit update later
    # put in try/finally block because otherwise if exception thrown in the following code
    # (which can happen) then you no longer get updates of the NmrResidue table

    with progressManager("connecting %s to %s" % (droppedNmrResidue.pid, nmrResidue.pid)):
      nmrResidue._startCommandEchoBlock("connecting %s to %s" % (droppedNmrResidue.pid, nmrResidue.pid))
      try:

        self.nmrResidueTable.setUpdateSilence(True)
        matchNmrResidue = None
        try:                          # display popup warning
          if data['shiftLeftMouse']:
            # leftShift drag; connect to previous
            nmrResidue.connectPrevious(droppedNmrResidue)
            matchNmrResidue = droppedNmrResidue.getOffsetNmrResidue(offset=-1)
            if matchNmrResidue is None:
              # Non -1 residue - stay with current
              getLogger().info("NmrResidue %s has no i-1 residue to display" % droppedNmrResidue)
              matchNmrResidue = nmrResidue
          else:
            nmrResidue.connectNext(droppedNmrResidue)
            matchNmrResidue = droppedNmrResidue
        except Exception as es:
          showWarning('Connect NmrResidue', str(es))
        finally:
          self.nmrResidueTable.setUpdateSilence(False)

        # update the NmrResidueTable
        self.nmrResidueTable.displayTableForNmrChain(droppedNmrResidue.nmrChain)
        if matchNmrResidue:
          self.navigateToNmrResidue(matchNmrResidue)

      except Exception as es:
        getLogger().warning(str(es))
        # raise es
      finally:
        nmrResidue._endCommandEchoBlock()


  def _centreStripForNmrResidue(self, nmrResidue, strip):
    """
    Centre y-axis of strip based on chemical shifts of from NmrResidue.nmrAtoms
    """
    if not nmrResidue:
      getLogger().warning('No NmrResidue specified')
      return

    if not strip:
      getLogger().warning('No Strip specified')
      return

    yShifts = matchAxesAndNmrAtoms(strip, nmrResidue.nmrAtoms)[strip.axisOrder[1]]
    yShiftValues = [x.value for x in yShifts]
    if yShiftValues:
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
      getLogger().warn('No assignment matrix specified')
      return

    # Assignment score has format {score: nmrResidue} where score is a float
    # assignMatrix[0] is a dict {score: nmrResidue} assignMatrix[1] is a concurrent list of scores
    numberOfMatches = int(self.numberOfMatchesWidget.getText())
    assignmentScores = sorted(list(assignMatrix.keys()))[:numberOfMatches]
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
        nmrResiduePid = nmrAtomPairs[ii][0].nmrResidue.pid
        strip.setStripLabelText(nmrResiduePid)
        strip.showStripLabel()

      self._centreStripForNmrResidue(assignMatrix[assignmentScores[0]], module.strips[0])

  def _closeModule(self):
    """
    Re-implementation of the closeModule method of the CcpnModule class required 
    """
    # TODO: use proper subclassing
    for notifier in self._stripNotifiers:
      notifier.unRegister()
    self._stripNotifiers = []
    super(BackboneAssignmentModule, self)._closeModule()

def nmrAtomsFromResidue(nmrResidue):
  """
  Retrieve a list of nmrAtoms from nmrResidue
  """
  # nmrResidue = nmrResidue.mainNmrResidue
  nmrResidues = []
  previousNmrResidue = nmrResidue.previousNmrResidue
  if previousNmrResidue:
    nmrResidues.append(previousNmrResidue)
  nmrResidues.append(nmrResidue)
  nextNmrResidue = nmrResidue.nextNmrResidue
  if nextNmrResidue:
    nmrResidues.append(nextNmrResidue)

  nmrAtoms = []
  for nr in nmrResidues:
    nmrAtoms.extend(nr.nmrAtoms)

  return nmrAtoms

def markNmrAtoms(mainWindow, nmrAtoms:typing.List[NmrAtom]):

  # get the display
  # displays = self._getDisplays()

  # application = mainWindow.application
  # project = mainWindow.application.project
  # current = mainWindow.application.current

  displays = [dp for dp in mainWindow.spectrumDisplays]

  if len(displays) == 0:
    getLogger().warning('No Spectrum Displays')
    showWarning('markNmrAtoms', 'No spectrum Displays')
    return

  # mainWindow.clearMarks()     # clear the marks for the minute

  for display in displays:
    strips = display.strips

    for strip in strips:
      # assume that this returns list of nmrAtoms in the display

      shiftDict = matchAxesAndNmrAtoms(strip, nmrAtoms)
      # atomPositions = shiftDict[strip.axisOrder[2]]
      atomPositions = [[x.value for x in shiftDict[axisCode]] for axisCode in strip.axisOrder]
      positions = []
      for atomPos in atomPositions:
        if atomPos:
          if len(atomPos) < 2:
            positions.append(atomPos[0])
          else:
            positions.append(max(atomPos) - min(atomPos) / 2)
        else:
          positions.append('')
        # navigateToPositionInStrip(strip, positions, widths=widths) # don't need to change display yet

        strip.spectrumDisplay.mainWindow.markPositions(list(shiftDict.keys()),
                                                       list(shiftDict.values()))



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
