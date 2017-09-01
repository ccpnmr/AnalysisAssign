"""
Module to assign peaks
Responds to current.peaks

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
__dateModified__ = "$dateModified: 2017-07-07 16:32:46 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b2 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import typing
import numpy as np
from functools import partial

from PyQt4 import QtGui

from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.Peak import Peak
from ccpn.core.lib import CcpnSorting
from ccpn.core.lib.AssignmentLib import ATOM_NAMES, nmrAtomsForPeaks, peaksAreOnLine, sameAxisCodes
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.Table import ObjectTable, Column
from ccpn.util.Logging import getLogger
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpnmodel.ccpncore.lib.Constants import  defaultNmrChainCode


logger = getLogger()


class PeakAssigner(CcpnModule):
  """Module for assignment of nmrAtoms to the different axes of a peak.
  Module responds to current.peak
  """

  # overide in specific module implementations
  includeSettingsWidget = True
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsPosition = 'top'
  className = 'PeakAssigner'

  def __init__(self, mainWindow,  name="Peak Assigner"):

    CcpnModule.__init__(self, mainWindow=mainWindow, name=name)

    # Derive application, project, and current from mainWindow
    self.mainWindow = mainWindow
    self.application = mainWindow.application
    self.project = mainWindow.application.project
    self.current = mainWindow.application.current
    self.colourScheme = self.application.colourScheme

    self.listWidgets = []
    self.objectTables = []
    self.labels = []
    self.assignmentWidgets = []
    self.chainPulldowns = []
    self.seqCodePulldowns = []
    self.resTypePulldowns = []
    self.atomTypePulldowns = []

    # settings
    doubleToleranceCheckboxLabel = Label(self.settingsWidget, text="Double Tolerances ", grid=(0,0))
    self.doubleToleranceCheckbox = CheckBox(self.settingsWidget, checked=False,
                                            callback=self._updateInterface,
                                            grid=(0,1))

    intraCheckboxLabel = Label(self.settingsWidget, text="Only Intra-residual ", grid=(0,2))
    self.intraCheckbox = CheckBox(self.settingsWidget, checked=False,
                                   callback=self._updateInterface,
                                   grid=(0,3))

    multiCheckboxLabel = Label(self.settingsWidget, text="Allow multiple peaks ", grid=(0,4))
    self.multiCheckbox = CheckBox(self.settingsWidget, checked=True,
                                   callback=self._updateInterface,
                                   grid=(0,5))

    expCheckBoxLabel = Label(self.settingsWidget, "Filter By Experiment", grid=(0,6))
    self.expCheckBox = CheckBox(self.settingsWidget, checked=True,
                                callback=self._updateInterface,
                                grid=(0,7))

    # Main content widgets
    self.peakLabel = Label(self.mainWidget, text='Peak:', bold=True, grid=(0,0), vAlign='center', margins=[2,5,2,5])
    self.selectionFrame = Frame(self.mainWidget, showBorder=True, fShape='noFrame', grid=(1, 0), vAlign='top')
    self.selectionLayout = QtGui.QGridLayout()
    self.selectionLayout.setSpacing(0)
    self.selectionLayout.setContentsMargins(0, 0, 0, 0)
    self.selectionFrame.setLayout(self.selectionLayout)

    # respond to peaks
    self.current.registerNotify(self._updateInterface, 'peaks')
    self._updateInterface()

    self.closeModule = self._closeModule

  def __del__(self):
    self.current.unRegisterNotify(self._updateInterface, 'peaks')

  def _createEmptyNmrAtomsTable(self, dim:int):
    """Create an empty table for the specified peak dimension to contain possible Nmr Atoms that
    can be assigned to that peak dimension.
    """
    columns = [Column('NmrAtom', lambda nmrAtom: str(nmrAtom.id)),
               Column('Shift', lambda nmrAtom: self._getShift(nmrAtom)),
               Column('Delta', lambda nmrAtom: self.getDeltaShift(nmrAtom, dim))]

    objectTable = ObjectTable(self, columns,
                              actionCallback=partial(self._assignNmrAtomToDim, dim),
                              objects=[], autoResize=True)

    self.objectTables.append(objectTable)

  def _createEmptyListWidget(self, dim:int):
    """
    Creates an empty ListWidget to contain the dimensionNmrAtoms assigned to a peak dimension.
    """
    listWidget = ListWidget(self, callback=partial(self._updateAssignmentWidget, dim),
                            rightMouseCallback=self._updateNmrAtomsFromListWidgets)
    listWidget.setFixedWidth(145)   # ejb - a bit wider to see the longer atom pids
    listWidget.setFixedHeight(100)
    self.listWidgets.append(listWidget)

  def _createEmptyWidgetLabel(self, dim:int):
    """
    Creates an empty Label to contain peak dimension position.
    """
    label = Label(self, text='', margins=[2,3,2,1])
    self.labels.append(label)

  def _createAssignmentWidget(self, dim:int):
    """
    Creates an assignment widget consisting of three PulldownLists.
    """
    newAssignmentWidget = QtGui.QWidget()
    newLayout = QtGui.QGridLayout()
    chainLabel = Label(self, 'Chain', hAlign='c')
    seqCodeLabel = Label(self, 'Sequence', hAlign='c')
    atomTypeLabel = Label(self, 'Atom', hAlign='c')
    chainPulldown = self._createChainPulldown(dim)
    seqCodePulldown = self._createSeqCodePulldown(dim)
    atomTypePulldown = self._createAtomTypePulldown(dim)
    applyButton = Button(self, 'New', callback=partial(self._createNewNmrAtom, dim))
    self.reassignButton = Button(self, 'Assign', callback=partial(self._setAssignment, dim))
    newLayout.addWidget(chainLabel, 0, 0)
    newLayout.addWidget(chainPulldown, 0, 1)
    newLayout.addWidget(seqCodeLabel, 1, 0)
    newLayout.addWidget(seqCodePulldown, 1, 1)
    newLayout.addWidget(atomTypeLabel, 2, 0)
    newLayout.addWidget(atomTypePulldown, 2, 1)
    newLayout.addWidget(applyButton, 3, 0, 1, 1)
    newLayout.addWidget(self.reassignButton, 3, 1, 1, 1)
    newAssignmentWidget.setLayout(newLayout)
    self.assignmentWidgets.append(newAssignmentWidget)

  def _setAssignment(self, dim:int):
    """
    Assigns dimensionNmrAtoms to peak dimension when called using Assign Button in assignment widget.
    """
    # FIXME Potential Bug: no error checks for dim. It can give easily an IndexError


    nmrChain = self.project.fetchNmrChain(self.chainPulldowns[dim].currentText())
    nmrResidue = nmrChain.fetchNmrResidue(self.seqCodePulldowns[dim].currentText())
    nmrAtom = nmrResidue.fetchNmrAtom(self.atomTypePulldowns[dim].currentText())
    for peak in self.current.peaks:
      # dimNmrAtoms = peak.dimensionNmrAtoms[dim]

      dimNmrAtoms = list(peak.dimensionNmrAtoms[dim])   # ejb - changed to list
      currentItem = self.listWidgets[dim].currentItem()
      if not currentItem:
        self.listWidgets[dim].addItem(nmrAtom.pid)
        currentItem = self.listWidgets[dim].item(self.listWidgets[dim].count()-1)
      currentObject = self.project.getByPid(currentItem.text())

      try:
        toAssign = dimNmrAtoms.index(currentObject)                   # error here..

        dimNmrAtoms[toAssign] = nmrAtom
        allAtoms = list(peak.dimensionNmrAtoms)
        allAtoms[dim] = dimNmrAtoms
        peak.dimensionNmrAtoms = allAtoms
      except Exception as es:
        showWarning(str(self.windowTitle()), str(es))

    self._updateInterface()

  def _createChainPulldown(self, dim:int) -> PulldownList:
    """
    Creates a PulldownList for selection of NmrChains.
    """
    pulldownList = PulldownList(self)
    pulldownList.setEditable(True)
    pulldownList.lineEdit().editingFinished.connect(partial(self._addItemToPulldown, pulldownList))
    self.chainPulldowns.append(pulldownList)
    return pulldownList

  def _createSeqCodePulldown(self, dim:int) -> PulldownList:
    """
    Creates a PulldownList for selection of NmrResidue Sequence codes.
    """
    pulldownList = PulldownList(self)
    pulldownList.setEditable(True)
    self.seqCodePulldowns.append(pulldownList)
    return pulldownList

  def _createAtomTypePulldown(self, dim:int) -> PulldownList:
    """
    Creates a PulldownList for selection of atom types.
    """
    pulldownList = PulldownList(self)
    pulldownList.setEditable(True)
    self.atomTypePulldowns.append(pulldownList)
    return pulldownList

  def _createEnoughTablesAndLists(self):
    '''Makes sure there are enough tables for the amount
       of dimensions of the currently selected peak(s).
       This method only runs when all peaks have the same
       amount of dimensions as is guaranteed by running
       _peaksAreCompatible.py

    '''

    Ndimensions = len(self.current.peak.position)


    # Create extra tables if needed.
    for dim in range(len(self.objectTables), Ndimensions):
      self._createEmptyNmrAtomsTable(dim)

    for dim in range(len(self.listWidgets), Ndimensions):
      self._createEmptyListWidget(dim)

    for dim in range(len(self.labels), Ndimensions):
      self._createEmptyWidgetLabel(dim)

    for dim in range(len(self.assignmentWidgets), Ndimensions):
      self._createAssignmentWidget(dim)

    self.widgetItems = list(zip(self.labels[:Ndimensions], self.listWidgets[:Ndimensions],
                    self.assignmentWidgets[:Ndimensions], self.objectTables[:Ndimensions]))

    for pair in self.widgetItems:
      widget = QtGui.QWidget(self)
      layout = QtGui.QGridLayout()
      #layout.setSpacing(10)
      #layout.setMargin(5)
      #layout.setContentsMargins(4, 4, 4, 4)
      layout.setSpacing(2)
      layout.setMargin(1)
      layout.setContentsMargins(2, 1, 2, 1)
      layout.addWidget(pair[0], 0, 0, 1, 1)
      layout.addWidget(pair[1], 1, 0, 2, 1)
      layout.addWidget(pair[2], 1, 1, 2, 1)
      layout.addWidget(pair[3], 3, 0, 1, 2)
      pair[2].setStyleSheet("PulldownList {border: 0px solid;}")
      pair[2].setStyleSheet("border: 0px solid")
      pair[3].setStyleSheet("color: black; border: 0px solid;")
      widget.setLayout(layout)
      self.widgets.append(widget)
      self.selectionLayout.addWidget(widget, 0, self.widgetItems.index(pair))
    #
    self._updateLayout(self.selectionLayout, Ndimensions)

  # Update functions

  def _updateInterface(self, peaks:typing.List[Peak]=None):
    """Updates the whole module, including recalculation
       of which nmrAtoms fit to the peaks.

    """
    self._emptyAllTablesAndLists()
    if not self.current.peaks or not self._peaksAreCompatible():
      return
    self.peakLabel.setText('Peak: %s' % self.current.peak.id)
    self._createEnoughTablesAndLists()
    self._updateTables()
    self._updateAssignedNmrAtomsListwidgets()
    self._updateWidgetLabels()

  def _updateWidgetLabels(self):

    Ndimensions = len(self.current.peak.position)

    for dim, label in zip(range(Ndimensions), self.labels):
      positions = [peak.position[dim] for peak in self.current.peaks]
      # FIXME Bug: could be ZERO DIVISION!
      avgPos = round(sum(positions)/len(positions), 3)
      axisCode = self.current.peak.peakList.spectrum.axisCodes[dim]
      text = 'Axis "%s": %.3f' % (axisCode, avgPos)
      label.setText(text)

  def _updateTables(self):
    '''Updates the tables indicating the different assignment
       possibilities of the peak dimensions.

    '''
    peaks = self.current.peaks
    doubleTolerance = self.doubleToleranceCheckbox.isChecked()
    intraResidual = self.intraCheckbox.isChecked()
    nmrAtomsForTables = nmrAtomsForPeaks(peaks, self.project.nmrAtoms,
                                             doubleTolerance=doubleTolerance,
                                             intraResidual=intraResidual)

    Ndimensions = len(nmrAtomsForTables)
    for dim, objectTable, nmrAtoms in zip(range(Ndimensions),
                                          self.objectTables,
                                          nmrAtomsForTables):
      if peaksAreOnLine(peaks, dim):
        objectTable.setObjects(nmrAtomsForTables[dim])
        objectTable.show()
      else:
        objectTable.setObjects([NOL])

  def _updateAssignedNmrAtomsListwidgets(self):
    '''Update the listWidget showing which nmrAtoms
       are assigned to which peak dimensions. If multiple
       peaks are selected, only the assignment that they
       have in common are shown. Maybe this should be all
       assignments. You can see that at the peak annotation
       though.
    '''

    Ndimensions = len(self.current.peak.position)

    if self.current.peaks:
      for dim, listWidget in zip(range(Ndimensions), self.listWidgets):

        ll = [set(peak.dimensionNmrAtoms[dim]) for peak in self.current.peaks]
        self.nmrAtoms = list(sorted(set.intersection(*ll)))
        listWidget.addItems([str(a.pid) for a in self.nmrAtoms])

  def _updateNmrAtomsFromListWidgets(self):

    assignmentArray = [0] * len(self.listWidgets)
    for listWidget in self.listWidgets:
      assignments = [self.project.getByPid(listWidget.item(i).text()) for i in range(listWidget.count())]
      index = self.listWidgets.index(listWidget)
      assignmentArray[index] = assignments

    print(assignmentArray, 'assignmentArray')
    for peak in self.current.peaks:
      peak.dimensionNmrAtoms = assignmentArray

  def _updateLayout(self, layout:QtGui.QLayout, ndim:int):
    """
    Remove excess assignment widgets if number of dimensions is less than number of assignment
    widgets displayed.
    """

    rowCount = layout.rowCount()
    colCount = layout.columnCount()

    for r in range(rowCount):
      for m in range(ndim, colCount):
        item = layout.itemAtPosition(r, m)
        if item:
          if item.widget():
            item.widget().hide()
        layout.removeItem(item)


  def _updateAssignmentWidget(self, dim:int, item:object):
    """
    Update all information in assignment widget when NmrAtom is selected in list widget of that
    assignment widget.
    """
    nmrAtom = self.project.getByPid(item.text())
    # self.project._appBase.current.nmrAtom = nmrAtom
    chain = nmrAtom.nmrResidue.nmrChain
    sequenceCode = nmrAtom.nmrResidue.sequenceCode
    self.chainPulldowns[dim].setData([chain.id for chain in self.project.nmrChains])
    self.chainPulldowns[dim].setIndex(self.chainPulldowns[dim].texts.index(chain.id))
    sequenceCodes = [nmrResidue.sequenceCode for nmrResidue in self.project.nmrResidues]
    self.seqCodePulldowns[dim].setData(sorted(sequenceCodes, key=CcpnSorting.stringSortKey))
    self.seqCodePulldowns[dim].setIndex(self.seqCodePulldowns[dim].texts.index(sequenceCode))
    atomPrefix = self.current.peak.peakList.spectrum.isotopeCodes[dim][-1]
    atomNames = [atomName for atomName in ATOM_NAMES if atomName[0] == atomPrefix] + [nmrAtom.name]
    self.atomTypePulldowns[dim].setData(atomNames)
    self.atomTypePulldowns[dim].setIndex(self.atomTypePulldowns[dim].texts.index(nmrAtom.name))


  def _setResidueType(self, dim:int, index:int):
    """
    Set residue type in assignment widget based on chain and sequence code.
    """
    sequenceCode = self.seqCodePulldowns[dim].texts[index]
    nmrChain = self.project.fetchNmrChain(self.chainPulldowns[dim].currentText())
    residueType = nmrChain.fetchNmrResidue(sequenceCode).residueType
    self.resTypePulldowns[dim].setIndex(self.resTypePulldowns[dim].texts.index(residueType.upper()))



  def _addItemToPulldown(self, pulldown:object):
    """
    Generic function to add items to pulldown list if text in pulldown list widget is changed
    """
    if pulldown.lineEdit().isModified():
      text = pulldown.lineEdit().text()
      if text not in pulldown.texts:
        pulldown.addItem(text)


  def getDeltaShift(self, nmrAtom:NmrAtom, dim:int) -> float:
    """
    Calculation of delta shift to add to the table.
    """
    if not self.current.peaks:
      return ''

    if nmrAtom is NOL:
      return ''

    deltas = []
    for peak in self.current.peaks:
      shiftList = peak.peakList.spectrum.chemicalShiftList
      if shiftList:
        shift = shiftList.getChemicalShift(nmrAtom.id)
        if shift:
          position = peak.position[dim]
          deltas.append(abs(shift.value-position))
    # average = sum(deltas)/len(deltas) #Bug: ZERO DIVISION!

    if len(deltas)>0:
      return '%6.3f' % np.mean(deltas)
    else:
      return ''

  def _getShift(self, nmrAtom:NmrAtom) -> float:
    """
    Calculation of chemical shift value to add to the table.
    """
    if not self.current.peaks:
      return ''

    if nmrAtom is NOL:
      return ''
    for peak in self.current.peaks:
      shiftList = peak.peakList.spectrum.chemicalShiftList
      if shiftList:
        shift = shiftList.getChemicalShift(nmrAtom.id)
        if shift:
          return '%8.3f' % shift.value


  def _peaksAreCompatible(self) -> bool:
    """
    If multiple peaks are selected, a check is performed
    to determine whether assignment of corresponding
    dimensions of a peak allowed.
    """

    if len(self.current.peaks) == 1:
      return True
    if not self.multiCheckbox.isChecked():
      self.project._logger.warning("Multiple peaks selected, not allowed.")
      return False
    dimensionalities = set([len(peak.position) for peak in self.current.peaks])
    if len(dimensionalities) > 1:
      self.project._logger.warning('Not all peaks have the same number of dimensions.')
      return False
    for dim in range(len(self.current.peak.position)):
      if not sameAxisCodes(self.current.peaks, dim):
        self.project._logger.warning('''The combination of axiscodes is different for multiple
                 selected peaks.''')
        return False
    return True


  def _emptyAllTablesAndLists(self):
    """
    Quick erase of all present information in ListWidgets and ObjectTables.
    """
    self.peakLabel.setText('Peak: <None>')
    for label in self.labels:
      label.setText('')
    for objectTable in self.objectTables:
      objectTable.setObjects([])
    for listWidget in self.listWidgets:
      listWidget.clear()


  def _createNewNmrAtom(self, dim):
    isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[dim]
    nmrAtom = self.project.fetchNmrChain(shortName=defaultNmrChainCode
                                           ).newNmrResidue().newNmrAtom(isotopeCode=isotopeCode)


    for peak in self.current.peaks:
      if nmrAtom not in peak.dimensionNmrAtoms[dim]:
        # newAssignments = peak.dimensionNmrAtoms[dim] + [nmrAtom]

        newAssignments = list(peak.dimensionNmrAtoms[dim]) + [nmrAtom]    # ejb - changed to list
        axisCode = peak.peakList.spectrum.axisCodes[dim]
        peak.assignDimension(axisCode, newAssignments)
    self.listWidgets[dim].addItem(nmrAtom.pid)
    self._updateTables()



  def _assignNmrAtomToDim(self, dim:int, row:int=None, col:int=None, obj:object=None):
    '''Assign the nmrAtom that is double clicked to the
       the corresponding dimension of the selected
       peaks.

    '''
    #### Should be
    objectTable = self.objectTables[dim]
    nmrAtom = objectTable.getCurrentObject()

    if nmrAtom is NOL:
      return

    for peak in self.current.peaks:
      #### Should be simplified with function in Peak class
      if nmrAtom not in peak.dimensionNmrAtoms[dim]:
        newAssignments = list(peak.dimensionNmrAtoms[dim]) + [nmrAtom]
        axisCode = peak.peakList.spectrum.axisCodes[dim]
        peak.assignDimension(axisCode, newAssignments)

    self.listWidgets[dim].addItem(nmrAtom.pid)

  def _closeModule(self):
    """
    Re-implementation of closeModule function from CcpnModule to unregister notification on current.peaks
    """
    self.project._appBase.current.unRegisterNotify(self._updateInterface, 'peaks')
    super(PeakAssigner, self)._closeModule()

class NotOnLine(object):
  """
  Small 'fake' object to get a message the user in the assignment
  Table that a specific dimension can not be assigned in one go
  since the frequencies of the peaks in this dimension are not on
  one line (i.e. the C frequencies of the CA and CB in a strip for
  instance).
  """

  def __init__(self):
    self.pid = 'Multiple selected peaks not on line.'
    self.id = 'Multiple selected peaks not on line.'

NOL = NotOnLine()
