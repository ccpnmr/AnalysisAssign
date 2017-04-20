"""
Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
Use settings to define the spectral displays, the active spectra and the tolerances for peak picking

This module closely works with the Atom Selector module

First version by SS
Refactored by GWV to be responsible to active state on start; proper callbacks 
and to include "Restricted pick and assign" button.

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
__dateModified__ = "$dateModified: 2017-04-07 11:40:22 +0100 (Fri, April 07, 2017) $"
__version__ = "$Revision: 3.0.b1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"

__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================
from PyQt4 import QtGui

from ccpn.ui.gui.lib import PeakList
from ccpn.ui.gui.lib import Strip
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTable, NmrResidueTableModule
from ccpn.ui.gui.widgets.Base import Base
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.DoubleSpinbox import DoubleSpinbox
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.ScrollArea import ScrollArea
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget, DoubleSpinBoxCompoundWidget

from ccpn.util.Logging import getLogger
logger = getLogger()


class PickAndAssignModule(NmrResidueTableModule):
  """
  Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
  Use settings to define the spectral displays, the active spectra and the tolerances for peak picking

  This module closely works with the Atom Selector module
  """
  includeSettingsWidget = True
  maxSettingsState = 2
  settingsMinimumSizes = (500, 200)

  def __init__(self, parent=None, project=None, name='Pick And Assign', **kw):

    super(PickAndAssignModule, self).__init__(parent=parent, name=name)
    # project, current, application and mainWindow are inherited from CcpnModule

    # Main widget
    self.restrictedPickButton = Button(self.nmrResidueTable._widget, text='Restricted Pick', grid=(0, 2),
                                       callback=self.restrictedPick)
    self.restrictedPickButton.setMinimumSize(120, 30)
    self.assignSelectedButton = Button(self.nmrResidueTable._widget, text='Assign Selected', grid=(0, 3),
                                       callback=self.assignSelected)
    self.assignSelectedButton.setMinimumSize(120, 30)
    self.restrictedPickAndAssignButton = Button(self.nmrResidueTable._widget, text='Restricted Pick and Assign', grid=(0, 4),
                                                callback=self.restrictedPickAndAssign)
    self.restrictedPickAndAssignButton.setMinimumSize(160, 30)

    # Settings widget

    # change some of the defaults setting inherited from NmrResidueTableModule
    self.sequentialStripsWidget.checkBox.setChecked(False)
    self.displaysWidget.addPulldownItem(0)

    # create row's of spectrum information
    self._spectraWidget = Frame(self.settingsWidget, grid=(0, 1), gridSpan=(4,1))
    self.spectrumLabel = Label(self._spectraWidget, 'Spectrum', bold=True, grid=(0,0), hAlign='left')
    self.useLabel = Label(self._spectraWidget, 'Use?', bold=True, grid=(0, 1), hAlign='left')

    self._spectraWidgets = {}  # spectrum.pid, frame dict to show/hide
    spectra = list(set([sp for dp in self.mainWindow.spectrumDisplays for sp in dp.strips[0].spectra]))
    for row, spectrum in enumerate(self.project.spectra):
      f = _SpectrumRow(self._spectraWidget, spectrum,
                       grid=(row+1,0), gridSpan=(1,1+len(spectrum.axisCodes)), vAlign='top')
      self._spectraWidgets[spectrum.pid] = f

  def _closeModule(self):
    """
    Unregister notifiers and close module.
    """
    #self._unRegisterNotifiers()
    self.close()

  def assignSelected(self):
    "Assign current.peaks on the bases of nmrAtoms of current.nmrResidue"

    if self.current.nmrResidue is None:
      logger.error('Undefined nmrResidue; select one first before proceeding')
      return

    if len(self.current.peaks) == 0:
      logger.error('Undefined peak(s); select one or more before proceeding')
      return

    self.project._appBase._startCommandBlock('application.pickAndAssignModule.assignSelected()')
    try:
      shiftDict = {}
      for atom in self.current.nmrResidue.nmrAtoms:
        shiftDict[atom.isotopeCode] = []

      for peak in self.current.peaks:
        shiftList = peak.peakList.spectrum.chemicalShiftList
        spectrum = peak.peakList.spectrum
        for nmrAtom in self.current.nmrResidue.nmrAtoms:
          if nmrAtom.isotopeCode in shiftDict.keys():
            shiftDict[nmrAtom.isotopeCode].append((nmrAtom, shiftList.getChemicalShift(nmrAtom.id).value))
        for ii, isotopeCode in enumerate(spectrum.isotopeCodes):
          if isotopeCode in shiftDict.keys():
            for shift in shiftDict[isotopeCode]:
              sValue = shift[1]
              pValue = peak.position[ii]
              if abs(sValue-pValue) <= spectrum.assignmentTolerances[ii]:
                peak.assignDimension(spectrum.axisCodes[ii], [shift[0]])
      self.current.peaks = []
      # update the NmrResidue table
      self.nmrResidueTable._update(self.current.nmrResidue.nmrChain)

    finally:
      self.project._endCommandEchoBlock()

  def restrictedPick(self, nmrResidue=None):
    """
    Routine refactored in revision 9381.
 
    Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
    spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
    centre points with tolerances and the y as the long axis to pick the whole region.
    """
    if not nmrResidue:
      nmrResidue = self.current.nmrResidue

    if nmrResidue is None:
      logger.error('Undefined nmrResidue; select one first before proceeding')
      return

    self.project._appBase._startCommandBlock('application.pickAndAssignModule.restrictedPick(nmrResidue)',
                                             nmrResidue=nmrResidue)
    try:
      peaks = []
      for module in self.project.spectrumDisplays:
        if len(module.axisCodes) > 2:
          for spectrumView in module.strips[0].spectrumViews:
            visiblePeakListViews = [peakListView for peakListView in spectrumView.peakListViews
                                    if peakListView.isVisible()]
            if len(visiblePeakListViews) == 0:
              continue
            else:
              peakList, pks = PeakList.restrictedPick(peakListView=visiblePeakListViews[0],
                                                        axisCodes=module.axisCodes[0::2], nmrResidue=nmrResidue)
              peaks = peaks + pks
      self.current.peaks = peaks
      # update the NmrResidue table
      self.nmrResidueTable._update(nmrResidue.nmrChain)

    finally:
      self.project._appBase._endCommandBlock()

  def restrictedPickAndAssign(self, nmrResidue=None):
    """
    Functionality for beta2 to include the Assign part
     
    Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
    spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
    centre points with tolerances and the y as the long axis to pick the whole region.
    
    Calls assignSelected to assign
    """
    if not nmrResidue:
      nmrResidue = self.current.nmrResidue
    elif not self.current.nmrResidue:
      print('No current nmrResidue')
      return

    self.project._appBase._startCommandBlock('application.pickAndAssignModule.restrictedPickAndAssign(nmrResidue)',
                                              nmrResidue=nmrResidue)
    try:
      peaks = []
      for module in self.project.spectrumDisplays:
        if len(module.axisCodes) > 2:
          for spectrumView in module.strips[0].spectrumViews:
            visiblePeakListViews = [peakListView for peakListView in spectrumView.peakListViews
                                    if peakListView.isVisible()]
            if len(visiblePeakListViews) == 0:
              continue
            else:
              peakList, pks = PeakList.restrictedPick(peakListView=visiblePeakListViews[0],
                                                      axisCodes=module.axisCodes[0::2], nmrResidue=nmrResidue)
              peaks = peaks + pks
      self.current.peaks = peaks
      self.assignSelected()
      # update the NmrResidue table
      self.nmrResidueTable._update(nmrResidue.nmrChain)

    finally:
      self.project._endCommandEchoBlock()

  def goToPositionInModules(self, nmrResidue=None, row=None, col=None):
    "Go to the positions defined my NmrAtoms of nmrResidue in the active displays"

    activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()
    self.project._appBase._startCommandBlock('application.pickAndAssignModule.goToPositionInModules(nmrResidue)', nmrResidue=nmrResidue)
    try:
      if nmrResidue is not None:
        if self.project._appBase.ui.mainWindow is not None:
          mainWindow = self.project._appBase.ui.mainWindow
        else:
          mainWindow = self.project._appBase._mainWindow
        mainWindow.clearMarks()
        for display in activeDisplays:
          strip = display.strips[0]
          n = len(strip.axisCodes)
          #Strip.navigateToNmrAtomsInStrip(strip=strip, nmrAtoms=nmrResidue.nmrAtoms, widths=['default', 'full', ''])
          if n == 2:
            widths = ['default', 'default']
          else:
            widths = ['default', 'full'] + (n-2)*['']
          Strip.navigateToNmrAtomsInStrip(strip=strip, nmrAtoms=nmrResidue.nmrAtoms, widths=widths, markPositions=(n==2))
        self.current.nmrResidue = nmrResidue
    finally:
      self.project._endCommandEchoBlock()


class _SpectrumRow(Frame):
  "Class to make a spectrum row"

  def __init__(self, parent, spectrum, **kwds):

    super(_SpectrumRow, self).__init__(parent, **kwds)

    col = 0
    self.checkbox = CheckBoxCompoundWidget(self, grid=(0, col), gridSpan=(1,1), hAlign='left',
                                           checked=True, labelText=spectrum.pid,
                                           minimumWidths = [100,20] )

    self.spinBoxes = []
    for ii, axisCode in enumerate(spectrum.axisCodes):
      decimals, step = (2, 0.01) if axisCode[0:1] == 'H' else (1, 0.1)
      col += 1; ds = DoubleSpinBoxCompoundWidget(
                                   self, grid=(0, col), gridSpan=(1,1), hAlign='left',
                                   minimumWidths=[30, 50], maximumWidths=[30, 50],
                                   labelText = axisCode,
                                   value = spectrum.assignmentTolerances[ii],
                                   decimals=decimals, step=step, range=(0, None))
      self.spinBoxes.append(ds)


class SpectrumSelectionWidget(QtGui.QWidget, Base):
  "Class to make a widget with spectral settings"
  def __init__(self, parent, project, getDisplays, **kw):

    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, **kw)
    self.project = project
    self._getDisplays = getDisplays # function to get the displays


    self.spectrumLabel = Label(self, 'Spectrum')
    self.layout().addWidget(self.spectrumLabel, 0, 0)
    self.useLabel = Label(self, 'Use?', grid=(0, 1), hAlign='c')
    self.refreshBox = CheckBox(self, grid=(0, 4))
    self.checkBoxLabel = Label(self, 'Auto Refresh', grid=(0, 5))
    self.update = self._update
    self.setObjectName('spectrumSelectionWidget')
    self.update()
    # self.setStyleSheet()

  def getActiveDisplays(self):
    if self.displayList.count() == 1 and self.displayList.item(0).text() == '<All>':
      return self.project.spectrumDisplays
    else:
      return [self.project.getByPid(self.displayList.item(ii).text()) for ii in range(self.displayList.count())]

  def _update(self):
    rowCount = self.layout().rowCount()
    colCount = self.layout().columnCount()

    for r in range(1, rowCount):
      for m in range(colCount):
        item = self.layout().itemAtPosition(r, m)
        if item:
          if item.widget():
            item.widget().hide()
        self.layout().removeItem(item)
    # if self.displayList.count() == 1 and self.displayList.item(0).text() == '<All>':
#    activeDisplays = self.getActiveDisplays()
    activeDisplays = self._getDisplays()
    spectra = set([spectrumView.spectrum for spectrumDisplay in activeDisplays
                   for spectrumView in spectrumDisplay.spectrumViews])


    for ii, spectrum in enumerate(spectra):
      for tol in spectrum.assignmentTolerances:
        if tol is None:
          index = spectrum.assignmentTolerances.index(tol)
          tolerance = spectrum.spectralWidths[index]/spectrum.pointCounts[index]
          spectrumTolerances = list(spectrum.assignmentTolerances)
          spectrumTolerances[index] = tolerance
          spectrum.assignmentTolerances = spectrumTolerances

      spectrumLabel1 = Label(self, spectrum.pid, grid=(ii+1, 0), vAlign='t')
      spectrumCheckBox1 = CheckBox(self, grid=(ii+1, 1), hAlign='c', vAlign='t')
      spectrumCheckBox1.setChecked(True)
      spectrumTol1 = Label(self, spectrum.axisCodes[0], grid=(ii+1, 2), vAlign='t')
      spectrumTol1Value = DoubleSpinbox(self, grid=(ii+1, 3), vAlign='t')
      spectrumTol1Value.setDecimals(3)
      spectrumTol1Value.setSingleStep(0.001)
      spectrumTol1Value.setValue(spectrum.assignmentTolerances[0])

      spectrumTol2 = Label(self, spectrum.axisCodes[1], grid=(ii+1, 4), vAlign='t')
      spectrumTol2Value = DoubleSpinbox(self, grid=(ii+1, 5), vAlign='t')
      spectrumTol2Value.setDecimals(3)
      spectrumTol2Value.setSingleStep(0.001)
      spectrumTol2Value.setValue(spectrum.assignmentTolerances[1])

      if spectrum.dimensionCount > 2:
        for jj in range(spectrum.dimensionCount-2):
          spectrumTol3 = Label(self, spectrum.axisCodes[2+jj], grid=(ii+1+jj, 6), vAlign='t')
          spectrumTol3Value = DoubleSpinbox(self, grid=(ii+jj+1, 7), vAlign='t')
          spectrumTol3Value.setDecimals(3)
          spectrumTol3Value.setSingleStep(0.001)
          spectrumTol3Value.setValue(spectrum.assignmentTolerances[2+jj])




