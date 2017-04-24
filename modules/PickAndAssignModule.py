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

  def __init__(self, parent, application):

    super(PickAndAssignModule, self).__init__(parent=parent, application=application, name='Pick And Assign')

    # derive project, current, and mainWindow from application
    self.application = application

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
    for row, spectrum in enumerate(self.application.project.spectra):
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

    if self.application.current.nmrResidue is None:
      logger.error('Undefined nmrResidue; select one first before proceeding')
      return

    if len(self.application.current.peaks) == 0:
      logger.error('Undefined peak(s); select one or more before proceeding')
      return

    self.application.project._appBase._startCommandBlock('application.pickAndAssignModule.assignSelected()')
    try:
      shiftDict = {}
      for atom in self.application.current.nmrResidue.nmrAtoms:
        shiftDict[atom.isotopeCode] = []

      for peak in self.application.current.peaks:
        shiftList = peak.peakList.spectrum.chemicalShiftList
        spectrum = peak.peakList.spectrum
        for nmrAtom in self.application.current.nmrResidue.nmrAtoms:
          if nmrAtom.isotopeCode in shiftDict.keys():
            shiftDict[nmrAtom.isotopeCode].append((nmrAtom, shiftList.getChemicalShift(nmrAtom.id).value))
        for ii, isotopeCode in enumerate(spectrum.isotopeCodes):
          if isotopeCode in shiftDict.keys():
            for shift in shiftDict[isotopeCode]:
              sValue = shift[1]
              pValue = peak.position[ii]
              if abs(sValue-pValue) <= spectrum.assignmentTolerances[ii]:
                peak.assignDimension(spectrum.axisCodes[ii], [shift[0]])
      self.application.current.peaks = []
      # update the NmrResidue table
      self.nmrResidueTable._update(self.application.current.nmrResidue.nmrChain)

    finally:
      self.application._endCommandBlock()

  #TODO:GEERTEN: compact the two routines
  def restrictedPick(self, nmrResidue=None):
    """
    Routine refactored in revision 9381.
 
    Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
    spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
    centre points with tolerances and the y as the long axis to pick the whole region.
    """
    if not nmrResidue:
      nmrResidue = self.application.current.nmrResidue

    if nmrResidue is None:
      logger.error('Undefined nmrResidue; select one first before proceeding')
      return

    self.application._startCommandBlock('application.pickAndAssignModule.restrictedPick(nmrResidue)',
                                             nmrResidue=nmrResidue)
    try:
      peaks = []
      for module in self.application.project.spectrumDisplays:
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
      self.application.current.peaks = peaks
      # update the NmrResidue table
      self.nmrResidueTable._update(nmrResidue.nmrChain)

    finally:
      self.application._endCommandBlock()

  def restrictedPickAndAssign(self, nmrResidue=None):
    """
    Functionality for beta2 to include the Assign part
     
    Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
    spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
    centre points with tolerances and the y as the long axis to pick the whole region.
    
    Calls assignSelected to assign
    """
    if not nmrResidue:
      nmrResidue = self.application.current.nmrResidue
    elif not self.application.current.nmrResidue:
      print('No current nmrResidue')
      return

    self.application._startCommandBlock('application.pickAndAssignModule.restrictedPickAndAssign(nmrResidue)',
                                              nmrResidue=nmrResidue)
    try:
      peaks = []
      for module in self.application.project.spectrumDisplays:
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
      self.application.current.peaks = peaks
      self.assignSelected()
      # update the NmrResidue table
      self.nmrResidueTable._update(nmrResidue.nmrChain)

    finally:
      self.application._endCommandBlock()

  def goToPositionInModules(self, nmrResidue=None, row=None, col=None):
    "Go to the positions defined my NmrAtoms of nmrResidue in the active displays"

    activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()
    self.application._startCommandBlock('application.pickAndAssignModule.goToPositionInModules(nmrResidue)', nmrResidue=nmrResidue)
    try:
      if nmrResidue is not None:
        mainWindow = self.application.ui.mainWindow
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
        self.application.current.nmrResidue = nmrResidue
    finally:
      self.application._endCommandBlock()


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
                                   minimumWidths=(30, 50), maximumWidths=(30, 50),
                                   labelText = axisCode,
                                   value = spectrum.assignmentTolerances[ii],
                                   decimals=decimals, step=step, range=(0, None))
      self.spinBoxes.append(ds)
