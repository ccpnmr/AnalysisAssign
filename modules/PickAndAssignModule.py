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
from PyQt4 import QtGui

from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.Base import Base
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Module import CcpnModule
from ccpn.ui.gui.widgets.DoubleSpinbox import DoubleSpinbox
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.ScrollArea import ScrollArea

from ccpn.ui.gui.lib.PeakList import restrictedPick
from ccpn.ui.gui.lib.Strip import navigateToNmrAtomsInStrip

from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTable

class PickAndAssignModule(CcpnModule, Base):

  def __init__(self, parent=None, project=None, name='Pick And Assign', **kw):

    CcpnModule.__init__(self, parent=None, name=name)
    Base.__init__(self, **kw)
    self.project = project
    self.current = project._appBase.current
    self.nmrResidueTable = NmrResidueTable(self.widget1, project=project, callback=self._goToPositionInModules, grid=(0, 0), gridSpan=(1, 5), stretch=(1, 1))
    self.restrictedPickButton = Button(self.nmrResidueTable, text='Restricted Pick', callback=self._restrictedPick, grid=(0, 2))
    self.assignSelectedButton = Button(self.nmrResidueTable, text='Assign Selected', callback=self._assignSelected, grid=(0, 3))
    self.refreshButton = Button(self.nmrResidueTable, text='Refresh', callback=self._refresh, grid=(0, 4))
    self.settingsButton = Button(self.nmrResidueTable, icon='icons/applications-system', grid=(0, 5), hPolicy='fixed', toggle=True)
    self.settingsButton.toggled.connect(self._toggleWidget2)
    self.settingsButton.setChecked(False)
    displaysLabel = Label(self.widget2, 'Selected Displays', grid=(0, 0))
    self.displaysPulldown = PulldownList(self.widget2, grid=(1, 0), callback=self._updateListWidget)
    self.displaysPulldown.setData([sd.pid for sd in project.spectrumDisplays])
    self.displayList = ListWidget(self.widget2, grid=(0, 1), gridSpan=(4, 1))
    self.displayList.addItem('<All>')
    self.displayList.setFixedWidth(self.displaysPulldown.width())
    self.scrollArea = ScrollArea(self.widget2, grid=(0, 2), gridSpan=(4, 4))
    self.spectrumSelectionWidget = SpectrumSelectionWidget(self.scrollArea, project, self.displayList)
    self.scrollArea.setWidget(self.spectrumSelectionWidget)
    self.displayList.removeItem = self._removeListWidgetItem
    self.refreshButton.hide()

    self.project.registerNotifier('NmrResidue', 'rename', self._updateListWidget)
    self.project.registerNotifier('NmrResidue', 'rename', self._updateNmrResidueTable)

    self.closeModule = self._closeModule

  def _closeModule(self):
    self.project.registerNotifier('NmrResidue', 'rename', self._updateListWidget)
    self.project.registerNotifier('NmrResidue', 'rename', self._updateNmrResidueTable)
    self.close()

  def _updateListWidget(self, item):
    if self.displayList.count() == 1 and self.displayList.item(0).text() == '<All>':
      self.displayList.takeItem(0)
    self.displayList.addItem(self.project.getByPid(item).pid)
    self.spectrumSelectionWidget.update()

  def _updateNmrResidueTable(self, nmrResidue, oldPid):
    self.nmrResidueTable.nmrResidueTable.updateTable()
    self.nmrResidueTable.nmrResidueTable._updateSelectorContents()


  def _removeListWidgetItem(self):
    self.displayList.takeItem(self.displayList.currentRow())

    if self.displayList.count() == 0:
      self.displayList.addItem('<All>')
    self.spectrumSelectionWidget.update()

  def _toggleWidget2(self):
    if self.settingsButton.isChecked():
      self.widget2.show()
    else:
      self.widget2.hide()

  def _refresh(self):
    pass

  def _assignSelected(self):

    self.project._appBase._startCommandBlock('application.pickAndAssignModule._assignSelected()')
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
    finally:
      self.project._appBase._endCommandBlock()


  def _restrictedPick(self, nmrResidue=None):
    """
    Routine refactored in revision 9381.
    Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
    spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
    centre points with tolerances and the y as the long axis to pick the whole region.
    """
    if not nmrResidue:
      nmrResidue = self.current.nmrResidue
    elif not self.current.nmrResidue:
      print('No current nmrResidue')
      return

    self.project._appBase._startCommandBlock('application.pickAndAssignModule._restrictedPick(nmrResidue)', nmrResidue=nmrResidue)
    try:
      for module in self.project.spectrumDisplays:
        if len(module.axisCodes) > 2:
          for spectrumView in module.strips[0].spectrumViews:
            visiblePeakListViews = [peakListView for peakListView in spectrumView.peakListViews
                                    if peakListView.isVisible()]
            if len(visiblePeakListViews) == 0:
              return
            else:
              peakList, peaks = restrictedPick(peakListView=visiblePeakListViews[0],
                                               axisCodes=module.axisCodes[0::2], nmrResidue=nmrResidue)
              # Lines below here to be removed when a notifier handles display of newly picked peaks
              if len(peaks) > 0:
                for strip in module.strips:
                  strip.showPeaks(peakList)
                for peakListView in module.peakListViews:
                  peakItems = [peakListView.peakItems[peak] for peak in peaks if peak in peakListView.peakItems.keys()]
                  for peakItem in peakItems:
                    peakItem.isSelected = True
    finally:
      self.project._appBase._endCommandBlock()


  def _goToPositionInModules(self, nmrResidue=None, row=None, col=None):

    self.project._appBase._startCommandBlock('application.pickAndAssignModule._goToPositionInModules(nmrResidue)', nmrResidue=nmrResidue)
    try:
      if self.project._appBase.ui.mainWindow is not None:
        mainWindow = self.project._appBase.ui.mainWindow
      else:
        mainWindow = self.project._appBase._mainWindow
      mainWindow.clearMarks()
      for displayPid in self.displayList.items():
        strip = self.project.getByPid(displayPid).strips[0]
        navigateToNmrAtomsInStrip(strip=strip, nmrAtoms=nmrResidue.nmrAtoms, widths=['default', 'full', ''])
      self.current.nmrResidue = nmrResidue
    finally:
      self.project._appBase._endCommandBlock()

class SpectrumSelectionWidget(QtGui.QWidget, Base):

  def __init__(self, parent, project, displayList, **kw):
    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, **kw)
    self.project = project
    self.spectrumLabel = Label(self, 'Spectrum')
    self.layout().addWidget(self.spectrumLabel, 0, 0)
    self.useLabel = Label(self, 'Use?', grid=(0, 1), hAlign='c')
    self.refreshBox = CheckBox(self, grid=(0, 4))
    self.checkBoxLabel = Label(self, 'Auto Refresh', grid=(0, 5))
    self.displayList = displayList
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
    activeDisplays = self.getActiveDisplays()
    spectra = set([spectrumView.spectrum for spectrumDisplay in activeDisplays
                   for spectrumView in spectrumDisplay.spectrumViews])


    for ii, spectrum in enumerate(spectra):
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




