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
__version__ = "$Revision: 3.0.b3 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Geerten Vuister $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

from PyQt5 import QtGui, QtWidgets

from ccpn.ui.gui.lib import PeakList
from ccpn.ui.gui.lib import Strip
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTable, NmrResidueTableModule
from ccpn.ui.gui.widgets.Base import Base
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.ButtonList import ButtonList
from ccpn.ui.gui.widgets.CheckBoxes import CheckBoxes
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget
from ccpn.ui.gui.widgets.DoubleSpinbox import DoubleSpinbox
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.ScrollArea import ScrollArea
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget, DoubleSpinBoxCompoundWidget
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.NmrResidue import NmrResidue
from ccpn.util.Logging import getLogger
from ccpn.ui.gui.guiSettings import getColours, DIVIDER
from ccpn.ui.gui.widgets.HLine import HLine


logger = getLogger()


class PickAndAssignModule(NmrResidueTableModule):
    """
    Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
    Use settings to define the spectral displays, the active spectra and the tolerances for peak picking
  
    This module closely works with the Atom Selector module
    """
    className = 'PickAndAssignModule'

    includeSettingsWidget = True
    maxSettingsState = 2
    settingsPosition = 'top'
    settingsMinimumSizes = (500, 200)

    def __init__(self, mainWindow, name='Pick and Assign'):

        super(PickAndAssignModule, self).__init__(mainWindow=mainWindow, name=name)  # ejb ='Pick And Assign')

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        # Main widget
        self.restrictedPickButton = Button(text='Restricted\nPick', callback=self.restrictedPick,
                                           setLayout=True, spacing=(0, 0))
        self.nmrResidueTable.addWidgetToPos(self.restrictedPickButton, row=1, col=2)

        self.assignSelectedButton = Button(text='Assign\nSelected', callback=self.assignSelected,
                                           setLayout=True, spacing=(0, 0))
        self.nmrResidueTable.addWidgetToPos(self.assignSelectedButton, row=1, col=3)

        self.restrictedPickAndAssignButton = Button(text='Restricted\nPick and Assign',
                                                    setLayout=True, spacing=(0, 0),
                                                    callback=self.restrictedPickAndAssign)

        self.restrictedPickButton.setEnabled(False)
        self.assignSelectedButton.setEnabled(False)
        self.restrictedPickAndAssignButton.setEnabled(False)
        self.nmrResidueTable.addWidgetToPos(self.restrictedPickAndAssignButton, row=1, col=4)

        # self._spacer = Frame(None, setLayout=True)
        # self._spacer.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        # self.nmrResidueTable.addWidgetToPos(self._spacer, row=1, col=6)

        # change some of the defaults setting inherited from NmrResidueTableModule
        self.sequentialStripsWidget.checkBox.setChecked(False)
        self.displaysWidget.addPulldownItem(0)  # select the <all> option

        # create row's of spectrum information
        self._spectraWidget = Frame(parent=self.settingsWidget,
                                    setLayout=True, showBorder=True, hPolicy='minimal',
                                    grid=(0, 1), gridSpan=(4, 1), vAlign='top', hAlign='left')

        # modifier for atomCode
        spectraRow = 0
        self.atomCodeFrame = Frame(self._spectraWidget, setLayout=True, showBorder=False, fShape='noFrame',
                                   grid=(spectraRow, 0), gridSpan=(1, 4), vAlign='top', hAlign='left')
        self.axisCodeLabel = Label(self.atomCodeFrame, 'Restricted Axes:', grid=(0, 0))
        self.axisCodeOptions = CheckBoxes(self.atomCodeFrame, selectedInd=0, texts=['C'],
                                          callback=self._changeAxisCode, grid=(0, 1))
        # self.nmrResidueTable.addWidgetToPos(self.atomCodeFrame, row=1, col=5)

        spectraRow += 1
        HLine(self._spectraWidget, grid=(spectraRow, 0), gridSpan=(1, 4),
              colour=getColours()[DIVIDER], height=15)

        spectraRow += 1
        Label(self._spectraWidget, 'Spectrum', grid=(spectraRow, 0))
        Label(self._spectraWidget, 'Tolerance', grid=(spectraRow, 1))
        Label(self._spectraWidget, 'Tolerance', grid=(spectraRow, 2))
        Label(self._spectraWidget, 'Tolerance', grid=(spectraRow, 3))
        self.spectraStartRow = spectraRow + 1

        self._spectraWidgets = {}  # spectrum.pid, frame dict to show/hide
        for row, spectrum in enumerate(self.application.project.spectra):

            spectraRow += 1
            # f = _SpectrumRow(parent=self._spectraWidget, setLayout=True, spectrum=spectrum,
            #                  grid=(self.spectraStartRow + row, 0), gridSpan=(1, 1 + len(spectrum.axisCodes)), vAlign='top')
            f = _SpectrumRow(parent=self._spectraWidget, row=spectraRow, col=0, setLayout=True, spectrum=spectrum)

            self._spectraWidgets[spectrum.pid] = f

        # add a spacer in the bottom-right corner to stop everything moving
        rows = self.settingsWidget.layout().rowCount()
        cols = self.settingsWidget.layout().columnCount()
        Spacer(self.settingsWidget, 5, 5,
               QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding,
               grid=(rows, cols), gridSpan=(1, 1))

        self.nmrResidueTable._setWidgetHeight(45)

        # need to feedback to current.nmrResidueTable
        self._selectOnTableCurrentNmrResiduesNotifier = None
        self._registerNotifiers()

        # these need to change whenever different spectrumDisplays are selected
        self._setAxisCodes()
        self.axisCodeOptions.selectAll()

        # just clear the 'C' - assume that this is generally the second checkBox
        self.axisCodeOptions.clearIndex(1)

    def _registerNotifiers(self):
        """
        set up the notifiers
        """
        self._selectOnTableCurrentNmrResiduesNotifier = Notifier(self.current,
                                                                 [Notifier.CURRENT],
                                                                 targetName=NmrResidue._pluralLinkName,
                                                                 callback=self._selectionCallback)

    def _unRegisterNotifiers(self):
        """
        clean up the notifiers
        """
        if self._selectOnTableCurrentNmrResiduesNotifier is not None:
            self._selectOnTableCurrentNmrResiduesNotifier.unRegister()

    def _selectionCallback(self, data):
        """
        enable/disable the pick buttons
        """
        selected = data[Notifier.OBJECT].nmrResidue

        if selected:
            self.restrictedPickButton.setEnabled(True)
            self.assignSelectedButton.setEnabled(True)
            self.restrictedPickAndAssignButton.setEnabled(True)
        else:
            self.restrictedPickButton.setEnabled(False)
            self.assignSelectedButton.setEnabled(False)
            self.restrictedPickAndAssignButton.setEnabled(False)

    def _closeModule(self):
        """
        Unregister notifiers and close module.
        """
        self._unRegisterNotifiers()
        super(PickAndAssignModule, self)._closeModule()

    def assignSelected(self):
        "Assign current.peaks on the bases of nmrAtoms of current.nmrResidue"

        if self.application.current.nmrResidue is None:
            logger.error('Undefined nmrResidue; select one first before proceeding')
            return

        if len(self.application.current.peaks) == 0:
            logger.error('Undefined peak(s); select one or more before proceeding')
            return

        self.application._startCommandBlock('application.pickAndAssignModule.assignSelected()')
        try:
            lastNmrResidue = self.application.current.nmrResidue

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
                            if abs(sValue - pValue) <= spectrum.assignmentTolerances[ii]:
                                peak.assignDimension(spectrum.axisCodes[ii], [shift[0]])

            # self.application.current.peaks = []
            # update the NmrResidue table
            self.nmrResidueTable._update(self.application.current.nmrResidue.nmrChain)

            # reset to the last selected nmrResidue - stops other tables messing up
            self.application.current.nmrResidue = lastNmrResidue

        finally:
            self.application._endCommandBlock()

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

        currentAxisCodes = self.axisCodeOptions.getSelectedText()
        currentAxisCodeIndexes = self.axisCodeOptions.getSelectedIndexes()

        self.application._startCommandBlock('application.pickAndAssignModule.restrictedPick(nmrResidue)',
                                            nmrResidue=nmrResidue)
        try:
            peaks = []
            displays = self._getDisplays()
            validPeakListViews = {}

            # loop through all the selected displays/spectrumViews/peakListViews that are visible
            for dp in displays:
                if dp.strips:
                    for sv in dp.strips[0].spectrumViews:
                        for plv in sv.peakListViews:
                            if plv.isVisible() and sv.isVisible():
                                if plv.peakList not in validPeakListViews:
                                    validPeakListViews[plv.peakList] = (sv.spectrum, plv)
                                else:

                                    # skip for now, only one valid peakListView needed per peakList
                                    # validPeakListViews[plv.peakList] += (plv,)
                                    pass

            for pk, specAndView in validPeakListViews.items():
                spectrum, peakListView = specAndView

                axisCodes = [spectrum.axisCodes[self.spectrumIndex[spectrum].index(ii)]
                             for ii in currentAxisCodeIndexes if ii in self.spectrumIndex[spectrum]]

                peakList, pks = PeakList.restrictedPick(peakListView=peakListView,
                                                        axisCodes=axisCodes, nmrResidue=nmrResidue)
                peaks = peaks + pks

            # for module in self.application.project.spectrumDisplays:
            #     if len(module.axisCodes) >= 2:
            #         for spectrumView in module.strips[0].spectrumViews:
            #
            #             visiblePeakListViews = [peakListView for peakListView in spectrumView.peakListViews
            #                                     if peakListView.isVisible()]
            #
            #             # if len(visiblePeakListViews) == 0:
            #             #     continue
            #             # else:
            #             #     peakList, pks = PeakList.restrictedPick(peakListView=visiblePeakListViews[0],
            #             #                                             axisCodes=module.axisCodes[0::2], nmrResidue=nmrResidue)
            #             #     peaks = peaks + pks
            #
            #             # if len(visiblePeakListViews) == 0:
            #             #     spectrum = spectrumView.spectrum
            #             #
            #             #     axisCodes = [axis for ]

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
            self.restrictedPick(nmrResidue)

            # if peaks have been selected then assign them
            if self.application.current.peaks:
                self.assignSelected()

                # notifier for other modules
                nmrResidue._finaliseAction('change')

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
                        widths = ['default', 'full'] + (n - 2) * ['']

                    Strip.navigateToNmrAtomsInStrip(strip=strip,
                                                    nmrAtoms=nmrResidue.nmrAtoms,
                                                    widths=strip._getCurrentZoomRatio(strip.viewBox.viewRange()),
                                                    markPositions=(n == 2))
                self.application.current.nmrResidue = nmrResidue
        finally:
            self.application._endCommandBlock()

    def _changeAxisCode(self):
        pass

    def _setAxisCodes(self):

        import difflib
        from ccpn.util.Common import _axisCodeMapIndices, axisCodeMapping

        spectra = self.application.project.spectra

        if spectra:

            maxLen = 0
            refAxisCodes = None
            for spectrum in spectra:
                if len(spectrum.axisCodes) > maxLen:
                    maxLen = len(spectrum.axisCodes)
                    refAxisCodes = list(spectrum.axisCodes)

            if not maxLen:
                return

            axisCodes = [[] for ii in range(maxLen)]
            axisLabels = [set() for ii in range(maxLen)]

            mappings = {}
            for spectrum in spectra:
                matchAxisCodes = spectrum.axisCodes

                mapping = axisCodeMapping(refAxisCodes, matchAxisCodes)
                for k, v in mapping.items():
                    if v not in mappings:
                        mappings[v] = set([k])
                    else:
                        mappings[v].add(k)

                mapping = axisCodeMapping(matchAxisCodes, refAxisCodes)
                for k, v in mapping.items():
                    if v not in mappings:
                        mappings[v] = set([k])
                    else:
                        mappings[v].add(k)

                # example of mappings dict
                # ('Hn', 'C', 'Nh')
                # {'Hn': {'Hn'}, 'Nh': {'Nh'}, 'C': {'C'}}
                # {'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'C'}}
                # {'CA': {'C'}, 'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'CA', 'C'}}
                # {'CA': {'C'}, 'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'CA', 'C'}}

                self.spectrumIndex = {}
                # go through the spectra
                for spectrum in spectra:
                    self.spectrumIndex[spectrum] = [0 for ii in range(len(spectrum.axisCodes))]

                    # get the spectrum dimension axisCode, nd see if is already there
                    for spectrumDim, spectrumAxis in enumerate(spectrum.axisCodes):

                        if spectrumAxis in refAxisCodes:
                            self.spectrumIndex[spectrum][spectrumDim] = refAxisCodes.index(spectrumAxis)
                            axisLabels[self.spectrumIndex[spectrum][spectrumDim]].add(spectrumAxis)

                        else:
                            # if the axisCode is not in the reference list then find the mapping from the dict
                            for k, v in mappings.items():
                                if spectrumAxis in v:
                                    # refAxisCodes[dim] = k
                                    self.spectrumIndex[spectrum][spectrumDim] = refAxisCodes.index(k)
                                    axisLabels[refAxisCodes.index(k)].add(spectrumAxis)

            axisLabels = [', '.join(ax) for ax in axisLabels]
            self.axisCodeOptions.setCheckBoxes(texts=axisLabels, tipTexts=axisLabels)

    def _getValidAxisCode(self, numChars=1):
        """Get the valid axis code from the buttons, numChars is included as this may be needed for DNA/RNA
        """
        code = self.axisCodeOptions.getSelectedText()
        return code[0:numChars]


class _SpectrumRow(Frame):
    "Class to make a spectrum row"

    def __init__(self, parent, spectrum, row=0, col=0, **kwds):
        super(_SpectrumRow, self).__init__(parent, **kwds)

        # col = 0
        # self.checkbox = CheckBoxCompoundWidget(self, grid=(0, col), gridSpan=(1, 1), hAlign='left',
        #                                        checked=True, labelText=spectrum.pid,
        #                                        fixedWidths=[100, 50])

        self.checkbox = Label(parent, spectrum.pid, grid=(row, col), gridSpan=(1, 1), hAlign='left')

        self.spinBoxes = []
        for ii, axisCode in enumerate(spectrum.axisCodes):
            decimals, step = (2, 0.01) if axisCode[0:1] == 'H' else (1, 0.1)
            col += 1
            ds = DoubleSpinBoxCompoundWidget(
                    parent, grid=(row, col), gridSpan=(1, 1), hAlign='left',
                    fixedWidths=(30, 50),
                    labelText=axisCode,
                    value=spectrum.assignmentTolerances[ii],
                    decimals=decimals, step=step, range=(0, None))
            ds.setObjectName(str(spectrum.pid + axisCode))
            self.spinBoxes.append(ds)
