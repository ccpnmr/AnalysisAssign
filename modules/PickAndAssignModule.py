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
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:21 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b5 $"
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
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTable, NmrResidueTableModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.CheckBoxes import CheckBoxes
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget, DoubleSpinBoxCompoundWidget
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.NmrResidue import NmrResidue
from ccpn.util.Logging import getLogger
from ccpn.ui.gui.guiSettings import getColours, DIVIDER
from ccpn.ui.gui.widgets.HLine import HLine
from ccpn.core.lib.ContextManagers import undoBlock

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

    includePeakLists = False
    includeNmrChains = False
    includeSpectrumTable = True

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

        self.restrictedPickButton.setEnabled(True)
        self.assignSelectedButton.setEnabled(True)
        self.restrictedPickAndAssignButton.setEnabled(True)
        self.nmrResidueTable.addWidgetToPos(self.restrictedPickAndAssignButton, row=1, col=4)

        # change some of the defaults setting inherited from NmrResidueTableModule
        self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.setChecked(False)

        if self.nmrResidueTableSettings.displaysWidget:
            self.nmrResidueTableSettings.displaysWidget.addPulldownItem(0)  # select the <all> option

        self.nmrResidueTableSettings.setLabelText('Navigate to\nDisplay(s):')
        self.nmrResidueTable._setWidgetHeight(45)

        # need to feedback to current.nmrResidueTable
        self._selectOnTableCurrentNmrResiduesNotifier = None
        self._registerNotifiers()

        # # these need to change whenever different spectrumDisplays are selected
        # self._setAxisCodes()
        if self.nmrResidueTableSettings.axisCodeOptions:
            self.nmrResidueTableSettings.axisCodeOptions.selectAll()

            # just clear the 'C' axes - this is the usual configuration
            for ii, box in enumerate(self.nmrResidueTableSettings.axisCodeOptions.checkBoxes):
                if box.text().upper().startswith('C'):
                    self.nmrResidueTableSettings.axisCodeOptions.clearIndex(ii)

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
        super()._closeModule()

    def assignSelected(self):
        "Assign current.peaks on the bases of nmrAtoms of current.nmrResidue"

        if self.application.current.nmrResidue is None:
            logger.error('Undefined nmrResidue; select one first before proceeding')
            return

        if len(self.application.current.peaks) == 0:
            logger.error('Undefined peak(s); select one or more before proceeding')
            return

        with undoBlock():
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

                    pValue = peak.position[ii]
                    if isotopeCode in shiftDict.keys():

                        shiftList = set()
                        for shift in shiftDict[isotopeCode]:
                            sValue = shift[1]
                            # pValue = peak.position[ii]
                            if abs(sValue - pValue) <= spectrum.assignmentTolerances[ii]:
                                # peak.assignDimension(spectrum.axisCodes[ii], [shift[0]])
                                shiftList.add(shift[0])

                        if shiftList:
                            peak.assignDimension(spectrum.axisCodes[ii], list(shiftList))

            # self.application.current.peaks = []
            # update the NmrResidue table
            self.nmrResidueTable._update(self.application.current.nmrResidue.nmrChain)

            # reset to the last selected nmrResidue - stops other tables messing up
            self.application.current.nmrResidue = lastNmrResidue

    #TODO:GEERTEN: compact the two routines
    def restrictedPick(self, nmrResidue=None):
        """
        Routine refactored in revision 9381.
     
        Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
        spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
        centre points with tolerances and the y as the long axis to pick the whole region.
        """
        nmrResidue = self.project.getByPid(nmrResidue) if isinstance(nmrResidue, str) else nmrResidue

        if not nmrResidue:
            nmrResidue = self.application.current.nmrResidue

        if nmrResidue is None:
            logger.error('Undefined nmrResidue; select one first before proceeding')
            return

        currentAxisCodes = self.nmrResidueTableSettings.axisCodeOptions.getSelectedText()
        currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptions.getSelectedIndexes()

        with undoBlock():

            peaks = []
            # displays = self._getDisplays()

            gid = self.nmrResidueTableSettings.spectrumDisplayPulldown.getText()
            displays = [self.application.getByGid(gid)]

            validPeakListViews = {}

            # loop through all the selected displays/spectrumViews/peakListViews that are visible
            for dp in displays:

                # ignore undefined displays
                if not dp:
                    continue

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

                axisCodes = [spectrum.axisCodes[self.nmrResidueTableSettings.spectrumIndex[spectrum].index(ii)]
                             for ii in currentAxisCodeIndexes if ii in self.nmrResidueTableSettings.spectrumIndex[spectrum]]

                peakList, pks = PeakList.restrictedPick(peakListView=peakListView,
                                                        axisCodes=axisCodes, nmrResidue=nmrResidue)
                if pks:
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

    # from ccpn.util.decorators import profile
    # @profile
    def restrictedPickAndAssign(self, nmrResidue=None):
        """
        Functionality for beta2 to include the Assign part
         
        Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
        spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
        centre points with tolerances and the y as the long axis to pick the whole region.
        
        Calls assignSelected to assign
        """
        nmrResidue = self.project.getByPid(nmrResidue) if isinstance(nmrResidue, str) else nmrResidue

        if not nmrResidue:
            nmrResidue = self.application.current.nmrResidue
        elif not self.application.current.nmrResidue:
            print('No current nmrResidue')
            return

        with undoBlock():

            self.restrictedPick(nmrResidue)

            # if peaks have been selected then assign them
            if self.application.current.peaks:
                self.assignSelected()

                # notifier for other modules
                nmrResidue._finaliseAction('change')

    def goToPositionInModules(self, nmrResidue=None, row=None, col=None):
        "Go to the positions defined my NmrAtoms of nmrResidue in the active displays"

        nmrResidue = self.project.getByPid(nmrResidue) if isinstance(nmrResidue, str) else nmrResidue

        activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

        with undoBlock():

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
                                                    widths=strip._getCurrentZoomRatio(strip.viewRange()),
                                                    markPositions=(n == 2))
                self.application.current.nmrResidue = nmrResidue

