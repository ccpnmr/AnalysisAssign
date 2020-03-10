"""Module Documentation here

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
__version__ = "$Revision: 3.0.0 $"
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
from ccpn.core.NmrChain import NmrChain
from ccpn.ui.gui.lib.SpectrumDisplay import makeStripPlot

from ccpn.ui.gui.lib.Strip import matchAxesAndNmrAtoms
from ccpn.ui.gui.lib.Strip import navigateToNmrResidueInDisplay

from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTableModule

from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.CompoundWidgets import ListCompoundWidget, PulldownListCompoundWidget
from ccpn.ui.gui.widgets.MessageDialog import showWarning, progressManager, showYesNo
from ccpn.ui.gui.widgets.PulldownListsForObjects import ChemicalShiftListPulldown
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.lib.GuiNotifier import GuiNotifier
from ccpn.ui.gui.widgets.DropBase import DropBase

from ccpn.util.Logging import getLogger
from ccpn.core.NmrAtom import NmrAtom
from ccpn.ui.gui.widgets.PlaneToolbar import STRIPLABEL_CONNECTDIR, STRIPLABEL_CONNECTNONE, \
    STRIPCONNECT_LEFT, STRIPCONNECT_RIGHT
from ccpn.core.lib.ContextManagers import undoBlock


ALL = '<all>'
MINMATCHES = 1
MAXMATCHES = 7
DEFAULTMATCHES = 2
STRIPBACKBONE = 'backboneAssignment'
MARKCONNECTED = False


class BackboneAssignmentModule(NmrResidueTableModule):
    className = 'BackboneAssignmentModule'

    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'left'
    settingsMinimumSizes = (500, 200)

    includeDisplaySettings = True
    activePulldownClass = NmrChain

    def __init__(self, mainWindow=None, name='Backbone Assignment'):

        super(BackboneAssignmentModule, self).__init__(mainWindow=mainWindow, name=name)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        if mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current
            self.nmrChains = self.application.project.nmrChains

        self.matchCheckBoxWidget = CheckBox(self.nmrResidueTable._widget,
                                            grid=(1, 2), checked=True, text='Find matches')

        ### Settings ###

        # change some of the defaults setting inherited from NmrResidueTableModule
        self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.setChecked(True)
        if self.nmrResidueTableSettings.displaysWidget:
            self.nmrResidueTableSettings.displaysWidget.addPulldownItem(0)

        colWidth0 = 180
        colWidth = 200  # for labels of the compound widgets
        colWidth2 = 120  # for the numberOfMatchesWidget

        row = self.nmrResidueTableSettings.maxRows  ## Number of widgets of NmrResidueTable - add extra widgets below
        col = 0

        # Number of matches to show
        row += 1
        self.numberOfMinusMatchesWidget = PulldownListCompoundWidget(self.nmrResidueTableSettings,
                                                                     grid=(row, col), vAlign='top', hAlign='left',
                                                                     fixedWidths=(colWidth0, colWidth0, None),
                                                                     orientation='left',
                                                                     labelText="i-1 Matches to show:",
                                                                     texts=[str(tt) for tt in range(MINMATCHES, MAXMATCHES)],
                                                                     default=DEFAULTMATCHES
                                                                     )
        row += 1
        self.numberOfPlusMatchesWidget = PulldownListCompoundWidget(self.nmrResidueTableSettings,
                                                                    grid=(row, col), vAlign='top', hAlign='left',
                                                                    fixedWidths=(colWidth0, colWidth0, None),
                                                                    orientation='left',
                                                                    labelText="i+1 Matches to show:",
                                                                    texts=[str(tt) for tt in range(MINMATCHES, MAXMATCHES)],
                                                                    default=DEFAULTMATCHES
                                                                    )

        # Match module selection
        # row += 1
        # # cannot set a notifier for displays, as these are not (yet?) implemented
        # self.matchWidget = ListCompoundWidget(self.nmrResidueTableSettings,
        #                                       grid=(row, col), vAlign='top', hAlign='left',
        #                                       fixedWidths=(colWidth0, colWidth0, colWidth0),
        #                                       orientation='left',
        #                                       labelText="Match module(s):",
        #                                       texts=[display.pid for display in self.mainWindow.spectrumDisplays]
        #                                       )
        # self.matchWidget.setPreSelect(self._fillDisplayWidget)
        # self.matchWidget.setFixedHeights((None, None, 40))

        # new match module pulldown list
        row += 1
        self.matchWidget = PulldownListCompoundWidget(self.nmrResidueTableSettings, labelText="Match module:",
                                                      fixedWidths=(colWidth0, colWidth0, None), grid=(row, col), gridSpan=(1, 2))
        self.matchWidget.setPreSelect(self._fillMatchWidget)
        self._fillMatchWidget()
        self.matchWidget.pulldownList.setIndex(0)

        # new target module pulldown list
        row += 1
        self.targetWidget = PulldownListCompoundWidget(self.nmrResidueTableSettings, labelText="Target module:",
                                                       fixedWidths=(colWidth0, colWidth0, None), grid=(row, col), gridSpan=(1, 2))
        self.targetWidget.setPreSelect(self._fillTargetWidget)
        self._fillTargetWidget()
        self.targetWidget.pulldownList.setIndex(0)

        # Chemical shift list selection
        row += 1
        self.shiftListWidget = ChemicalShiftListPulldown(self.nmrResidueTableSettings, self.mainWindow,
                                                         grid=(row, col), vAlign='top', hAlign='left',
                                                         fixedWidths=(colWidth0, colWidth0, None),
                                                         callback=self._setupShiftDicts, default=None
                                                         )
        self._setupShiftDicts()
        self._spacer = Spacer(self.settingsWidget, 5, 5,
                              QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding,
                              grid=(row + 20, 10), gridSpan=(1, 1))

        # for compatibility with previous implementation
        #self.moduleList = self.matchWidget.listWidget

        self._stripNotifiers = []  # list to store GuiNotifiers for strips
        self.nmrResidueTable.multiSelect = True
        self.nmrResidueTable.setSelectionMode(self.nmrResidueTable.SingleSelection)
        self.nmrResidueTable.setStyleSheet('''
                    NmrResidueTable {border: 1px solid  #a9a9a9;
                    border-bottom-right-radius: 2px;
                    border-bottom-left-radius: 2px;}
                    ''')

        corner = QtGui.QWidget()
        corner.setStyleSheet('''
            border-top: 1px solid #a9a9a9;
            border-left: 1px solid #a9a9a9;
        ''')
        self.nmrResidueTable.setCornerWidget(corner)

        #self.nmrResidueTable._setWidgetHeight(48)

        self.mainWidget.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        self.layout.setContentsMargins(0, 1, 0, 0)

    def _fillMatchWidget(self):
        ll = ['> select-to-add <'] + [display.pid for display in self.mainWindow.spectrumDisplays]
        thisText = self.matchWidget.getText()
        self.matchWidget.pulldownList.setData(texts=ll)
        if thisText:
            self.matchWidget.select(thisText, True)

    def _fillTargetWidget(self):
        ll = ['> select-to-add <'] + [display.pid for display in self.mainWindow.spectrumDisplays]
        thisText = self.targetWidget.getText()
        self.targetWidget.pulldownList.setData(texts=ll)
        if thisText:
            self.targetWidget.select(thisText, True)

    def _getDisplays(self):
        "return list of displays to navigate"
        displays = []

        if self.nmrResidueTableSettings.displaysWidget:
            dGids = self.nmrResidueTableSettings.displaysWidget.getTexts()  # gid's of displays
            if len(dGids) == 0: return displays

            matchGids = self.matchWidget.getText()  # gid of the match module
            targetGids = None  #.targetWidget.getText()         # gid of the targets module - don't discard for the minute

            if ALL in dGids:
                displays = [dp for dp in self.application.ui.mainWindow.spectrumDisplays if dp.pid not in (matchGids, targetGids)]
            else:
                displays = [self.application.getByGid(gid) for gid in dGids if (gid != ALL and gid not in (matchGids, targetGids))]

        return displays

    def _getMatchDisplays(self):
        """return list of displays to display matches
        """
        mGids = self.matchWidget.getText()  # gid of the match displays
        displays = [self.application.getByGid(gid) for gid in (mGids,)]
        return displays

    def _getTargetDisplays(self):
        """return list of displays to display targets
        """
        mGids = self.targetWidget.getText()  # gid of the match displays
        displays = [self.application.getByGid(gid) for gid in (mGids,)]
        return displays

    def navigateToNmrResidueCallBack(self, data):
        """
        Navigate in selected displays to nmrResidue; skip if none defined
        """
        from ccpn.core.lib.CallBack import CallBack

        # nmrResidue = data[CallBack.OBJECT]
        # if not nmrResidue:
        #     return
        # if isinstance(nmrResidue, (tuple, list)):
        #     nmrResidue = nmrResidue[0]

        nmrResidue = data[CallBack.ROWOBJECT]  # the item clicked, not everything selected
        row = data[CallBack.ROW]
        col = data[CallBack.COL]
        self.navigateToNmrResidue(nmrResidue, row=row, col=col)

    def navigateToNmrResidue(self, nmrResidue, row=None, col=None):
        """
        Navigate in selected displays to nmrResidue; skip if no displays defined
        If matchCheckbox is checked, also call findAndDisplayMatches
        """
        displays = self._getDisplays()
        if len(displays) == 0 and self.nmrResidueTableSettings.displaysWidget:
            getLogger().warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        matchIndex = self.matchWidget.getIndex()
        targetIndex = self.targetWidget.getIndex()
        if self.matchCheckBoxWidget.isChecked() and matchIndex == 0:
            getLogger().warning('Undefined match module; select in settings first or unselect "Find matches"')
            showWarning('startAssignment', 'Undefined match module;\nselect in settings first or unselect "Find matches"')
            return

        if self.matchCheckBoxWidget.isChecked() and targetIndex == 0:
            getLogger().warning('Undefined target module; select in settings first or unselect "Find matches"')
            showWarning('startAssignment', 'Undefined target module;\nselect in settings first or unselect "Find matches"')
            return

        if (matchIndex == targetIndex) and matchIndex != 0:
            getLogger().warning('Match module and Target module cannot be the same')
            showWarning('startAssignment', 'Match module and Target module cannot be the same')
            return

        with undoBlock():
            # self.application._startCommandBlock(
            #         'BackboneAssignmentModule.navigateToNmrResidue(project.getByPid(%r))' % nmrResidue.pid)
            # try:

            # optionally clear the marks
            if self.nmrResidueTableSettings.autoClearMarksWidget.checkBox.isChecked():
                self.mainWindow.clearMarks()

            # clear any notifiers of previous strips
            for notifier in self._stripNotifiers:
                notifier.unRegister()
                del (notifier)
            self._stripNotifiers = []

            nr = nmrResidue.mainNmrResidue
            targetDisplays = self._getTargetDisplays()

            # navigate the displays
            for display in displays:

                display.showAllStripHeaders()  # tag all headers with backboneAssignment module as handler

                if len(display.strips) > 0:

                    # if contains 2D's (e.g. a hsqc) then keep zoom
                    if display.spectrumViews[0].spectrum.dimensionCount <= 2:
                        newWidths = []  #_getCurrentZoomRatio(display.strips[0].viewBox.viewRange())
                    else:
                        # set the width in case of nD (n>2)
                        _widths = {'H': 2.5, 'C': 1.0, 'N': 1.0}
                        _ac = display.strips[0].axisCodes[0]
                        _w = _widths.setdefault(_ac[0], 1.0)
                        newWidths = [_w, 'full']

                    strips = navigateToNmrResidueInDisplay(nr, display, stripIndex=0,
                                                           widths=newWidths,
                                                           showSequentialResidues=(len(display.axisCodes) > 2) and
                                                                                  self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.isChecked(),
                                                           markPositions=False,  #self.nmrResidueTableSettings.markPositionsWidget.checkBox.isChecked()
                                                           showDropHeaders=display in targetDisplays,
                                                           )

                    # activate a callback notifiers; allow dropping onto the NmrResidueLabel
                    for st, strip in enumerate(strips):
                        if strip is not None:
                            # NB connections are made as connectPrevious / connectNext to passed-in NmrResidue
                            # It follows that it IS the mainNMr Residue that should be passed in here
                            # Note, though, that you get the same connections WHICHEVER strip you drop on
                            # label = strip.header.getLabel(position='l')
                            # notifier = GuiNotifier(label,  #getStripLabel(),
                            #                        [GuiNotifier.DROPEVENT], [DropBase.TEXT],
                            #                        self._processDroppedNmrResidueLabel,
                            #                        toLabel=strip.header.getLabel(position='c'),
                            #                        plusChain=False)
                            # self._stripNotifiers.append(notifier)
                            # label = strip.header.getLabel(position='r')
                            # notifier = GuiNotifier(label,  #getStripLabel(),
                            #                        [GuiNotifier.DROPEVENT], [DropBase.TEXT],
                            #                        self._processDroppedNmrResidueLabel,
                            #                        toLabel=strip.header.getLabel(position='c'),
                            #                        plusChain=True)
                            # self._stripNotifiers.append(notifier)

                            strip.header.handle = STRIPBACKBONE
                            strip.header.headerVisible = True

                        strip.spectrumDisplay.setColumnStretches(True)

                    # layout.setColumnStretch(col, colStr
                    # strips[0].spectrumDisplay.stripFrame.setStretch(1,1)

            # ejb
            # if 'i-1' residue, take CA CB, and take H, N from the 'i' residue (.mainNmrResidue)
            # check if contains '-1' in pid, is this robust? no :)

            if self.nmrResidueTableSettings.markPositionsWidget.checkBox.isChecked():
                if nmrResidue.relativeOffset is not None and nmrResidue.relativeOffset != 0:
                    # -1, +1 residue so need to split the CA, CB from the N, H
                    nmrAtomsOffset = nmrAtomsFromResidue(nmrResidue)
                    nmrAtomsCentre = nmrAtomsFromResidue(nmrResidue.mainNmrResidue)

                    nmrAtoms = []
                    # this should check the experiment type and choose the correct atoms from there
                    for naOffset in nmrAtomsOffset:
                        if naOffset.name.startswith('CA') or naOffset.name.startswith('CB'):
                            nmrAtoms.append(naOffset)
                    for naCentre in nmrAtomsCentre:
                        if naCentre.name.startswith('N') or naCentre.name.startswith('H'):
                            nmrAtoms.append(naCentre)

                    markNmrAtoms(mainWindow=self.mainWindow, nmrAtoms=nmrAtoms)
                else:
                    if MARKCONNECTED:
                        nmrAtoms = nmrAtomsFromResidue(nmrResidue.mainNmrResidue)
                    else:
                        nmrAtoms = nmrResidue.mainNmrResidue.nmrAtoms
                    markNmrAtoms(mainWindow=self.mainWindow, nmrAtoms=nmrAtoms)

            if self.matchCheckBoxWidget.isChecked():
                self.findAndDisplayMatches(nmrResidue)

            # # update current (should trigger SequenceGraph)
            # self.application.current.nmrChain = nmrResidue.nmrChain
            # self.application.current.nmrResidue = nmrResidue

        # finally:
        #     self.application._endCommandBlock()

    def findAndDisplayMatches(self, nmrResidue):
        "Find and displays the matches to nmrResidue"

        # If NmrResidue is a -1 offset NmrResidue, set queryShifts as value from self.interShifts dictionary
        # Set matchShifts as self.intraShifts
        # if nmrResidue.sequenceCode.endswith('-1'):
        if nmrResidue.relativeOffset == -1:
            # direction = '-1'
            # iNmrResidue = nmrResidue.mainNmrResidue

            if nmrResidue not in self.interShifts:
                queryShifts = []
            else:
                queryShifts = [shift for shift in self.interShifts[nmrResidue]
                               if shift.nmrAtom.isotopeCode == '13C']
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

            if nmrResidue not in self.intraShifts:
                queryShifts = []
            else:
                queryShifts = [shift for shift in self.intraShifts[nmrResidue]
                               if shift.nmrAtom.isotopeCode == '13C']
            matchShifts = self.interShifts

        if True:
            queryShifts = [shift for shift in self.allShifts[nmrResidue]
                           if shift.nmrAtom.isotopeCode == '13C']
            assignMatrix = getNmrResidueMatches(queryShifts, self.allShifts, 'averageQScore')
            print('assignMatrix',assignMatrix)
        else:
            assignMatrix = getNmrResidueMatches(queryShifts, matchShifts, 'averageQScore')

        if not assignMatrix.values():
            getLogger().info('No matches found for NmrResidue: %s' % nmrResidue.pid)
            return
        self._createMatchStrips(assignMatrix)

    def _processDroppedNmrResidrueLabel(self, data, toLabel=None, plusChain=None):
        if toLabel and toLabel.obj:
            self._processDroppedNmrResidue(data, toLabel.obj, plusChain)

    def _processDroppedNmrResidue(self, data, nmrResidue, plusChain=None):
        """Process the dropped NmrResidue id"""

        droppedNmrResidue = None
        if DropBase.TEXT in data and len(data[DropBase.TEXT]) > 0:
            droppedNmrResidue = self.application.project.getByPid(data[DropBase.TEXT])
        if droppedNmrResidue is None:
            showWarning(str(self.windowTitle()), 'Backbone assignment: invalid dropped item')
            getLogger().warning('Backbone assignment: invalid "pid" of dropped item')
            return

        if not isinstance(droppedNmrResidue, NmrResidue):
            showWarning(str(self.windowTitle()), 'Backbone assignment: item is not an nmrResidue')
            getLogger().warning('Backbone assignment: item is not an nmrResidue')
            return

        getLogger().debug('nmrResidue:%s, droppedNmrResidue:%s', nmrResidue, droppedNmrResidue)
        if droppedNmrResidue == nmrResidue:
            return

        allNmrResidues = nmrResidue._getAllConnectedList()

        isPlus = data[STRIPLABEL_CONNECTDIR] if STRIPLABEL_CONNECTDIR in data else STRIPLABEL_CONNECTNONE
        if isPlus == STRIPCONNECT_RIGHT:
            data['shiftLeftMouse'] = False
        elif isPlus == STRIPCONNECT_LEFT:
            data['shiftLeftMouse'] = True
        else:
            data['shiftLeftMouse'] = None

        index = allNmrResidues.index(nmrResidue)
        lenNmr = len(allNmrResidues) - 1

        if data['shiftLeftMouse'] and not plusChain and index == 0:
            # okay to connect to left
            okay = True

        elif not data['shiftLeftMouse'] and plusChain and index == lenNmr:
            # okay to connect to right
            okay = True

        elif data['shiftLeftMouse'] and plusChain and index == lenNmr:
            # check connecting i-1 nmrResidue to the right
            yesNo = showYesNo(str(self.windowTitle()), "Trying to connect 'i-1' nmrResidue to end of chain.\n\n"
                                                       "Do you want to continue?")
            getLogger().warning("Trying to connect 'i-1' nmrResidue to end of chain")
            if not yesNo:
                return

            # force the connection to the start of chain
            data['shiftLeftMouse'] = False
            okay = True

        elif not data['shiftLeftMouse'] and not plusChain and index == 0:
            # check connecting i+1 nmrResidue to the left
            yesNo = showYesNo(str(self.windowTitle()), "Trying to connect 'i+1' nmrResidue to start of chain.\n\n"
                                                       "Do you want to continue?")
            getLogger().warning("Trying to connect 'i+1' nmrResidue to start of chain")
            if not yesNo:
                return

            # force the connection to the end of chain
            data['shiftLeftMouse'] = True
            okay = True

        else:
            # connecting to the middle of a stretch - may do disconnect later
            showWarning(str(self.windowTitle()), "Illegal connection, cannot connect to the middle of a chain")
            getLogger().warning("Illegal connection, cannot connect to the middle of a chain")
            return

        # silence the update of the nmrResidueTable as we will to an explicit update later
        # put in try/finally block because otherwise if exception thrown in the following code
        # (which can happen) then you no longer get updates of the NmrResidue table

        if data['shiftLeftMouse']:
            progressText = "connecting  %s  >  %s" % (droppedNmrResidue.pid, nmrResidue.pid)
        else:
            progressText = "connecting  %s  <  %s" % (nmrResidue.pid, droppedNmrResidue.pid)

        with undoBlock():
            with progressManager(self.mainWindow, progressText):

                try:
                    matchNmrResidue = None
                    try:  # display popup warning
                        if data['shiftLeftMouse']:
                            # leftShift drag; connect to previous

                            if not nmrResidue.residue and not droppedNmrResidue.residue:
                                nmrResidue.connectPrevious(droppedNmrResidue)

                            # SPECIAL CASES
                            elif nmrResidue.residue and not droppedNmrResidue.residue:
                                # connected an unassigned nmrChain to the current assigned chain
                                if droppedNmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    droppedNmrResidue.nmrChain.assignSingleResidue(droppedNmrResidue, nmrResidue.residue.previousResidue)
                                else:
                                    nRes = nmrResidue.residue
                                    for ii in range(len(droppedNmrResidue.nmrChain.mainNmrResidues)):
                                        nRes = nRes.previousResidue
                                    droppedNmrResidue.nmrChain.assignConnectedResidues(nRes)

                            elif not nmrResidue.residue and droppedNmrResidue.residue:
                                # connected an assigned chain to the current unassigned nmrChain

                                if nmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    nmrResidue.nmrChain.assignSingleResidue(nmrResidue, droppedNmrResidue.residue.nextResidue)
                                else:
                                    nmrResidue.nmrChain.assignConnectedResidues(droppedNmrResidue.residue.nexResidue)

                            matchNmrResidue = droppedNmrResidue.getOffsetNmrResidue(offset=-1)
                            if matchNmrResidue is None:
                                # Non -1 residue - stay with current
                                getLogger().info("NmrResidue %s has no i-1 residue to display" % droppedNmrResidue)
                                matchNmrResidue = nmrResidue
                        else:

                            if not nmrResidue.residue and not droppedNmrResidue.residue:
                                nmrResidue.connectNext(droppedNmrResidue)

                            # SPECIAL CASES
                            elif nmrResidue.residue and not droppedNmrResidue.residue:
                                # connected an unassigned nmrChain to the current assigned chain
                                if droppedNmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    droppedNmrResidue.nmrChain.assignSingleResidue(droppedNmrResidue, nmrResidue.residue.nextResidue)
                                else:
                                    droppedNmrResidue.nmrChain.assignConnectedResidues(nmrResidue.residue.nextResidue)

                            elif not nmrResidue.residue and droppedNmrResidue.residue:
                                # connected an assigned chain to the current unassigned nmrChain

                                if nmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    nmrResidue.nmrChain.assignSingleResidue(nmrResidue, droppedNmrResidue.residue.previousResidue)
                                else:
                                    dropRes = droppedNmrResidue.residue
                                    for ii in range(len(nmrResidue.nmrChain.mainNmrResidues)):
                                        dropRes = dropRes.previousResidue
                                    nmrResidue.nmrChain.assignConnectedResidues(dropRes)

                            matchNmrResidue = droppedNmrResidue

                    except Exception as es:
                        showWarning('Connect NmrResidue', str(es))
                    finally:
                        if matchNmrResidue:
                            self.navigateToNmrResidue(matchNmrResidue)

                    # # update the NmrResidueTable
                    # getLogger().info('>>>DISPLAYTABLE', droppedNmrResidue.nmrChain, self.project.nmrChains)
                    # self.nmrResidueTable.displayTableForNmrChain(droppedNmrResidue.nmrChain)

                    # from ccpn.ui.gui.lib.OpenGL.CcpnOpenGL import GLNotifier
                    #
                    # GLSignals = GLNotifier(parent=self)
                    # GLSignals.emitEvent(triggers=[GLNotifier.GLMARKS])

                except Exception as es:
                    getLogger().warning(str(es))
                    # raise es
                finally:
                    pass

        # update the NmrResidueTable - outside of the undoBlock for notifiers to catch up
        self.nmrResidueTable.displayTableForNmrChain(droppedNmrResidue.nmrChain)

        from ccpn.ui.gui.lib.OpenGL.CcpnOpenGL import GLNotifier

        GLSignals = GLNotifier(parent=self)
        GLSignals.emitEvent(triggers=[GLNotifier.GLMARKS])

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
            yPosition = (max(yShiftValues) + min(yShiftValues)) / 2
            yWidth = max(yShiftValues) - min(yShiftValues) + 10
            strip.orderedAxes[1].position = yPosition
            strip.orderedAxes[1].width = yWidth

            try:
                axisCode = strip.axisCodes[1]
                strip._CcpnGLWidget.setAxisPosition(axisCode=axisCode, position=yPosition, update=False)
                strip._CcpnGLWidget.setAxisWidth(axisCode=axisCode, width=yWidth, update=False)
                strip._CcpnGLWidget._rescaleAllAxis()



            except Exception as es:
                getLogger().debugGL('OpenGL widget not instantiated')

    def _centreCcpnStripsForNmrResidue(self, nmrResidue, strips):
        """
        Centre y-axis of strip based on chemical shifts of from NmrResidue.nmrAtoms
        """
        if not nmrResidue:
            getLogger().warning('No NmrResidue specified')
            return

        if not strips:
            getLogger().warning('No Strip specified')
            return

        yShifts = matchAxesAndNmrAtoms(strips[0], nmrResidue.nmrAtoms)[strips[0].axisOrder[1]]
        yShiftValues = [x.value for x in yShifts]
        if yShiftValues:

            _minPpmWidths = {'H': 0.5, 'C': 8.0, 'N': 2.0}              # based on standard ratios

            yPosition = (max(yShiftValues) + min(yShiftValues)) / 2
            yWidth = max(yShiftValues) - min(yShiftValues)

            # original strips match axes
            strips[0].orderedAxes[1].position = yPosition
            strips[0].orderedAxes[1].width = yWidth

            try:
                axisCode = strips[0].axisCodes[1]

                for strip in strips:
                    # adjust the position of the strip to be clear of the new headers

                    minPpm = 1.0 if axisCode[0] not in _minPpmWidths else _minPpmWidths[axisCode[0]]

                    yPixel = max(yWidth, minPpm) / strip._CcpnGLWidget.height()
                    yPos = yPosition - (40 * yPixel)
                    yW = max(yWidth, minPpm) + (140 * yPixel)

                    strip._CcpnGLWidget.setAxisPosition(axisCode=axisCode, position=yPos, update=False)
                    strip._CcpnGLWidget.setAxisWidth(axisCode=axisCode, width=yW, update=False)
                    strip._CcpnGLWidget._scaleToYAxis()

                from ccpn.ui.gui.lib.OpenGL.CcpnOpenGL import GLNotifier

                GLSignals = GLNotifier(parent=self)
                GLSignals.emitPaintEvent()

            except Exception as es:
                getLogger().debugGL('OpenGL widget not instantiated')

    def _setupShiftDicts(self, *args):
        """
        Creates two ordered dictionaries for the inter residue and intra residue CA and CB shifts for
        all NmrResidues in the project.
        """
        self.intraShifts = OrderedDict()
        self.interShifts = OrderedDict()
        self.allShifts = OrderedDict()
        chemicalShiftList = self.application.project.getByPid(self.shiftListWidget.pulldownList.currentText())

        if chemicalShiftList:
            for nmrResidue in self.application.project.nmrResidues:
                nmrAtoms = [nmrAtom for nmrAtom in nmrResidue.nmrAtoms]
                shifts = [chemicalShiftList.getChemicalShift(atom.id) for atom in nmrAtoms]
                if nmrResidue.sequenceCode.endswith('-1'):
                    self.interShifts[nmrResidue] = shifts
                else:
                    self.intraShifts[nmrResidue] = shifts
                self.allShifts[nmrResidue] = shifts

    def _createMatchStrips(self, assignMatrix: typing.Tuple[typing.Dict[NmrResidue, typing.List[ChemicalShift]], typing.List[float]]):
        """
        Creates strips in match module corresponding to the best assignment possibilities
        in the assignMatrix.
        """
        if not assignMatrix:
            getLogger().warn('No assignment matrix specified')
            return

        # Assignment score has format {score: nmrResidue} where score is a float
        # assignMatrix[0] is a dict {score: nmrResidue} assignMatrix[1] is a concurrent list of scores
        # numberOfMatches = int(self.numberOfMatchesWidget.getText())
        assignmentScores = sorted(list(assignMatrix.keys()))[:MAXMATCHES]
        nmrAtomPairs = []
        scoreAssignment = []
        scoreLabelling = []

        matchDirection = 0
        for assignmentScore in assignmentScores:
            matchResidue = assignMatrix[assignmentScore]
            if matchResidue.sequenceCode.endswith('-1'):
                iNmrResidue = matchResidue.mainNmrResidue

                # this is where the nmrResidue can be dropped in the existing nmrChain
                scoreLabelling.append('[ i+1 ]')
                matchDirection = 1
            else:
                iNmrResidue = matchResidue
                scoreLabelling.append('[ i-1 ]')
                matchDirection = -1

            scoreAssignment.append('[ %i' % int(100 - min(1000 * assignmentScore, 100)) + '% ]')

            nmrAtomPairs.append((iNmrResidue.fetchNmrAtom(name='N'), iNmrResidue.fetchNmrAtom(name='H')))

        if matchDirection == 1:
            numberOfMatches = int(self.numberOfPlusMatchesWidget.getText())
        else:
            numberOfMatches = int(self.numberOfMinusMatchesWidget.getText())

        nmrAtomPairs = nmrAtomPairs[:numberOfMatches]
        scoreAssignment = scoreAssignment[:numberOfMatches]
        scoreLabelling = scoreLabelling[:numberOfMatches]

        for module in self._getMatchDisplays():

            # skip of the module if not defined - possibly in the case that spectrumDisplays have been closed
            if not module:
                continue

            makeStripPlot(module, nmrAtomPairs)

            for ii, strip in enumerate(module.strips):
                nmrResiduePid = nmrAtomPairs[ii][0].nmrResidue.pid

                # strip.setStripLabelText(nmrResiduePid)
                # strip.showStripLabel()
                # strip.setStripLabelisPlus(True if scoreLabelling[ii].startswith('i+1') else False)
                # strip.setStripResidueIdText(scoreLabelling[ii])
                # strip.showStripResidueId()
                # strip.setStripResidueDirText(scoreAssignment[ii])
                # strip.showStripResidueDir()

                strip.header.reset()
                strip.header.setLabelText(position='l', text=scoreLabelling[ii])
                strip.header.setLabelText(position='c', text=nmrResiduePid)

                # TODO:ED need to improve this
                # strip.header.setLabelConnectDir(position='c', connectDir=STRIPCONNECT_LEFT if scoreLabelling[ii].startswith('i-1') else STRIPCONNECT_RIGHT)
                strip.header.setLabelConnectDir(position='c', connectDir=STRIPCONNECT_LEFT if 'i-1' in scoreLabelling[ii] else STRIPCONNECT_RIGHT)
                strip.header.setLabelText(position='r', text=scoreAssignment[ii])

                # disable dropping onto these labels
                strip.header.setLabelObject(position='l', obj=None)
                strip.header.setLabelObject(position='c', obj=None)
                strip.header.setLabelObject(position='r', obj=None)

                strip.header.handle = STRIPBACKBONE
                strip.header.headerVisible = True

            # self._centreStripForNmrResidue(assignMatrix[assignmentScores[0]], module.strips[0])
            self._centreCcpnStripsForNmrResidue(assignMatrix[assignmentScores[0]], module.strips)
            module.setColumnStretches(stretchValue=True)

    def _closeModule(self):
        """
        Re-implementation of the closeModule method of the CcpnModule class required
        """
        # TODO: use proper subclassing

        for display in self._getDisplays() + self._getMatchDisplays():
            if display:
                display.hideAllStripHeaders(handle=STRIPBACKBONE)

        for notifier in self._stripNotifiers:
            if notifier:
                notifier.unRegister()

        self._stripNotifiers = []
        super()._closeModule()


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


def nmrAtomsFromOffsets(nmrResidue):
    """
    Retrieve a list of nmrAtoms from nmrResidue
    """
    # nmrResidue = nmrResidue.mainNmrResidue
    nmrResidues = []
    nmrResidues.append(nmrResidue)
    if nmrResidue.offsetNmrResidues:
        nmrResidues.extend(nmrResidue.offsetNmrResidues)

    nmrAtoms = []
    for nr in nmrResidues:
        nmrAtoms.extend(nr.nmrAtoms)

    return nmrAtoms


def markNmrAtoms(mainWindow, nmrAtoms: typing.List[NmrAtom]):
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

        if strips:
            strip = strips[0]

            # for strip in strips:
            # assume that this returns list of nmrAtoms in the display

            shiftDict = matchAxesAndNmrAtoms(strip, nmrAtoms)
            # atomPositions = shiftDict[strip.axisOrder[2]]
            # atomPositions = [[x.value for x in shiftDict[axisCode]] for axisCode in strip.axisOrder]
            # positions = []

            # for atomPos in atomPositions:
            #     if atomPos:
            #         if len(atomPos) < 2:
            #             positions.append(atomPos[0])
            #         else:
            #             positions.append(max(atomPos) - min(atomPos) / 2)
            #     else:
            #         positions.append('')
            # navigateToPositionInStrip(strip, positions, widths=widths) # don't need to change display yet

            mainWindow.markPositions(list(shiftDict.keys()),
                                     list(shiftDict.values()))


#=====  Just some code to 'save' =====
# def hasNmrResidue(nmrChain, residueCode):
#     "Simple function to check if sequenCode is found within the nmrResidues of nmrChain"
#     resCodes = [res.sequenceCode for res in nmrChain.nmrResidues]
#     return (residueCode in resCodes)
#
#
# def endOfchain(nmrResidue):
#     # changes to end of connected chain; not a good idea
#     if nmrResidue.nmrChain.isConnected:
#         if nmrResidue.sequenceCode.endswith('-1'):
#             nmrResidue = nmrResidue.nmrChain.mainNmrResidues[0].getOffsetNmrResidue(-1)
#         else:
#             nmrResidue = nmrResidue.nmrChain.mainNmrResidues[-1]
#     return nmrResidue
#
#
# def getPids(fromObject, attributeName):
#     "Get a list of pids fromObject.attributeName or None on error"
#     if not hasattr(fromObject, attributeName): return None
#     return [obj.pid for obj in getattr(fromObject, attributeName)]
#
#
#===== end code save =====


if __name__ == '__main__':
    from ccpn.ui.gui.widgets.Application import TestApplication


    app = TestApplication()

    popup = BackboneAssignmentModule()

    popup.show()
    popup.raise_()
    app.start()
