"""This file contains AssignmentInspectorModule class

modified by Geerten 1-9/12/2016:
- intialisation with 'empty' settings possible,
- now responsive to current.nmrResidues
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
__author__ = "$Author: geertenv $"
__date__ = "$Date: 2016-07-09 14:17:30 +0100 (Sat, 09 Jul 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

from PyQt5 import QtCore, QtWidgets
from contextlib import contextmanager
from ccpn.util.OrderedSet import OrderedSet
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
# from ccpn.ui.gui.widgets.Table import ObjectTable, Column
from ccpn.ui.gui.widgets.GuiTable import GuiTable
from ccpn.ui.gui.widgets.Column import ColumnClass, Column
from ccpn.util.Logging import getLogger
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget
from ccpn.ui.gui.widgets.CompoundWidgets import ListCompoundWidget
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.core.lib.peakUtils import getPeakPosition, getPeakAnnotation
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.NmrAtom import NmrAtom, NmrResidue
from ccpn.core.Peak import Peak
from ccpn.ui.gui.modules.ChemicalShiftTable import ChemicalShiftTable
from ccpn.ui.gui.widgets.Splitter import Splitter
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.core.lib.CallBack import CallBack
from ccpn.ui.gui.lib.Strip import navigateToNmrAtomsInStrip, \
    _getCurrentZoomRatio, navigateToNmrResidueInDisplay
from ccpn.core.PeakList import PeakList
from ccpn.util.Common import makeIterableList


logger = getLogger()
ALL = '<all>'
NMRRESIDUES = 'nmrResidues'
NMRATOMS = 'nmrAtoms'


class _emptyObject():
    # small object to facilitate passing data to simulated event
    def __init__(self):
        pass


class AssignmentInspectorModule(CcpnModule):
    """
    This Module allows inspection of the NmrAtoms of a selected NmrResidue
    It responds to current.nmrResidues, taking the last added residue to this list
    The NmrAtom listWidget allows for selection of the nmrAtom; subsequently its assignedPeaks
    are displayed.
    """

    # overide in specific module implementations
    className = 'AssignmentInspectorModule'
    attributeName = 'peaks'

    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'top'

    def __init__(self, mainWindow, name='Assignment Inspector', chemicalShiftList=None):
        super().__init__(mainWindow=mainWindow, name=name)  # gwv

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        # self.sampledDims = {} #GWV: not sure what this is supposed to do
        # self.ids = []  # list of currently displayed NmrAtom ids + <all>
        #
        # policies = dict(vAlign='top')

        # settings window

        # self.splitter = Splitter(self.mainWidget, QtCore.Qt.Vertical)
        self.splitter = Splitter(self.mainWidget, horizontal=False)
        self._chemicalShiftFrame = Frame(self.splitter, setLayout=True)
        self._assignmentFrame = Frame(self.splitter, setLayout=True)
        # self._chemicalShiftFrame = Frame(self.mainWidget, setLayout=True, grid=(0,0))
        # self._assignmentFrame = Frame(self.mainWidget, setLayout=True, grid=(1,0))
        self.mainWidget.getLayout().addWidget(self.splitter)

        self.splitter.setStretchFactor(0, 3)
        self.splitter.setStretchFactor(1, 2)
        self.splitter.setChildrenCollapsible(False)
        self._assignmentFrame.setMinimumHeight(100)

        self._AIwidget = Widget(self.settingsWidget, setLayout=True,
                                grid=(0, 0), vAlign='top', hAlign='left')

        # cannot set a notifier for displays, as these are not (yet?) implemented and the Notifier routines
        # underpinning the addNotifier call do not allow for it either
        colwidth = 140
        self.displaysWidget = ListCompoundWidget(self._AIwidget,
                                                 grid=(0, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                                                 vPolicy='minimal',
                                                 #minimumWidths=(colwidth, 0, 0),
                                                 fixedWidths=(colwidth, 2 * colwidth, None),
                                                 orientation='left',
                                                 labelText='Display(s):',
                                                 tipText='SpectrumDisplay modules to respond to double-click',
                                                 texts=[ALL] + [display.pid for display in self.application.ui.mainWindow.spectrumDisplays],
                                                 defaults=[ALL]
                                                 )
        self.displaysWidget.setPreSelect(self._fillDisplayWidget)
        self.displaysWidget.setFixedHeights((None, None, 40))

        self.sequentialStripsWidget = CheckBoxCompoundWidget(
                self._AIwidget,
                grid=(1, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Show sequential strips:',
                checked=False
                )

        self.markPositionsWidget = CheckBoxCompoundWidget(
                self._AIwidget,
                grid=(2, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Mark positions:',
                checked=True
                )
        self.autoClearMarksWidget = CheckBoxCompoundWidget(
                self._AIwidget,
                grid=(3, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Auto clear marks:',
                checked=True
                )
        self.showNmrAtomListWidget = CheckBoxCompoundWidget(
                self._AIwidget,
                grid=(4, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Show nmrAtom list:',
                checked=True,
                callback=self._setNmrAtomListVisible,
                )

        # main window

        # AssignedPeaksTable need to be intialised before chemicalShiftTable, as the callback of the latter requires
        # the former to be present
        self.assignedPeaksTable = AssignmentInspectorTable(parent=self._assignmentFrame,
                                                           mainWindow=self.mainWindow,
                                                           moduleParent=self,
                                                           setLayout=True,
                                                           selectionCallback=self._setCurrentPeak,
                                                           actionCallback=self._navigateToPeak,
                                                           grid=(0, 0),
                                                           hiddenColumns=['Pid'])
        self.chemicalShiftTable = ChemicalShiftTable(parent=self._chemicalShiftFrame,
                                                     mainWindow=self.mainWindow,
                                                     moduleParent=self.assignedPeaksTable,  # just to give a unique id
                                                     setLayout=True,
                                                     actionCallback=self.navigateToNmrResidueCallBack,
                                                     selectionCallback=self._selectionCallback,
                                                     grid=(0, 0),
                                                     hiddenColumns=['Pid', 'Shift list peaks', 'All peaks'])

        # settingsWidget
        if chemicalShiftList is not None:
            self.chemicalShiftTable.selectChemicalShiftList(chemicalShiftList)

        # install the event filter to handle maximising from floated dock
        self.installMaximiseEventHandler(self._maximise, self._closeModule)
        self._setNmrAtomListVisible(False)

        self._registerNotifiers()

    def _fillDisplayWidget(self):
        ll = ['> select-to-add <'] + [ALL] + [display.pid for display in self.mainWindow.spectrumDisplays]
        self.displaysWidget.pulldownList.setData(texts=ll)

    def _getDisplays(self):
        """
        Return list of displays to navigate - if needed
        """
        displays = []
        # check for valid displays
        gids = self.displaysWidget.getTexts()
        if len(gids) == 0: return displays
        if ALL in gids:
            displays = self.mainWindow.spectrumDisplays
        else:
            displays = [self.application.getByGid(gid) for gid in gids if gid != ALL]
        return displays

    def _maximise(self):
        """
        refresh the table on a maximise event
        """
        self._refreshTable()

    def _registerNotifiers(self):
        """Set up the notifiers
        """
        self.setNotifier(self.current, [Notifier.CURRENT], targetName=NmrResidue._pluralLinkName,
                         callback=self._highlightNmrResidues)
        self.setNotifier(self.current, [Notifier.CURRENT], targetName=NmrAtom._pluralLinkName,
                         callback=self._highlightNmrAtoms)

    def _closeModule(self):
        """
        CCPN-INTERNAL: used to close the module
        """
        self.assignedPeaksTable._close()
        self.chemicalShiftTable._close()
        super()._closeModule()

    def close(self):
        """
        Close the table from the commandline
        """
        self._closeModule()

    def _setNmrAtomListVisible(self, visible=None):
        """change the visibility of the nmrAtom list in the peak tables widget
        """
        if visible is None:
            visible = self.showNmrAtomListWidget.get()
        else:
            if not isinstance(visible, bool):
                raise TypeError('nmrList visibility must be True/False')

            self.showNmrAtomListWidget.set(visible)

        self.assignedPeaksTable.nmrAtomListFrame.setVisible(visible)

    def _setCurrentPeak(self, data):
        """
        PeakTable select callback
        """
        peaks = data[CallBack.OBJECT]
        # multiselection allowed, set current to all selected peaks
        if peaks:
            self.application.current.peaks = peaks

    @contextmanager
    def _notifierBlanking(self):
        """block nmrAtom and nmrResidue notifiers for this module only
        """
        self.setBlankingAllNotifiers(True)
        try:
            # transfer control to the calling function
            yield

        except Exception as es:
            getLogger().warning('Error in AssignmentInspectorModule', str(es))

        finally:
            # reset notifiers
            self.setBlankingAllNotifiers(False)

    def navigateToNmrResidueCallBack(self, data):
        """Navigate in selected displays to nmrResidue; skip if none defined
        """
        # handle a single chemicalShift - SHOULD always contain an object
        objs = data[CallBack.OBJECT]
        if not objs:
            return
        if isinstance(objs, (tuple, list)):
            chemicalShift = objs[0]
        else:
            chemicalShift = objs

        nmrResidue = chemicalShift.nmrAtom.nmrResidue

        getLogger().debug('nmrResidue=%s' % (nmrResidue.id))

        displays = self._getDisplays()
        if len(displays) == 0:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        from ccpn.core.lib.ContextManagers import undoBlock

        with undoBlock():
            # optionally clear the marks
            if self.autoClearMarksWidget.checkBox.isChecked():
                self.mainWindow.clearMarks()

            # navigate the displays
            for display in displays:
                if len(display.strips) > 0:
                    newWidths = []  #_getCurrentZoomRatio(display.strips[0].viewBox.viewRange())
                    navigateToNmrResidueInDisplay(nmrResidue, display, stripIndex=0,
                                                  widths=newWidths,  #['full'] * len(display.strips[0].axisCodes),
                                                  showSequentialResidues=(len(display.axisCodes) > 2) and
                                                                         self.sequentialStripsWidget.checkBox.isChecked(),
                                                  markPositions=self.markPositionsWidget.checkBox.isChecked()
                                                  )

    def _actionCallback(self, data):
        # def _selectionCallback(self, data):
        """
        Notifier single-click action on item in table
        Highlight nmrAtoms from chemicalShifts
        """
        # multiselection table will return a list of objects
        objs = data[CallBack.OBJECT]
        if not objs:
            return
        if isinstance(objs, (tuple, list)):
            objList = objs
        else:
            objList = (objs,)

        if objList:
            getLogger().debug('AssignmentInspector_ChemicalShift>>> action', objList)
            nmrResidues = [cs.nmrAtom.nmrResidue for cs in objList]

            if nmrResidues:
                # navigate to nmrResidues in displays

                pass

                # with self._notifierBlanking():
                #     # SHOULD be only 1, but multi-selection may give more
                #     self.current.nmrAtoms = [cs.nmrAtom for cs in objList]
                #     self.current.nmrResidues = nmrResidues
                #
                #     # don't need to update the selection on the chemicalShiftTable
                #
                #     self.assignedPeaksTable._updateModuleCallback({NMRRESIDUES: nmrResidues,
                #                                                    NMRATOMS   : self.current.nmrAtoms},
                #                                                   updateFromNmrResidues=False)

    def _selectionCallback(self, data):
        # def _actionCallback(self, data):
        """
        Notifier Callback for double-clicking a row in the table
        Highlight all nmrAtoms belonging to the same nmrResidue as nmrAtom in checmicalShifts
        """

        objList = data[CallBack.OBJECT]

        if objList:
            getLogger().debug('AssignmentInspector_ChemicalShift>>> action', objList)

            nmrResidues = [cs.nmrAtom.nmrResidue for cs in objList]

            if nmrResidues:
                with self._notifierBlanking():

                    nmrAtoms = OrderedSet()
                    for nmrRes in nmrResidues:
                        for nmrAtom in nmrRes.nmrAtoms:
                            nmrAtoms.add(nmrAtom)

                    self.current.nmrAtoms = tuple(nmrAtoms)
                    self.current.nmrResidues = nmrResidues

                    self._highlightChemicalShifts(nmrResidues)

                    self.assignedPeaksTable._updateModuleCallback({NMRRESIDUES: nmrResidues,
                                                                   NMRATOMS   : tuple(nmrAtoms)},
                                                                  updateFromNmrResidues=True)

    def _highlightChemicalShifts(self, nmrResidues):
        """
        Highlight chemical shifts in the table
        """
        if self.chemicalShiftTable._dataFrameObject:
            getLogger().debug('_highlightChemicalShifts ', nmrResidues)

            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects

            residues = set(nmrResidues)
            # highlightList = [cs for cs in chemicalShifts if cs.nmrAtom.nmrResidue in residues]
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in residues]
            self.chemicalShiftTable._highLightObjs(highlightList)

    def _navigateToPeak(self, data):
        """
        PeakTable double-click callback; navigate in to peak in current.strip
        """
        displays = self._getDisplays()
        markPositions = self.markPositionsWidget.checkBox.isChecked()

        if len(displays) == 0:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        # get the first object from the callback
        objs = data[CallBack.OBJECT]
        if not objs:
            return
        if isinstance(objs, (tuple, list)):
            peak = objs[0]
        else:
            peak = objs

        if peak:
            self.current.peak = peak

            from ccpn.core.lib.ContextManagers import undoBlock

            with undoBlock():
                # optionally clear the marks
                if self.autoClearMarksWidget.checkBox.isChecked():
                    self.application.ui.mainWindow.clearMarks()

                # navigate the displays
                for display in displays:
                    for strip in display.strips:

                        validPeakListViews = [pp.peakList for pp in strip.peakListViews if isinstance(pp.peakList, PeakList)]

                        if peak.peakList in validPeakListViews:
                            widths = None
                            if peak.peakList.spectrum.dimensionCount <= 2:
                                widths = _getCurrentZoomRatio(strip.viewRange())

                            navigateToNmrAtomsInStrip(strip, makeIterableList(peak.assignedNmrAtoms),
                                                      widths=widths, markPositions=markPositions,
                                                      setNmrResidueLabel=False)

    def _highlightNmrResidues(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._dataFrameObject:
            getLogger().debug('_highlightNmrResidues ', objList)

            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects
            # peaks = self.assignedPeaksTable._dataFrameObject._objects

            nmrResidues = set(objList.nmrResidues)  #        set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in nmrResidues]

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self.assignedPeaksTable._updateModuleCallback({NMRRESIDUES: nmrResidues},
                                                          updateFromNmrResidues=True)

    def _highlightNmrAtoms(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._dataFrameObject:
            getLogger().debug('_highlightNmrAtoms ', objList)

            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects
            # peaks = self.assignedPeaksTable._dataFrameObject._objects

            nmrResidues = set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in nmrResidues]
            # print ('>>>', highlightList)

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self.assignedPeaksTable._updateModuleCallback({NMRRESIDUES: nmrResidues,
                                                           NMRATOMS   : self.current.nmrAtoms},
                                                          updateFromNmrResidues=False)

        return

        # if objList:
        #     residues = [cs.nmrAtom.nmrResidue for cs in objList]
        #
        #     if residues:
        #         self.current.nmrAtoms = [cs.nmrAtom for cs in objList]
        #         self.current.nmrResidues = [cs.nmrAtom.nmrResidue for cs in objList]
        #
        #         self.assignedPeaksTable._updateModuleCallback({'value': residues})
        #
        # getLogger().debug('AssignmentInspector>>> highlight nmrAtoms', objList)


class AssignmentInspectorTable(GuiTable):
    """
    Class to present a NmrResidue Table and a NmrChain pulldown list, wrapped in a Widget
    """
    className = 'AssignmentInspectorTable'
    attributeName = 'chemicalShifts'

    OBJECT = 'object'
    TABLE = 'table'

    def __init__(self, parent=None, mainWindow=None, moduleParent=None, actionCallback=None, selectionCallback=None,
                 checkBoxCallback=None, nmrChain=None, multiSelect=False,
                 **kwds):
        """
        Initialise the widgets for the module.
        """
        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        if mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current
        else:
            self.application = None
            self.project = None
            self.current = None

        self.sampledDims = {}  #GWV: not sure what this is supposed to do
        self.ids = []  # list of currently displayed NmrAtom ids + <all>

        # main window

        # Frame-1: NmrAtoms list, hidden when table first opens
        self._nmrAtomListFrameWidth = 130
        self.nmrAtomListFrame = Frame(parent, grid=(0, 0), gridSpan=(1, 1), setLayout=True)  # ejb

        self.nmrAtomListFrame.setFixedWidth(self._nmrAtomListFrameWidth)
        self.nmrAtomLabel = Label(self.nmrAtomListFrame, 'NmrAtom(s):', bold=True,
                                  grid=(0, 0), gridSpan=(1, 1), vAlign='center', margins=[2, 5, 2, 5])

        self.attachedNmrAtomsList = ListWidget(self.nmrAtomListFrame,
                                               contextMenu=False,
                                               grid=(1, 0), gridSpan=(1, 1)
                                               )
        self.attachedNmrAtomsList.itemSelectionChanged.connect(self._updatePeakTableCallback)
        self.attachedNmrAtomsList.setDragEnabled(False)

        self.attachedNmrAtomsList.setFixedWidth(self._nmrAtomListFrameWidth - 2)
        self.attachedNmrAtomsList.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.MinimumExpanding)
        self.nmrAtomListFrame.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)

        # Frame-2: peaks
        self.frame2 = Frame(parent, grid=(0, 1), gridSpan=(1, 1), setLayout=True)  # ejb
        self.peaksLabel = Label(self.frame2, 'Peaks assigned to NmrAtom(s):', bold=True,
                                grid=(0, 0), gridSpan=(1, 4), margins=[2, 5, 2, 5])

        self.peaksLabel.setFixedHeight(24)
        self.peaksLabel.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        self.frame2.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)

        # initialise the currently attached dataFrame
        self._hiddenColumns = ['Pid']
        self.dataFrameObject = None

        # initialise the table and put in the second column (first hidden at start)
        super().__init__(parent=self.frame2,
                         mainWindow=self.mainWindow,
                         dataFrameObject=None,
                         setLayout=True,
                         autoResize=True, multiSelect=True,
                         selectionCallback=selectionCallback,
                         actionCallback=actionCallback,
                         grid=(1, 0), gridSpan=(1, 4),
                         enableDelete=False, enableSearch=False
                         )
        self.moduleParent = moduleParent

        # self._registerNotifiers() - removed for testing

        self._peakList = None
        # update if current.nmrResidue is defined
        if self.application.current.nmrResidue is not None and self.application.current.chemicalShiftList is not None:
            self._updateModuleCallback({NMRRESIDUES: [self.application.current.nmrResidue]},
                                       updateFromNmrResidues=True)

        # set the required table notifiers
        self.setTableNotifiers(tableClass=None,
                               rowClass=Peak,
                               cellClassNames=None,
                               tableName='peakList', rowName='peak',
                               changeFunc=self._refreshTable,
                               className=self.attributeName,
                               updateFunc=self._refreshTable,
                               tableSelection='_peakList',
                               pullDownWidget=None,  #self.ncWidget
                               callBackClass=Peak,
                               selectCurrentCallBack=self._selectOnTableCurrentPeaksNotifierCallback,
                               moduleParent=moduleParent)

    def _updateModuleCallback(self, data: dict, updateFromNmrResidues=True):
        """
        Callback function: Module responsive to nmrResidues; updates the list widget with nmrAtoms and updates peakTable if
        current.nmrAtom belongs to nmrResidue
        """
        if data:

            # update list from nmrResidues, show all nmrAtoms connected
            nmrResidues = data[NMRRESIDUES] if NMRRESIDUES in data else []

            if updateFromNmrResidues:
                # generate nmrAtom list from nmrResidues
                nmrAtoms = OrderedSet()
                for nmrRes in nmrResidues:
                    for nmrAtom in nmrRes.nmrAtoms:
                        nmrAtoms.add(nmrAtom)

            else:
                # get from the data dict
                nmrAtoms = data[NMRATOMS] if NMRATOMS in data else []

            # there is currently a hidden list widget containing the nmrAtom ids
            self.attachedNmrAtomsList.clear()
            self.ids = [atm.id for atm in nmrAtoms]
            self.attachedNmrAtomsList.addItems(self.ids)

            # populate peak table with the correct peaks
            self._updatePeakTable(nmrAtoms, messageAll=updateFromNmrResidues)

            # self.attachedNmrAtomsList.clear()
            #
            # if nmrResidues is not None and len(nmrResidues) > 0 \
            #         and nmrResidues[-1] and len(nmrResidues[-1].nmrAtoms) > 0:
            #
            #     # get the pids and append <all>
            #     self.ids = [atm.id for atm in nmrResidues[-1].nmrAtoms] + [ALL]
            #     self.attachedNmrAtomsList.addItems(self.ids)
            #
            #     # # clear and fill the peak table
            #     # self.assignedPeaksTable.setObjects([])
            #     if self.application.current.nmrAtom is not None and self.application.current.nmrAtom.id in self.ids:
            #         logger.debug('UPDATING selection')
            #
            #         self._updatePeakTable(self.application.current.nmrAtom.id)
            #     else:
            #         logger.debug('UPDATING All')
            #
            #         self._updatePeakTable(ALL)
            #
            #     # new to populate table
            # else:
            #     logger.debug('No valid nmrAtom/nmrResidue defined')

    def _selectOnTableCurrentPeaksNotifierCallback(self, data):
        """
        Callback from a notifier to highlight the peaks on the peak table
        :param data:
        """

        currentPeaks = data['value']
        self._selectOnTableCurrentPeaks(currentPeaks)

    def _selectOnTableCurrentPeaks(self, currentPeaks):
        """
        Highlight the list of peaks on the table
        :param currentPeaks:
        """

        self.highlightObjects(currentPeaks)
        # if len(currentPeaks) > 0:
        #     self._highLightObjs(currentPeaks)
        # else:
        #     self.clearSelection()

    def _updatePeakTableCallback(self, data=None):
        """
        Update the peakTable using item.text (which contains a NmrAtom pid or <all>)
        """
        # something has been selected, so get all selected items
        numTexts = self.attachedNmrAtomsList.count()
        selectedTexts = self.attachedNmrAtomsList.getSelectedTexts()
        nmrAtoms = [self.project.getByPid('NA:' + id) for id in selectedTexts]

        # populate the table with valid nmrAtoms
        self._updatePeakTable([atm for atm in nmrAtoms if atm is not None],
                              messageAll=True if numTexts == len(nmrAtoms) else False)

    # @contextmanager
    # def _projectBlanking(self):
    #
    #     self.project.blankNotification()
    #     objs = self.getSelectedObjects()
    #
    #     try:
    #         # transfer control to the calling function
    #         yield
    #
    #     except Exception as es:
    #         getLogger().warning('Error in AssignmentInspectorModule', str(es))
    #
    #     finally:
    #         # populate from the Pandas dataFrame inside the dataFrameObject
    #         self.setTableFromDataFrameObject(dataFrameObject=self._dataFrameObject)
    #         self._highLightObjs(objs)
    #         self.project.unblankNotification()

    def _updatePeakTable(self, nmrAtoms, messageAll=True):
        """
        Update peak table depending on value of id;
        clears peakTable if pid is None
        """
        if not nmrAtoms:
            # get all items from the table
            nmrAtoms = [self.project.getByPid('NA:' + id) for id in self.attachedNmrAtomsList.getTexts()]

            # # populate the table with valid nmrAtoms
            # self._updatePeakTable([atm for atm in nmrAtoms if atm is not None], messageAll=True)

        self._peakList = _emptyObject()
        self._peakList.peaks = list(set([pk for nmrAtom in nmrAtoms for pk in nmrAtom.assignedPeaks]))

        # with self._projectBlanking():
        #     self._dataFrameObject = self.getDataFrameFromList(table=self,
        #                                                       buildList=self._peakList.peaks,
        #                                                       colDefs=self.getColumns(),
        #                                                       hiddenColumns=self._hiddenColumns)

        self.populateTable(rowObjects=self._peakList.peaks,
                           columnDefs=self.getColumns())

        ids = [atm.id for atm in nmrAtoms]

        if messageAll:
            self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % ALL)
        else:
            atomList = ', '.join([str(id) for id in ids])
            self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % atomList)  # nmrAtom.id)

        # if messageAll:
        #
        #     self._peakList = _emptyObject()
        #     self._peakList.peaks = list(set([pk for nmrAtom in nmrAtoms for pk in nmrAtom.assignedPeaks]))
        #
        #     with self._projectBlanking():
        #         self._dataFrameObject = self.getDataFrameFromList(table=self,
        #                                                           buildList=self._peakList.peaks,
        #                                                           colDefs=self.getColumns(),
        #                                                           hiddenColumns=self._hiddenColumns)
        #
        #     # self.assignedPeaksTable.setObjects(peaks)
        #     # highlight current.nmrAtom in the list widget
        #     self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
        #     self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % ALL)
        # else:
        #     pid = 'NA:' + id
        #     nmrAtom = self.application.project.getByPid(pid)
        #
        #     if nmrAtom is not None:
        #         with self._projectBlanking():
        #             self._dataFrameObject = self.getDataFrameFromList(table=self,
        #                                                               buildList=nmrAtom.assignedPeaks,
        #                                                               colDefs=self.getColumns(),
        #                                                               hiddenColumns=self._hiddenColumns)
        #
        #         # self.assignedPeaksTable.setObjects(nmrAtom.assignedPeaks)
        #         # highlight current.nmrAtom in the list widget
        #
        #         self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
        #         atomList = ', '.join([str(id) for id in self.ids if id != ALL])
        #         self.peaksLabel.setText('Peaks assigned to NmrAtom: %s' % id)  # nmrAtom.id)

    def getColumns(self):
        "get columns for initialisation of table"
        columns = ColumnClass([('Peak', lambda pk: pk.id, '', None),
                               ('Pid', lambda pk: pk.pid, 'Pid of peak', None),
                               ('_object', lambda pk: pk, 'Object', None),
                               ('serial', lambda pk: pk.serial, '', None)])
        tipTexts = []
        # get the maxmimum number of dimensions from all spectra in the project
        numDim = max([sp.dimensionCount for sp in self.application.project.spectra] + [1])

        for i in range(numDim):
            j = i + 1
            c = Column('Assign F%d' % j,
                       lambda pk, dim=i: getPeakAnnotation(pk, dim),
                       'NmrAtom assignments of peak in dimension %d' % j,
                       None)
            columns._columns.append(c)

            # columns.append(c)
            # tipTexts.append('NmrAtom assignments of peak in dimension %d' % j)

            sampledDim = self.sampledDims.get(i)
            if sampledDim:
                text = 'Sampled\n%s' % sampledDim.conditionVaried
                tipText = 'Value of sampled plane'
                unit = sampledDim

            else:
                text = 'Pos F%d' % j
                tipText = 'Peak position in dimension %d' % j
                unit = 'ppm'

            c = Column(text,
                       lambda pk, dim=i, unit=unit: getPeakPosition(pk, dim, unit),
                       tipText,
                       None)
            columns._columns.append(c)

            # columns.append(c)
            # tipTexts.append(tipText)
        columns._columns.append(Column('height', lambda pk: pk.height, '', None))
        columns._columns.append(Column('volume', lambda pk: pk.volume, '', None))
        return columns

    # def _setCurrentPeak(self, data):
    #   """
    #   PeakTable select callback
    #   """
    #   from ccpn.core.lib.CallBack import CallBack
    #
    #   peak = data[CallBack.OBJECT]
    #   # multiselection not allowed, sot only return the first object in list
    #   if peak:
    #     self.application.current.peaks = peak

    # def _navigateToPeak(self, data):
    #   """
    #   PeakTable double-click callback; navigate in to peak in current.strip
    #   """
    #   displays = self._getDisplays()
    #   if len(displays) == 0:
    #     logger.warning('Undefined display module(s); select in settings first')
    #     showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
    #     return
    #
    #   peak = data[CallBack.OBJECT]
    #   if peak:
    #     self.current.peak = peak
    #
    #     self.application._startCommandBlock('%s.navigateToPositionInStrip(project.getByPid(%r))' %
    #         (self.className, peak.position))
    #     try:
    #       # optionally clear the marks
    #       if self.autoClearMarksWidget.checkBox.isChecked():
    #           self.application.ui.mainWindow.clearMarks()
    #
    #       # navigate the displays
    #       for display in displays:
    #         for strip in display.strips:
    #
    #           validPeakListViews = [pp.peakList for pp in strip.peakListViews if isinstance(pp.peakList, PeakList)]
    #
    #           if peak.peakList in validPeakListViews:
    #             widths = None
    #             if peak.peakList.spectrum.dimensionCount <= 2:
    #               widths = _getCurrentZoomRatio(strip.viewRange())
    #
    #             navigateToPositionInStrip(strip=strip, positions=peak.position, widths=widths)
    #
    #     finally:
    #         self.application._endCommandBlock()

    # peak = data[CallBack.OBJECT]
    #
    # from ccpn.ui.gui.lib.Strip import navigateToPositionInStrip
    # #print('>peakTableDoubleClick>', peak)
    # if peak is not None and self.application.current.strip is not None:
    #   self.application.current.peak = peak
    #   navigateToPositionInStrip(strip=self.application.current.strip, positions=peak.position)

    def _getPeakHeight(self, peak):
        """
        Returns the height of the specified peak as formatted string or 'None' if undefined
        """
        if peak.height:
            return '%7.2E' % float(peak.height)
        else:
            return '%s' % None

    # def _getSearchWidget(self):
    #   """
    #   CCPN-INTERNAL: used to get searchWidget
    #   """
    #   return self.searchWidget

    # def _close(self):
    #     """
    #     Cleanup the notifiers when the window is closed
    #     """
    #     self.clearTableNotifiers()

    def _refreshTable(self, *args):
        self.update()
