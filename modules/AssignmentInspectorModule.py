"""This file contains AssignmentInspectorModule class

modified by Geerten 1-9/12/2016:
- intialisation with 'empty' settings possible,
- now responsive to current.nmrResidues
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
__version__ = "$Revision: 3.0.b4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: geertenv $"
__date__ = "$Date: 2016-07-09 14:17:30 +0100 (Sat, 09 Jul 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

from PyQt5 import QtCore, QtWidgets
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
# from ccpn.ui.gui.widgets.Table import ObjectTable, Column
from ccpn.ui.gui.widgets.QuickTable import QuickTable
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
from ccpn.ui.gui.lib.Strip import navigateToPositionInStrip, navigateToNmrAtomsInStrip, _getCurrentZoomRatio
from ccpn.core.PeakList import PeakList
from ccpn.util.Common import makeIterableList


logger = getLogger()
ALL = '<all>'


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
    settingsPosition = 'left'

    def __init__(self, mainWindow, name='Assignment Inspector', chemicalShiftList=None):
        # CcpnModule.__init__(self, parent=mainWindow.moduleArea, name=name)
        # CcpnModule.__init__(self, mainWindow=mainWindow, name=name)  # ejb
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
                                                     moduleParent=self,
                                                     setLayout=True,
                                                     actionCallback=self._actionCallback,
                                                     selectionCallback=self._selectionCallback,
                                                     grid=(0, 0),
                                                     hiddenColumns=['Pid', 'Shift list peaks', 'All peaks'])

        self._selectCurrentNmrAtomsNotifier = Notifier(self.current, [Notifier.CURRENT], targetName=NmrAtom._pluralLinkName,
                                                       callback=self._highlightNmrAtoms)

        # settingsWidget
        if chemicalShiftList is not None:
            self.chemicalShiftTable.selectChemicalShiftList(chemicalShiftList)

        # install the event filter to handle maximising from floated dock
        self.installMaximiseEventHandler(self._maximise, self._closeModule)

        self._registerNotifiers()

    def _fillDisplayWidget(self):
        list = ['> select-to-add <'] + [ALL] + [display.pid for display in self.mainWindow.spectrumDisplays]
        self.displaysWidget.pulldownList.setData(texts=list)

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

    # def _registerNotifiers(self):
    #   # self.application.current.registerNotify(self._updateModuleCallback, 'nmrResidues')
    #   # self.project.registerNotifier('NmrAtom', 'change', self._refreshTable)   # just refresh the table
    #   # self.project.registerNotifier('Peak', 'change', self._refreshTable, onceOnly=True)
    #
    #   self._updateNotifier = Notifier(self.current
    #                                   , triggers=[Notifier.CURRENT]
    #                                   , targetName='nmrResidues'
    #                                   , callback=self._updateModuleCallback)
    #   self._nmrAtomNotifier = Notifier(self.project
    #                                    , triggers=[Notifier.CHANGE]
    #                                    , targetName='NmrAtom'
    #                                    , callback=self._refreshTable)
    #   self._peakNotifier = Notifier(self.project
    #                                 , triggers=[Notifier.CHANGE]
    #                                 , targetName='Peak'
    #                                 , callback=self._refreshTable)
    #
    # def _unRegisterNotifiers(self):
    #   # self.application.current.unRegisterNotify(self._updateModuleCallback, 'nmrResidues')
    #   # self.project.unregisterNotifier('NmrAtom', 'change', self.assignedPeaksTable.update)   # just refresh the table
    #   # self.project.unRegisterNotifier('Peak', 'change', self.assignedPeaksTable.update)
    #
    #   if self._updateNotifier:
    #     self._updateNotifier.unRegister()
    #   if self._nmrAtomNotifier:
    #     self._nmrAtomNotifier.unRegister()
    #   if self._peakNotifier:
    #     self._peakNotifier.unRegister()

    def _registerNotifiers(self):
        """Set up the notifiers
        """
        self._selectOnTableCurrentNmrResiduesNotifier = Notifier(self.current,
                                                                 [Notifier.CURRENT],
                                                                 targetName=NmrResidue._pluralLinkName,
                                                                 callback=self._highlightNmrResidues)

    def _unRegisterNotifiers(self):
        """Clean up the notifiers
        """
        if self._selectOnTableCurrentNmrResiduesNotifier is not None:
            self._selectOnTableCurrentNmrResiduesNotifier.unRegister()

    def _closeModule(self):
        """
        CCPN-INTERNAL: used to close the module
        """
        if self._selectCurrentNmrAtomsNotifier:
            self._selectCurrentNmrAtomsNotifier.unRegister()
        self.assignedPeaksTable._close()
        self.chemicalShiftTable._close()
        self._unRegisterNotifiers()
        super(AssignmentInspectorModule, self)._closeModule()

    def close(self):
        """
        Close the table from the commandline
        """
        self._closeModule()

    def _setCurrentPeak(self, data):
        """
        PeakTable select callback
        """
        from ccpn.core.lib.CallBack import CallBack

        peak = data[CallBack.OBJECT]
        # multiselection not allowed, sot only return the first object in list
        if peak:
            self.application.current.peaks = peak

    def _actionCallback(self, data):
        """
        Notifier DoubleClick action on item in table
        """
        obj = data[CallBack.OBJECT]

        getLogger().debug('AssignmentInspector_ChemicalShift>>> action', obj)

    def _selectionCallback(self, data):
        """
        Notifier Callback for selecting a row in the table
        """
        objList = data[CallBack.OBJECT]

        if objList:
            residues = [cs.nmrAtom.nmrResidue for cs in objList]

            if residues:
                self.current.nmrAtoms = [cs.nmrAtom for cs in objList]
                self.current.nmrResidues = [cs.nmrAtom.nmrResidue for cs in objList]

                self.assignedPeaksTable._updateModuleCallback({'value': residues})

        getLogger().debug('AssignmentInspector_ChemicalShift>>> selection', objList)

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

        peak = data[CallBack.OBJECT]
        if peak:
            self.current.peak = peak

            self.application._startCommandBlock('%s.navigateToPositionInStrip(project.getByPid(%r))' %
                                                (self.className, peak.position))
            try:
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

                            # navigateToPositionInStrip(strip=strip, positions=peak.position, widths=widths)
                            navigateToNmrAtomsInStrip(strip, makeIterableList(peak.assignedNmrAtoms),
                                                      widths=widths, markPositions=markPositions,
                                                      setNmrResidueLabel=False)
            finally:
                self.application._endCommandBlock()

    def _highlightNmrResidues(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._dataFrameObject:
            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects
            # peaks = self.assignedPeaksTable._dataFrameObject._objects

            residues = set(objList.nmrResidues)  #        set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom.nmrResidue in residues]
            # print ('>>>', highlightList)

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self.assignedPeaksTable._updateModuleCallback({'value': list(residues)})

    def _highlightNmrAtoms(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._dataFrameObject:
            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects
            # peaks = self.assignedPeaksTable._dataFrameObject._objects

            residues = set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom.nmrResidue in residues]
            # print ('>>>', highlightList)

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self.assignedPeaksTable._updateModuleCallback({'value': list(residues)})

        return

        if objList:
            residues = [cs.nmrAtom.nmrResidue for cs in objList]

            if residues:
                self.current.nmrAtoms = [cs.nmrAtom for cs in objList]
                self.current.nmrResidues = [cs.nmrAtom.nmrResidue for cs in objList]

                self.assignedPeaksTable._updateModuleCallback({'value': residues})

        getLogger().debug('AssignmentInspector>>> highlight nmrAtoms', objList)


class AssignmentInspectorTable(QuickTable):
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
        # Frame-1: NmrAtoms
        width = 130
        self.frame1 = Frame(None, setLayout=True)  # ejb
        self.frame1.setFixedWidth(width)
        self.nmrAtomLabel = Label(self.frame1, 'NmrAtom(s):', bold=True,
                                  grid=(0, 0), gridSpan=(1, 1), vAlign='center', margins=[2, 5, 2, 5])

        self.attachedNmrAtomsList = ListWidget(self.frame1,
                                               callback=self._updatePeakTableCallback, contextMenu=False,
                                               grid=(1, 0), gridSpan=(1, 1)
                                               )
        self.attachedNmrAtomsList.setFixedWidth(width - 2)

        self.frame1.hide()

        # Frame-2: peaks
        self.frame2 = Frame(parent, grid=(0, 0), gridSpan=(1, 1), setLayout=True)  # ejb
        self.peaksLabel = Label(self.frame2, 'Peaks assigned to NmrAtom(s):', bold=True,
                                grid=(0, 0), gridSpan=(1, 1), margins=[2, 5, 2, 5])

        self.frame2.setFixedHeight(24)
        self.peaksLabel.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        self.frame2.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Minimum)
        # initialise the currently attached dataFrame
        self._hiddenColumns = ['Pid']
        self.dataFrameObject = None

        # self.assignedPeaksTable = ObjectTable(self.frame2, self.getColumns(),
        #                                       selectionCallback=self._setCurrentPeak,
        #                                       actionCallback=self._navigateToPeak,
        #                                       objects=[], autoResize=True,
        #                                       grid=(1, 0), gridSpan=(1, 5), **policies
        #                                       )

        QuickTable.__init__(self, parent=parent,
                            mainWindow=self.mainWindow,
                            dataFrameObject=None,
                            setLayout=True,
                            autoResize=True, multiSelect=True,
                            selectionCallback=selectionCallback,
                            actionCallback=actionCallback,
                            grid=(3, 0), gridSpan=(1, 6),
                            enableDelete=False, enableSearch=False
                            )

        # self._assignmentFrame.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        #self.attachedNmrAtomsList.setFixedHeight(200)
        #self.assignedPeaksTable.setFixedHeight(200)

        # self._registerNotifiers()

        self._peakList = None
        # update if current.nmrResidue is defined
        if self.application.current.nmrResidue is not None and self.application.current.chemicalShiftList is not None:
            self._updateModuleCallback({'value': [self.application.current.nmrResidue]})

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
                               callBackClass=NmrResidue,
                               selectCurrentCallBack=None)  #self._updateModuleCallback)   #self._selectOnTableCurrentNmrResiduesNotifierCallback)

    def _updateModuleCallback(self, data: dict):
        """
        Callback function: Module responsive to nmrResidues; updates the list widget with nmrAtoms and updates peakTable if
        current.nmrAtom belongs to nmrResidue
        """
        if data and 'value' in data:
            nmrResidues = data['value']
            self.attachedNmrAtomsList.clear()

            if nmrResidues is not None and len(nmrResidues) > 0 \
                    and nmrResidues[-1] and len(nmrResidues[-1].nmrAtoms) > 0:

                # get the pids and append <all>
                self.ids = [atm.id for atm in nmrResidues[-1].nmrAtoms] + [ALL]
                self.attachedNmrAtomsList.addItems(self.ids)

                # # clear and fill the peak table
                # self.assignedPeaksTable.setObjects([])
                if self.application.current.nmrAtom is not None and self.application.current.nmrAtom.id in self.ids:
                    self._updatePeakTable(self.application.current.nmrAtom.id)
                else:
                    self._updatePeakTable(ALL)

                # new to populate table
            else:
                logger.debug('No valid nmrAtom/nmrResidue defined')

    def _updatePeakTableCallback(self, item):
        """
        Update the peakTable using item.text (which contains a NmrAtom pid or <all>)
        """
        if item is not None:
            text = item.text()
            self._updatePeakTable(text)
        else:
            logger.error('No valid item selected')


    class emptyObject():
        def __init__(self):
            pass


    def _updatePeakTable(self, id):
        """
        Update peak table depending on value of id;
        clears peakTable if pid is None
        """
        if id is None:
            # self.assignedPeaksTable.setObjects([])
            self.assignedPeaksTable.clearTable()
            return

        if id == ALL:
            self._peakList = self.emptyObject()

            self._peakList.peaks = list(set([pk for nmrAtom in self.application.current.nmrResidue.nmrAtoms for pk in nmrAtom.assignedPeaks]))

            self.project.blankNotification()
            objs = self.getSelectedObjects()
            self._dataFrameObject = self.getDataFrameFromList(table=self,
                                                              buildList=self._peakList.peaks,
                                                              colDefs=self.getColumns(),
                                                              hiddenColumns=self._hiddenColumns)

            # populate from the Pandas dataFrame inside the dataFrameObject
            self.setTableFromDataFrameObject(dataFrameObject=self._dataFrameObject)
            self._highLightObjs(objs)
            self.project.unblankNotification()

            # self.assignedPeaksTable.setObjects(peaks)
            # highlight current.nmrAtom in the list widget
            self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
            self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % ALL)
        else:
            pid = 'NA:' + id
            nmrAtom = self.application.project.getByPid(pid)
            #print('>>', pid, nmrAtom)
            if nmrAtom is not None:
                self.project.blankNotification()
                objs = self.getSelectedObjects()
                self._dataFrameObject = self.getDataFrameFromList(table=self,
                                                                  buildList=nmrAtom.assignedPeaks,
                                                                  colDefs=self.getColumns(),
                                                                  hiddenColumns=self._hiddenColumns)

                # populate from the Pandas dataFrame inside the dataFrameObject
                self.setTableFromDataFrameObject(dataFrameObject=self._dataFrameObject)
                self._highLightObjs(objs)
                self.project.unblankNotification()

                # self.assignedPeaksTable.setObjects(nmrAtom.assignedPeaks)
                # highlight current.nmrAtom in the list widget

                self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
                atomList = ', '.join([str(id) for id in self.ids if id != ALL])
                self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % atomList)  # nmrAtom.id)

    def getColumns(self):
        "get columns for initialisation of table"
        columns = ColumnClass([('Peak', lambda pk: pk.id, '', None)
                                  , ('Pid', lambda pk: pk.pid, 'Pid of peak', None)
                                  , ('_object', lambda pk: pk, 'Object', None)
                                  , ('serial', lambda pk: pk.serial, '', None)])
        tipTexts = []
        # get the maxmimum number of dimensions from all spectra in the project
        numDim = max([sp.dimensionCount for sp in self.application.project.spectra] + [1])

        for i in range(numDim):
            j = i + 1
            c = Column('Assign F%d' % j
                       , lambda pk, dim=i: getPeakAnnotation(pk, dim)
                       , 'NmrAtom assignments of peak in dimension %d' % j
                       , None)
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

            c = Column(text
                       , lambda pk, dim=i, unit=unit: getPeakPosition(pk, dim, unit)
                       , tipText
                       , None)
            columns._columns.append(c)

            # columns.append(c)
            # tipTexts.append(tipText)

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

    def _close(self):
        """
        Cleanup the notifiers when the window is closed
        """
        self.clearTableNotifiers()

    def _refreshTable(self, *args):
        self.update()
