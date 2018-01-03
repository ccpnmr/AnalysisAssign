"""This file contains AssignmentInspectorModule class

derived from ModifyAssignmentModule by Simon;
extensively modified by Geerten 1-9/12/2016:
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
__version__ = "$Revision: 3.0.b2 $"
#=========================================================================================
# Created
#=========================================================================================

__author__ = "$Author: geertenv $"
__date__ = "$Date: 2016-07-09 14:17:30 +0100 (Sat, 09 Jul 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

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
  includeSettingsWidget = True
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  Position = 'top'

  def __init__(self, mainWindow, name='Assignment Inspector'):
    # CcpnModule.__init__(self, parent=mainWindow.moduleArea, name=name)
    CcpnModule.__init__(self, mainWindow=mainWindow, name=name)    # ejb

    # Derive application, project, and current from mainWindow
    self.mainWindow = mainWindow
    self.application = mainWindow.application
    self.project = mainWindow.application.project
    self.current = mainWindow.application.current

    self.sampledDims = {} #GWV: not sure what this is supposed to do
    self.ids = []  # list of currently displayed NmrAtom ids + <all>

    policies = dict(vAlign='top')

    # settings window

    self._AIwidget = Widget(self.settingsWidget, setLayout=True,
                             grid=(0,0), vAlign='top', hAlign='left')

    # cannot set a notifier for displays, as these are not (yet?) implemented and the Notifier routines
    # underpinning the addNotifier call do not allow for it either
    colwidth = 140
    self.displaysWidget = ListCompoundWidget(self._AIwidget,
                                             grid=(0,0), vAlign='top', stretch=(0,0), hAlign='left',
                                             vPolicy='minimal',
                                             #minimumWidths=(colwidth, 0, 0),
                                             fixedWidths=(colwidth, 2*colwidth, None),
                                             orientation = 'left',
                                             labelText='Display(s):',
                                             tipText = 'SpectrumDisplay modules to respond to double-click',
                                             texts=[ALL] + [display.pid for display in self.application.ui.mainWindow.spectrumDisplays]
                                             )
    self.displaysWidget.setFixedHeigths((None, None, 40))

    self.sequentialStripsWidget = CheckBoxCompoundWidget(
                                             self._AIwidget,
                                             grid=(1,0), vAlign='top', stretch=(0,0), hAlign='left',
                                             #minimumWidths=(colwidth, 0),
                                             fixedWidths=(colwidth, 30),
                                             orientation = 'left',
                                             labelText = 'Show sequential strips:',
                                             checked = False
                                            )

    self.markPositionsWidget = CheckBoxCompoundWidget(
                                             self._AIwidget,
                                             grid=(2,0), vAlign='top', stretch=(0,0), hAlign='left',
                                             #minimumWidths=(colwidth, 0),
                                             fixedWidths=(colwidth, 30),
                                             orientation = 'left',
                                             labelText = 'Mark positions:',
                                             checked = True
                                            )
    self.autoClearMarksWidget = CheckBoxCompoundWidget(
                                             self._AIwidget,
                                             grid=(3,0), vAlign='top', stretch=(0,0), hAlign='left',
                                             #minimumWidths=(colwidth, 0),
                                             fixedWidths=(colwidth, 30),
                                             orientation = 'left',
                                             labelText = 'Auto clear marks:',
                                             checked = True
                                            )

    # main window
    # Frame-1: NmrAtoms
    width = 130
    self.frame1 = Frame(self.mainWidget, grid=(0,0), **policies, fShape='styledPanel', fShadow='plain', setLayout=True) # ejb
    self.frame1.setFixedWidth(width)
    self.nmrAtomLabel = Label(self.frame1, 'NmrAtom(s):', bold=True,
                              grid=(0, 0), gridSpan=(1, 1), vAlign='center', margins=[2,5,2,5])

    self.attachedNmrAtomsList = ListWidget(self.frame1,
                                           callback=self._updatePeakTableCallback, contextMenu=False,
                                           grid=(1, 0), gridSpan=(1, 1), **policies
                                           )
    self.attachedNmrAtomsList.setFixedWidth(width-2)


    # Frame-2: peaks
    self.frame2 = Frame(self.mainWidget, grid=(0,1), gridSpan=(1,5), **policies, fShape='styledPanel', fShadow='plain', setLayout=True) # ejb
    self.peaksLabel = Label(self.frame2, 'Peaks assigned to NmrAtom(s):', bold=True,
                            grid=(0, 0), gridSpan=(1, 1), vAlign='center', margins=[2,5,2,5])

    # initialise the currently attached dataFrame
    self._hiddenColumns = ['Pid']
    self.dataFrameObject = None

    # self.assignedPeaksTable = ObjectTable(self.frame2, self.getColumns(),
    #                                       selectionCallback=self._setCurrentPeak,
    #                                       actionCallback=self._navigateToPeak,
    #                                       objects=[], autoResize=True,
    #                                       grid=(1, 0), gridSpan=(1, 5), **policies
    #                                       )

    self.assignedPeaksTable = QuickTable(parent=self.frame2
                                         , mainWindow=self.mainWindow
                                         , dataFrameObject=None
                                         , setLayout=True
                                         , autoResize=True, multiSelect=False
                                         , selectionCallback=self._setCurrentPeak
                                         , actionCallback=self._navigateToPeak
                                         , grid=(1, 0), gridSpan=(1, 5), **policies
                                         , enableDelete=False)

    #self.attachedNmrAtomsList.setFixedHeight(200)
    #self.assignedPeaksTable.setFixedHeight(200)

    self._registerNotifiers()

    # update if current.nmrResidue is defined
    if self.application.current.nmrResidue is not None:
      self._updateModuleCallback([self.application.current.nmrResidue])

    # set the required table notifiers
    # self.setTableNotifiers(tableClass=NmrChain
    #                        , rowClass=NmrResidue
    #                        , cellClassNames=(NmrAtom, 'nmrAtom')
    #                        , tableName='nmrChain', rowName='nmrResidue'
    #                        , changeFunc=self.displayTableForNmrChain
    #                        , className=self.attributeName
    #                        , updateFunc=self._update
    #                        , tableSelection='nmrChain'
    #                        , pullDownWidget=self.ncWidget
    #                        , selectCurrentCallBack=self._selectOnTableCurrentNmrResiduesNotifierCallback)

    # install the event filter to handle maximising from floated dock
    self.installMaximiseEventHandler(self._maximise)

  def _maximise(self):
    """
    refresh the table on a maximise event
    """
    self._refreshTable()

  def _registerNotifiers(self):
    # self.application.current.registerNotify(self._updateModuleCallback, 'nmrResidues')
    # self.project.registerNotifier('NmrAtom', 'change', self._refreshTable)   # just refresh the table
    # self.project.registerNotifier('Peak', 'change', self._refreshTable, onceOnly=True)

    self._updateNotifier = Notifier(self.current
                                    , triggers=[Notifier.CURRENT]
                                    , targetName='nmrResidues'
                                    , callback=self._updateModuleCallback)
    self._nmrAtomNotifier = Notifier(self.project
                                     , triggers=[Notifier.CHANGE]
                                     , targetName='NmrAtom'
                                     , callback=self._refreshTable)
    self._peakNotifier = Notifier(self.project
                                  , triggers=[Notifier.CHANGE]
                                  , targetName='Peak'
                                  , callback=self._refreshTable)

  def _unregisterNotifiers(self):
    # self.application.current.unRegisterNotify(self._updateModuleCallback, 'nmrResidues')
    # self.project.unregisterNotifier('NmrAtom', 'change', self.assignedPeaksTable.update)   # just refresh the table
    # self.project.unRegisterNotifier('Peak', 'change', self.assignedPeaksTable.update)

    if self._updateNotifier:
      self._updateNotifier.unRegister()
    if self._nmrAtomNotifier:
      self._nmrAtomNotifier.unRegister()
    if self._peakNotifier:
      self._peakNotifier.unRegister()

  def _refreshTable(self, *args):
    self.assignedPeaksTable.update()

  def _closeModule(self):
    """
    CCPN-INTERNAL: used to close the module
    """
    # self._unregisterNotifiers()
    self.assignedPeaksTable.clearTableNotifiers()
    super(AssignmentInspectorModule, self)._closeModule()

  def close(self):
    """
    Close the table from the commandline
    """
    self._closeModule()

  def _updateModuleCallback(self, data:dict):
    """
    Callback function: Module responsive to nmrResidues; updates the list widget with nmrAtoms and updates peakTable if
    current.nmrAtom belongs to nmrResidue
    """
    if data and 'value' in data:
      nmrResidues = data['value']
      self.attachedNmrAtomsList.clear()
      if nmrResidues is not None and len(nmrResidues) > 0 and len(nmrResidues[-1].nmrAtoms) > 0:
        # get the pids and append <all>
        self.ids = [atm.id for atm in nmrResidues[-1].nmrAtoms] + [ALL]
        self.attachedNmrAtomsList.addItems(self.ids)

        # # clear and fill the peak table
        # self.assignedPeaksTable.setObjects([])
        # if self.application.current.nmrAtom is not None and self.application.current.nmrAtom.id in self.ids:
        #   self._updatePeakTable(self.application.current.nmrAtom.id)
        # else:
        #   self._updatePeakTable(ALL)

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
      peaks = list(set([pk for nmrAtom in self.application.current.nmrResidue.nmrAtoms for pk in nmrAtom.assignedPeaks]))

      self.project.blankNotification()
      self._dataFrameObject = self.assignedPeaksTable.getDataFrameFromList(table=self.assignedPeaksTable
                                                                           , buildList=peaks
                                                                           , colDefs=self.getColumns()
                                                                           , hiddenColumns=self._hiddenColumns)

      # populate from the Pandas dataFrame inside the dataFrameObject
      self.assignedPeaksTable.setTableFromDataFrameObject(dataFrameObject=self._dataFrameObject)
      self.project.unblankNotification()

      # self.assignedPeaksTable.setObjects(peaks)
      # highlight current.nmrAtom in the list widget
      self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
      self.peaksLabel.setText('Assigned peaks of NmrAtoms(s): %s' % ALL)
    else:
      pid = 'NA:'+ id
      nmrAtom = self.application.project.getByPid(pid)
      #print('>>', pid, nmrAtom)
      if nmrAtom is not None:

        self.project.blankNotification()
        self._dataFrameObject = self.assignedPeaksTable.getDataFrameFromList(table=self.assignedPeaksTable
                                                          , buildList=nmrAtom.assignedPeaks
                                                          , colDefs=self.getColumns()
                                                          , hiddenColumns=self._hiddenColumns)

        # populate from the Pandas dataFrame inside the dataFrameObject
        self.assignedPeaksTable.setTableFromDataFrameObject(dataFrameObject=self._dataFrameObject)
        self.project.unblankNotification()

        # self.assignedPeaksTable.setObjects(nmrAtom.assignedPeaks)
        # highlight current.nmrAtom in the list widget

        self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
        self.peaksLabel.setText('Assigned peaks of NmrAtom: %s' % nmrAtom.id)

  def getColumns(self):
    "get columns for intialisation of table"
    columns = ColumnClass([('Peak', lambda pk: pk.serial, '', None)
                            , ('Pid', lambda pk: pk.pid, 'Pid of peak', None)
                            , ('id', lambda pk: pk.serial, '', None)])
    tipTexts = []
    # get the maxmimum number of dimensions from all spectra in the project
    numDim = max([sp.dimensionCount for sp in self.application.project.spectra] + [1])

    for i in range(numDim):
      j = i + 1
      c = Column('Assign F%d' % j
                 , lambda pk, dim=i:getPeakAnnotation(pk, dim)
                 , 'NmrAtom assignments of peak in dimension %d' % j
                 , None)
      columns._columns.append(c)

      # columns.append(c)
      # tipTexts.append('NmrAtom assignments of peak in dimension %d' % j)

      sampledDim = self.sampledDims.get(i)
      if sampledDim:
        text = 'Sampled\n%s' % sampledDim.conditionVaried
        tipText='Value of sampled plane'
        unit = sampledDim

      else:
        text = 'Pos F%d' % j
        tipText='Peak position in dimension %d' % j
        unit = 'ppm'

      c = Column(text
                 , lambda pk, dim=i, unit=unit:getPeakPosition(pk, dim, unit)
                 , tipText
                 , None)
      columns._columns.append(c)

      # columns.append(c)
      # tipTexts.append(tipText)

    return columns

  # def _setCurrentPeak(self, peak, row, col):
  def _setCurrentPeak(self, data):
    """
    PeakTable select callback
    """
    peak = data[Notifier.OBJECT]
    # multiselection not allowed, sot only return the first object in list
    if peak:
      self.application.current.peak = peak[0]

  # def _navigateToPeak(self, peak, row, col):
  def _navigateToPeak(self, data):
    """
    PeakTable double-click callback; navigate in to peak in current.strip
    """
    peak = data[Notifier.OBJECT]

    from ccpn.ui.gui.lib.Strip import navigateToPositionInStrip
    #print('>peakTableDoubleClick>', peak)
    if peak is not None and self.application.current.strip is not None:
      self.application.current.peak = peak
      navigateToPositionInStrip(strip=self.application.current.strip, positions=peak.position)

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
