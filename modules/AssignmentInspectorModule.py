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
__dateModified__ = "$dateModified: 2017-04-07 11:40:21 +0100 (Fri, April 07, 2017) $"
__version__ = "$Revision: 3.0.b1 $"
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
from ccpn.ui.gui.widgets.Table import ObjectTable, Column
from ccpn.util.Logging import getLogger
from core.lib.peakUtils import getPeakPosition, getPeakAnnotation

logger = getLogger()


class AssignmentInspectorModule(CcpnModule):
  """
  This Module allows inspection of the NmrAtoms of a selected NmrResidue
  It responds to current.nmrResidues, taking the last added residue to this list
  The NmrAtom listWidget allows for selection of the nmrAtom; subsequently its assignedPeaks
  are displayed.

  """

  ALL = '<all>'

  # overide in specific module implementations
  className = 'AssignmentInspectorModule'
  includeSettingsWidget = False
  maxSettingsState = 3  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsOnTop = True

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
    self.assignedPeaksTable = ObjectTable(self.frame2, self.getColumns(),
                                          selectionCallback=self._setCurrentPeak,
                                          actionCallback=self._navigateToPeak,
                                          objects=[], autoResize=True,
                                          grid=(1, 0), gridSpan=(1, 5), **policies
                                          )
    #self.attachedNmrAtomsList.setFixedHeight(200)
    #self.assignedPeaksTable.setFixedHeight(200)

    self.application.current.registerNotify(self._updateModuleCallback, 'nmrResidues')
    # update if current.nmrResidue is defined
    if self.application.current.nmrResidue is not None:
      self._updateModuleCallback([self.application.current.nmrResidue])

  def _closeModule(self):
    self.application.current.unRegisterNotify(self._updateModuleCallback, 'nmrResidues')
    super(AssignmentInspectorModule, self)._closeModule()

  def _updateModuleCallback(self, nmrResidues):
    """
    Callback function: Module responsive to nmrResidues; updates the list widget with nmrAtoms and updates peakTable if
    current.nmrAtom belongs to nmrResidue
    """
    self.attachedNmrAtomsList.clear()
    if nmrResidues is not None and len(nmrResidues) > 0 and len(nmrResidues[-1].nmrAtoms) > 0:
      # get the pids and append <all>
      self.ids = [atm.id for atm in nmrResidues[-1].nmrAtoms] + [self.ALL]
      self.attachedNmrAtomsList.addItems(self.ids)
      # clear and fill the peak table
      self.assignedPeaksTable.setObjects([])
      if self.application.current.nmrAtom is not None and self.application.current.nmrAtom.id in self.ids:
        self._updatePeakTable(self.application.current.nmrAtom.id)
      else:
        self._updatePeakTable(self.ALL)
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
      self.assignedPeaksTable.setObjects([])
      return

    if id == self.ALL:
      peaks = list(set([pk for nmrAtom in self.application.current.nmrResidue.nmrAtoms for pk in nmrAtom.assignedPeaks]))
      self.assignedPeaksTable.setObjects(peaks)
      # highlight current.nmrAtom in the list widget
      self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
      self.peaksLabel.setText('Assigned peaks of NmrAtoms(s): %s' % self.ALL)
    else:
      pid = 'NA:'+ id
      nmrAtom = self.application.project.getByPid(pid)
      #print('>>', pid, nmrAtom)
      if nmrAtom is not None:
        self.assignedPeaksTable.setObjects(nmrAtom.assignedPeaks)
        # highlight current.nmrAtom in the list widget
        self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
        self.peaksLabel.setText('Assigned peaks of NmrAtom: %s' % nmrAtom.id)

  def getColumns(self):
    "get collumns for intialisation of table"
    columns = [Column('Peak', 'id')]
    tipTexts = []
    # get the maxmimum number of dimensions from all spectra in the project
    numDim = max([sp.dimensionCount for sp in self.application.project.spectra] + [1])

    for i in range(numDim):
      j = i + 1
      c = Column('Assign F%d' % j, lambda pk, dim=i:getPeakAnnotation(pk, dim))
      columns.append(c)
      tipTexts.append('NmrAtom assignments of peak in dimension %d' % j)

      sampledDim = self.sampledDims.get(i)
      if sampledDim:
        text = 'Sampled\n%s' % sampledDim.conditionVaried
        tipText='Value of sampled plane'
        unit = sampledDim

      else:
        text = 'Pos F%d' % j
        tipText='Peak position in dimension %d' % j
        unit = 'ppm'
      c = Column(text, lambda pk, dim=i, unit=unit:getPeakPosition(pk, dim, unit))
      columns.append(c)
      tipTexts.append(tipText)

    return columns

  def _setCurrentPeak(self, peak, row, col):
    """
    PeakTable select callback
    """
    if peak is not None:
      self.application.current.peak = peak

  def _navigateToPeak(self, peak, row, col):
    """
    PeakTable double-click callback; navigate in to peak in current.strip
    """
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
