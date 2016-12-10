"""This file contains ModifyAssignmentModule class

intial version by Simon;
extensively modified by Geerten 1-9/12/2016:
- intialisation with 'empty' settings possible,
- now responsive to current.nmrResidues
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date: 2016-07-09 14:17:30 +0100 (Sat, 09 Jul 2016) $"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon P Skinner, Geerten W Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/license"
              "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/license"
                " or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author: geertenv $"
__date__ = "$Date: 2016-07-09 14:17:30 +0100 (Sat, 09 Jul 2016) $"
__version__ = "$Revision: 9605 $"

#=========================================================================================
# Start of code
#=========================================================================================


from ccpn.ui.gui.widgets.Base import Base
from ccpn.ui.gui.widgets.Module import CcpnModule
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.Table import ObjectTable, Column
from ccpn.ui.gui.modules.PeakTable import PeakListSimple
from ccpn.util.Logging import getLogger

logger = getLogger()

from ccpn.ui.gui.modules.peakUtils import getPeakPosition, getPeakAnnotation

class ModifyAssignmentModule(CcpnModule):
  """
  This Module allows inspection of the NmrAtoms of a selected NmrResidue
  It responds to current.nmrResidues, taking the last added residue to this list
  The NmrAtom listWidget allows for selection of the nmrAtom; subsequently its assignedPeaks
  are displayed.

  """
  ALL = '<all>'

  def __init__(self, parent=None, project=None, **kw):

    CcpnModule.__init__(self, parent=parent, name='Modify Assignment')
    #Base.__init__(self, **kw)

    self.project = project
    self.current = project._appBase.current
    self.sampledDims = {}
    self.pids = []  # list of currently displayed NmrAtom pids + <all>

    #assignedPeaksLabel = Label(self, '', grid=(0, 0))
    #self.mainWidget.layout().setSpacing(0)
    #self.mainWidget.layout().setContentsMargins(0,0,0,0)

    gridLine = 0
    self.nmrAtomLabel = Label(self.mainWidget, 'NmrAtom(s):', grid=(gridLine, 0), gridSpan=(1, 1), vAlign='top')
    self.peaksLabel = Label(self.mainWidget, 'Assigned peaks of selected NmrAtom:', grid=(gridLine, 1), gridSpan=(1, 1))
    gridLine += 1

    self.attachedNmrAtomsList = ListWidget(self.mainWidget, grid=(gridLine, 0), gridSpan=(1, 1), vAlign='top',
                                           callback=self._updatePeakTableCallback, contextMenu=False)
    self.attachedNmrAtomsList.setFixedWidth(120)
    #gridLine += 1

    self.assignedPeaksTable = ObjectTable(self.mainWidget, self.getColumns(), selectionCallback=self.setCurrentPeak,
                              objects=[], grid=(gridLine, 1), gridSpan=(1, 5))

    self.current.registerNotify(self._updateModuleCallback, 'nmrResidues')
    # update if current.nmrResidue is defined
    if self.current.nmrResidue is not None:
      self._updateModuleCallback([self.current.nmrResidue])

  def _updateModuleCallback(self, nmrResidues):
    """
    Callback function: Module responsive to nmrResidues; updates the list widget with nmrAtoms and updates peakTable if
    current.nmrAtom belongs to nmrResidue
    """
    self.attachedNmrAtomsList.clear()
    if nmrResidues is not None and len(nmrResidues) > 0 and len(nmrResidues[-1].nmrAtoms) > 0:
      # get the pids and append <all>
      self.pids = [atm.pid for atm in nmrResidues[-1].nmrAtoms] + [self.ALL]
      self.attachedNmrAtomsList.addItems(self.pids)
      # clear and fill the peak table
      self.assignedPeaksTable.setObjects([])
      if self.current.nmrAtom is not None and self.current.nmrAtom.pid in self.pids:
        self._updatePeakTable(self.current.nmrAtom.pid)
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

  def _updatePeakTable(self, pid):
    """
    Update peak table depending on value of pid;
    clears peakTable if pid is None
    """
    if pid is None:
      self.assignedPeaksTable.setObjects([])
      return

    if pid == self.ALL:
      peaks = list(set([pk for nmrAtom in self.current.nmrResidue.nmrAtoms for pk in nmrAtom.assignedPeaks]))
      self.assignedPeaksTable.setObjects(peaks)
      # highlight current.nmrAtom in the list widget
      self.attachedNmrAtomsList.setCurrentRow(self.pids.index(pid))
    else:
      nmrAtom = self.project.getByPid(pid)
      if nmrAtom is not None:
        self.assignedPeaksTable.setObjects(nmrAtom.assignedPeaks)
        # highlight current.nmrAtom in the list widget
        self.attachedNmrAtomsList.setCurrentRow(self.pids.index(pid))

  def getColumns(self):

    columns = [Column('Peak', 'id')]
    tipTexts = []
    # get the maxmimum number of dimensions from all spectra in the project
    numDim = max([sp.dimensionCount for sp in self.project.spectra] + [1])

    for i in range(numDim):
      j = i + 1
      c = Column('Assign F%d' % j, lambda pk, dim=i:getPeakAnnotation(pk, dim))
      columns.append(c)
      tipTexts.append('NmrAtom assignments of peak in dimension %d' % j)

    for i in range(numDim):
      j = i + 1

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

  def setCurrentPeak(self, peak, row, col):
    if peak is not None:
      self.current.peak = peak

  def _getPeakHeight(self, peak):
    """
    Returns the height of the specified peak.
    """
    if peak.height:
      return '%7.2E' % float(peak.height*peak.peakList.spectrum.scale)

  def __del__(self):
    self.current.unRegisterNotify(self.updateModule, 'nmrAtoms')