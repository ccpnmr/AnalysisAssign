"""This file contains ModifyAssignmentModule class

intial version by Simon;
extensively modified by Geerten 1-7/12/2016:
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

class ModifyAssignmentModule(CcpnModule, Base):
  """
  This Module allows inspection of the NmrAtoms of a given NmrResidue
  It responds to current.nmrResidues, taking the last added residue to this list
  The NmrAtom listWidget allows for selection of the nmrAtom; subsequently its assignedPeaks
  are displayed.

  """
  def __init__(self, parent=None, project=None, **kw):

    CcpnModule.__init__(self, parent=parent, name='Modify Assignment')
    Base.__init__(self, **kw)

    self.project = project
    self.current = project._appBase.current
    self.sampledDims = {}

    #assignedPeaksLabel = Label(self, '', grid=(0, 0))
    gridLine = 0
    self.nmrAtomLabel = Label(self, 'NmrAtom(s):', grid=(gridLine, 0), gridSpan=(1, 1), vAlign='top')
    self.peaksLabel = Label(self, 'Attached peaks of selected NmrAtom:', grid=(gridLine, 1), gridSpan=(1, 1))
    gridLine += 1

    self.attachedNmrAtomsList = ListWidget(self, grid=(gridLine, 0), gridSpan=(1, 1), hPolicy='fixed', vAlign='top',
                                           callback=self.updatePeakTable, contextMenu=False)
    #gridLine += 1

    self.assignedPeaksTable = ObjectTable(self, self.getColumns(), selectionCallback=self.setCurrentPeak,
                              objects=[], grid=(gridLine, 1), gridSpan=(1, 5))

    self.current.registerNotify(self.updateModule, 'nmrResidues')
    # update if current.nmrResidue is defined
    if self.current.nmrResidue is not None:
      self.updateModule([self.current.nmrResidue])

  def updateModule(self, nmrResidues):
    self.attachedNmrAtomsList.clear()
    if nmrResidues is not None and len(nmrResidues) > 0 and len(nmrResidues[-1].nmrAtoms) > 0:
#      nmrAtom = nmrAtoms[0]
#      self.attachedNmrAtomsList.addItems(list(set([x.id for peak in self.current.nmrAtom.assignedPeaks
#                                                   for dim in peak.dimensionNmrAtoms for x in dim])))
      self.attachedNmrAtomsList.addItems([atm.pid for atm in nmrResidues[-1].nmrAtoms])
      self.assignedPeaksTable.setObjects([])
    else:
      logger.error('No valid nmrAtom/nmrResidue defined')

  def getColumns(self):

    columns = [Column('Id', 'id')]
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
    self.current.peak = peak

  def _getPeakHeight(self, peak):
    """
    Returns the height of the specified peak.
    """
    if peak.height:
      return '%7.2E' % float(peak.height*peak.peakList.spectrum.scale)


  def updatePeakTable(self, item):
    """
    Update the peakTable using item (NmrAtom pid)
    """
    if item is not None:
      self.assignedPeaksTable.setObjects(self.project.getByPid(item.text()).assignedPeaks)
    else:
      logger.error('No valid nmrAtom selected')

  def __del__(self):
    self.current.unRegisterNotify(self.updateModule, 'nmrAtoms')