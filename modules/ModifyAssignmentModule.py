"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - : 2014-06-04 18:13:10 +0100 (Wed, 04 Jun 2014) $"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon P Skinner, Geerten W Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/license"
               "or ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/license"
                 " or ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = ": rhfogh $"
__date__ = ": 2014-06-04 18:13:10 +0100 (Wed, 04 Jun 2014) $"
__version__ = ": 7686 $"

#=========================================================================================
# Start of code
#=========================================================================================


from PyQt4 import QtGui

from ccpn.ui.gui.widgets.Base import Base
from ccpn.ui.gui.widgets.Module import CcpnModule
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.Table import ObjectTable, Column
from ccpn.ui.gui.modules.PeakTable import PeakListSimple

from ccpn.ui.gui.modules.peakUtils import getPeakPosition, getPeakAnnotation

class ModifyAssignmentModule(CcpnModule, Base):

  def __init__(self, parent=None, project=None, nmrAtom=None, **kw):

    CcpnModule.__init__(self, parent=None, name='Modify Assignment')
    Base.__init__(self, **kw)
    self.project = project
    assignedPeaksLabel = Label(self, '', grid=(0, 0))
    nmrAtomLabel = Label(self, 'Current NmrAtom', grid=(1, 0))
    self.currentNmrAtomLabel = Label(self, nmrAtom.id, grid=(1, 1))
    attachedNmrAtomsLabel = Label(self, '', grid=(2, 0))
    self.attachedNmrAtomsList = ListWidget(self, grid=(3, 0), gridSpan=(1, 2), callback=self.updatePeakTable)
    assignedPeaksLabel = Label(self, '', grid=(4, 0))
    self.current = project._appBase.current
    self.sampledDims = {}
    self.current.nmrAtom = nmrAtom
    self.assignedPeaksTable = ObjectTable(self, self.getColumns(), actionCallback=self.setCurrentPeak,
                              objects=[], grid=(5, 0), gridSpan=(1, 2))

    self.updateModule([nmrAtom])
    self.current.registerNotify(self.updateModule, 'nmrAtoms')




  def updateModule(self, nmrAtoms):
    self.attachedNmrAtomsList.clear()
    nmrAtom = nmrAtoms[0]
    self.attachedNmrAtomsList.addItems(list(set([x.id for peak in nmrAtom.assignedPeaks
                                                 for dim in peak.dimensionNmrAtoms for x in dim])))

  def getColumns(self):

    columns = [Column('Id', 'id')]
    tipTexts = []
    numDim = max([len(pk.position) for pk in self.current.nmrAtom.assignedPeaks])
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
    self.assignedPeaksTable.setObjects(self.project.getByPid('NA:'+item.text()).assignedPeaks)


  def __del__(self):
    self.current.unRegisterNotify(self.updateModule, 'nmrAtoms')