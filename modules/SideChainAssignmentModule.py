
from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule
from ccpn.ui.gui.lib.SpectrumDisplay import makeStripPlot, makeStripPlotFromSingles

from ccpnmodel.ccpncore.lib.spectrum import Spectrum as SpectrumLib

class SideChainAssignmentModule(PickAndAssignModule):

  def __init__(self, parent=None, project=None):

    PickAndAssignModule.__init__(self, parent, project, name='Sidechain Assignment')

    self.refreshButton.show()
    self.spectrumSelectionWidget.refreshBox.setCallback(self.mediateRefresh)
    self.nmrResidueTable.nmrResidueTable.setTableCallback(self._startAssignment)
    self.mode = 'singles'

  def mediateRefresh(self):
    print('box is %s' %str(self.spectrumSelectionWidget.refreshBox.isChecked()))
    if self.spectrumSelectionWidget.refreshBox.isChecked():
      self.__adminsterNotifiers()
    elif not self.spectrumSelectionWidget.refreshBox.isChecked():
      self.__unRegisterNotifiers()


  def updateModules(self, nmrAtom):
    if not nmrAtom.nmrResidue is self.current.nmrResidue:
      return
    else:
      print('nmrAtom: %s is changing' % nmrAtom.pid)


  def __unRegisterNotifiers(self):
    self.project.unRegisterNotifier('NmrAtom', 'rename', self.updateModules)
    self.project.unRegisterNotifier('NmrAtom', 'create', self.updateModules)
    self.project.unRegisterNotifier('NmrAtom', 'modify', self.updateModules)
    self.project.unRegisterNotifier('NmrAtom', 'delete', self.updateModules)

  def __adminsterNotifiers(self):
    self.project.registerNotifier('NmrAtom', 'rename', self.updateModules)
    self.project.registerNotifier('NmrAtom', 'create', self.updateModules)
    self.project.registerNotifier('NmrAtom', 'modify', self.updateModules)
    self.project.registerNotifier('NmrAtom', 'delete', self.updateModules)

  def _closeModule(self):
    self.__unRegisterNotifiers()
    self.close()


  def _startAssignment(self):
    if self.mode == 'singles':
      self._startAssignmentFromSingles()
    elif self.mode == 'pairs':
      self._startAssignmentFromPairs()


  def _startAssignmentFromPairs(self):
    from ccpn.core.lib.AssignmentLib import getBoundNmrAtomPairs
    activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

    for display in activeDisplays:
      axisCodes = display.strips[0].axisCodes
      nmrAtomPairs = getBoundNmrAtomPairs(self.current.nmrResidue.nmrAtoms, axisCodes[2][0])
      displayIsotopeCodes = [SpectrumLib.name2IsotopeCode(code) for code in axisCodes]
      print(nmrAtomPairs, '1')
      for nmrAtomPair in nmrAtomPairs:
        if not all(x.isotopeCode in displayIsotopeCodes for x in nmrAtomPair):
          nmrAtomPairs.remove(nmrAtomPair)
        # if nmrAtomPair[0].isotopeCode == nmrAtomPair[1].isotopeCode:
        #   if nmrAtomPair[0].isotopeCode == displayIsotopeCodes[2]:
        #     nmrAtomPairs.remove(nmrAtomPair)
      print(nmrAtomPairs, '2')
      makeStripPlot(display, nmrAtomPairs)

  def _startAssignmentFromSingles(self):
    from ccpn.core.lib.AssignmentLib import getBoundNmrAtomPairs
    activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

    for display in activeDisplays:
      axisCodes = display.strips[0].axisCodes
      nmrAtoms = set()
      displayIsotopeCodes = [SpectrumLib.name2IsotopeCode(code) for code in axisCodes]

      for nmrAtom in self.current.nmrResidue.nmrAtoms:
        if nmrAtom.isotopeCode in displayIsotopeCodes:
          nmrAtoms.add(nmrAtom)

      print(nmrAtoms)
      makeStripPlotFromSingles(display, list(nmrAtoms))
      #
      # nmrAtomPairs = getBoundNmrAtomPairs(self.current.nmrResidue.nmrAtoms, axisCodes[2][0])
      # displayIsotopeCodes = [SpectrumLib.name2IsotopeCode(code) for code in axisCodes]
      # print(nmrAtomPairs, '1')
      # for nmrAtomPair in nmrAtomPairs:
      #   if not all(x.isotopeCode in displayIsotopeCodes for x in nmrAtomPair):
      #     nmrAtomPairs.remove(nmrAtomPair)
      #   # if nmrAtomPair[0].isotopeCode == nmrAtomPair[1].isotopeCode:
      #   #   if nmrAtomPair[0].isotopeCode == displayIsotopeCodes[2]:
      #   #     nmrAtomPairs.remove(nmrAtomPair)
      # print(nmrAtomPairs, '2')
      # makeStripPlot(display, nmrAtomPairs)




