
from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule
from ccpn.ui.gui.lib.SpectrumDisplay import makeStripPlot, makeStripPlotFromSingles

from ccpn.ui.gui.lib.Strip import matchAxesAndNmrAtoms
from ccpn.ui.gui.lib.Window import markPositions


from ccpnmodel.ccpncore.lib.spectrum import Spectrum as SpectrumLib

class SideChainAssignmentModule(PickAndAssignModule):

  def __init__(self, parent=None, project=None):

    PickAndAssignModule.__init__(self, parent, project, name='Sidechain Assignment')

    self.refreshButton.show()
    self.refreshButton.setCallback(self._startAssignment)
    self.spectrumSelectionWidget.refreshBox.setCallback(self.mediateRefresh)
    self.nmrResidueTable.nmrResidueTable.setTableCallback(self._startAssignment)
    self.mode = 'pairs'

  def mediateRefresh(self):
    print('box is %s' % str(self.spectrumSelectionWidget.refreshBox.isChecked()))
    if self.spectrumSelectionWidget.refreshBox.isChecked():
      self.__adminsterNotifiers()
    elif not self.spectrumSelectionWidget.refreshBox.isChecked():
      self.__unRegisterNotifiers()


  def updateModules(self, nmrAtom):
    if not nmrAtom.nmrResidue is self.current.nmrResidue:
      return
    else:
      print('nmrAtom: %s is changing' % nmrAtom.pid)

      self._startAssignment()


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
    self.project._appBase.ui.mainWindow.clearMarks()
    if self.mode == 'singles':
      self._startAssignmentFromSingles()
    elif self.mode == 'pairs':
      self._startAssignmentFromPairs()


  def _startAssignmentFromPairs(self):
    from ccpn.core.lib.AssignmentLib import getBoundNmrAtomPairs
    activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

    for display in activeDisplays:
      axisCodes = display.strips[0].axisCodes
      nmrAtomPairs = getBoundNmrAtomPairs(self.current.nmrResidue.nmrAtoms, axisCodes[-1][0])
      displayIsotopeCodes = [SpectrumLib.name2IsotopeCode(code) for code in axisCodes]
      pairsToRemove = []
      for nmrAtomPair in nmrAtomPairs:
        pairIsotopeCodes = [nap.isotopeCode for nap in nmrAtomPair]
        nmrAtoms = set()
        if displayIsotopeCodes[1] in pairIsotopeCodes and displayIsotopeCodes[0] not in pairIsotopeCodes:
          pairsToRemove.append(nmrAtomPair)
          nmrAtoms.add(nmrAtomPair[0])
          nmrAtoms.add(nmrAtomPair[1])
        if not all(x.isotopeCode in displayIsotopeCodes for x in nmrAtomPair):
          pairsToRemove.append(nmrAtomPair)
          nmrAtoms.add(nmrAtomPair[0])
          nmrAtoms.add(nmrAtomPair[1])
        elif nmrAtomPair[0].isotopeCode == nmrAtomPair[1].isotopeCode and not self.has_duplicates(displayIsotopeCodes):
          pairsToRemove.append(nmrAtomPair)
          nmrAtoms.add(nmrAtomPair[0])
          nmrAtoms.add(nmrAtomPair[1])
        if len(displayIsotopeCodes) > 2:
          if nmrAtomPair[0].isotopeCode == nmrAtomPair[1].isotopeCode and displayIsotopeCodes[0] != displayIsotopeCodes[2]:
            if displayIsotopeCodes.count(nmrAtomPair[0].isotopeCode) != 2:
              nmrAtoms.add(nmrAtomPair[0])
              nmrAtoms.add(nmrAtomPair[1])
              pairsToRemove.append(nmrAtomPair)
      for pair in pairsToRemove:
        if pair in nmrAtomPairs:
          nmrAtomPairs.remove(pair)
      if len(nmrAtomPairs) > 1:
        sortedNmrAtomPairs = self.sortNmrAtomPairs(nmrAtomPairs)

      else:
        sortedNmrAtomPairs = nmrAtomPairs

      if len(display.strips[0].axisCodes) > 2:
        makeStripPlot(display, sortedNmrAtomPairs, autoWidth=False)
      nmrAtoms = [x for x in nmrAtomPairs for x in x]
      axisCodePositionDict = matchAxesAndNmrAtoms(display.strips[0], nmrAtoms)
      markPositions(self.project, list(axisCodePositionDict.keys()), list(axisCodePositionDict.values()))


  def has_duplicates(self, seq):
    return any(seq.count(x) > 1 for x in seq)

  def sortNmrAtomPairs(self, nmrAtomPairs):
    order = ['C', 'CA', 'CB', 'CG', 'CG1', 'CG2', 'CD1', 'CD2', 'CE', 'CZ', 'N', 'ND', 'NE', 'NZ', 'NH',
             'H', 'HA', 'HB', 'HG', 'HD',  'HE',  'HZ', 'HH']
    ordering = []
    for p in nmrAtomPairs:
      if p[0].name[:len(p[0].name)] in order:
          ordering.append((order.index(p[0].name[:len(p[0].name)]), p))

    if len(nmrAtomPairs) > 1:
      sortedNmrAtomPairs = [x[1] for x in sorted(ordering, key=lambda x: x[0])]
    else:
      sortedNmrAtomPairs = nmrAtomPairs
    return sortedNmrAtomPairs

  def _startAssignmentFromSingles(self):
    activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

    for display in activeDisplays:
      axisCodes = display.strips[0].axisCodes
      nmrAtoms = set()
      displayIsotopeCodes = [SpectrumLib.name2IsotopeCode(code) for code in axisCodes]

      for nmrAtom in self.current.nmrResidue.nmrAtoms:
        if nmrAtom.isotopeCode in displayIsotopeCodes and nmrAtom.isotopeCode == displayIsotopeCodes[2]:
          nmrAtoms.add(nmrAtom)

      makeStripPlotFromSingles(display, list(nmrAtoms))
      axisCodePositionDict = matchAxesAndNmrAtoms(display.strips[0], nmrAtoms)
      markPositions(self.project, list(axisCodePositionDict.keys()), list(axisCodePositionDict.values()))


