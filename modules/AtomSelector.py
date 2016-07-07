"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date$"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon P Skinner, Geerten W Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/license"
              "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/license"
                " or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Revision$"

#=========================================================================================
# Start of code
#=========================================================================================

import os
import typing
from functools import partial

from PyQt4 import QtCore

from ccpn.core.Peak import Peak
from ccpn.core.lib.Assignment import isInterOnlyExpt, getNmrAtomPrediction, CCP_CODES
from ccpn.ui.gui.lib.PeakAssignment import peaksAreOnLine
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Module import CcpnModule
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.RadioButton import RadioButton
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.util import Path
from ccpnmodel.ccpncore.lib.assignment.ChemicalShift import PROTEIN_ATOM_NAMES, ALL_ATOMS_SORTED
from ccpnmodel.ccpncore.lib.spectrum import Spectrum as spectrumLib


class AtomSelector(CcpnModule):
  """
  Module to be used with PickAndAssignModule for prediction of nmrAtom names and assignment of nmrAtoms
  to peak dimensions
  """
  def __init__(self, parent, project=None):
    CcpnModule.__init__(self, name='Atom Selector')
    self.orientation = 'vertical'
    self.moveLabel=False
    self.pickAndAssignWidget = Widget(self)
    # self.pickAndAssignWidget.setMaximumSize(500, 300)
    self.parent = parent

    if self.parent is not None:
      self.pythonConsole = self.parent.pythonConsole
    else:
      self.pythonConsole = self._parent.pythonConsole

    self.current = self.parent._appBase.current
    self.project = project
    self.current.registerNotify(self._predictAssignments, 'peaks')
    # self.current.registerNotify(self.predictAssignments, 'peaks')
    nmrResidueLabel = Label(self, 'Current NmrResidue', grid=(0, 0))
    gridLine = -1
    # gridLine 0
    gridLine += 1
    self.molTypeLabel = Label(self, 'Molecule Type', grid=(gridLine, 0))
    self.molTypePulldown = PulldownList(self, grid=(gridLine, 1))
    self.molTypePulldown.setData(['protein', 'DNA', 'RNA', 'carbohydrate', 'other'])
    self.radioButton1 = RadioButton(self, grid=(gridLine, 2), hAlign='r', callback=self._createBackBoneButtons)
    self.radioButton1.setChecked(True)
    self.label1 = Label(self, 'Backbone', grid=(gridLine, 3), hAlign='l')
    self.radioButton2 = RadioButton(self, grid=(gridLine, 4), hAlign='r', callback=self._createSideChainButtons)
    self.label2 = Label(self, 'Side chain', grid=(gridLine, 5), hAlign='l')
    # gridLine 1
    gridLine += 1
    nmrResidueLabel = Label(self, 'Current NmrResidue:', grid=(gridLine, 0))
    self.currentNmrResidueLabel = Label(self, grid=(gridLine, 1))
    # gridLine 2
    gridLine += 1
    self.layout.addWidget(self.pickAndAssignWidget, gridLine, 0, 3, 8)
    self.current.registerNotify(self._updateWidget, 'nmrResidues')
    self.closeModule = self._closeModule
    self.buttons = {}
    
  def _closeModule(self):
    self.current.unRegisterNotify(self._predictAssignments, 'peaks')
    self.current.unRegisterNotify(self._updateWidget, 'nmrResidues')
    self.close()

  def _updateWidget(self, nmrResidues=None):
    self.currentNmrResidueLabel.setText(nmrResidues[0].id)
    if self.radioButton1.isChecked():
      self._createBackBoneButtons()
    elif self.radioButton2.isChecked():
      self._createSideChainButtons()
    else:
      return


  def _createBackBoneButtons(self):
    self._cleanupPickAndAssignWidget()
    headerLabel = Label(self, text='i-1')
    self.pickAndAssignWidget.layout().addWidget(headerLabel, 0, 0)
    headerLabel2 = Label(self.pickAndAssignWidget, text='i', grid=(0, 1))
    headerLabel3 = Label(self.pickAndAssignWidget, text='i+1', grid=(0, 2))
    self.buttons = {}
    atoms = ['H', 'N', 'CA', 'CB', 'CO', 'HA', 'HB']
    for ii, atom in enumerate(atoms):
      self.buttons[atom] = []
      button1 = Button(self.pickAndAssignWidget, text=atom+' [i-1]', grid=(1+ii, 0), callback=partial(self._pickAndAssign, '-1', atom))
      button2 = Button(self.pickAndAssignWidget, text=atom+' [i]', grid=(1+ii, 1), callback=partial(self._pickAndAssign, '', atom))
      button3 = Button(self.pickAndAssignWidget, text=atom+' [i+1]', grid=(1+ii, 2), callback=partial(self._pickAndAssign, '+1', atom))
      self.buttons[atom].append(button1)
      self.buttons[atom].append(button2)
      self.buttons[atom].append(button3)

    for buttons in self.buttons.values():
      for button in buttons:
        button.clicked.connect(self._returnButtonsToNormal)



  def _createSideChainButtons(self):
    self._cleanupPickAndAssignWidget()
    self.hCheckBox = CheckBox(self.pickAndAssignWidget, hAlign='r', callback=self._toggleBox)
    self.hCheckBox.setChecked(True)
    self.pickAndAssignWidget.layout().addWidget(self.hCheckBox, 0, 0, QtCore.Qt.AlignRight)
    self.hLabel = Label(self.pickAndAssignWidget, 'H', grid=(0, 1), hAlign='l')
    self.cCheckBox = CheckBox(self.pickAndAssignWidget, grid=(0, 2), hAlign='r', callback=self._toggleBox)
    self.cCheckBox.setChecked(True)
    self.cLabel = Label(self.pickAndAssignWidget, 'C', grid=(0, 3), hAlign='l')
    self.nCheckBox = CheckBox(self.pickAndAssignWidget, grid=(0, 4), hAlign='r', callback=self._toggleBox)
    self.nLabel = Label(self.pickAndAssignWidget, 'N', grid=(0, 5), hAlign='l')
    self.nCheckBox.setChecked(True)
    self.otherCheckBox = CheckBox(self.pickAndAssignWidget, grid=(0, 6), hAlign='r', callback=self._toggleBox)
    self.otherLabel = Label(self.pickAndAssignWidget, 'Other', grid=(0, 7), hAlign='l')

    self._updateLayout()


  def _toggleBox(self):
    self._updateLayout()

  def _getAtomsForButtons(self, atomList, atomName):
    [atomList.remove(atom) for atom in sorted(atomList) if atom[0] == atomName]

  def _updateLayout(self):

    # group atoms in useful categories based on usage

    alphaAtoms = [x for x in ALL_ATOMS_SORTED['alphas']]
    betaAtoms = [x for x in ALL_ATOMS_SORTED['betas']]
    gammaAtoms = [x for x in ALL_ATOMS_SORTED['gammas']]
    moreGammaAtoms = [x for x in ALL_ATOMS_SORTED['moreGammas']]
    deltaAtoms = [x for x in ALL_ATOMS_SORTED['deltas']]
    moreDeltaAtoms = [x for x in ALL_ATOMS_SORTED['moreDeltas']]
    epsilonAtoms = [x for x in ALL_ATOMS_SORTED['epsilons']]
    moreEpsilonAtoms = [x for x in ALL_ATOMS_SORTED['moreEpsilons']]
    zetaAtoms = [x for x in ALL_ATOMS_SORTED['zetas']]
    etaAtoms = [x for x in ALL_ATOMS_SORTED['etas']]
    moreEtaAtoms = [x for x in ALL_ATOMS_SORTED['moreEtas']]
    atomButtonList = [alphaAtoms, betaAtoms, gammaAtoms, moreGammaAtoms, deltaAtoms, moreDeltaAtoms,
                      epsilonAtoms, moreEpsilonAtoms, zetaAtoms, etaAtoms, moreEtaAtoms]
    # Activate button for Carbons
    if not self.cCheckBox.isChecked():
      [self._getAtomsForButtons(atomList, 'C') for atomList in atomButtonList]

    if not self.hCheckBox.isChecked():
      [self._getAtomsForButtons(atomList, 'H') for atomList in atomButtonList]

    if not self.nCheckBox.isChecked():
      [self._getAtomsForButtons(atomList, 'N') for atomList in atomButtonList]

    rowCount = self.pickAndAssignWidget.layout().rowCount()
    colCount = self.pickAndAssignWidget.layout().columnCount()

    for r in range(1, rowCount):
      for m in range(colCount):
        item = self.pickAndAssignWidget.layout().itemAtPosition(r, m)
        if item:
          if item.widget():
            item.widget().hide()
        self.pickAndAssignWidget.layout().removeItem(item)

    if self.current.nmrResidue:
      self.currentNmrResidueLabel.setText(self.current.nmrResidue.id)
      if self.current.nmrResidue.residueType == '':
        for ii, atomList in enumerate(atomButtonList):

          for jj, atom in enumerate(atomList):
            button = Button(self.pickAndAssignWidget, text=atom, grid=(ii+1, jj), hAlign='t', callback=partial(self._pickAndAssign, '0', atom))
            button.setMinimumSize(45, 20)


      else:
        self.buttons = {}
        residueType = self.current.nmrResidue.residueType.upper()
        residueAtoms = PROTEIN_ATOM_NAMES[residueType]
        residueAlphas = [atom for atom in alphaAtoms if atom in residueAtoms]
        residueBetas = [atom for atom in betaAtoms if atom in residueAtoms]
        residueGammas = [atom for atom in gammaAtoms if atom in residueAtoms]
        residueMoreGammas = [atom for atom in moreGammaAtoms if atom in residueAtoms]
        residueDeltas = [atom for atom in deltaAtoms if atom in residueAtoms]
        residueMoreDeltas = [atom for atom in moreDeltaAtoms if atom in residueAtoms]
        residueEpsilons = [atom for atom in epsilonAtoms if atom in residueAtoms]
        residueMoreEpsilons = [atom for atom in moreEpsilonAtoms if atom in residueAtoms]
        residueZetas = [atom for atom in zetaAtoms if atom in residueAtoms]
        residueEtas = [atom for atom in etaAtoms if atom in residueAtoms]
        residueMoreEtas = [atom for atom in moreEtaAtoms if atom in residueAtoms]
        atomButtonList2 = [residueAlphas, residueBetas, residueGammas, residueMoreGammas, residueDeltas, residueMoreDeltas,
                      residueEpsilons, residueMoreEpsilons, residueZetas, residueEtas, residueMoreEtas]
        for ii, atomList in enumerate(atomButtonList2):
          for jj, atom in enumerate(atomList):
            self.buttons[atom] = []
            button = Button(self.pickAndAssignWidget, text=atom, grid=(ii+1, jj), hAlign='t', callback=partial(self._pickAndAssign, '0', atom))
            button.setMinimumSize(45, 20)
            self.buttons[atom].append(button)

  def _showMoreAtomButtons(self, buttons, moreButton):
    if moreButton.isChecked():
      [button.show() for button in buttons]
    else:
      [button.hide() for button in buttons]

  def _cleanupPickAndAssignWidget(self):

    layout = self.pickAndAssignWidget.layout()
    for r in range(layout.rowCount()):
      for c in range(layout.columnCount()):
        item = layout.itemAtPosition(r, c)
        if item:
          if item.widget():
            item.widget().hide()
        layout.removeItem(item)


  def _pickAndAssign(self, position:int, atomType:str):
    """
    Takes a position either -1, 0 or +1 and an atom type, fetches an NmrAtom with name corresponding
    to the atom type and the position and assigns it to correct dimension of current.peaks
    """
    if not self.current.nmrResidue:
      return
    name = atomType
    if position == '-1' and '-1' not in self.current.nmrResidue.sequenceCode:
      r = self.current.nmrResidue.previousNmrResidue
      if not r:
        r = self.current.nmrResidue.nmrChain.fetchNmrResidue(sequenceCode=self.current.nmrResidue.sequenceCode+'-1')
    else:
      r = self.current.nmrResidue

    newNmrAtom = r.fetchNmrAtom(name=name)
    # Replaced  by in-command echo
    # self.pythonConsole.writeConsoleCommand(
    #   "nmrAtom = nmrResidue.fetchNmrAtom(name='%s')" % name, nmrResidue=r.pid
    # )
    for peak in self.current.peaks:
      for strip in self.project.strips:
        for peakListView in strip.peakListViews:
          if peak in peakListView.peakItems.keys():
            # NBNB TBD FIXME consider using new _displayOrderSpectrumDimensionIndices
            # function to map display axes to spectrum axes
            axisCode = spectrumLib.axisCodeMatch(strip.axisCodes[1], peak.peakList.spectrum.axisCodes)
            index = peak.peakList.spectrum.axisCodes.index(axisCode)
            nmrAtoms = peak.dimensionNmrAtoms[index] + [newNmrAtom]
            peak.assignDimension(axisCode, nmrAtoms)
    self._returnButtonsToNormal()


  def _returnButtonsToNormal(self):
    """
    Returns all buttons in Atom Selector to original colours and style.
    """
    if self.parent._appBase.preferences.general.colourScheme == 'dark':
      styleSheet = open(os.path.join(Path.getPathToImport('ccpn.ui.gui.widgets'),
                                     'DarkStyleSheet.qss')).read()
      self.setStyleSheet('''ModuleLabel  {
                                        background-color: #BEC4F3;
                                        color: #122043;
                                        border: 1px solid #00092D;
                                       }''')
      for buttons in self.buttons.values():
        for button in buttons:
          button.setStyleSheet(styleSheet)
      self.setStyleSheet(styleSheet)

    elif self.parent._appBase.preferences.general.colourScheme == 'light':
      styleSheet = open(os.path.join(Path.getPathToImport('ccpn.ui.gui.widgets'),
                                     'LightStyleSheet.qss')).read()
      for buttons in self.buttons.values():
        for button in buttons:
          button.setStyleSheet(styleSheet)
      self.setStyleSheet(styleSheet)


  def _predictAssignments(self, peaks:typing.List[Peak]):
    """
    Predicts atom type for selected peaks and highlights the relevant buttons with confidence of
    that assignment prediction, green is very confident, orange is less confident.
    """
    self._returnButtonsToNormal()
    if not self.current.nmrResidue or None in peaks:
      return

    else:

      if len(peaks) == 0:
        for buttons in self.buttons.values():
          for button in buttons:
            button.clicked.connect(self._returnButtonsToNormal)

      else:

        if peaksAreOnLine(peaks, 1):
          # RHF 15/12/2015 bug fix. From usage elsewhere this ought to use types, not name.
          # NBNB TBD CHECK THIS
          # experiments = [peak.peakList.spectrum.experimentName for peak in peaks]
          types = set(peak.peakList.spectrum.experimentType for peak in peaks)
          anyInterOnlyExperiments = any(isInterOnlyExpt(x) for x in types)

          for peak in peaks:
            peakListViews = [peakListView for peakListView in self.project.peakListViews if peakListView.peakList == peak.peakList]
            axisCodes = [peakListView.spectrumView._parent.axisOrder for peakListView in peakListViews]
            mapping = spectrumLib._axisCodeMapIndices(axisCodes[0], peak.peakList.spectrum.axisCodes)
            isotopeCode = peak.peakList.spectrum.isotopeCodes[mapping[1]]
            if self.radioButton1.isChecked():
              predictedAtomTypes = [getNmrAtomPrediction(ccpCode, peak.position[mapping[1]], isotopeCode, strict=True) for ccpCode in CCP_CODES]
              refinedPreds = [[type[0][0][1], type[0][1]] for type in predictedAtomTypes if len(type) > 0]
              atomPredictions = set()
              for pred in refinedPreds:
                if pred[1] > 90:
                  atomPredictions.add(pred[0])
              for atomPred in atomPredictions:
                if atomPred == 'CB':
                  # if(any(isInterOnlyExpt(experiment) for experiment in experiments)):
                  if anyInterOnlyExperiments:
                    self.buttons['CB'][0].setStyleSheet('background-color: green')
                  else:
                    # self.cbButton1.setStyleSheet('background-color: orange')
                    self.buttons['CB'][0].setStyleSheet('background-color: green')
                    self.buttons['CB'][1].setStyleSheet('background-color: green')
                if atomPred == 'CA':
                  # if(any(isInterOnlyExpt(experiment) for experiment in experiments)):
                  if anyInterOnlyExperiments:
                    self.buttons['CA'][0].setStyleSheet('background-color: green')
                  else:
                    # self.caButton1.setStyleSheet('background-color: orange')
                    self.buttons['CA'][0].setStyleSheet('background-color: green')
                    self.buttons['CA'][1].setStyleSheet('background-color: green')
            elif self.radioButton2.isChecked():
              if self.current.nmrResidue.residueType == '':
                predictedAtomTypes = [getNmrAtomPrediction(ccpCode, peak.position[mapping[1]], isotopeCode) for ccpCode in CCP_CODES]
              else:
                predictedAtomTypes = getNmrAtomPrediction(self.current.nmrResidue.residueType.title(), peak.position[mapping[1]], isotopeCode)
              for type in predictedAtomTypes:
                for atomType, buttons in self.buttons.items():
                  if type[0][1] == atomType:
                    for button in buttons:
                      if type[1] > 85:
                        button.setStyleSheet('background-color: green')
                      elif 50 < type[1] < 85:
                        button.setStyleSheet('background-color: orange')
                      if type[1] < 50:
                        button.setStyleSheet('background-color: red')



