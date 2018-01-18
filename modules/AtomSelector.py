"""
AtomSelector module: 

used for backbone and sidechain assignments to set the nmrAtom of a (number of) selected
peaks to a specific nucleus. Thus far, it set the 'y-axis' dimension of the peak (to be
made flexible in the settings tab). Options other than 'protein' are not yet implemented.

Original by SS
First rework by GWV
"""

#TODO:GEERTEN Needs complete refactoring:
"""
- buttons are destroyed (!?) and create with every refresh; 
- atom type prediction in sidechain mode (when type of nmrResidue is not yet known) gives
  the result of the last residue (i.e. Tyr), rather then an average
- Should retain predictions of selected peaks when switching backbone/sidechain
- assignSelected should be refactored to be proper for all general cases

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
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

# import os
import typing
from functools import partial

from PyQt4 import QtCore, QtGui

from ccpn.core.Peak import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.lib import Pid
from ccpn.core.lib.AssignmentLib import isInterOnlyExpt, getNmrAtomPrediction, CCP_CODES
from ccpn.core.lib.AssignmentLib import peaksAreOnLine
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.RadioButton import RadioButton
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.MessageDialog import showYesNo
from ccpnmodel.ccpncore.lib.assignment.ChemicalShift import PROTEIN_ATOM_NAMES, ALL_ATOMS_SORTED

from ccpn.util.Logging import getLogger
from ccpn.ui.gui.widgets.MessageDialog import showWarning

logger = getLogger()


class AtomSelectorModule(CcpnModule):
  """
  Module to be used with PickAndAssignModule for prediction of nmrAtom names and assignment of nmrAtoms
  to peak dimensions
  Responds to current.nmrResidue and current.peaks
  """
  includeSettingsWidget = True
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsPosition = 'top'

  def __init__(self, mainWindow=None, name='Atom Selector', nmrAtom=None):
    CcpnModule.__init__(self, mainWindow=mainWindow, name=name)

    # Derive application, project, and current from mainWindow
    self.mainWindow = mainWindow
    self.application = mainWindow.application
    self.project = mainWindow.application.project
    self.current = mainWindow.application.current

    self.current.registerNotify(self._predictAssignments, 'peaks')
    self.current.registerNotify(self._nmrResidueCallBack, 'nmrResidues')

    # Settings Widget
    self.molTypeLabel = Label(self.settingsWidget, 'Molecule Type', grid=(0, 0))
    self.molTypePulldown = PulldownList(self.settingsWidget, grid=(0, 1), texts=['protein', 'DNA', 'RNA', 'carbohydrate', 'other'])
    self._spacer = Spacer(self.settingsWidget, 15, 5
                         , QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed
                         , grid=(0,2), gridSpan=(1,1))
    self.radioButton1 = RadioButton(self.settingsWidget, grid=(0, 3), hAlign='r', callback=self._createBackBoneButtons)
    self.radioButton1.setChecked(True)
    self.label1 = Label(self.settingsWidget, 'Backbone', grid=(0, 4), hAlign='l')
    self.radioButton2 = RadioButton(self.settingsWidget, grid=(0, 5), hAlign='r', callback=self._createSideChainButtons)
    self.label2 = Label(self.settingsWidget, 'Side chain', grid=(0, 6), hAlign='l')
    self._spacer = Spacer(self.settingsWidget, 5, 5
                         , QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding
                         , grid=(3,12), gridSpan=(1,1))

    # modifiers for sidechain
    self.offsetLabel = Label(self.settingsWidget, 'Display for:   offset =', grid=(1,0))
    self.offsetSelector = PulldownList(self.settingsWidget, grid=(1, 1), texts = ['0', '-1', '+1'],
                                       callback= self._offsetPullDownCallback)
    self.hCheckBox = CheckBox(self.settingsWidget, checked=True, grid=(1,2), hAlign='r', callback=self._toggleBox)
    self.hLabel = Label(self.settingsWidget, 'H', grid=(1,3), hAlign='l')
    self.cCheckBox = CheckBox(self.settingsWidget, checked=True, grid=(1,4), hAlign='r', callback=self._toggleBox)
    self.cLabel = Label(self.settingsWidget, 'C', grid=(1,5), hAlign='l')
    self.nCheckBox = CheckBox(self.settingsWidget, checked=True, grid=(1, 6), hAlign='r', callback=self._toggleBox)
    self.nLabel = Label(self.settingsWidget, 'N', grid=(1, 7), hAlign='l')
    self.otherCheckBox = CheckBox(self.settingsWidget, checked=False, grid=(1, 8), hAlign='r', callback=self._toggleBox)
    self.otherLabel = Label(self.settingsWidget, 'Other', grid=(1, 9), hAlign='l')
    # just store all these widget to be able to toggle them
    self._sidechainModifiers = [self.offsetLabel, self.offsetSelector,
                                self.hCheckBox, self.hLabel,
                                self.cCheckBox, self.cLabel,
                                self.nCheckBox, self.nLabel,
                                self.otherCheckBox, self.otherLabel]
    for w in self._sidechainModifiers:  w.hide()

    # Main widget
    gridLine = -1
    # gridline 0
    gridLine += 1
    nmrResidueLabel = Label(self.mainWidget, 'Current NmrResidue:', grid=(gridLine, 0))
    self.currentNmrResidueLabel = Label(self.mainWidget, grid=(gridLine, 1))

    # gridLine 1
    gridLine += 1
    self.pickAndAssignWidget = Widget(self.mainWidget, setLayout=True, grid=(gridLine, 0), gridSpan=(11,10), vAlign='top')

    self.buttons = {}

    self._updateWidget()
    self._predictAssignments(self.current.peaks)

  def _closeModule(self):
    self.current.unRegisterNotify(self._predictAssignments, 'peaks')
    self.current.unRegisterNotify(self._nmrResidueCallBack, 'nmrResidues')
    super(AtomSelectorModule, self)._closeModule()

  def _nmrResidueCallBack(self, nmrResidues=None):
    "Callback if current.nmrResidue changes"
    if nmrResidues is not None and self.current.nmrResidue:
      self._updateWidget()
      self._predictAssignments(self.current.peaks)

  def _offsetPullDownCallback(self, tmp=None):
    "Callback if offset pullDown changes"
    if self.current.nmrResidue:
      self._updateWidget()
      self._predictAssignments(self.current.peaks)

  def _updateWidget(self):
    "Update the widget to reflect the proper state"
    if self.current.nmrResidue is not None:
      self.currentNmrResidueLabel.setText(self.current.nmrResidue.id)
    else:
      self.currentNmrResidueLabel.setText('<not-defined>')

    if self.radioButton1.isChecked():
      for w in self._sidechainModifiers: w.hide()
      self._createBackBoneButtons()
    elif self.radioButton2.isChecked():
      for w in self._sidechainModifiers: w.show()
      self._createSideChainButtons()
    return

  def _createBackBoneButtons(self):
    self._cleanupPickAndAssignWidget()
    for w in self._sidechainModifiers: w.hide()

    # wb104 27 Jun 2017: changed _cleanupPickAndAssignWidget so removes widgets in reverse order
    # not sure if there was anything else leading to this cludge (and one below) though
    # cludge: don't know why I have to do this for the button to appear: TODO: fix this
    ###_Label = Label(self, text='')
    ###self.pickAndAssignWidget.layout().addWidget(_Label, 0, 0)
    self.buttons = {}
    atoms = ['H', 'N', 'CA', 'CB', 'CO', 'HA', 'HB']
    for ii, atom in enumerate(atoms):
      self.buttons[atom] = []
      for jj, offset in enumerate(['-1', '0', '+1']):
        btext = self.atomLabel(atom, offset)
        button = Button(self.pickAndAssignWidget, text=btext, grid=(ii, jj)
                        , callback=partial(self.assignSelected, offset, atom))
        button.setMinimumSize(45, 24)

        self.buttons[atom].append(button)

    # for buttons in self.buttons.values():
    #   for button in buttons:
    #     button.clicked.connect(self._returnButtonsToNormal)

    self._predictAssignments(self.current.peaks)

  def _createSideChainButtons(self):
    self._cleanupPickAndAssignWidget()
    for w in self._sidechainModifiers: w.show()

    # see comment about cludge above
    # cludge: don't know why I have to do this for the button to appear: TODO: fix this
    ###_label= Label(self.pickAndAssignWidget, '',  hAlign='l')
    ###self.pickAndAssignWidget.layout().addWidget(_label, 0, 0, QtCore.Qt.AlignRight)

    self._updateLayout()
    self._predictAssignments(self.current.peaks)

  def _toggleBox(self):
    if self.radioButton1.isChecked():
      for w in self._sidechainModifiers: w.hide()
    elif self.radioButton2.isChecked():
      for w in self._sidechainModifiers: w.show()
    self._updateLayout()
    self._predictAssignments(self.current.peaks)

  def _getAtomsForButtons(self, atomList, atomName):
    [atomList.remove(atom) for atom in sorted(atomList) if atom[0] == atomName]

  def _getAtomButtonList(self, residueType=None):

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

    if residueType:
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
      residueAtomButtonList = [residueAlphas, residueBetas, residueGammas, residueMoreGammas,
                               residueDeltas, residueMoreDeltas, residueEpsilons,
                               residueMoreEpsilons, residueZetas, residueEtas, residueMoreEtas]
      return residueAtomButtonList

    return atomButtonList

  def _updateLayout(self):

    # group atoms in useful categories based on usage
    atomButtonList = self._getAtomButtonList()

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
        self.buttons = {}
        for ii, atomList in enumerate(atomButtonList):

          for jj, atom in enumerate(atomList):
            self.buttons[atom] = []
            offset = self.offsetSelector.currentText()
            btext = self.atomLabel(atom, offset)
            button = Button(self.pickAndAssignWidget, text=btext, grid=(ii, jj), hAlign='t',
                            callback=partial(self.assignSelected, offset, atom))
            button.setMinimumSize(45, 24)

            self.buttons[atom].append(button)

      else:
        self.buttons = {}
        if self.offsetSelector.currentText() == '-1':
          nmrResidue = self.current.nmrResidue.previousNmrResidue
        elif self.offsetSelector.currentText() == '+1':
          nmrResidue = self.current.nmrResidue.nextNmrResidue
        else:
          nmrResidue = self.current.nmrResidue
        residueType = nmrResidue.residueType.upper()
        atomButtonList2 = self._getAtomButtonList(residueType)
        for ii, atomList in enumerate(atomButtonList2):
          for jj, atom in enumerate(atomList):
            self.buttons[atom] = []
            button = Button(self.pickAndAssignWidget, text=atom, grid=(ii, jj), hAlign='t',
                    callback=partial(self.assignSelected, self.offsetSelector.currentText(), atom))
            button.setMinimumSize(45, 24)
            self.buttons[atom].append(button)

  def _showMoreAtomButtons(self, buttons, moreButton):
    if moreButton.isChecked():
      [button.show() for button in buttons]
    else:
      [button.hide() for button in buttons]

  def _cleanupPickAndAssignWidget(self):

    layout = self.pickAndAssignWidget.layout()
    for r in range(layout.rowCount()-1,-1,-1):
      for c in range(layout.columnCount()-1,-1,-1):
        item = layout.itemAtPosition(r, c)
        if item:
          if item.widget():
            item.widget().hide()
        layout.removeItem(item)

  def atomLabel(self, atom, offset):
    return str(atom + ' [i]' if offset == '0' else atom + ' [i' + offset + ']')

  def checkAssignedAtoms(self, nmrResidue, atoms, predictAtoms):
    """
    Check if the i-1, i, i+1 nmrAtoms for the current residue exist
    :param nmrResidue:
    :return foundAtoms - dict containing True for each found nmrAtom:
    """
    foundAtoms = {}
    # atoms = ['H', 'N', 'CA', 'CB', 'CO', 'HA', 'HB']
    for ii, atom in enumerate(atoms):
      for jj, offset in enumerate(['-1', '0', '+1']):
        bText = self.atomLabel(atom, offset)

        if self._checkAssignedAtom(nmrResidue, offset, atom):
          foundAtoms[bText] = True
          if atom in self.buttons.keys():
            for button in self.buttons[atom]:
              if button.getText() == bText:

                # colour the button if the atom exists
                if bText in predictAtoms:
                  score = predictAtoms[bText]
                  
                  if score >= 85:
                    button.setStyleSheet('background-color: mediumseagreen')
                  elif 50 < score < 85:
                    button.setStyleSheet('background-color: lightsalmon')
                  if score < 50:
                    button.setStyleSheet('background-color: mediumvioletred')

                else:
                  button.setStyleSheet('background-color: cornflowerblue')

    return foundAtoms

  def _getNmrResidue(self, nmrChain, sequenceCode: typing.Union[int, str] = None,
                       residueType: str = None) -> typing.Optional[NmrResidue]:
    partialId = '%s.%s.' % (nmrChain.id, str(sequenceCode).translate(Pid.remapSeparators))
    ll = self.project.getObjectsByPartialId(className='NmrResidue', idStartsWith=partialId)
    if ll:
      return ll[0]
    else:
      return None

  def _getCorrectResidue(self, nmrResidue, offset:int, atomType:str):
    name = atomType
    r = None
    if offset == '-1' and '-1' not in self.current.nmrResidue.sequenceCode:
      r = nmrResidue.previousNmrResidue
      if not r:
        r = self._getNmrResidue(nmrResidue.nmrChain, sequenceCode=nmrResidue.sequenceCode + '-1')
    elif offset == '+1' and '+1' not in self.current.nmrResidue.sequenceCode:
      r = nmrResidue.nextNmrResidue
      if not r:
        r = self._getNmrResidue(nmrResidue.nmrChain, sequenceCode=nmrResidue.sequenceCode + '+1')
    else:
      r = nmrResidue

    return r

  def _checkAssignedAtom(self, nmrResidue, offset:int, atomType:str):
    r = self._getCorrectResidue(nmrResidue=nmrResidue, offset=offset, atomType=atomType)
    if r:
      return r.getNmrAtom(atomType.translate(Pid.remapSeparators))
    else:
      return None

  def assignSelected(self, offset:int, atomType:str):
    """
    Takes a position either -1, 0 or +1 and an atom type, fetches an NmrAtom with name corresponding
    to the atom type and the position and assigns it to correct dimension of current.peaks
    """

    if not self.current.nmrResidue:
      return

    # self.project._appBase._startCommandBlock('application.atomSelector.assignSelected(atomType={!r}, offset={})'.format(atomType, offset))
    self.application._startCommandBlock('application.atomSelector.assignSelected(atomType={!r}, offset={})'.format(atomType, offset))
    try:
      name = atomType

      # search for and create the nmrResidue if it doesn't exists
      if offset == '-1' and '-1' not in self.current.nmrResidue.sequenceCode:
        r = self.current.nmrResidue.previousNmrResidue
        if not r:
          r = self.current.nmrResidue.nmrChain.fetchNmrResidue(sequenceCode=self.current.nmrResidue.sequenceCode+'-1')
      elif offset == '+1' and '+1' not in self.current.nmrResidue.sequenceCode:
        r = self.current.nmrResidue.nextNmrResidue
        if not r:
          r = self.current.nmrResidue.nmrChain.fetchNmrResidue(sequenceCode=self.current.nmrResidue.sequenceCode + '+1')
      else:
        r = self.current.nmrResidue

      # check whether the nmrAtom already exists
      newNmrAtom = r.getNmrAtom(name.translate(Pid.remapSeparators))

      if newNmrAtom:
        if showYesNo('Atom Selector', 'nmrAtom %s exists, do you want to delete?' % newNmrAtom):
          newNmrAtom.delete()
      else:
        newNmrAtom = r.fetchNmrAtom(name=name)
        for peak in self.current.peaks:
          for strip in self.project.strips:
            for peakListView in strip.peakListViews:
              if peak in peakListView.peakItems.keys():
                spectrumIndices = peakListView.spectrumView._displayOrderSpectrumDimensionIndices
                index = spectrumIndices[1]
                axisCode = peak.axisCodes[index]

                currentList = list(peak.dimensionNmrAtoms[index])
                if newNmrAtom not in currentList:
                  nmrAtoms = list(peak.dimensionNmrAtoms[index]) + [newNmrAtom]
                  peak.assignDimension(axisCode, nmrAtoms)
                else:
                  getLogger().warning('Error adding new nmrAtom to %s' % self.current.nmrResidue)

      # self.current.peaks = []
      # flag a change on this peak

    except Exception as es:
      showWarning(str(self.windowTitle()), str(es))
    finally:
      # self.project._endCommandEchoBlock()
      self.application._endCommandBlock()
      # self._returnButtonsToNormal()
      self._predictAssignments(self.current.peaks)

  def _returnButtonsToNormal(self):
    """
    Returns all buttons in Atom Selector to original colours and style.
    """
    if self.application.colourScheme == 'dark':
      backgroundColour1 = '#535a83'
      backgroundColour2 = '#e4e15b'
    else:
      backgroundColour1 = '#bd8413'
      backgroundColour2 = '#fdfdfc'

    for buttons in self.buttons.values():
      for button in buttons:
        button.setStyleSheet(
          """Dock QPushButton { background-color: %s }
             Dock QPushButton::hover { background-color: %s}""" % (backgroundColour1, backgroundColour2)
        )

  def _predictAssignments(self, peaks:typing.List[Peak]):
    """
    Predicts atom type for selected peaks and highlights the relevant buttons with confidence of
    that assignment prediction, green is very confident, orange is less confident.
    """
    self._returnButtonsToNormal()
    if self.current.nmrResidue is None or len(peaks) == 0:
      return

    # check if peaks coincide
    for dim in range(peaks[0].peakList.spectrum.dimensionCount):
      if not peaksAreOnLine(peaks, dim):
        logger.debug('dimension %s: peaksAreonLine=False' % dim)
        return

    types = set(peak.peakList.spectrum.experimentType for peak in peaks)
    anyInterOnlyExperiments = any(isInterOnlyExpt(x) for x in types)

    logger.debug('peaks=%s' % (peaks,))
    logger.debug('types=%s, anyInterOnlyExperiments=%s' % (types, anyInterOnlyExperiments))

    peak = peaks[0]
    peakListViews = [peakListView for peakListView in self.project.peakListViews if peakListView.peakList == peak.peakList]

    #TODO: AtomSelector crashes if there is nothing in the view
    if peakListViews:
      spectrumIndices = peakListViews[0].spectrumView._displayOrderSpectrumDimensionIndices
      isotopeCode = peak.peakList.spectrum.isotopeCodes[spectrumIndices[1]]

      # backbone
      if self.radioButton1.isChecked():
        predictedAtomTypes = [getNmrAtomPrediction(ccpCode, peak.position[spectrumIndices[1]], isotopeCode, strict=True)
                              for ccpCode in CCP_CODES]
        refinedPreds = [(type[0][0][1], type[0][1]) for type in predictedAtomTypes if len(type) > 0]
        atomPredictions = set()
        for atomPred, score in refinedPreds:
          if score > 90:
            atomPredictions.add(atomPred)

        # list containing those atoms that exist - used for colouring in 'checkAssignedAtoms'
        foundPredictList = {}
        for atomPred in atomPredictions:
          if atomPred == 'CB' and self.buttons['CB']:
            if anyInterOnlyExperiments:
              self.buttons['CB'][0].setStyleSheet('background-color: green')
              foundPredictList[self.atomLabel('CB', '-1')] = 100
            else:
              self.buttons['CB'][0].setStyleSheet('background-color: green')
              self.buttons['CB'][1].setStyleSheet('background-color: green')
              foundPredictList[self.atomLabel('CB', '-1')] = 100
              foundPredictList[self.atomLabel('CB', '0')] = 100
          if atomPred == 'CA' and self.buttons['CA']:
            if anyInterOnlyExperiments:
              self.buttons['CA'][0].setStyleSheet('background-color: green')
              foundPredictList[self.atomLabel('CA', '-1')] = 100
            else:
              self.buttons['CA'][0].setStyleSheet('background-color: green')
              self.buttons['CA'][1].setStyleSheet('background-color: green')
              foundPredictList[self.atomLabel('CA', '-1')] = 100
              foundPredictList[self.atomLabel('CA', '0')] = 100

        # new routine to colour any existing atoms
        foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, ['H', 'N', 'CA', 'CB', 'CO', 'HA', 'HB']
                                             , foundPredictList)

      # sidechain is checked
      elif self.radioButton2.isChecked():
        foundPredictList = {}

        if self.current.nmrResidue.residueType == '':
          # In this case, we loop over all CCP_CODES (i.e. residue types)
          predictedAtomTypes = []
          for residueType in CCP_CODES:
            for type, score in getNmrAtomPrediction(residueType, peak.position[spectrumIndices[1]], isotopeCode):
              if len(type) > 0 and score>50:
                predictedAtomTypes.append((type,score))

        else:
          if self.offsetSelector.currentText() == '-1':
            nmrResidue = self.current.nmrResidue.previousNmrResidue
          elif self.offsetSelector.currentText() == '+1':
            nmrResidue = self.current.nmrResidue.nextNmrResidue
          else:
            nmrResidue = self.current.nmrResidue
          predictedAtomTypes = getNmrAtomPrediction(nmrResidue.residueType.title(), peak.position[spectrumIndices[1]], isotopeCode)

        print('>predictAtomTypes>', predictedAtomTypes)
        # find the maximum of each atomType
        predictedDict = {}
        for type, score in predictedAtomTypes:
          if type[1] not in predictedDict:
            predictedDict[type[1]] = (type[0], score)
          else:
            if score > predictedDict[type[1]][1]:
              predictedDict[type[1]] = (type[0], score)
        print ('>>>predictedDict', predictedDict)

        for atomDictType in predictedDict.keys():
          bText = self.atomLabel(atomDictType, '0')
          for atomType, buttons in self.buttons.items():      # get the correct button list
            if atomDictType == atomType:
              for button in buttons:
                if bText == button.getText():
                  print('>type[1], atomType, button>', atomDictType, bText)
                  foundPredictList[self.atomLabel(atomDictType, '0')] = score

                  if score >= 85:
                    button.setStyleSheet('background-color: green')
                  elif 50 < score < 85:
                    button.setStyleSheet('background-color: orange')
                  if score < 50:
                    button.setStyleSheet('background-color: red')

        # new routine to colour any existing atoms
        atomButtonList = self._getAtomButtonList()
        atomButtonList = [x for i in atomButtonList for x in i]
        foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, atomButtonList
                                             , foundPredictList)

        #
        #
        # for type, score in predictedAtomTypes:
        #   #print('>type, score>', type, score)
        #   for atomType, buttons in self.buttons.items():
        #     if type[1] == atomType:
        #       for button in buttons:
        #         print('>type[1], atomType, button>', type[1], atomType,  button.getText())
        #         if score >= 85:
        #           button.setStyleSheet('background-color: green')
        #         elif 50 < score < 85:
        #           button.setStyleSheet('background-color: orange')
        #         if score < 50:
        #           button.setStyleSheet('background-color: red')
        #
        #
        #
