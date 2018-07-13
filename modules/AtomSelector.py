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
__version__ = "$Revision: 3.0.b3 $"
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
import copy
from functools import partial
from collections import OrderedDict
from PyQt5 import QtCore, QtGui, QtWidgets

from ccpn.core.Peak import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.lib import Pid
from ccpn.core.lib.AssignmentLib import isInterOnlyExpt, getNmrAtomPrediction, CCP_CODES
from ccpn.core.lib.AssignmentLib import peaksAreOnLine
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.RadioButton import RadioButton
from ccpn.ui.gui.widgets.RadioButtons import RadioButtons
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.MessageDialog import showYesNo
from ccpnmodel.ccpncore.lib.assignment.ChemicalShift import PROTEIN_ATOM_NAMES, ALL_ATOMS_SORTED
from ccpn.util.Logging import getLogger
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.lib.AssignmentLib import _assignNmrAtomsToPeaks


logger = getLogger()

# TODO:ED Add DNA, RNA structures to the list
# MOLECULE_TYPES = ['protein', 'DNA', 'RNA', 'carbohydrate', 'other']
MOLECULE_TYPES = ['protein']
ATOM_TYPES = ['H', 'N', 'CA', 'CB', 'CO', 'HA', 'HB']


class AtomSelectorModule(CcpnModule):
  """
  Module to be used with PickAndAssignModule for prediction of nmrAtom names and assignment of nmrAtoms
  to peak dimensions
  Responds to current.nmrResidue and current.peaks
  """
  className = 'AtomSelectorModule'

  includeSettingsWidget = True
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsPosition = 'top'

  def __init__(self, mainWindow=None, name='Atom Selector', nmrAtom=None):
    CcpnModule.__init__(self, mainWindow=mainWindow, name=name)

    # Derive application, project, and current from mainWindow
    self.mainWindow = mainWindow
    if mainWindow:
      self.application = mainWindow.application
      self.project = mainWindow.application.project
      self.current = mainWindow.application.current
      self._registerNotifiers()

    # Settings Widget
    self.molTypeLabel = Label(self.settingsWidget, 'Molecule Type', grid=(0, 0))
    self.molTypePulldown = PulldownList(self.settingsWidget, grid=(0, 1), texts=MOLECULE_TYPES , callback=self._changeMoleculeType)

    self.modeTypeLabel = Label(self.settingsWidget, 'Mode', grid=(1, 0))
    self.modeRadioButtons = RadioButtons(self.settingsWidget, texts=['Backbone','Side chain'], selectedInd=0, callback=self._createButtonsCallback, grid=(1,1))
    self.radioButton1, self.radioButton2 = self.modeRadioButtons.radioButtons

    # modifiers for sidechain
    self.offsetLabel = Label(self.settingsWidget, 'Offset', grid=(2,0))
    self.offsetSelector = PulldownList(self.settingsWidget, grid=(2, 1), texts = ['0', '-1', '+1'], callback= self._offsetPullDownCallback)

    self.atomTypeLabel = Label(self.settingsWidget, 'Atom Type', grid=(3, 0))
    self.atomOptions = RadioButtons(self.settingsWidget,selectedInd=1, texts=['H','C','N', 'Other'],callback=self._toggleBox, grid=(3,1))
    self.hCheckBox, self.cCheckBox, self.nCheckBox, self.otherCheckBox  = self.atomOptions.radioButtons

    self._sidechainModifiers = [self.offsetLabel, self.offsetSelector]

    self.settingsWidget.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
    self.settingsWidget.setContentsMargins(10,10,10,10)
    self.mainWidget.setContentsMargins(10, 10, 10, 10)

    for w in self._sidechainModifiers:  w.hide()

    # Main widget
    gridLine = 0
    self._residueFrame = Frame(self.mainWidget, setLayout=True, grid=(gridLine, 0), gridSpan=(1,1))
    self._nmrResidueLabel = Label(self._residueFrame, 'Current NmrResidue:', grid=(0, 0)
                                  , hPolicy='minimal')
    self.currentNmrResidueLabel = Label(self._residueFrame, grid=(0, 1), gridSpan=(1, 3)
                                        , hPolicy='minimalexpanding', hAlign='l')
    self._residueFrame.setFixedHeight(25)
    self.mainWidget.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Expanding)

    # gridLine 1
    gridLine += 1
    self.pickAndAssignWidget = Widget(self.mainWidget, setLayout=True, grid=(gridLine, 0), gridSpan=(1,1), vAlign='top')
    Spacer(self.mainWidget, 5, 5
           , QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding
           , grid=(gridLine, 0), gridSpan=(1,1))

    self.buttons = {}

    self._updateWidget()
    # self._predictAssignments(self.current.peaks)

  def _registerNotifiers(self):
    self._nmrAtomNotifier = Notifier(self.project
                                    , [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE]
                                    , NmrAtom.__name__
                                    , self._nmrResidueCallBack)
    self._peakNotifier = Notifier(self.current
                                 , [Notifier.CURRENT]
                                 , Peak._pluralLinkName
                                 , self._predictAssignmentsCallBack)
    self._nmrResidueNotifier = Notifier(self.current
                                 , [Notifier.CURRENT]
                                 , NmrResidue._pluralLinkName
                                 , self._nmrResidueCallBack)

  def _unRegisterNotifiers(self):
    """
    clean up the notifiers
    """
    if self._peakNotifier is not None:
      self._peakNotifier.unRegister()
    if self._nmrResidueNotifier is not None:
      self._nmrResidueNotifier.unRegister()
    if self._nmrAtomNotifier is not None:
      self._nmrAtomNotifier.unRegister()

  def _closeModule(self):
    self._unRegisterNotifiers()
    super(AtomSelectorModule, self)._closeModule()

  def _createButtonsCallback(self):
    if self.radioButton1.isChecked():
      self._createBackBoneButtons()
    if self.radioButton2.isChecked():
      self._createSideChainButtons()

  def _nmrResidueCallBack(self, nmrResidues=None):
    "Callback if current.nmrResidue changes"
    if nmrResidues is not None and self.current.nmrResidue:
      self._updateWidget()
      if self.current.peaks:
        # self._predictAssignments(self.current.peaks)
        self.pickAndAssignWidget.show()
      else:
        self.pickAndAssignWidget.hide()
    else:
      self.pickAndAssignWidget.hide()

  def _offsetPullDownCallback(self, tmp=None):
    "Callback if offset pullDown changes"
    if self.current.nmrResidue:
      self._updateWidget()
      self._predictAssignments(self.current.peaks)

  def _updateWidget(self):
    "Update the widget to reflect the proper state"
    try:
      if self.current.nmrResidue is not None:
        self.currentNmrResidueLabel.setText(self.current.nmrResidue.id)
        if self.radioButton1.isChecked():
          for w in self._sidechainModifiers:
            w.hide()
          self._createBackBoneButtons()
        elif self.radioButton2.isChecked():
          for w in self._sidechainModifiers:
            w.show()
          self._createSideChainButtons()
      else:
        self.currentNmrResidueLabel.setText('<not-defined>')
      self._setCheckedButtonOfAssignedAtoms(self.current.nmrResidue)
      return
    except:
      return

  def _removeOffsetFromButtonText(self, text:str):
    p = text.split(' ')
    if len(p)>0: return p[0]
    else: return text

  def _setCheckedButtonOfAssignedAtoms(self, nmrResidue, offset=None):
    '''setChecked the radioButton Of Assigned Nmr Atoms. "Check" as tick the button '''
    from ccpn.util.Common import makeIterableList

    if len(self.current.peaks) > 0:
      peaks = self.current.peaks
    else:
      return
    if not nmrResidue:
      return
    allButtons = [b for key in self.buttons.keys() for b in self.buttons[key]]

    for peak in peaks:
      buttonsToCheck = set()
      for assignedNmrAtom in  makeIterableList(peak.assignedNmrAtoms):
        if assignedNmrAtom in nmrResidue.nmrAtoms:
          for button in allButtons:
            if assignedNmrAtom:
              if offset =='0' or offset is None:
                if assignedNmrAtom.name == button.getText():
                  buttonsToCheck.add(button)
              else:
                btext = self.atomLabel(assignedNmrAtom.name, offset)
                if btext ==  button.getText():
                  buttonsToCheck.add(button)
      for b in allButtons:
        if b in list(buttonsToCheck):
          b.setChecked(True)
        else:
          b.setChecked(False)





  def _createBackBoneButtons(self):
    self._cleanupPickAndAssignWidget()
    for w in self._sidechainModifiers: w.hide()

    # wb104 27 Jun 2017: changed _cleanupPickAndAssignWidget so removes widgets in reverse order
    # not sure if there was anything else leading to this cludge (and one below) though
    # cludge: don't know why I have to do this for the button to appear: TODO: fix this
    ###_Label = Label(self, text='')
    ###self.pickAndAssignWidget.layout().addWidget(_Label, 0, 0)
    self.buttons = {}
    atoms = ATOM_TYPES

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
      for ii, atom in enumerate(atoms):
        self.buttons[atom] = []

        # skip of startswith these atomTypes
        if not self.cCheckBox.isChecked() and atom.startswith('C'):
          continue
        if not self.hCheckBox.isChecked() and atom.startswith('H'):
          continue
        if not self.nCheckBox.isChecked() and atom.startswith('N'):
          continue
        if not self.otherCheckBox.isChecked() and not atom.startswith('C')\
                                              and not atom.startswith('H') \
                                              and not atom.startswith('N'):
          continue

        for jj, offset in enumerate(['-1', '0', '+1']):
          btext = self.atomLabel(atom, offset)
          button = RadioButton(self.pickAndAssignWidget, text=btext, grid=(ii, jj)
                          , callback=partial(self.assignSelected, offset, atom))
          button.setMinimumSize(45, 24)
          self.buttons[atom].append(button)

    self._predictAssignments(self.current.peaks)

  def _createSideChainButtons(self):
    self._cleanupPickAndAssignWidget()
    for w in self._sidechainModifiers: w.show()

    # see comment about cludge above
    # cludge: don't know why I have to do this for the button to appear: TODO: fix this
    ###_label= Label(self.pickAndAssignWidget, '',  hAlign='l')
    ###self.pickAndAssignWidget.layout().addWidget(_label, 0, 0, QtCore.Qt.AlignRight)

    self._updateChainLayout()
    self._predictAssignments(self.current.peaks)

  def _toggleBox(self):
    if self.radioButton1.isChecked():
      for w in self._sidechainModifiers: w.hide()
    elif self.radioButton2.isChecked():
      for w in self._sidechainModifiers: w.show()
    self._updateWidget()
    # self._predictAssignments(self.current.peaks)

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

  def _getDnaRnaButtonList(self, atomList=None, residueType=None):
    residueAtomButtonList = copy.deepcopy(ALL_DNARNA_ATOMS_SORTED)

    if residueType and atomList:
      residueAtoms = atomList[residueType]

      for atomType in ALL_DNARNA_ATOMS_SORTED.keys():
        atomTypeList = ALL_DNARNA_ATOMS_SORTED[atomType]
        for atom in atomTypeList:
          if atom not in residueAtoms:
            residueAtomButtonList[atomType].remove(atom)

    return [residueAtomButtonList[atom] for atom in residueAtomButtonList.keys()]

  def _updateChainLayout(self):

    # group atoms in useful categories based on usage
    atomButtonList = self._getAtomButtonList()

    # testing DNA/RNA buttonlist
    # atomButtonList = self._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DT')

    # Activate button for Carbons
    if not self.cCheckBox.isChecked():
      [self._getAtomsForButtons(atomList, 'C') for atomList in atomButtonList]

    if not self.hCheckBox.isChecked():
      [self._getAtomsForButtons(atomList, 'H') for atomList in atomButtonList]

    if not self.nCheckBox.isChecked():
      [self._getAtomsForButtons(atomList, 'N') for atomList in atomButtonList]

    if not self.otherCheckBox.isChecked():
      for atomList in atomButtonList:
        [atomList.remove(atom) for atom in sorted(atomList) if not atom.startswith('C') \
        and not atom.startswith('H') \
        and not atom.startswith('N')]

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
            button = RadioButton(self.pickAndAssignWidget, text=btext, grid=(ii, jj), hAlign='t',
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
            button = RadioButton(self.pickAndAssignWidget, text=atom, grid=(ii, jj), hAlign='t',
                    callback=partial(self.assignSelected, self.offsetSelector.currentText(), atom))
            # button = Button(self.pickAndAssignWidget, text=atom, grid=(ii, jj), hAlign='t',
            #         callback=partial(self.assignSelected, self.offsetSelector.currentText(), atom))
            button.setMinimumSize(45, 24)
            # button.setAutoExclusive(True)

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

  def atomLabel(self, atom, offset, showAll=False):
    if showAll:
      return str(atom + ' [i]' if offset == '0' else atom + ' [i' + offset + ']')
    else:
      return str(atom if offset == '0' else atom + ' [i' + offset + ']')

  def checkAssignedAtoms(self, nmrResidue, atoms, predictAtoms, checkMode='backbone'):
    """
    Check if the i-1, i, i+1 nmrAtoms for the current residue exist
    :param nmrResidue:
    :return foundAtoms - dict containing True for each found nmrAtom:
    """
    foundAtoms = {}
    if checkMode == 'backbone':
      # all residues are displayed so use the central mainNmrResidue
      nmrResidue = self._getMainNmrResidue(nmrResidue)

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

                # else: # Users don't need to know/ or be notified here that the nmrAtom exist in the selected nmrResidue.
                #   button.setStyleSheet('background-color: cornflowerblue')


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
    if offset == '-1' and '-1' not in nmrResidue.sequenceCode:
      r = nmrResidue.previousNmrResidue
      if not r:
        r = self._getNmrResidue(nmrResidue.nmrChain, sequenceCode=nmrResidue.sequenceCode + '-1')
    elif offset == '+1' and '+1' not in nmrResidue.sequenceCode:
      r = nmrResidue.nextNmrResidue
      if not r:
        r = self._getNmrResidue(nmrResidue.nmrChain, sequenceCode=nmrResidue.sequenceCode + '+1')
    else:
      r = nmrResidue

    return r

  def _checkAssignedAtom(self, nmrResidue, offset:int, atomType:str):
    '''
    This checks only if an NMRAtom exists not if is assigned to a peak!
    :param nmrResidue:
    :param offset:
    :param atomType:
    :return:
    '''
    r = self._getCorrectResidue(nmrResidue=nmrResidue, offset=offset, atomType=atomType)
    if r:
      atom = r.getNmrAtom(atomType.translate(Pid.remapSeparators))
      if atom and not atom.isDeleted:
        return r
      else:
        return None
    else:
      return None

  def _getMainNmrResidue(self, nmrResidue):
    if nmrResidue:
      if nmrResidue.relativeOffset and nmrResidue.relativeOffset != 0:
        return nmrResidue.mainNmrResidue

    return nmrResidue

  def _assignDimension(self):
    """
    update the peak assignment and create a event to update the module
    """
    pass

  def deassignAtomFromSelectedPeaks(self, peaks, nmrAtom):

    if not peaks: return
    if not nmrAtom: return

    newAssignedAtoms = ()
    for peak in peaks:
      for subTuple in peak.assignedNmrAtoms:
        a = tuple(None if na is nmrAtom else na for na in subTuple)
        newAssignedAtoms += (a,)

      peak.assignedNmrAtoms = newAssignedAtoms



  def assignSelected(self, offset:int, atomType:str):
    """
    Takes a position either -1, 0 or +1 and an atom type, fetches an NmrAtom with name corresponding
    to the atom type and the position and assigns it to correct dimension of current.peaks
    """

    if not self.current.nmrResidue or not self.current.peaks:
      return

    assignResidue = self._getMainNmrResidue(self.current.nmrResidue)

    self.application._startCommandBlock('application.atomSelector.assignSelected(atomType={!r}, offset={})'.format(atomType, offset))
    try:
      name = atomType

      # search for and create the nmrResidue if it doesn't exists
      if offset == '-1' and '-1' not in assignResidue.sequenceCode:
        r = assignResidue.previousNmrResidue
        if not r:
          r = assignResidue.nmrChain.fetchNmrResidue(sequenceCode=assignResidue.sequenceCode+'-1')
      elif offset == '+1' and '+1' not in assignResidue.sequenceCode:
        r = assignResidue.nextNmrResidue
        if not r:
          r = assignResidue.nmrChain.fetchNmrResidue(sequenceCode=assignResidue.sequenceCode + '+1')
      else:
        r = assignResidue

      button = self.sender() #This mechanisms allows to uncheck a radiobutton and still keep the autoexclusionas
      button.setAutoExclusive(False)
      if not button.isChecked():
        # means we are unchecking. Therefore Needs to deassign that atom to the selected peak:
         nmrAtom = r.getNmrAtom(name.translate(Pid.remapSeparators))
         if nmrAtom:
           self.deassignAtomFromSelectedPeaks(self.current.peaks, nmrAtom)
           self._setCheckedButtonOfAssignedAtoms(r, offset=offset)
           return

      else:
        nmrAtom = r.fetchNmrAtom(name=name)
        _assignNmrAtomsToPeaks(strip=self.current.strip,
                               nmrAtoms=[nmrAtom], peaks=self.current.peaks)
        self._setCheckedButtonOfAssignedAtoms(r, offset=offset)

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

    for buttons in self.buttons.values():
      for button in buttons:
        button.setStyleSheet(
          """Dock QRadioButton { background-color: %s }
             Dock QRadioButton::hover { background-color: %s}""" % ('lightgrey', 'white'))

  def _predictAssignmentsCallBack(self, data):
    peaks = data[Notifier.VALUE]

    self._predictAssignments(peaks)
    self._setCheckedButtonOfAssignedAtoms(self.current.nmrResidue)

  def _predictAssignments(self, peaks:typing.List[Peak]):
    """
    Predicts atom type for selected peaks and highlights the relevant buttons with confidence of
    that assignment prediction, green is very confident, orange is less confident.
    """
    self._returnButtonsToNormal()
    if self.current.nmrResidue is None or len(peaks) == 0:
      self.pickAndAssignWidget.hide()
      return

    # make sure that the widget is visible
    self.pickAndAssignWidget.show()

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
              self.buttons['CB'][0].setStyleSheet('background-color: mediumseagreen')
              foundPredictList[self.atomLabel('CB', '-1')] = 100
            else:
              self.buttons['CB'][0].setStyleSheet('background-color: mediumseagreen')
              self.buttons['CB'][1].setStyleSheet('background-color: mediumseagreen')
              foundPredictList[self.atomLabel('CB', '-1')] = 100
              foundPredictList[self.atomLabel('CB', '0')] = 100
          if atomPred == 'CA' and self.buttons['CA']:
            if anyInterOnlyExperiments:
              self.buttons['CA'][0].setStyleSheet('background-color: mediumseagreen')
              foundPredictList[self.atomLabel('CA', '-1')] = 100
            else:
              self.buttons['CA'][0].setStyleSheet('background-color: mediumseagreen')
              self.buttons['CA'][1].setStyleSheet('background-color: mediumseagreen')
              foundPredictList[self.atomLabel('CA', '-1')] = 100
              foundPredictList[self.atomLabel('CA', '0')] = 100

        # new routine to colour any existing atoms
        # foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, ATOM_TYPES
        #                                      , foundPredictList, 'backbone')

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

        # print('>predictAtomTypes>', predictedAtomTypes)
        # find the maximum of each atomType
        predictedDict = {}
        for type, score in predictedAtomTypes:
          if type[1] not in predictedDict:
            predictedDict[type[1]] = (type[0], score)
          else:
            if score > predictedDict[type[1]][1]:
              predictedDict[type[1]] = (type[0], score)
        # print ('>>>predictedDict', predictedDict)

        for atomDictType in predictedDict.keys():
          bText = self.atomLabel(atomDictType, '0')
          for atomType, buttons in self.buttons.items():      # get the correct button list
            if atomDictType == atomType:
              for button in buttons:
                if bText == button.getText():
                  # print('>type[1], atomType, button>', atomDictType, bText)
                  foundPredictList[self.atomLabel(atomDictType, '0')] = score

                  if score >= 85:
                    button.setStyleSheet('background-color: green')
                  elif 50 < score < 85:
                    button.setStyleSheet('background-color: orange')
                  if score < 50:
                    button.setStyleSheet('background-color: red')

        # new routine to colour any existing atoms
        # atomButtonList = self._getAtomButtonList()
        # atomButtonList = [x for i in atomButtonList for x in i]
        # foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, atomButtonList
        #                                      , foundPredictList, 'sideChain')

  def _changeMoleculeType(self, data):
    """
    change the  available atomList depending on the moleculeType
    :param data - str from pullDown:
    """
    pass

  def getResidueTypes(self, moleculeType:str='protein'):
    """
    return a list of residue types assiciated with the moleculeType
    :param moleculeType - str ['protein', 'DNA', 'RNA', 'carbohydrate', 'other']
    :return list of str:
    """
    if moleculeType in MOLECULE_TYPES:
      if moleculeType == 'protein':
        return [atomName for atomName in PROTEIN_ATOM_NAMES.keys()]
    else:
      return None


DNA_ATOMS = """
Res    Name     Atom     Count        Min.        Max.       Avg.     Std Dev

DA          H2       H        950          2.60         8.64       7.49       0.80       
DA          H61      H        145          2.56        13.40       6.83       1.16       
DA          H62      H        144          5.31        12.90       6.92       1.11       
DA          H8       H        1097         2.35         9.10       8.00       0.83       
DA          H1'      H        1092         4.19         8.50       6.04       0.33       
DA          H2'      H        1083         0.84         4.82       2.66       0.37       
DA          H2''     H        1065         0.56         6.89       2.79       0.37       
DA          H3'      H        1031         4.00         6.04       4.94       0.19       
DA          H4'      H        819          2.08         9.13       4.34       0.34       
DA          H5'      H        449          2.79         5.95       4.06       0.30       
DA          H5''     H        406          2.96         8.41       4.07       0.45       
DA          C2       C        134        151.00       159.33     154.69       1.25       
DA          C4       C        2          150.30       150.60     150.45       0.21       
DA          C5       C        2          119.80       119.90     119.85       0.07       
DA          C6       C        2          157.00       157.40     157.20       0.28       
DA          C8       C        145        138.00       142.79     141.33       0.94       
DA          C1'      C        138         81.50        91.90      85.19       1.26       
DA          C2'      C        123         36.23        42.93      40.58       1.05       
DA          C3'      C        99          72.71        83.26      78.46       1.70       
DA          C4'      C        68          84.43        91.33      87.21       1.19       
DA          C5'      C        45          47.77        70.96      67.18       3.42       
DA          N1       N        3          223.40       227.00     225.20       1.80       
DA          N3       N        3          214.50       216.40     215.73       1.07       
DA          N6       N        9           77.40        81.70      80.22       1.29       
DA          N7       N        2          233.50       234.10     233.80       0.42       
DA          N9       N        1          170.20       170.20     170.20       0.00       
DA          P        P        239        -45.29        35.48      -2.32       5.04       

DC          H41      H        593         -0.26        12.08       7.40       1.15       
DC          H42      H        576         -0.26        10.62       7.48       1.07       
DC          H5       H        1263        -6.22         8.25       5.40       0.89       
DC          H6       H        1277       -14.72         8.31       7.36       0.98       
DC          H1'      H        1269         2.17         6.60       5.72       0.69       
DC          H2'      H        1251       -29.45         6.93       2.11       1.42       
DC          H2''     H        1235       -20.82         8.48       2.39       1.08       
DC          H3'      H        1177        -7.37         8.62       4.72       0.65       
DC          H4'      H        923        -38.59        28.54       4.18       2.40       
DC          H5'      H        458        -18.35        11.19       4.17       1.81       
DC          H5''     H        381        -41.21        15.96       3.98       3.50       
DC          C2       C        2          158.80       159.20     159.00       0.28       
DC          C4       C        2          167.80       168.00     167.90       0.14       
DC          C5       C        125         84.55        99.90      98.12       2.25       
DC          C6       C        129        133.73       144.70     142.76       1.63       
DC          C1'      C        134         77.90        88.53      86.94       1.65       
DC          C2'      C        105         30.62        42.27      40.34       1.45       
DC          C3'      C        101         66.41        80.30      76.57       2.75       
DC          C4'      C        52          77.80        88.39      85.90       1.93       
DC          C5'      C        33          59.30        71.72      66.28       2.63       
DC          N1       N        1          150.70       150.70     150.70       0.00       
DC          N3       N        1          196.30       196.30     196.30       0.00       
DC          N4       N        16          95.10       100.41      98.29       1.33       
DC          P        P        332        -14.20        52.72      -1.65       6.05       

DG          H1       H        1299        -2.91        13.94      11.62       2.06       
DG          H21      H        322        -26.86        18.44       7.46       3.00       
DG          H22      H        309        -26.86        18.44       6.36       3.02       
DG          H8       H        1975         6.07         9.31       7.78       0.25       
DG          H1'      H        1979         2.53        12.45       5.90       0.45       
DG          H2'      H        1932         1.43         8.23       2.67       0.41       
DG          H2''     H        1902         1.16        11.50       2.72       0.42       
DG          H3'      H        1863         0.00        13.24       4.97       0.68       
DG          H4'      H        1547         0.00        13.26       4.46       0.98       
DG          H5'      H        885          0.00         7.82       4.12       0.36       
DG          H5''     H        761          0.00         7.76       4.10       0.35       
DG          C2       C        2          156.50       156.70     156.60       0.14       
DG          C4       C        2          153.10       154.00     153.55       0.64       
DG          C5       C        2          117.40       117.60     117.50       0.14       
DG          C6       C        2          161.00       161.40     161.20       0.28       
DG          C8       C        182        134.38       142.93     138.20       1.62       
DG          C1'      C        166         74.90        91.82      85.07       2.02       
DG          C2'      C        128         32.64        49.00      40.23       2.00       
DG          C3'      C        124         67.72        81.41      77.69       2.27       
DG          C4'      C        77          75.27        90.64      86.92       2.19       
DG          C5'      C        57          58.50        71.42      67.06       2.21       
DG          N1       N        33         142.74       147.70     145.74       1.96       
DG          N2       N        5           75.10        76.10      75.52       0.38       
DG          N7       N        4          236.90       238.30     237.28       0.68       
DG          N9       N        1          168.50       168.50     168.50       0.00       
DG          P        P        406        -20.92        47.10      -1.39       6.34       

DT          H3       H        575          2.12        14.53      12.85       1.97       
DT          H6       H        1420         1.61         8.10       7.36       0.32       
DT          H7       H        8            1.46         1.79       1.65       0.13       
DT          H71      H        1249         0.18         7.53       1.60       0.41       
DT          H72      H        1248         0.18         7.53       1.60       0.41       
DT          H73      H        1248         0.18         7.53       1.60       0.41       
DT          H1'      H        1412         1.65         6.72       5.88       0.54       
DT          H2'      H        1402         0.50         6.21       2.16       0.47       
DT          H2''     H        1383         0.78         4.87       2.39       0.36       
DT          H3'      H        1337         3.79        14.11       4.80       0.47       
DT          H4'      H        997          2.00        14.15       4.16       0.55       
DT          H5'      H        608          2.54        14.20       4.20       1.56       
DT          H5''     H        547          0.98         7.19       3.97       0.40       
DT          C2       C        2          152.90       153.40     153.15       0.35       
DT          C4       C        2          168.50       169.10     168.80       0.42       
DT          C5       C        9           14.30       113.60      36.43      43.64       
DT          C6       C        229        135.80       143.60     139.24       0.92       
DT          C7       C        59          11.96        14.90      14.19       0.62       
DT          C1'      C        215         80.60        92.22      86.88       1.37       
DT          C2'      C        195         34.38        42.96      39.99       0.88       
DT          C3'      C        176         71.58        80.10      77.49       1.64       
DT          C4'      C        64          78.80        89.31      86.04       1.46       
DT          C5'      C        41          59.40        68.40      66.36       1.91       
DT          N1       N        1          142.90       142.90     142.90       0.00       
DT          N3       N        41         156.40       160.50     159.27       0.71       
DT          P        P        261         -5.65        47.79      -1.36       7.03       
"""

RNA_ATOMS = """
Res    Name     Atom     Count        Min.        Max.       Avg.     Std Dev

A          H2       H        1617         0.00         9.39       7.63       0.47       
A          H61      H        298          0.00        10.85       7.41       1.18       
A          H62      H        265          0.00         9.98       6.78       1.11       
A          H8       H        1592         7.00         9.20       8.01       0.27       
A          H1'      H        1616         4.56         7.16       5.88       0.21       
A          HO2'     H        44           4.27         7.75       5.75       1.18       
A          H2'      H        1414         3.61         5.91       4.59       0.22       
A          H3'      H        1213         3.91         5.45       4.63       0.21       
A          H4'      H        1042         0.00         5.29       4.46       0.25       
A          H5'      H        803          0.00         5.03       4.30       0.27       
A          H5''     H        786          0.00       153.90       4.37       5.35       
A          C2       C        1015         0.00       167.71     151.71       8.98       
A          C4       C        35           0.00       151.89     131.04      47.78       
A          C5       C        35           0.00       153.63     104.96      39.80       
A          C6       C        37           0.00       159.60     139.05      49.16       
A          C8       C        993          0.00       145.80     138.18       9.80       
A          C1'      C        944         84.23       130.39      91.64       3.57       
A          C2'      C        703          4.75        99.00      75.25       5.81       
A          C3'      C        640          0.00       101.40      73.68       4.97       
A          C4'      C        629          0.00        92.60      82.53       4.09       
A          C5'      C        540          0.00       109.20      66.22       6.79       
A          N1       N        253          0.00       229.70     216.11      30.52       
A          N3       N        203          0.00       231.08     210.51      30.33       
A          N6       N        157          0.00        97.01      79.33      13.30       
A          N7       N        151          0.00       237.10     224.40      37.21       
A          N9       N        189          0.00       176.80     166.48      24.64       
A          P        P        176         -5.96         2.43      -2.97       1.64       

C          H41      H        973          5.64        10.38       7.91       0.72       
C          H42      H        951          5.07         9.06       7.25       0.69       
C          H5       H        1865         4.20         7.09       5.48       0.29       
C          H6       H        1915         5.62         8.46       7.68       0.21       
C          H1'      H        1871         3.42         6.42       5.55       0.24       
C          HO2'     H        66           3.96         9.70       5.82       1.49       
C          H2'      H        1578         0.00         5.07       4.32       0.24       
C          H3'      H        1389         0.00         5.63       4.41       0.23       
C          H4'      H        1135         0.00        12.58       4.33       0.37       
C          H5'      H        850          0.00         7.17       4.26       0.38       
C          H5''     H        871          0.00         5.04       4.10       0.46       
C          C2       C        65           0.00       189.36     136.58      59.06       
C          C4       C        78           0.00       169.52     147.08      53.62       
C          C5       C        1101         0.00       141.05      97.28       4.33       
C          C6       C        1132        93.40       145.02     139.71       7.35       
C          C1'      C        1076        84.00       118.95      92.68       2.10       
C          C2'      C        839          0.00        99.40      75.76       5.32       
C          C3'      C        780          0.00       104.90      73.18       6.41       
C          C4'      C        735          0.00        92.60      81.94       5.86       
C          C5'      C        609          0.00       110.20      64.72       9.47       
C          N1       N        229          0.00       178.83     145.28      30.98       
C          N3       N        181          0.00       204.30     184.53      43.16       
C          N4       N        390          0.00       104.50      96.67      10.21       
C          P        P        186         -5.13         0.62      -3.15       1.65       

G          H1       H        1527         6.09        14.32      12.35       0.98       
G          H21      H        444          0.00         9.39       7.37       1.42       
G          H22      H        402          0.00         9.06       6.31       1.22       
G          H8       H        2293         0.00         8.63       7.60       0.39       
G          H1'      H        2257         3.41         7.62       5.66       0.33       
G          HO2'     H        53           4.40         7.10       5.77       1.16       
G          H2'      H        1931         3.26         6.30       4.58       0.25       
G          H3'      H        1724         3.80         5.78       4.56       0.26       
G          H4'      H        1405         0.00         5.11       4.43       0.22       
G          H5'      H        1187         2.89         5.42       4.28       0.23       
G          H5''     H        1157         2.55         5.11       4.18       0.23       
G          C2       C        59           0.00       162.50     131.48      56.36       
G          C4       C        44           0.00       158.62     120.73      61.96       
G          C5       C        78           0.00       167.24     106.95      40.80       
G          C6       C        86           0.00       162.74     142.40      49.29       
G          C8       C        1429         0.00       146.73     135.38       8.27       
G          C1'      C        1244        79.80        95.05      91.50       1.97       
G          C2'      C        975          0.00        99.40      75.27       4.88       
G          C3'      C        887          0.00       102.20      73.85       5.27       
G          C4'      C        855         72.08        96.10      82.65       2.63       
G          C5'      C        758         50.36       108.20      66.41       6.34       
G          N1       N        891          0.00       166.07     145.58      13.10       
G          N2       N        151          0.00        83.38      71.00      17.00       
G          N3       N        25           0.00       234.10      97.88      80.79       
G          N7       N        184          0.00       240.99     221.87      50.58       
G          N9       N        280          0.00       176.40     164.01      30.01       
G          P        P        208         -6.00         0.56      -2.86       1.59       

U          H3       H        898          1.20       154.54      13.24       4.89       
U          H5       H        1524         4.20         7.83       5.47       0.31       
U          H6       H        1590         5.81         8.52       7.75       0.20       
U          H1'      H        1559         3.69         6.59       5.61       0.26       
U          HO2'     H        45           4.08         9.74       6.06       1.51       
U          H2'      H        1357         0.00         6.63       4.37       0.27       
U          H3'      H        1117         0.00         7.83       4.48       0.25       
U          H4'      H        949          0.00         4.82       4.36       0.23       
U          H5'      H        697          0.00         4.85       4.24       0.30       
U          H5''     H        713          0.00         4.86       4.15       0.29       
U          C2       C        80           0.00       183.20     144.31      38.05       
U          C4       C        86           0.00       169.70     157.82      39.47       
U          C5       C        944          0.00       154.13     102.90       6.62       
U          C6       C        975          0.00       169.95     140.14       8.53       
U          C1'      C        990         80.62        96.30      91.90       2.14       
U          C2'      C        737          0.00        99.80      75.04       4.49       
U          C3'      C        656          0.00       102.50      73.74       5.43       
U          C4'      C        643          0.00        94.42      82.70       4.14       
U          C5'      C        525          0.00       106.60      65.07       5.67       
U          N1       N        209          0.00       192.14     147.96      25.60       
U          N3       N        553          0.00       167.30     158.99      16.67       
U          O4       O        5            0.00         0.00       0.00       0.00       
U          P        P        149         -5.30         1.58      -2.99       1.78       
"""

DNA_ATOM_NAMES = {
  'DT': ['H3', 'H6', 'H7', 'H71', 'H72', 'H73', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
         'C2', 'C4', 'C5', 'C6', 'C7', "C1'", "C2'", "C3'", "C4'", "C5'",
         'N1', 'N3',
         'P'],
  'DC': ['H41', 'H42', 'H5', 'H6', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
         'C2', 'C4', 'C5', 'C6', "C1'", "C2'", "C3'", "C4'", "C5'",
         'N1', 'N3', 'N4',
         'P'],
  'DA': ['H2', 'H61', 'H62', 'H8', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
         'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
         'N1', 'N3', 'N6', 'N7', 'N9', 'P'],
  'DG': ['H1', 'H21', 'H22', 'H8', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
         'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
         'N1', 'N2', 'N7', 'N9',
         'P']
}
RNA_ATOM_NAMES = {
  'G': ['H1', 'H21', 'H22', 'H8', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
        'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
        'N1', 'N2', 'N3', 'N7', 'N9',
        'P'],
  'U': ['H3', 'H5', 'H6', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
        'C2', 'C4', 'C5', 'C6', "C1'", "C2'", "C3'", "C4'", "C5'",
        'N1', 'N3', 'O4',
        'P'],
  'A': ['H2', 'H61', 'H62', 'H8', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
        'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
        'N1', 'N3', 'N6', 'N7', 'N9',
        'P'],
  'C': ['H41', 'H42', 'H5', 'H6', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
        'C2', 'C4', 'C5', 'C6', "C1'", "C2'", "C3'", "C4'", "C5'",
        'N1', 'N3', 'N4',
        'P']
}

ALL_DNARNA_ATOMS_SORTED = OrderedDict([
  ('O', ["HO2'"]),
  ('1', ["C1'", 'H1', "H1'", 'N1']),
  ('2', ['C2', "C2'", 'H2', "H2'", "H2''", 'H21', 'H22', 'N2']),
  ('3', ["C3'", 'H3', "H3'", 'N3']),
  ('4', ['C4', "C4'", "H4'", 'H41', 'H42', 'N4', 'O4']),
  ('5', ['C5', "C5'", 'H5', "H5'", "H5''"]),
  ('6', ['C6', 'H6', 'H61', 'H62', 'N6']),
  ('7', ['C7', 'H7', 'H71', 'H72', 'H73', 'N7']),
  ('8', ['C8', 'H8']),
  ('9', ['N9']),
  ('P', ['P'])
])

if __name__ == '__main__':
  from ccpn.ui.gui.widgets.Application import TestApplication
  from ccpn.ui.gui.widgets.TextEditor import TextEditor
  app = TestApplication()

  popup = AtomSelectorModule()

  textBox = TextEditor(popup.mainWidget, grid=(3,0), gridSpan=(1,1))

  all_atoms = dict()
  for atom in ['O', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'P']:
    all_atoms[atom] = []

  atomList = [DNA_ATOMS, RNA_ATOMS]
  startAtoms = ['DA', 'DC', 'DG', 'DT', 'A', 'G', 'C', 'U']

  for atomText in atomList:
    atoms = {}
    atomText = atomText.split('\n')
    for line in atomText:
      ll = line.split()
      if ll:
        if ll[0] in startAtoms:
          if ll[0] in atoms.keys():
            atoms[ll[0]].append(ll[1])
          else:
            atoms[ll[0]] = [ll[1]]

          if ll[1] != 'P' and ll[1] not in all_atoms[ll[1][1]]:
            all_atoms[ll[1][1]].append(ll[1])
            all_atoms[ll[1][1]].sort()

    textBox.append(str(atoms))
  all_atoms['P'] = ['P']
  textBox.append(str(all_atoms))

  popup.currentNmrResidueLabel.setText('NmrResidue here')

  print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DT'))
  print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DC'))
  print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DA'))
  print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DG'))

  print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'G'))
  print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'U'))
  print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'A'))
  print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'C'))

  popup.show()
  popup.raise_()
  app.start()

