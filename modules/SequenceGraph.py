"""Module Documentation here

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
__author__ = "$Author: skinnersp $"
__date__ = "$Date: 2016-05-23 10:02:47 +0100 (Thu, 26 May 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

import json
import typing

import numpy as np
from PyQt4 import QtGui, QtCore

from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.lib.AssignmentLib import getNmrResiduePrediction
from ccpn.core.lib.AssignmentLib import nmrAtomPairsByDimensionTransfer
from ccpn.ui.gui.guiSettings import textFont, textFontBold, textFontLarge
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Icon import Icon
from ccpn.ui.gui.widgets.ToolBar import ToolBar

from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget
from ccpn.ui.gui.widgets.PulldownListsForObjects import NmrChainPulldown
from ccpn.core.NmrChain import NmrChain

from ccpn.util.Constants import ccpnmrJsonData
from ccpn.util.Logging import getLogger

logger = getLogger()


class GuiNmrAtom(QtGui.QGraphicsTextItem):
  """
  A graphical object specifying the position and name of an atom when created by the Assigner.
  Can be linked to a Nmr Atom.
  """
  def __init__(self, project, text, pos=None, nmrAtom=None):

    super(GuiNmrAtom, self).__init__()

    self.setPlainText(text)
    self.setPos(QtCore.QPointF((pos[0]-self.boundingRect().x()), (pos[1]-self.boundingRect().y())))

    self.project = project
    self.current = project._appBase.current
    self.nmrAtom = nmrAtom
    ###if nmrAtom:
    ###  self.name = nmrAtom.name
    self.connectedAtoms = 0

    # wb104: not sure why below is needed rather than setFlags() but it is
    self.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
    ###self.setFlags(QtGui.QGraphicsItem.ItemIsSelectable | self.flags())
    #self.setFlag(QtGui.QGraphicsItem.ItemIsSelectable)

    if project._appBase.colourScheme == 'dark':
      colour1 = '#F7FFFF'
      colour2 = '#BEC4F3'
    elif project._appBase.colourScheme == 'light':
      colour1 = '#FDFDFC'
      colour2 = '#555D85'

    #self.setDefaultTextColor(QtGui.QColor(colour1))
    if self.isSelected:
      self.setDefaultTextColor(QtGui.QColor(colour2))
    else:
      self.setDefaultTextColor(QtGui.QColor(colour1))

  def mouseDoubleClickEvent(self, event):
    """
    CCPN INTERNAL - re-implementation of double click event
    """
    #print('>>doubleClickEvent')
    # if self.nmrAtom is Not None:
    #   self.current.nmrAtom = self.nmrAtom
    #   self.current.nmrResidue = self.nmrAtom.nmrResidue
    pass

  def mousePressEvent(self, event):
    """
    CCPN INTERNAL - re-implementation of mouse press event
    """
    ###print('>>pressEvent')
    if self.nmrAtom is not None:
      self.current.nmrAtom = self.nmrAtom
      self.current.nmrResidue = self.nmrAtom.nmrResidue

  def mouseReleaseEvent(self, event):
    """
    CCPN INTERNAL - re-implementation of mouse press event
    """
    #print('>>release Event')
    pass


class GuiNmrResidue(QtGui.QGraphicsTextItem):
  """
  Object linking residues displayed in Assigner and Nmr Residues. Contains functionality for drag and
  drop assignment in conjunction with the Sequence Module.
  """

  def __init__(self, parent, nmrResidue, caAtom):

    super(GuiNmrResidue, self).__init__()
    self.setPlainText(nmrResidue.id)
    project = nmrResidue.project
    self.setFont(textFont)
    if project._appBase.colourScheme == 'dark':
      self.setDefaultTextColor(QtGui.QColor('#F7FFFF'))
    elif project._appBase.colourScheme == 'light':
      self.setDefaultTextColor(QtGui.QColor('#555D85'))
    self.setPos(caAtom.x()-caAtom.boundingRect().width()/2, caAtom.y()+30)
    ###self.setFlags(QtGui.QGraphicsItem.ItemIsSelectable | self.flags())
    self.setFlag(QtGui.QGraphicsItem.ItemIsSelectable)
    self.parent = parent
    self.nmrResidue = nmrResidue
    self.mousePressEvent = self._mousePressEvent
    self.mouseMoveEvent = self._mouseMoveEvent

  def _update(self):
    self.setPlainText(self.nmrResidue.id)


  def _mouseMoveEvent(self, event):

    nmrItem = None
    if (event.buttons() == QtCore.Qt.LeftButton) and (event.modifiers() & QtCore.Qt.ShiftModifier):
        for item in self.parent.scene.items():
          if isinstance(item, GuiNmrResidue) and item.isSelected():
            nmrChain = item.nmrResidue.nmrChain
            nmrItem = item    #.parentWidget()

        if nmrItem:
          drag = QtGui.QDrag(event.widget())
          mimeData = QtCore.QMimeData()
          itemData = json.dumps({'pids': [nmrChain.pid]})
          mimeData.setData(ccpnmrJsonData, itemData)
          mimeData.setText(itemData)
          drag.setMimeData(mimeData)

          #TODO:ED get rid of _appBase
          dragLabel = QtGui.QLabel()
          dragLabel.setText(self.toPlainText())
          dragLabel.setFont(textFontLarge)
          if nmrItem.nmrResidue.project._appBase.colourScheme == 'dark':
            dragLabel.setStyleSheet('color : #F7FFFF')
          elif nmrItem.nmrResidue.project._appBase.colourScheme == 'light':
            dragLabel.setStyleSheet('color : #555D85')

          pixmap = QtGui.QPixmap.grabWidget(dragLabel)    # ejb -    this gets the whole window   event.widget())
          painter = QtGui.QPainter(pixmap)
          painter.setCompositionMode(painter.CompositionMode_DestinationIn)
          painter.fillRect(pixmap.rect(), QtGui.QColor(0, 0, 0, 240))
          painter.end()
          drag.setPixmap(pixmap)
          drag.setHotSpot(QtCore.QPoint(dragLabel.width() / 2, dragLabel.height() / 2))

          drag.start(QtCore.Qt.MoveAction)      # ejb - same as BackboneAssignment

          # if drag.exec_(QtCore.Qt.MoveAction | QtCore.Qt.CopyAction, QtCore.Qt.CopyAction) == QtCore.Qt.MoveAction:
          #   pass
          #   # self.close()
          # else:
          #   self.show()

  def _mousePressEvent(self, event):
    self.nmrResidue.project._appBase.current.nmrResidue = self.nmrResidue
    self.setSelected(True)



class AssignmentLine(QtGui.QGraphicsLineItem):
  """
  Object to create lines between GuiNmrAtoms with specific style, width, colour and displacement.
  """
  def __init__(self, x1, y1, x2, y2, colour, width, style=None):
    QtGui.QGraphicsLineItem.__init__(self)
    self.pen = QtGui.QPen()
    self.pen.setColor(QtGui.QColor(colour))
    self.pen.setCosmetic(True)
    self.pen.setWidth(width)
    ###if style and style == 'dash':
    if style == 'dash':
      self.pen.setStyle(QtCore.Qt.DotLine)
    self.setPen(self.pen)
    self.setLine(x1, y1, x2, y2)


class SequenceGraphModule(CcpnModule):
  """
  A module for the display of stretches of sequentially linked and assigned stretches of
  NmrResidues.
  """

  className = 'SequenceGraph'

  includeSettingsWidget = False
  maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
  settingsPosition = 'top'

  def __init__(self, mainWindow, name='Sequence Graph', nmrChain=None):

    CcpnModule.__init__(self, mainWindow=mainWindow, name=name)

    # Derive application, project, and current from mainWindow
    self.mainWindow = mainWindow
    self.application = mainWindow.application
    self.project = mainWindow.application.project
    self.current = mainWindow.application.current

    ###frame = Frame(parent=self.mainWidget)
    self._sequenceGraphScrollArea = QtGui.QScrollArea()
    self._sequenceGraphScrollArea.setWidgetResizable(True)
    self.scene = QtGui.QGraphicsScene(self)
    self.scrollContents = QtGui.QGraphicsView(self.scene, self)
    self.scrollContents.setRenderHints(QtGui.QPainter.Antialiasing)
    self.scrollContents.setInteractive(True)
    self.scrollContents.setGeometry(QtCore.QRect(0, 0, 380, 1000))
    ###self.horizontalLayout2 = QtGui.QHBoxLayout(self.scrollContents)
    self._sequenceGraphScrollArea.setWidget(self.scrollContents)
    self.mainWidget.getLayout().addWidget(self._sequenceGraphScrollArea, 2, 0, 1, 6)
    #frame.addWidget(self._sequenceGraphScrollArea, 4, 0, 1, 6)

    self.residueCount = 0

    #TODO:GEERTEN: StyleSheet
    if self.application.colourScheme == 'dark':
      self._lineColour = '#f7ffff'
    elif self.application.colourScheme == 'light':
      self._lineColour = ''  # TODO: check if correct

    """
    self.modeLabel = Label(self, 'Mode: ', grid=(0, 3))
    self.modePulldown = PulldownList(self, grid=(0, 4), gridSpan=(1, 1), callback=self.setMode)
    self.modePulldown.setData(['fragment', 'Assigned - backbone'])
    """
    self.nmrChainPulldown = NmrChainPulldown(self.mainWidget, self.project, grid=(0, 0), gridSpan=(1, 1),
                                             showSelectName=True,
                                             callback=self.setNmrChainDisplay)

    self.refreshCheckBox = CheckBoxCompoundWidget(self.mainWidget,
                                                  labelText='Auto refresh NmrChain:',
                                                  checked=True,
                                                  tipText='Update display when current.nmrChain changes',
                                                  grid=(0, 1), gridSpan=(1,1))

    self.assignmentsCheckBox = CheckBoxCompoundWidget(self.mainWidget,
                                                      labelText='Show peak assignments:',
                                                      checked=False,
                                                      tipText='Show peak assignments on display coloured by positiveContourColour',
                                                      callback=self._updateShownAssignments,
                                                      grid=(0, 2), gridSpan=(1,1))

    self.nmrResiduesCheckBox = CheckBoxCompoundWidget(self.mainWidget,
                                                      labelText='Show all NmrResidues:',
                                                      checked=False,
                                                      tipText='Show all the NmrResidues in the NmrChain',
                                                      callback=self._updateShownAssignments,
                                                      grid=(0, 3), gridSpan=(1,1))

    self.editingToolbar = ToolBar(self.mainWidget, grid=(0, 5), gridSpan=(1, 1), hAlign='right', iconSizes=(24,24))

    self.disconnectPreviousAction = self.editingToolbar.addAction("disconnectPrevious", self.disconnectPreviousNmrResidue)
    self.disconnectPreviousIcon = Icon('icons/previous')
    self.disconnectPreviousAction.setIcon(self.disconnectPreviousIcon)
    self.disconnectAction = self.editingToolbar.addAction("disconnect", self.disconnectNmrResidue)
    self.disconnectIcon = Icon('icons/minus')
    self.disconnectAction.setIcon(self.disconnectIcon)
    self.disconnectNextAction = self.editingToolbar.addAction("disconnectNext", self.disconnectNextNmrResidue)
    self.disconnectNextIcon = Icon('icons/next')
    self.disconnectNextAction.setIcon(self.disconnectNextIcon)
    #self.editingToolbar.hide()

    self.atomSpacing = 66
    self.guiResiduesShown = []
    self.predictedStretch = []
    self.direction = None
    self.selectedStretch = []
    self.scene.dragEnterEvent = self.dragEnterEvent
    self.guiNmrResidues = []
    self.guiNmrAtomDict = {}

    ###self.setMode('fragment')  # cannot be moved up!
    self._registerNotifiers()

    if nmrChain is not None:
      self.selectSequence(nmrChain)

  def selectSequence(self, nmrChain=None):
    """
    Manually select a Sequence from the pullDown
    """
    if nmrChain is None:
      logger.warning('select: No Sequence selected')
      raise ValueError('select: No Sequence selected')
    else:
      if not isinstance(nmrChain, NmrChain):
        logger.warning('select: Object is not of type Sequence')
        raise TypeError('select: Object is not of type Sequence')
      else:
        for widgetObj in self.nmrChainPulldown.textList:
          if nmrChain.pid == widgetObj:
            self.nmrChainPulldown.select(nmrChain.pid)
            self.setNmrChainDisplay(nmrChain)

  def _registerNotifiers(self):
    self.current.registerNotify(self._updateModule, 'nmrChains')

    self.project.registerNotifier('NmrResidue', 'rename', self._resetNmrResiduePidForAssigner)
#    self.project.registerNotifier('NmrChain', 'delete', self.removeNmrChainFromPulldown)
#    self.project.registerNotifier('NmrChain', 'create', self.addNmrChainToPulldown)
    self.project.registerNotifier('Peak', 'change', self._updateShownAssignments)

  def _unRegisterNotifiers(self):
    self.project.unRegisterNotifier('NmrResidue', 'rename', self._resetNmrResiduePidForAssigner)
#    self.project.unRegisterNotifier('NmrChain', 'delete', self.removeNmrChainFromPulldown)
#    self.project.unRegisterNotifier('NmrChain', 'create', self.addNmrChainToPulldown)
    self.project.unRegisterNotifier('Peak', 'change', self._updateShownAssignments)

  def _updateModule(self, nmrChains=None):
    """
    Update in reponse to change of current.nmrChains
    """
    #if nmrChains is None or len(nmrChains)==0: return
    nmrChain = self.current.nmrChain
    if not nmrChain:
      return

    if not self.refreshCheckBox.isChecked():
      return

    #self.sequenceGraph.clearAllItems()
    ###self.nmrChainPulldown.pulldownList.select(self.current.nmrChain.pid)
    self.nmrChainPulldown.select(nmrChain.pid)
    ###self.setNmrChainDisplay(nmrChain.pid)
    self.setNmrChainDisplay(nmrChain)

  """
  def setMode(self, mode):
    if self.project.nmrChains:
      self.editingToolbar.hide()
      if mode == 'fragment':
        self.editingToolbar.show()
        #self.nmrChainPulldown.setData([c.pid for c in self.project.nmrChains])
        #self.nmrChainLabel.setText('NmrChain: ')
      elif mode == 'Assigned - backbone':
        pass
        #self.nmrChainLabel.setText('Chain: ')
        #self.nmrChainPulldown.setData([self.project.getByPid('NC:%s' % chain.shortName).pid for chain in self.project.chains if self.project.getByPid('NC:%s' % chain.shortName)])
      self.modePulldown.select(mode)
      self.setNmrChainDisplay(self.nmrChainPulldown.getText())
    else:
      logger.warning('No valid NmrChain is selected.')
"""

  def setNmrChainDisplay(self, nmrChainOrPid):

    if isinstance(nmrChainOrPid, str):
      nmrChain = self.project.getByPid(nmrChainOrPid)
    else:
      nmrChain = nmrChainOrPid

    # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
    if not nmrChain:
      return

    ###self.project._appBase._startCommandBlock('application.sequenceGraph.setNmrChainDisplay({!r})'.format(nmrChainPid))
    self.project._appBase._startCommandBlock('application.sequenceGraph.setNmrChainDisplay({!r})'.format(nmrChain.pid))
    try:
      #self.current.nmrChain = self.project.getByPid(nmrChainPid)
      #if not self.current.nmrChain:
      #  logger.warning('No NmrChain selected.')
      #  return
      self.clearAllItems()

      ###nmrChain = self.project.getByPid(nmrChainPid)
      ###if self.modePulldown.currentText() == 'fragment':
      if True:
        """
        if nmrChain.isConnected:
          for nmrResidue in nmrChain.mainNmrResidues:
            self.addResidue(nmrResidue, '+1')
        elif self.current.nmrResidue is not None and self.current.nmrResidue in nmrChain.nmrResidues:
          self.addResidue(self.current.nmrResidue, '+1')
"""
        connectingLinesNeeded = set()
        if self.nmrResiduesCheckBox.isChecked():
          for nmrResidue in nmrChain.nmrResidues:
            if nmrResidue is nmrResidue.mainNmrResidue:
              self.addResidue(nmrResidue, '+1')
              if nmrResidue.nextNmrResidue:
                connectingLinesNeeded.add(len(self.guiResiduesShown)-1)
        else:
          nmrResidue = self.current.nmrResidue
          if nmrResidue in nmrChain.nmrResidues:
            while nmrResidue.previousNmrResidue: # go to start of connected stretch
              nmrResidue = nmrResidue.previousNmrResidue
          elif nmrChain.isConnected or nmrChain.chain: # either NC:# or NC:A type nmrChains but not NC:@
            nmrResidue = nmrChain.mainNmrResidues[0]

          while nmrResidue:  # add all of connected stretch
            self.addResidue(nmrResidue, '+1')
            nmrResidue = nmrResidue.nextNmrResidue

        if len(self.predictedStretch) > 2:
          self.predictSequencePosition(self.predictedStretch)

      ###elif self.modePulldown.currentText() == 'Assigned - backbone':
      ###  self._showBackboneAssignments(nmrChain)

      for ii, res in enumerate(self.guiResiduesShown[:-1]):
        if not self.nmrResiduesCheckBox.isChecked() or ii in connectingLinesNeeded:
          self._addConnectingLine(res['CO'], self.guiResiduesShown[ii + 1]['N'], self._lineColour, 1.0, 0)

      if self.assignmentsCheckBox.isChecked():
        self._getAssignmentsFromSpectra()

    finally:
      self.project._endCommandEchoBlock()

  def resetSequenceGraph(self):

    self.nmrChainPulldown.pulldownList.select('NC:@-')

  def disconnectPreviousNmrResidue(self):
    self.current.nmrResidue.disconnectPrevious()
    self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
    ###self.updateNmrResidueTable()

  def _closeModule(self):
    self._unRegisterNotifiers()
    #delattr(self.parent, 'sequenceGraph')
    super(SequenceGraphModule, self)._closeModule()

  def disconnectNextNmrResidue(self):
    self.current.nmrResidue.disconnectNext()
    self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
    #self.updateNmrResidueTable()

  def disconnectNmrResidue(self):
    self.current.nmrResidue.disconnect()
    self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
    #self.updateNmrResidueTable()

  def _resetNmrResiduePidForAssigner(self, nmrResidue, oldPid:str):
    """Reset pid for NmrResidue and all offset NmrResidues"""
    for nr in [nmrResidue] + list(nmrResidue.offsetNmrResidues):
      for guiNmrResidue in self.guiNmrResidues:
        if guiNmrResidue.nmrResidue is nr:
          guiNmrResidue._update()

  def clearAllItems(self):
    """
    Removes all displayed residues in the sequence graph and resets items count to zero.
    """
    for item in self.scene.items():
      self.scene.removeItem(item)

    self.residueCount = 0
    self.predictedStretch = []
    self.guiResiduesShown = []
    self.guiNmrResidues = []
    self.guiNmrAtomDict = {}
    self.scene.clear()

  def _assembleResidue(self, nmrResidue:NmrResidue, atoms:typing.Dict[str, GuiNmrAtom]):
    """
    Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
    creates a graphical representation of a residue in the assigner
    """

    for item in atoms.values():
      self.scene.addItem(item)

    nmrAtoms = [atom.name for atom in nmrResidue.nmrAtoms]
    if "CB" in list(atoms.keys()):
      self._addConnectingLine(atoms['CA'], atoms['CB'], self._lineColour, 1.0, 0)
    if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
      self._addConnectingLine(atoms['H'], atoms['N'], self._lineColour, 1.0, 0)
    if nmrResidue.residueType != 'PRO':
        self._addConnectingLine(atoms['H'], atoms['N'], self._lineColour, 1.0, 0)
    else:
      self.scene.removeItem(atoms['H'])
    # if not 'CB' in nmrAtoms:
    #   self.scene.removeItem(atoms['CB'])
    #   self.scene.removeItem(cbLine)

    self._addConnectingLine(atoms['N'], atoms['CA'], self._lineColour, 1.0, 0)
    self._addConnectingLine(atoms['CO'], atoms['CA'], self._lineColour, 1.0, 0)
    self.nmrResidueLabel = GuiNmrResidue(self, nmrResidue, atoms['CA'])
    self.guiNmrResidues.append(self.nmrResidueLabel)
    self.scene.addItem(self.nmrResidueLabel)
    self._addResiduePredictions(nmrResidue, atoms['CA'])

  """
  def addSideChainAtoms(self, nmrResidue, cbAtom, colour):
    residue = {}
    for k, v in ATOM_POSITION_DICT[nmrResidue.residueType].items():
      if k != 'boundAtoms':
        position = [cbAtom.x()+v[0], cbAtom.y()+v[1]]
        nmrAtom = nmrResidue.fetchNmrAtom(name=k)
        newAtom = self._createGuiNmrAtom(k, position, nmrAtom)
        self.scene.addItem(newAtom)
        residue[k] = newAtom
        self.guiNmrAtomDict[nmrAtom] = newAtom

    for boundAtomPair in ATOM_POSITION_DICT[nmrResidue.residueType]['boundAtoms']:
      atom1 = residue[boundAtomPair[0]]
      atom2 = residue[boundAtomPair[1]]
      newLine = AssignmentLine(atom1.x(), atom1.y(), atom2.x(), atom2.y(), colour, 1.0)
      self.scene.addItem(newLine)
"""

  def addResidue(self, nmrResidue:NmrResidue, direction:str, atomSpacing=None):
    """
    Takes an Nmr Residue and a direction, either '-1 or '+1', and adds a residue to the sequence graph
    corresponding to the Nmr Residue.
    Nmr Residue name displayed beneath CA of residue drawn and residue type predictions displayed
    beneath Nmr Residue name
    """
    atoms = {}
    if atomSpacing:
      self.atomSpacing = atomSpacing
    nmrAtoms = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]

    residueAtoms = {"H":  np.array([0, 0]),
                    "N":  np.array([0, -1*self.atomSpacing]),
                    "CA": np.array([self.atomSpacing, -1*self.atomSpacing]),
                    "CB": np.array([self.atomSpacing, -2*self.atomSpacing]),
                    "CO": np.array([2*self.atomSpacing, -1*self.atomSpacing])
                    }
    if nmrResidue.residueType == 'GLY':
      del residueAtoms['CB']
    if self.residueCount == 0:
      for k, v in residueAtoms.items():
        if k in nmrAtoms:
          nmrAtom = nmrResidue.fetchNmrAtom(name=k)
        else:
          nmrAtom = None
        atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)
      self.guiResiduesShown.append(atoms)
      self.predictedStretch.append(nmrResidue)

    else:
      for k, v in residueAtoms.items():
        if k in nmrAtoms:
          nmrAtom = nmrResidue.fetchNmrAtom(name=k)
        else:
          nmrAtom = None

        if direction == '-1':
          pos = np.array([self.guiResiduesShown[0]['H'].x()-3*self.atomSpacing, self.guiResiduesShown[0]['H'].y()])
          atoms[k] = self._createGuiNmrAtom(k, v+pos, nmrAtom)

        else:
          pos = np.array([self.guiResiduesShown[-1]['H'].x()+3*self.atomSpacing, self.guiResiduesShown[-1]['H'].y()])
          atoms[k] = self._createGuiNmrAtom(k, v+pos, nmrAtom)

      if direction == '-1':
        self.guiResiduesShown.insert(0, atoms)
        self.predictedStretch.insert(0, nmrResidue)
      else:
        self.guiResiduesShown.append(atoms)
        self.predictedStretch.append(nmrResidue)

    self._assembleResidue(nmrResidue, atoms)

    self.residueCount += 1

  def _addResiduePredictions(self, nmrResidue:NmrResidue, caAtom:GuiNmrAtom):
    """
    Gets predictions for residue type based on BMRB statistics and determines label positions
    based on caAtom position.
    """

    predictions = list(set(map(tuple, (getNmrResiduePrediction(nmrResidue, self.project.chemicalShiftLists[0])))))
    predictions.sort(key=lambda a: float(a[1][:-1]), reverse=True)
    for prediction in predictions:
      predictionLabel = QtGui.QGraphicsTextItem()
      predictionLabel.setPlainText(prediction[0]+' '+prediction[1])
      if self.project._appBase.colourScheme == 'dark':
        predictionLabel.setDefaultTextColor(QtGui.QColor('#F7FFFF'))
      elif self.project._appBase.colourScheme == 'light':
        predictionLabel.setDefaultTextColor(QtGui.QColor('#555D85'))
      predictionLabel.setFont(textFontBold)
      predictionLabel.setPos(caAtom.x()-caAtom.boundingRect().width()/2,
                             caAtom.y()+(30*(predictions.index(prediction)+2)))
      self.scene.addItem(predictionLabel)

  def predictSequencePosition(self, nmrResidues:list):
    """
    Predicts sequence position for Nmr residues displayed in the Assigner and highlights appropriate
    positions in the Sequence Module if it is displayed.
    """
    from ccpn.core.lib.AssignmentLib import getSpinSystemsLocation

    possibleMatches = getSpinSystemsLocation(self.project, nmrResidues,
                      self.project.chains[0], self.project.chemicalShiftLists[0])

    for possibleMatch in possibleMatches:
      if possibleMatch[0] > 1 and not len(possibleMatch[1]) < len(nmrResidues):
        if hasattr(self.project._appBase, 'sequenceModule'):
          self.project._appBase.sequenceModule._highlightPossibleStretches(possibleMatch[1])

  def _updateShownAssignments(self, peak=None):
    ###if self.current.nmrChain is not None:
    ###  self.setNmrChainDisplay(self.current.nmrChain.pid)

    nmrChainPid = self.nmrChainPulldown.getText()
    if nmrChainPid:
      self.setNmrChainDisplay(nmrChainPid)

  """
  def _showBackboneAssignments(self, nmrChain):
    self.project._startCommandEchoBlock('_showBackboneAssignments', nmrChain)
    try:

      for residue in nmrChain.chain.residues:
        if not residue.nmrResidue:
          newNmrResidue = nmrChain.fetchNmrResidue(sequenceCode=residue.sequenceCode, residueType=residue.residueType)
          for atom in residue.atoms:
            newNmrResidue.fetchNmrAtom(name=atom.name)
        self.addResidue(residue.nmrResidue, direction='+1')
      for ii, res in enumerate(self.guiResiduesShown):
        if ii % 10 == 0:
          if self.project._appBase.ui.mainWindow is not None:
            mainWindow = self.project._appBase.ui.mainWindow
          else:
            mainWindow = self.project._appBase._mainWindow
          mainWindow.pythonConsole.writeConsoleCommand('%s residues added' % str(ii))
        ###if ii+1 < len(self.guiResiduesShown)-1:
        if ii + 1 < len(self.guiResiduesShown):
            self._addConnectingLine(res['CO'], self.guiResiduesShown[ii+1]['N'], self._lineColour, 1.0, 0)

      self._getAssignmentsFromSpectra()
    finally:
      self.project._endCommandEchoBlock()
"""

  def _addConnectingLine(self, atom1:GuiNmrAtom, atom2:GuiNmrAtom, colour:str, width:float, displacement:float, style:str=None):
    """
    Adds a line between two GuiNmrAtoms using the width, colour, displacement and style specified.
    """
    if atom1.y() > atom2.y():
      y1 = atom1.y() - (atom1.boundingRect().height()*.05)-displacement
      y2 = atom2.y() + (atom2.boundingRect().height())-displacement

    elif atom1.y() < atom2.y():
      y1 = atom1.y() + (atom1.boundingRect().height())-displacement
      y2 = atom2.y() - (atom2.boundingRect().height()*0.08)-displacement

    else:
      y1 = atom1.y() + (atom1.boundingRect().height()*0.5)
      y2 = atom2.y() + (atom2.boundingRect().height()*0.5)

    if atom1.x() > atom2.x():
      x1 = atom1.x()
      x2 = atom2.x() + atom2.boundingRect().width()

    elif atom1.x() < atom2.x():
      x1 = atom1.x() + atom1.boundingRect().width()
      x2 = atom2.x()

    else:
      x1 = atom1.x() + (atom1.boundingRect().width()/2)+displacement
      x2 = atom2.x() + (atom1.boundingRect().width()/2)+displacement
      y1 += displacement
      y2 += displacement

    newLine = AssignmentLine(x1, y1, x2, y2, colour, width, style)
    self.scene.addItem(newLine)
    return newLine

  def _createGuiNmrAtom(self, atomType:str, position:tuple, nmrAtom:NmrAtom=None) -> GuiNmrAtom:
    """
    Creates a GuiNmrAtom specified by the atomType and graphical position supplied.
    GuiNmrAtom can be linked to an NmrAtom by supplying it to the function.
    """
    atom = GuiNmrAtom(self.project, text=atomType, pos=position, nmrAtom=nmrAtom)
    self.guiNmrAtomDict[nmrAtom] = atom
    return atom

  def _getAssignmentsFromSpectra(self):
    for spectrum in self.project.spectra:
      connections = [x for y in list(nmrAtomPairsByDimensionTransfer(spectrum.peakLists).values())
                     for x in y]
      for ii, connection in enumerate(connections):
        # nmrAtomPair = [self.project._data2Obj.get(connection[0]).nmrAtom,
        #                self.project._data2Obj.get(connection[1]).nmrAtom]
        # sorting makes sure drawing is done properly
        guiNmrAtomPair = [self.guiNmrAtomDict.get(a) for a in sorted(connection, reverse=True)]
        if None not in guiNmrAtomPair:
          displacement = min(guiNmrAtomPair[0].connectedAtoms, guiNmrAtomPair[1].connectedAtoms)
          self._addConnectingLine(guiNmrAtomPair[0], guiNmrAtomPair[1], spectrum.positiveContourColour, 2.0, displacement)
          guiNmrAtomPair[0].connectedAtoms += 1.0
          guiNmrAtomPair[1].connectedAtoms += 1.0

import math
atomSpacing = 66
cos36 = math.cos(math.pi/5)
sin36 = math.sin(math.pi/5)
tan36 = math.tan(math.pi/5)

cos54 = math.cos(3*math.pi/10)
sin54 = math.sin(3*math.pi/10)

cos60 = math.cos(math.pi/3)
sin60 = math.sin(math.pi/3)
sin72 = math.sin(2*math.pi/5)
cos72 = math.cos(2*math.pi/5)

ATOM_POSITION_DICT = {

  'ALA': {'HB%': [0.0, -0.75*atomSpacing],
          'boundAtoms': ['']},
  'CYS': {'SG':  [0.0, -1*atomSpacing], 'HG': [0, -1.75*atomSpacing]},
  'ASP': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'CG': [0, -1*atomSpacing]},
  'ASN': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'CG': [0, -1*atomSpacing], 'ND2': [0, -2*atomSpacing],
          'HD2x': [atomSpacing*-0.75, -2*atomSpacing-(0.75*atomSpacing*cos60)],
          'HD2y': [atomSpacing*+0.75, -2*atomSpacing-(0.75*atomSpacing*cos60)],
          },
  'GLU': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'HGx': [atomSpacing*-0.75, -1*atomSpacing], 'HGy': [atomSpacing*0.75, -1*atomSpacing],
          'CG':  [0, -1*atomSpacing], 'CD': [0, -2*atomSpacing]
          },
  'GLN': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'HGx': [atomSpacing*-0.75, -1*atomSpacing], 'HGy': [atomSpacing*0.75, -1*atomSpacing],
          'CG':  [0, -1*atomSpacing], 'CD': [0, -2*atomSpacing], 'NE2': [0, -3*atomSpacing],
          'HD2x': [atomSpacing*-0.75, -3*atomSpacing-(0.75*atomSpacing*cos60)],
          'HD2y': [atomSpacing*+0.75, -3*atomSpacing-(0.75*atomSpacing*cos60)],
          },
  'PHE': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'CG':  [0, -1*atomSpacing], 'CD1': [-1*atomSpacing, (-1-cos60)*atomSpacing],
          'CD2': [1*atomSpacing, (-1-cos60)*atomSpacing],
          'CE1': [-1*atomSpacing, (-2-cos60)*atomSpacing],
          'CE2': [1*atomSpacing, (-2-cos60)*atomSpacing],
          'HD1': [-1.75*atomSpacing, (-1-cos60)*atomSpacing],
          'HD2': [1.75*atomSpacing, (-1-cos60)*atomSpacing],
          'HE1': [-1.75*atomSpacing, (-2-cos60)*atomSpacing],
          'HE2': [1.75*atomSpacing, (-2-cos60)*atomSpacing],
          'CZ': [0, (-2-cos60-sin60)*atomSpacing], 'HZ': [0, (-2-cos60-sin60-0.75)*atomSpacing]
          },
  'TYR': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'CG':  [0, -1*atomSpacing], 'CD1': [-1*atomSpacing, (-1-cos60)*atomSpacing],
          'CD2': [1*atomSpacing, (-1-cos60)*atomSpacing],
          'CE1': [-1*atomSpacing, (-2-cos60)*atomSpacing],
          'CE2': [1*atomSpacing, (-2-cos60)*atomSpacing],
          'HD1': [-1.75*atomSpacing, (-1-cos60)*atomSpacing],
          'HD2': [1.75*atomSpacing, (-1-cos60)*atomSpacing],
          'HE1': [-1.75*atomSpacing, (-2-cos60)*atomSpacing],
          'HE2': [1.75*atomSpacing, (-2-cos60)*atomSpacing],
          'CZ': [0, (-2-cos60-sin60)*atomSpacing], 'HH': [0, (-2-cos60-sin60-0.75)*atomSpacing]
          },
  'SER': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'HG': [0, -1*atomSpacing]
          },
  'THR': {'HG1': [atomSpacing*-0.75, 0.0], 'HB': [atomSpacing*0.75, 0.0],
          'CG2': [0, -1*atomSpacing], 'HG2%': [0, -1.75*atomSpacing]
          },
  'MET': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'HGx': [atomSpacing*-0.75, -1*atomSpacing], 'HGy': [atomSpacing*0.75, -1*atomSpacing],
          'CG':  [0, -1*atomSpacing], 'SD': [0, -2*atomSpacing], 'CE': [0, -3*atomSpacing],
          'HE%': [0, -3.75*atomSpacing]
          },
  'ARG': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'HGx': [atomSpacing*-0.75, -1*atomSpacing], 'HGy': [atomSpacing*0.75, -1*atomSpacing],
          'CG':  [0, -1*atomSpacing], 'CD': [0, -2*atomSpacing], 'NE': [0, -3*atomSpacing],
          'CZ': [0, -4*atomSpacing], 'NH1': [atomSpacing*-1, -4*atomSpacing-(0.75*atomSpacing*cos60)],
          'NH2': [atomSpacing*+1, -4*atomSpacing-(0.75*atomSpacing*cos60)],
          },
  'VAL': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'CGx': [-1*atomSpacing, -1*(cos60*atomSpacing)],
          'CGy': [1*atomSpacing, -1*(cos60*atomSpacing)],
          'HGx%': [atomSpacing*-1, -1*(cos60*atomSpacing)-(0.75*atomSpacing)],
          'HGy%': [atomSpacing*+1, -1*(cos60*atomSpacing)-(0.75*atomSpacing)]
          },
  'LEU': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'HGx': [atomSpacing*-0.75, -1*atomSpacing], 'HGy': [atomSpacing*0.75, -1*atomSpacing],
          'CG':  [0, -1*atomSpacing],
          'CDx': [-1*atomSpacing, (-1-cos60)*atomSpacing],
          'CDy': [1*atomSpacing, (-1-cos60)*atomSpacing],
          'HDx%': [atomSpacing*-1, ((-1-cos60)*atomSpacing)-(0.75*atomSpacing)],
          'HDy%': [atomSpacing*+1, ((-1-cos60)*atomSpacing)-(0.75*atomSpacing)]
          },
  'ILE': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'CG1': [-1*atomSpacing, -1*(cos60*atomSpacing)],
          'CG2%': [1*atomSpacing, -1*(cos60*atomSpacing)],
          'HG1x': [atomSpacing*-1.75, -1*(cos60*atomSpacing)],
          'HG1y': [atomSpacing*-0.25, -1*(cos60*atomSpacing)],
          'HG2%': [1*atomSpacing, -1*(cos60*atomSpacing)-(0.75*atomSpacing)],
          'CD1%': [-1*atomSpacing, -1*(cos60*atomSpacing)-atomSpacing],
          'HD1%': [-1*atomSpacing, -1*(cos60*atomSpacing)-(1.75*atomSpacing)],
          },
  'LYS': {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
          'HGx': [atomSpacing*-0.75, -1*atomSpacing], 'HGy': [atomSpacing*0.75, -1*atomSpacing],
          'CG':  [0, -1*atomSpacing], 'CD': [0, -2*atomSpacing],
          'HDx': [atomSpacing*-0.75, -2*atomSpacing], 'HDy': [atomSpacing*0.75, -2*atomSpacing],
          'HEx': [atomSpacing*-0.75, -3*atomSpacing], 'HEy': [atomSpacing*0.75, -3*atomSpacing],
          'CE': [0, -3*atomSpacing],
          'NZ': [0, -4*atomSpacing], 'HZ%': [0, -4.75*atomSpacing],
          },
  'HIS':  {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
           'CG':  [0, -1*atomSpacing], 'ND1': [-1*atomSpacing, -1*(atomSpacing+(atomSpacing/(2*tan36)))],
           'CD2': [atomSpacing, -1*(atomSpacing+(atomSpacing/(2*tan36)))],
           'NE2': [atomSpacing/2, -1*(atomSpacing+(atomSpacing/(2*sin36))+(atomSpacing/(2*tan36)))],
           'CD1': [-0.5*atomSpacing, -1*(atomSpacing+(atomSpacing/(2*sin36))+(atomSpacing/(2*tan36)))],
          },


  'TRP':  {'HBx': [atomSpacing*-0.75, 0.0], 'HBy': [atomSpacing*0.75, 0.0],
           'CG':  [0, -1*atomSpacing], 'CD1': [atomSpacing, -1*atomSpacing],
           'NE1': [atomSpacing+(atomSpacing*cos72), -1*(atomSpacing+(atomSpacing*sin72))],
           'CE2': [atomSpacing+(atomSpacing*cos72)-(atomSpacing*sin54),
                   -1*(atomSpacing+(atomSpacing*sin72)+(atomSpacing*cos54))],
           'CD2': [-1*(atomSpacing*cos72), -1*(atomSpacing+(atomSpacing*sin72))],
           'CE3': [atomSpacing+(atomSpacing*cos72)-(atomSpacing*sin54)-(2*(atomSpacing*sin60)),
                   -1*(atomSpacing+(atomSpacing*sin72)+(atomSpacing*cos54))],
           'CZ2': [atomSpacing+(atomSpacing*cos72)-(atomSpacing*sin54),
                   -1*(2*atomSpacing+(atomSpacing*sin72)+(atomSpacing*cos54))],
           'CZ3': [atomSpacing+(atomSpacing*cos72)-(atomSpacing*sin54)-(2*(atomSpacing*sin60)),
                   -1*(2*atomSpacing+(atomSpacing*sin72)+(atomSpacing*cos54))],
           'CH2': [-1*(atomSpacing*cos72), -1*(2*atomSpacing+(atomSpacing*sin72)+(atomSpacing*cos54)+(atomSpacing*cos60))],

           'boundAtoms': [['CG', 'CD1'], ['CG', 'CD2'], ['CD2', 'CE3'], ['CD2', 'CE2'],
                          ['CD1', 'NE1'], ['CE2', 'CZ2'], ['CE3', 'CZ3'], ['CZ3', 'CH2'],
                          ['CZ2', 'CH2'], ['NE1', 'CE2']]
          },

  'PRO':  {
           'CB': [atomSpacing*cos72, -1*(atomSpacing*sin72)+atomSpacing],
           'CG': [-0.5*atomSpacing, -1*atomSpacing/(2*tan36)],
           'CD': [-1*(atomSpacing+(atomSpacing*cos72)), -1*(atomSpacing*sin72)+atomSpacing],
          }
}

# if residueType == 'ALA':
#       hb = self._createGuiNmrAtom('HB%', (cbAtom.x(), cbAtom.y()-self.atomSpacing))
#       self.scene.addItem(hb)
#       self._addConnectingLine(hb, cbAtom, 'white', 1.0, 0.0)
#
#     if residueType == 'CYS':
#       sg = self._createGuiNmrAtom('SG', (cbAtom.x(), cbAtom.y()-self.atomSpacing))
#       hg = self._createGuiNmrAtom('HG', (cbAtom.x(), cbAtom.y()-(2*self.atomSpacing)))
#       self.scene.addItem(sg)
#       self.scene.addItem(hg)
#       self._addConnectingLine(sg, cbAtom, 'white', 1.0, 0.0)
#       self._addConnectingLine(sg, hg, 'white', 1.0, 0.0)
#
#     if residueType == 'ASP':
#       hb2 = self._createGuiNmrAtom('HBx', (cbAtom.x()-(self.atomSpacing*0.75), cbAtom.y()))
#       hb3 = self._createGuiNmrAtom('HBy', (cbAtom.x()+(self.atomSpacing*0.75), cbAtom.y()))
#       cg = self._createGuiNmrAtom('CG', (cbAtom.x(), cbAtom.y()-self.atomSpacing))
#       self.scene.addItem(hb2)
#       self.scene.addItem(hb3)
#       self.scene.addItem(cg)
#       self._addConnectingLine(hb2, cbAtom, 'white', 1.0, 0.0)
#       self._addConnectingLine(hb3, cbAtom, 'white', 1.0, 0.0)
#       self._addConnectingLine(cg, cbAtom, 'white', 1.0, 0.0)
#
#     if residueType == 'GLU':
#       hb2 = self._createGuiNmrAtom('HBx', (cbAtom.x()-(self.atomSpacing*0.75), cbAtom.y()))
#       hb3 = self._createGuiNmrAtom('HBy', (cbAtom.x()+(self.atomSpacing*0.75), cbAtom.y()))
#       cg = self._createGuiNmrAtom('CG', (cbAtom.x(), cbAtom.y()-self.atomSpacing))
#       hg2 = self._createGuiNmrAtom('HBx', (cg.x()-(self.atomSpacing*0.75), cbAtom.y()))
#       hg3 = self._createGuiNmrAtom('HBy', (cg.x()+(self.atomSpacing*0.75), cbAtom.y()))
#       cd = self._createGuiNmrAtom('CD', (cbAtom.x(), cbAtom.y()-(2*self.atomSpacing)))
#       self.scene.addItem(hb2)
#       self.scene.addItem(hb3)
#       self.scene.addItem(cg)
#       self.scene.addItem(hg2)
#       self.scene.addItem(hg3)
#       self.scene.addItem(cd)
#       self._addConnectingLine(hb2, cbAtom, 'white', 1.0, 0.0)
#       self._addConnectingLine(hb3, cbAtom, 'white', 1.0, 0.0)
#       self._addConnectingLine(cg, hg2, 'white', 1.0, 0.0)
#       self._addConnectingLine(cg, hg3, 'white', 1.0, 0.0)
#       self._addConnectingLine(cg, cd, 'white', 1.0, 0.0)
