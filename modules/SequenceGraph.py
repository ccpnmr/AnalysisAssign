"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2020"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2020-03-17 00:13:56 +0000 (Tue, March 17, 2020) $"
__version__ = "$Revision: 3.0.1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2016-05-23 10:02:47 +0100 (Thu, 26 May 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

import json
import typing
import numpy as np
from functools import partial
from PyQt5 import QtGui, QtWidgets, QtCore
from collections import OrderedDict
from ccpn.util.OrderedSet import OrderedSet
from contextlib import contextmanager
from ccpn.core.lib.Pid import Pid
from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.Peak import Peak
from ccpn.core.Spectrum import Spectrum
from ccpn.core.lib.AssignmentLib import getNmrResiduePrediction
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.lib.CallBack import CallBack
from ccpn.ui.gui.lib.Strip import navigateToNmrResidueInDisplay, _getCurrentZoomRatio
from ccpn.ui.gui.lib.mouseEvents import makeDragEvent
from ccpn.ui.gui.widgets.Widget import Widget
# from ccpn.ui.gui.guiSettings import textFontSmall, textFontSmallBold, textFont
from ccpn.ui.gui.guiSettings import getColours
from ccpn.ui.gui.guiSettings import GUINMRATOM_NOTSELECTED, GUINMRATOM_SELECTED, \
    GUINMRRESIDUE, SEQUENCEGRAPHMODULE_LINE, SEQUENCEGRAPHMODULE_TEXT
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Menu import Menu
from ccpn.ui.gui.widgets.Icon import Icon
from ccpn.ui.gui.widgets.ToolBar import ToolBar
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget, ListCompoundWidget
from ccpn.ui.gui.widgets.PulldownListsForObjects import NmrChainPulldown
from ccpn.core.NmrChain import NmrChain
from ccpn.util.Common import makeIterableList
from ccpn.util.Constants import ccpnmrJsonData
from ccpn.util.Logging import getLogger
from ccpn.ui.gui.widgets.MessageDialog import showWarning, progressManager
from ccpn.ui.gui.widgets.Splitter import Splitter
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.modules.SequenceModule import SequenceModule
from ccpn.ui.gui.widgets.SettingsWidgets import SequenceGraphSettings
from ccpn.core.lib.AssignmentLib import getSpinSystemsLocation
from ccpn.core.lib.ContextManagers import notificationEchoBlocking, catchExceptions, undoBlock
from ccpnc.clibrary import Clibrary


_getNmrIndex = Clibrary.getNmrResidueIndex

logger = getLogger()
ALL = '<all>'


#==========================================================================================
# GuiNmrAtom
#==========================================================================================

class GuiNmrAtom(QtWidgets.QGraphicsTextItem):
    """
    A graphical object specifying the position and name of an atom when created by the Assigner.
    Can be linked to a Nmr Atom.
    """

    def __init__(self, mainWindow, text, pos=None, nmrAtom=None):
        super().__init__()

        self.setPlainText(text)
        self.setPos(QtCore.QPointF((pos[0] - self.boundingRect().x()), (pos[1] - self.boundingRect().y())))

        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current
        self.nmrAtom = nmrAtom

        self.connectedAtoms = 0
        self.connectedList = {}  # maintain connectivity between guiNmrAtoms
        # so that lines do not overlap

        self.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)

        # set the highlight colour for dragging to chain
        self.colours = getColours()
        if self.isSelected:
            self.setDefaultTextColor(QtGui.QColor(self.colours[GUINMRATOM_SELECTED]))
        else:
            self.setDefaultTextColor(QtGui.QColor(self.colours[GUINMRATOM_NOTSELECTED]))

    def mouseDoubleClickEvent(self, event):
        """CCPN INTERNAL - re-implementation of double click event
        """
        pass

    def mousePressEvent(self, event):
        """CCPN INTERNAL - re-implementation of mouse press event
        """
        if self.nmrAtom is not None:
            self.current.nmrAtom = self.nmrAtom
            self.current.nmrResidue = self.nmrAtom.nmrResidue
            event.accept()

    def _raiseContextMenu(self, event: QtGui.QMouseEvent):
        """Creates and raises a context menu enabling items to be disconnected
        """
        from ccpn.ui.gui.widgets.Menu import Menu
        from functools import partial

        contextMenu = Menu('', event.widget(), isFloatWidget=True)
        contextMenu.addAction('deassign all Peaks', partial(self._deassignAllPeaksFromNmrAtom))
        cursor = QtGui.QCursor()
        contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
        contextMenu.exec()

    def addConnectedList(self, connectedAtom):
        """maintain number of links between adjacent nmrAtoms.
        """
        keyVal = connectedAtom
        if keyVal in self.connectedList:
            self.connectedList[keyVal] += 1
        else:
            self.connectedList[keyVal] = 1

    def removeConnectedList(self, connectedAtom):
        """maintain number of links between adjacent nmrAtoms.
        """
        keyVal = connectedAtom
        if keyVal in self.connectedList:
            self.connectedList[keyVal] -= 1
        else:
            raise RuntimeError('Connection does not exist')

    def getConnectedList(self, connectedAtom):
        """Retrieve number of links between adjacent nmrAtoms.
        """
        keyVal = connectedAtom
        if keyVal in self.connectedList:
            return self.connectedList[keyVal]
        else:
            return 0

    def clearConnectedList(self):
        """Clear all connections for this guiNmrAtom but do not delete.
        """
        for keyVal in self.connectedList:
            keyVal.connectedList[self] = 0
            self.connectedList[keyVal] = 0

    def deleteConnectedList(self):
        """Delete all connections for this guiNmrAtom.
        """
        for keyVal in self.connectedList:
            del keyVal.connectedList[self]
        self.connectedList = {}


#==========================================================================================
# GuiNmrResidue
#==========================================================================================

class GuiNmrResidue(QtWidgets.QGraphicsTextItem):
    """
    Object linking residues displayed in Assigner and Nmr Residues. Contains functionality for drag and
    drop assignment in conjunction with the Sequence Module.
    """

    def __init__(self, parent, nmrResidue, caAtom):

        super().__init__()
        self.setPlainText(nmrResidue.id)

        self.mainWindow = parent.mainWindow
        self.application = self.mainWindow.application
        self.project = self.mainWindow.project
        self.current = self.mainWindow.application.current

        self.setFont(self.mainWindow.application._fontSettings.textFontSmall)
        self.colours = getColours()
        self.setDefaultTextColor(QtGui.QColor(self.colours[GUINMRRESIDUE]))

        self.setPos(caAtom.x() - caAtom.boundingRect().width() / 2, caAtom.y() + 30)

        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable)
        self._parent = parent
        self.nmrResidue = nmrResidue

        self.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)

    def _update(self):
        self.setPlainText(self.nmrResidue.id)

    def _mouseMoveEvent(self, event):
        """create a drag item if left button pressed
        """
        if (event.buttons() == QtCore.Qt.LeftButton):

            nmrItem = self

            if nmrItem:
                makeDragEvent(event.widget(), {'pids': [nmrItem.nmrResidue.pid]}, self.toPlainText(), action=QtCore.Qt.MoveAction)

    def _mousePressEvent(self, event):
        self.current.nmrResidue = self.nmrResidue
        self.setSelected(True)

    def _mouseDoubleClickEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self._parent.showNmrResidue(self)
            event.accept()
        else:
            super().mouseDoubleClickEvent(event)

    # def _raiseLineMenu(self, scene, event):
    #     """Creates and raises a context menu enabling items to be disconnected
    #     """
    #     from ccpn.ui.gui.widgets.Menu import Menu
    #     from functools import partial
    #
    #     cursor = QtGui.QCursor()
    #     contextMenu = Menu('', event.widget(), isFloatWidget=True)
    #     pressed = self.scene.mouseGrabberItem()
    #
    #     if self.selectedLine:
    #         thisLine = self.selectedLine
    #         contextMenu.addAction('deassign nmrAtoms from Peak: %s' % str(thisLine._peak.id))
    #         contextMenu.addSeparator()
    #         if thisLine._peak:
    #
    #             # skip if nothing connected
    #             if not thisLine._peak.assignedNmrAtoms:
    #                 return
    #
    #             # add the nmrAtoms to the menu
    #             for nmrAtomList in thisLine._peak.assignedNmrAtoms:
    #                 for nmrAtom in nmrAtomList:
    #                     if nmrAtom:
    #                         contextMenu.addAction(nmrAtom.id, partial(self._deassignPeak, thisLine._peak, nmrAtom))


#==========================================================================================
# Assignment line
#==========================================================================================

class AssignmentLine(QtWidgets.QGraphicsLineItem):
    """
    Object to create lines between GuiNmrAtoms with specific style, width, colour and displacement.
    Displacement allows multiplet peak lines to be shown as adjacent lines between the same nmrAtom.
    """

    def __init__(self, x1, y1, x2, y2, colour, width,
                 parent=None, style=None, peak=None, guiAtom1: GuiNmrAtom = None, guiAtom2: GuiNmrAtom = None, displacement=None):
        super().__init__()

        # set the pen colour and style
        self.pen = QtGui.QPen()
        self.pen.setColor(QtGui.QColor(colour))
        self.pen.setCosmetic(True)
        self.pen.setWidth(width)

        if style == 'dash':
            self.pen.setStyle(QtCore.Qt.DotLine)
        self.setPen(self.pen)
        self.setLine(x1, y1, x2, y2)
        self._parent = parent

        # store the peak and the guiNmrAtoms that the line connects
        self._peak = peak
        self.guiAtom1 = guiAtom1
        self.guiAtom2 = guiAtom2
        self.displacement = displacement

        # enable hovering so the current line can be set
        self.setAcceptedMouseButtons(QtCore.Qt.RightButton)
        self.setAcceptHoverEvents(True)

    def updateEndPoints(self):
        """Update the endPoints of the line to point. Co-ordinates are relative to the group to
        which the graphicsItem belongs, in this case the guiNmrResidue group. GuiNmrResidue group is the top level relative to the scene.
        """
        guiAtom1 = self.guiAtom1
        guiAtom2 = self.guiAtom2
        residue1 = guiAtom1.guiNmrResidueGroup
        residue2 = guiAtom2.guiNmrResidueGroup

        rx1 = residue1.x()
        ry1 = residue1.y()
        rx2 = residue2.x()
        ry2 = residue2.y()

        atom1Rect = guiAtom1.boundingRect()
        atom2Rect = guiAtom2.boundingRect()
        w1 = atom1Rect.width()
        h1 = atom1Rect.height()
        w2 = atom2Rect.width()
        h2 = atom2Rect.height()

        x1 = guiAtom1.x()  # + rx1
        y1 = guiAtom1.y()  # + ry1
        x2 = guiAtom2.x() + (rx2 - rx1)
        y2 = guiAtom2.y() + (ry2 - ry1)

        dx = x2 - x1
        dy = y2 - y1
        length = 2.0 * pow(dx * dx + dy * dy, 0.5)
        if self.displacement is not None:
            count = (guiAtom1.connectedList[guiAtom2] - 1) // 2
            disp = 6.0 * (self.displacement - count) / length
        else:
            disp = 0.0
        offsetX = dy * disp
        offsetY = -dx * disp
        kx1 = (w1 * dx) / length  # shorten the lines along length
        ky1 = (h1 * dy) / length
        kx2 = (w2 * dx) / length
        ky2 = (h2 * dy) / length

        xOff1 = w1 / 2.0  # offset to centre of bounding box
        yOff1 = h1 / 2.0
        xOff2 = w2 / 2.0
        yOff2 = h2 / 2.0

        x1 += xOff1 + kx1 + offsetX
        y1 += yOff1 + ky1 + offsetY
        x2 += xOff2 - kx2 + offsetX
        y2 += yOff2 - ky2 + offsetY

        self.setLine(x1, y1, x2, y2)

    def paint(self, painter, option, widget):
        """Automatically update the end-points of the assignment lines to point to the correct guiNmrAtoms
        """
        self.updateEndPoints()
        super().paint(painter, option, widget)

    def hoverEnterEvent(self, event):
        self._parent.selectedLine = self

    def hoverLeaveEvent(self, event):
        self._parent.selectedLine = None

    def mousePressEvent(self, event):
        """Required to respond to mouse press events.
        """
        pass


#==========================================================================================
# GuiNmrResidueGroup
#==========================================================================================

class GuiNmrResidueGroup(QtWidgets.QGraphicsItemGroup):
    """
    Group item to group all nmrAtoms/connecting lines of nmrResidue
    """

    def __init__(self, parent, nmrResidue, caAtom, pos):
        super().__init__()

        self.mainWindow = parent.mainWindow
        self.application = self.mainWindow.application
        self.project = self.mainWindow.project
        self.current = self.mainWindow.application.current
        self.nmrResidue = nmrResidue
        self._parent = parent
        self.setPos(QtCore.QPointF(pos, 0.0))
        self.crossChainCount = None
        self.crossChainResidue = None

        self.nmrResidueLabel = GuiNmrResidue(parent, nmrResidue, caAtom)
        self.addToGroup(self.nmrResidueLabel)

    def mousePressEvent(self, event):
        self.nmrResidueLabel._mousePressEvent(event)

    def mouseMoveEvent(self, event):
        self.nmrResidueLabel._mouseMoveEvent(event)

    def mouseDoubleClickEvent(self, event):
        self.nmrResidueLabel._mouseDoubleClickEvent(event)


#==========================================================================================
# NmrResidueList
#==========================================================================================

class NmrResidueList(object):
    """
    A Class to hold the information about each gui object in the scene
    """

    def __init__(self, mainWindow, settingsWidget, lineColour, textColour, atomSpacing, scene, module):
        """Initialise the new object.
        """
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.project
        self.current = mainWindow.application.current

        self.reset()
        self._SGwidget = settingsWidget
        self._lineColour = lineColour
        self._textColour = textColour
        self._atomSpacing = atomSpacing
        self._scene = scene
        self._module = module
        self.nmrChain = None

    def reset(self):
        self.residueCount = 0
        self.direction = None
        self.selectedStretch = []
        self.selectedLine = None

        # change to an orderedDict so that more nmrChains can be seen in the future
        self.nmrChains = OrderedDict()  # referenced by nmrChain

        # store all visible gui items
        self.guiNmrResidues = OrderedDict()  # referenced by nmrResidue
        self.guiNmrAtoms = OrderedDict()  # referenced by nmrAtom
        self.guiGhostNmrResidues = OrderedDict()  # referenced by nmrResidue

        self.guiNmrAtomsFromNmrResidue = OrderedDict()  # referenced by nmrResidue -> list(guiNmrAtoms)

        self.ghostList = {}

        self.connectingLines = {}  # referenced by peak?
        self.assignmentLines = {}

        self.nmrChain = None  # current active nmrChain

    def size(self, nmrChainId):
        """return the number of elements in the list nmrChain.
        """
        return len(self.nmrChains[nmrChainId]) if nmrChainId in self.nmrChains else None

    #==========================================================================================

    def getIndexNmrResidue(self, nmrResidue):
        """get the index in the nmrResidueList of the required nmrResidue.
        """
        for id, nmrCh in self.nmrChains.items():
            if nmrResidue in nmrCh:
                return nmrCh.index(nmrResidue)

        # resList = [item[0] for item in self.allNmrResidues]
        # if nmrResidue in resList:
        #     return resList.index(nmrResidue)

    def deleteNmrResidue(self, nmrResidue):
        """get the index in the nmrResidueList of the required nmrResidue.
        """
        for id, nmrCh in self.nmrChains.items():
            if nmrResidue in nmrCh:
                nmrCh.remove(nmrResidue)
                print('>>>removing from nmrChains')

    # def deleteNmrResidueIndex(self, index):
    #     """insert into the list as a tuple (obj, dict).
    #     """
    #     return self.allNmrResidues.pop(index)
    #
    # def _appendNmrResiduePair(self, nmrResidue, guiAtoms):
    #     """insert into the list as a tuple (obj, dict).
    #     """
    #     self.allNmrResidues.append((nmrResidue, guiAtoms))

    def _insertNmrResiduePair(self, nmrChainId, index, nmrResidue, guiAtoms):
        """insert into the list as a tuple (obj, dict).
        """
        if nmrChainId not in self.nmrChains:
            self.nmrChains[nmrChainId] = [nmrResidue]
        else:
            self.nmrChains[nmrChainId].insert(index, nmrResidue)

        # self.allNmrResidues.insert(index, (nmrResidue, guiAtoms))
        # self.allNmrResidues[index] = (nmrResidue, guiAtoms)

    # def _getNmrResiduePair(self, index):
    #     """insert into the list as a tuple.
    #     """
    #     val = self.allNmrResidues[index]
    #     return val[0], val[1]
    #
    # def _setNmrResiduePair(self, index, nmrResidue, guiAtoms):
    #     """insert into the list as a tuple.
    #     """
    #     self.allNmrResidues[index] = (nmrResidue, guiAtoms)

    # def _updateAllNmrResidueSize(self):
    #     """Change the size of this list based on whether mainNmrResidues has changed
    #     """
    #     if self.nmrChain:
    #         mainNmrResidues = self.nmrChain.mainNmrResidues
    #         if self.nmrChain not in self.allNmrResidues:
    #             self.allNmrResidues[self.nmrChain] = [None] * len(mainNmrResidues)
    #
    #         else:
    #             resList = self.allNmrResidues[self.nmrChain]
    #             lenAll = len(resList)
    #             if lenAll == 0:
    #                 # create a new list
    #                 self.allNmrResidues[self.nmrChain] = [None] * len(mainNmrResidues)
    #
    #             else:
    #                 # update the list size copying the original data into it
    #                 ii = _getNmrIndex(resList[0])
    #                 newList = [None] * len(mainNmrResidues)
    #                 for res in self.allNmrResidues[self.nmrChain]:
    #                     index = _getNmrIndex(res)
    #                     if index and index < len(mainNmrResidues):
    #                         newList[_getNmrIndex(res)] = res

    #==========================================================================================

    def addNmrResidue(self, nmrChainId, nmrResidue, index=0, _insertNmrRes=True):
        """Add a new nmrResidue at the required position.
        """
        self._addNmrResidue(nmrChainId, nmrResidue, nmrResidueIndex=index, lineList=self.connectingLines, _insertNmrRes=_insertNmrRes)

    def _addNmrResidue(self, nmrChainId, nmrResidue, nmrResidueIndex=0, lineList=None, _insertNmrRes=True):
        """Takes an Nmr Residue, and adds a residue to the sequence graph
        corresponding to the Nmr Residue at the required index.
        Nmr Residue name displayed beneath CA of residue drawn and residue type predictions displayed
        beneath Nmr Residue name
        """

        # check whether the guiNmrResidue already exists, and skip
        if nmrResidue in self.guiNmrResidues:
            return

        guiAtoms = {}
        # mainNmrResidues = nmrResidue.nmrChain.mainNmrResidues

        if atomSpacing:
            self.atomSpacing = atomSpacing
        nmrAtomNames = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
        backboneAtoms = DEFAULT_RESIDUE_ATOMS.copy()

        if nmrResidue.residueType == 'GLY':
            # GLY doesn't have CB
            del backboneAtoms['CB']

        # add extra guiAtoms to the list based on the backbone atoms dict above - to self.guiNmrAtomDict[nmrAtoms]
        self.addBackboneAtoms(nmrResidue, backboneAtoms, nmrAtomNames, guiAtoms)

        # add the sideChain atoms to self.guiNmrAtomDict[nmrAtoms]
        if self._SGwidget.checkBoxes['showSideChain']['checkBox'].isChecked():
            if 'CB' in backboneAtoms and nmrResidue.residueType:
                cbAtom = guiAtoms['CB']
                self.addSideChainAtoms(nmrResidue, cbAtom, nmrAtomNames, guiAtoms)

        # store the reference to all gui NmrAtoms from nmrResidue
        self.guiNmrAtomsFromNmrResidue[nmrResidue] = guiAtoms

        # self._updateAllNmrResidueSize()

        # # add the new nmrResidue to the current list
        # if not self.allNmrResidues:
        #
        #     # append to the list self.allNmrResidues
        #     self.allNmrResidues = [None] * len(mainNmrResidues)
        #     self.allNmrResidues[mainNmrResidues.index(nmrResidue)] = (nmrResidue, guiAtoms)
        #     # self._appendNmrResiduePair(nmrResidue, guiAtoms)
        #
        # else:

        # insert into the correct nmrChain
        if _insertNmrRes:
            self._insertNmrResiduePair(nmrChainId, nmrResidueIndex, nmrResidue, guiAtoms)

        # index = mainNmrResidues.index(nmrResidue)
        # self.allNmrResidues[nmrResidue][index] = (nmrResidue, guiAtoms)

        # compile a gui group for this nmrResidue - self.guiNmrResidues[nmrResidue
        newGuiResidueGroup = self._assembleGroupResidue(nmrResidue, guiAtoms)

        return newGuiResidueGroup

    def _addGhostResidue(self, nmrResidueCon1: NmrResidue,
                         guiRef: GuiNmrAtom,
                         nmrResidueCon0: NmrResidue,
                         name1: str, name0: str,
                         offsetAdjust,
                         atomSpacing=None, lineList=None):
        """Takes an Nmr Residue and a direction, either '-1 or '+1', and adds a residue to the sequence graph
        corresponding to the Nmr Residue.
        Nmr Residue name displayed beneath CA of residue drawn and residue type predictions displayed
        beneath Nmr Residue name
        """

        # need to keep a list of the atoms that have been added so don't repeat
        count = 0
        if nmrResidueCon0 in self.ghostList:
            count = len(self.ghostList[nmrResidueCon0])
            if nmrResidueCon1 in self.ghostList[nmrResidueCon0]:
                # already exists in the dict so exit
                return
        else:
            self.ghostList[nmrResidueCon0] = ()

        nmrResidue = nmrResidueCon1
        atoms = {}
        if atomSpacing:
            self.atomSpacing = atomSpacing
        nmrAtoms = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
        residueAtoms = DEFAULT_RESIDUE_ATOMS.copy()

        if nmrResidue.residueType == 'GLY':
            del residueAtoms['CB']

        for k, v in residueAtoms.items():
            if k in nmrAtoms:
                nmrAtom = nmrResidue.fetchNmrAtom(name=k)
            else:
                nmrAtom = None
            atoms[k] = self._createGhostGuiNmrAtom(k, v, nmrAtom)

        newGuiResidueGroup = self._assembleGhostResidue(nmrResidue, atoms, lineList=lineList)
        newGuiResidueGroup.crossChainCount = count
        newGuiResidueGroup.crossChainResidue = nmrResidueCon0

        self.ghostList[nmrResidueCon0] += (nmrResidueCon1,)
        return atoms

    #==========================================================================================

    def _createGuiNmrAtom(self, atomType: str, position: tuple, nmrAtom: NmrAtom = None) -> GuiNmrAtom:
        """Creates a GuiNmrAtom specified by the atomType and graphical position supplied.
        GuiNmrAtom can be linked to an NmrAtom by supplying it to the function.
        """
        guiAtom = GuiNmrAtom(self.mainWindow, text=atomType, pos=position, nmrAtom=nmrAtom)
        if nmrAtom:
            # only add to the dict if the nmrAtom exists
            self.guiNmrAtoms[nmrAtom] = guiAtom
        return guiAtom

    def _createGhostGuiNmrAtom(self, atomType: str, position: tuple, nmrAtom: NmrAtom = None) -> GuiNmrAtom:
        """Creates a GuiNmrAtom specified by the atomType and graphical position supplied.
        GuiNmrAtom can be linked to an NmrAtom by supplying it to the function.
        """
        guiAtom = GuiNmrAtom(self.mainWindow, text=atomType, pos=position, nmrAtom=nmrAtom)
        if nmrAtom:
            # only add to the dict if the nmrAtom exists
            self.guiNmrAtoms[nmrAtom] = guiAtom
        return guiAtom

    #==========================================================================================

    def _nmrAtomID(self, nmrAtom):
        """Create an id string to retrieve guiNmrAtoms.
        """
        return nmrAtom.nmrResidue.id + nmrAtom.name

    def addBackboneAtoms(self, nmrResidue, backboneAtoms, atomNames, guiAtoms):
        """add the backbone atoms.
        """
        for k, v in backboneAtoms.items():
            if k in atomNames:
                nmrAtom = nmrResidue.fetchNmrAtom(name=k)
            else:
                nmrAtom = None
            guiAtoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)

    def addSideChainAtoms(self, nmrResidue, cbAtom, atomNames, guiAtoms):
        """Add the sideChain atoms above the backbone line.
        """
        # residue = {}
        for k, v in ATOM_POSITION_DICT[nmrResidue.residueType].items():
            if k != 'boundAtoms':
                position = [cbAtom.x() + v[0], cbAtom.y() + v[1]]
                if k in atomNames:
                    nmrAtom = nmrResidue.fetchNmrAtom(name=k)
                else:
                    nmrAtom = None
                newAtom = self._createGuiNmrAtom(k, position, nmrAtom)

                # self.scene.addItem(newAtom)
                # residue[k] = newAtom
                guiAtoms[k] = newAtom

    #==========================================================================================

    def _assembleGroupResidue(self, nmrResidue: NmrResidue, guiAtoms: typing.Dict[str, GuiNmrAtom]):
        """Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
        creates a graphical representation of a residue in the assigner
        """
        guiResidueGroup = GuiNmrResidueGroup(self._module, nmrResidue, guiAtoms['CA'], 0)

        # insert into the gui list and add to the scene
        self.guiNmrResidues[nmrResidue] = guiResidueGroup
        self._scene.addItem(guiResidueGroup)

        # add the atoms to the group and set the reverse link
        for item in guiAtoms.values():
            guiResidueGroup.addToGroup(item)
            item.guiNmrResidueGroup = guiResidueGroup

        # add the backbone lines - sidechain will be added in the future
        if "CB" in guiAtoms:
            self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['CA'], guiAtoms['CB'],
                                           self._lineColour, 1.0, lineList=self.connectingLines, lineId=nmrResidue)

        if "H" in guiAtoms and nmrResidue.residueType != 'PRO':
            self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['H'], guiAtoms['N'],
                                           self._lineColour, 1.0, lineList=self.connectingLines, lineId=nmrResidue)

        self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['N'], guiAtoms['CA'],
                                       self._lineColour, 1.0, lineList=self.connectingLines, lineId=nmrResidue)
        self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['C'], guiAtoms['CA'],
                                       self._lineColour, 1.0, lineList=self.connectingLines, lineId=nmrResidue)

        self._addGroupResiduePredictions(guiResidueGroup, nmrResidue, guiAtoms['CA'])

        return guiResidueGroup

    def _assembleGhostResidue(self, nmrResidue: NmrResidue, guiAtoms: typing.Dict[str, GuiNmrAtom], lineList=None):
        """Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
        creates a graphical representation of a residue in the assigner
        """
        guiResidueGroup = GuiNmrResidueGroup(self._module, nmrResidue, guiAtoms['CA'], 0)
        self.guiGhostNmrResidues[nmrResidue] = guiResidueGroup
        self._scene.addItem(guiResidueGroup)

        # add the atoms to the group and set the reverse link
        for item in guiAtoms.values():
            guiResidueGroup.addToGroup(item)
            item.guiNmrResidueGroup = guiResidueGroup

        # add the backbone lines - sidechain will be added in the future
        if "CB" in list(guiAtoms.keys()):
            self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['CA'], guiAtoms['CB'],
                                           self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        if "H" in list(guiAtoms.keys()) and nmrResidue.residueType != 'PRO':
            self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['H'], guiAtoms['N'],
                                           self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['N'], guiAtoms['CA'],
                                       self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
        self._addConnectingLineToGroup(guiResidueGroup, guiAtoms['C'], guiAtoms['CA'],
                                       self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        return guiResidueGroup

    #==========================================================================================

    def _addConnectingLineToGroup(self, group: GuiNmrResidueGroup, guiAtom1: GuiNmrAtom, guiAtom2: GuiNmrAtom,
                                  colour: str, width: float, displacement: float = None, style: str = None,
                                  peak: Peak = None, lineList=None, lineId=typing.Any):
        """Adds a line between two GuiNmrAtoms using the width, colour, displacement and style specified.
        """
        itemKey = id(lineId)
        if itemKey not in lineList:
            lineList[itemKey] = []

        for line in lineList[itemKey]:
            if (line.guiAtom1 == guiAtom1 and line.guiAtom2 == guiAtom2 and line._peak == peak) or \
                    (line.guiAtom2 == guiAtom1 and line.guiAtom1 == guiAtom2 and line._peak == peak):
                break
        else:
            # add a newLine only if it doesn't already exists - may cause some gaps in the displacements
            newLine = AssignmentLine(0, 0, 0, 0, colour, width,
                                     parent=self, style=style, peak=peak,
                                     guiAtom1=guiAtom1, guiAtom2=guiAtom2, displacement=displacement)

            newLine.setParentItem(group)

            lineList[itemKey].append(newLine)
            return newLine

        return None

    def addConnectionsBetweenGroups(self, nmrChainId):
        """Add the connections between the groups.
        """
        # connectsNeeded = self._module.nmrResiduesCheckBox.isChecked()     # why was this here?

        # mainNmrResidues = self.nmrChain.mainNmrResidues
        # get the list of nmrResidues in the required nmrChain referenced by nmrChainId
        mainNmrResidues = self.nmrChains[nmrChainId] if nmrChainId in self.nmrChains else []  #[resPair[0] for resPair in self.nmrChains[nmrChainId]]

        # iterate through the adjacent pairs
        for prevRes, thisRes in zip(mainNmrResidues[:-1], mainNmrResidues[1:]):

            # add the connection the minus residue and point to the right - may need to change for +1 residues
            # prevRes, prevGuiAtoms = prev
            # thisRes, thisGuiAtoms = this
            prevGuiAtoms = self.guiNmrAtomsFromNmrResidue[prevRes]
            thisGuiAtoms = self.guiNmrAtomsFromNmrResidue[thisRes]

            if (prevRes.nextNmrResidue and prevRes.nextNmrResidue is thisRes):  # and connectsNeeded:
                # connect from this 'N' to the previous 'C'
                self._addConnectingLineToGroup(self.guiNmrResidues[prevRes],
                                               prevGuiAtoms['C'], thisGuiAtoms['N'],
                                               self._lineColour, 1.0, lineList=self.connectingLines,
                                               lineId=thisRes)

    #==========================================================================================

    def _addGroupResiduePredictions(self, guiResidueGroup: GuiNmrResidueGroup, nmrResidue: NmrResidue, caAtom: GuiNmrAtom):
        """Gets predictions for residue type based on BMRB statistics and determines label positions
        based on caAtom position.
        """
        predictions = list(set(map(tuple, (getNmrResiduePrediction(nmrResidue, self.project.chemicalShiftLists[0])))))
        predictions.sort(key=lambda a: float(a[1][:-1]), reverse=True)
        for prediction in predictions:
            predictionLabel = QtWidgets.QGraphicsTextItem()
            predictionLabel.setPlainText(prediction[0] + ' ' + prediction[1])
            predictionLabel.setDefaultTextColor(QtGui.QColor(self._textColour))
            predictionLabel.setFont(self.mainWindow.application._fontSettings.textFontSmallBold)
            predictionLabel.setPos(caAtom.x() - caAtom.boundingRect().width() / 2,
                                   caAtom.y() + (30 * (predictions.index(prediction) + 2)))

            guiResidueGroup.addToGroup(predictionLabel)

    #==========================================================================================

    def updateMainChainPositions(self, nmrChainId):
        """Update the positions of the group in the scene.
        May need to set vertical positions when more groups are allowed
        """
        if nmrChainId not in self.nmrChains:
            return

        # get the list of nmrResidues in the required nmrChain referenced by nmrChainId
        # mainNmrResidues = self.nmrChains[nmrChainId]  #[resPair[0] for resPair in self.nmrChains[nmrChainId]]

        mainNmrResidues = [nmrResidue for nmrResidue in self.nmrChains[nmrChainId] if nmrResidue.nmrChain.pid == nmrChainId and
                           not nmrResidue._flaggedForDelete]

        for ii, nmrResidue in enumerate(mainNmrResidues):
            # nmrResidue, guiAtoms = item
            # guiAtoms = self.guiNmrAtomsFromNmrResidue[nmrResidue]

            if nmrResidue in self.guiNmrResidues:
                guiItem = self.guiNmrResidues[nmrResidue]
                guiItem.setPos(QtCore.QPointF(ii * self.atomSpacing * 3.0, 0.0))

    def updateConnectedChainPositions(self, nmrChainId):
        """Update the positions of the groups in the scene.
        """

        # update crossChainResidue positions
        for res in self.guiGhostNmrResidues.values():
            if res.crossChainResidue and res.crossChainResidue in self.guiNmrResidues:
                link = self.guiNmrResidues[res.crossChainResidue]
                count = res.crossChainCount

                newPosx = link.x()
                newPosy = link.y()
                res.setPos(QtCore.QPointF(newPosx + (count * 0.5 - 1.0) * self.atomSpacing,
                                          newPosy + (count * 2.5 + 5.0) * self.atomSpacing))

    #==========================================================================================

    def _addAllPeakAssignments(self, nmrChainId):
        """Add all the peak assignments to the scene.
        """
        if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():

            # # create a set of sets ordered by spectra for active lines
            # self.LOCALinterResidueAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
            # self.LOCALinterChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
            # self.LOCALcrossChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())

            # mainNmrResidues = self.nmrChain.mainNmrResidues
            # get the list of nmrResidues in the required nmrChain referenced by nmrChainId
            mainNmrResidues = self.nmrChains[nmrChainId] if nmrChainId in self.nmrChains else []  #[resPair[0] for resPair in self.nmrChains[nmrChainId]]

            for nmrResidue in mainNmrResidues:
                # guiRes = self.guiNmrAtomsFromNmrResidue[nmrResidue]

                # add the internally connected Lines
                internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue)

                self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)

                # add assignment lines and create ghost nmrResidues if required
                self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                            self.assignmentLines, self.connectingLines)

            # self._addPeakAssignmentLinesToGroup(self.LOCALinterResidueAtomPairing, self.assignmentLines)
            # self._addPeakAssignmentLinesToGroup(self.LOCALinterChainAtomPairing, self.assignmentLines)
            # self._addPeakAssignmentLinesToAdjacentGroup(self._module.nmrChain, self.LOCALcrossChainAtomPairing,
            #                                             self.assignmentLines, self.connectingLines)

    #==========================================================================================

    def _addPeakAssignmentLinesToGroup(self, assignments, lineList):
        """Add the local peak assignments to each guiNmrResidue.
        """
        for specAssignments in assignments.values():
            for nmrAtomPair in specAssignments:

                guiNmrAtomPair = (self.guiNmrAtoms.get(nmrAtomPair[0]),
                                  self.guiNmrAtoms.get(nmrAtomPair[1]),
                                  nmrAtomPair[2]
                                  )

                # skip if not defined
                if None in guiNmrAtomPair:
                    continue

                # get the peak and the spectrum
                peak = guiNmrAtomPair[2]
                spectrum = peak.peakList.spectrum
                displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])

                # add the internal line to the guiNmrResidueGroup, should now move when group is moved
                guiNmrResidue = self.guiNmrResidues[guiNmrAtomPair[0].nmrAtom.nmrResidue]

                self._addConnectingLineToGroup(guiNmrResidue,
                                               guiNmrAtomPair[0],
                                               guiNmrAtomPair[1],
                                               spectrum.positiveContourColour,
                                               2.0, displacement=displacement,
                                               peak=peak, lineList=lineList, lineId=peak)

                # update displacements for both guiNmrAtoms
                guiNmrAtomPair[0].addConnectedList(guiNmrAtomPair[1])
                guiNmrAtomPair[1].addConnectedList(guiNmrAtomPair[0])

    def _addPeakAssignmentLinesToAdjacentGroup(self, nmrResidue, assignments, peaklineList, connectingLineList):
        """Add the peak assignments to nmrResidues in the same chain.
        """
        for specAssignments in assignments.values():
            for nmrAtomPair in specAssignments:

                guiNmrAtomPair = (self.guiNmrAtoms.get(nmrAtomPair[0]),
                                  self.guiNmrAtoms.get(nmrAtomPair[1]),
                                  nmrAtomPair[2]
                                  )

                # get the peak and the spectrum
                peak = guiNmrAtomPair[2]
                spectrum = peak.peakList.spectrum

                if guiNmrAtomPair[0] is None and guiNmrAtomPair[1] is None:
                    continue

                if guiNmrAtomPair[0] is None:
                    if nmrAtomPair[1].nmrResidue.nmrChain is not nmrResidue.nmrChain:
                        continue

                    newGhostResidue = self._addGhostResidue(nmrAtomPair[0].nmrResidue,
                                                            guiNmrAtomPair[1],
                                                            nmrAtomPair[1].nmrResidue,
                                                            nmrAtomPair[0].name,
                                                            nmrAtomPair[1].name,
                                                            True,
                                                            lineList=connectingLineList)
                    guiNmrAtomPair = (self.guiNmrAtoms.get(nmrAtomPair[0]),
                                      self.guiNmrAtoms.get(nmrAtomPair[1]),
                                      nmrAtomPair[2]
                                      )

                    group = self.guiNmrResidues[nmrAtomPair[1].nmrResidue]
                    displacement = guiNmrAtomPair[1].getConnectedList(guiNmrAtomPair[0])
                    self._addConnectingLineToGroup(group,
                                                   guiNmrAtomPair[1],
                                                   guiNmrAtomPair[0],
                                                   spectrum.positiveContourColour,
                                                   2.0, displacement=displacement,
                                                   peak=peak, lineList=peaklineList, lineId=nmrResidue)

                elif guiNmrAtomPair[1] is None:
                    if nmrAtomPair[0].nmrResidue.nmrChain is not nmrResidue.nmrChain:
                        continue

                    newGhostResidue = self._addGhostResidue(nmrAtomPair[1].nmrResidue,
                                                            guiNmrAtomPair[0],
                                                            nmrAtomPair[0].nmrResidue,
                                                            nmrAtomPair[1].name,
                                                            nmrAtomPair[0].name,
                                                            True,
                                                            lineList=connectingLineList)
                    guiNmrAtomPair = (self.guiNmrAtoms.get(nmrAtomPair[0]),
                                      self.guiNmrAtoms.get(nmrAtomPair[1]),
                                      nmrAtomPair[2]
                                      )

                    group = self.guiNmrResidues[nmrAtomPair[0].nmrResidue]
                    displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
                    self._addConnectingLineToGroup(group,
                                                   guiNmrAtomPair[0],
                                                   guiNmrAtomPair[1],
                                                   spectrum.positiveContourColour,
                                                   2.0, displacement=displacement,
                                                   peak=peak, lineList=peaklineList, lineId=nmrResidue)

                else:
                    if nmrAtomPair[0].nmrResidue.nmrChain is nmrResidue.nmrChain:
                        group = self.guiNmrResidues[nmrAtomPair[0].nmrResidue]
                        displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
                        self._addConnectingLineToGroup(group,
                                                       guiNmrAtomPair[0],
                                                       guiNmrAtomPair[1],
                                                       spectrum.positiveContourColour,
                                                       2.0, displacement=displacement,
                                                       peak=peak, lineList=peaklineList, lineId=nmrResidue)

                    elif nmrAtomPair[1].nmrResidue.nmrChain is nmrResidue.nmrChain:
                        group = self.guiNmrResidues[nmrAtomPair[1].nmrResidue]
                        displacement = guiNmrAtomPair[1].getConnectedList(guiNmrAtomPair[0])
                        self._addConnectingLineToGroup(group,
                                                       guiNmrAtomPair[1],
                                                       guiNmrAtomPair[0],
                                                       spectrum.positiveContourColour,
                                                       2.0, displacement=displacement,
                                                       peak=peak, lineList=peaklineList, lineId=nmrResidue)
                    else:
                        continue

                guiNmrAtomPair[0].addConnectedList(guiNmrAtomPair[1])
                guiNmrAtomPair[1].addConnectedList(guiNmrAtomPair[0])

    def _getPeakAssignmentsForResidue(self, nmrResidue, nmrAtomIncludeList=None):
        """Get the list of peak assignments from the nmrAtoms
        interResidueAtomPairing is the linking within the same nmrResidue
        interChainAtomPairing is the linking within the same chain but to different nmrResidues
        crossChainAtomPairing is the linking to different chains
        """

        # create a set of sets ordered by spectra
        interResidueAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
        interChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
        crossChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
        # emptyAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())

        nmrChain = nmrResidue.nmrChain

        for nmrAtom in nmrResidue.nmrAtoms:

            if nmrAtom._flaggedForDelete or nmrAtom.isDeleted:
                continue

            for peak in nmrAtom.assignedPeaks:

                # ignore peaks that are due for delete (can probably also use the notifier list)
                if peak._flaggedForDelete or peak.isDeleted:
                    continue

                spec = peak.peakList.spectrum
                for assignment in peak.assignments:

                    # find the mainNmrResidue for -1 and +1 connections
                    newCon = list(assignment)
                    for conNum in range(len(assignment)):

                        # assignments could be None
                        if assignment[conNum] and assignment[conNum].nmrResidue.relativeOffset == -1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:

                            # this is a minus residue so find connected, have to traverse to the previousNmrResidue
                            # will it always exist?
                            conName = assignment[conNum].name
                            preN = assignment[conNum].nmrResidue.mainNmrResidue.previousNmrResidue
                            if preN:
                                newConSwap = [nmrA for nmrA in preN.nmrAtoms if nmrA.name == conName]
                                if newConSwap:
                                    newCon[conNum] = newConSwap[0]
                            else:
                                newCon[conNum] = None  # not connected so skip

                        elif assignment[conNum] and assignment[conNum].nmrResidue.relativeOffset == +1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:

                            # this is a plus residue so find connected, have to traverse to the nextNmrResidue
                            # will it always exist?
                            conName = assignment[conNum].name
                            preN = assignment[conNum].nmrResidue.mainNmrResidue.nextNmrResidue
                            if preN:
                                newConSwap = [nmrA for nmrA in preN.nmrAtoms if nmrA.name == conName]
                                if newConSwap:
                                    newCon[conNum] = newConSwap[0]
                            else:
                                newCon[conNum] = None  # not connected so skip

                    assignment = newCon

                    # only get the assignments a-b if a and b are defined in the spectrum magnetisationTransfers list
                    for mag in self._module.magnetisationTransfers[spec]:
                        nmrAtom0 = assignment[mag[0] - 1]
                        nmrAtom1 = assignment[mag[1] - 1]
                        nmrAtom0 = nmrAtom0 if nmrAtom0 and not (nmrAtom0.isDeleted or nmrAtom0._flaggedForDelete) else None
                        nmrAtom1 = nmrAtom1 if nmrAtom1 and not (nmrAtom1.isDeleted or nmrAtom1._flaggedForDelete) else None

                        if not None in (nmrAtom0, nmrAtom1):

                            # ignore nmrAtoms that are not in the include list (if specified)
                            if nmrAtomIncludeList is not None and not (nmrAtom0 in nmrAtomIncludeList or nmrAtom1 in nmrAtomIncludeList):
                                continue

                            if (nmrAtom0.nmrResidue is nmrResidue) and (nmrAtom1.nmrResidue is nmrResidue):

                                # interResidueAtomPairing
                                if (nmrAtom1, nmrAtom0, peak) not in interResidueAtomPairing[spec]:
                                    interResidueAtomPairing[spec].add((nmrAtom0, nmrAtom1, peak))

                            elif (nmrAtom0.nmrResidue.nmrChain is nmrChain) and (nmrAtom1.nmrResidue.nmrChain is nmrChain):

                                # connections within the same chain
                                if (nmrAtom1, nmrAtom0, peak) not in interChainAtomPairing[spec]:
                                    interChainAtomPairing[spec].add((nmrAtom0, nmrAtom1, peak))

                            # elif (nmrAtom0.nmrResidue.nmrChain is nmrChain) and (nmrAtom1.nmrResidue.nmrChain is not nmrChain):
                            else:

                                # connections to a dif
                                if (nmrAtom1, nmrAtom0, peak) not in crossChainAtomPairing[spec]:
                                    crossChainAtomPairing[spec].add((nmrAtom0, nmrAtom1, peak))

        return interResidueAtomPairing, interChainAtomPairing, crossChainAtomPairing
        # return emptyAtomPairing, emptyAtomPairing, crossChainAtomPairing

    #==========================================================================================
    # spectrum update

    def _removeLinesFromScene(self, lineDist):
        """Remove all the lines in the group from the scene.
        """
        for lineList in lineDist.values():
            for line in lineList:
                if line in self._scene.items():
                    self._scene.removeItem(line)
        lineDist.clear()

    def removeAssignmentLinesFromScene(self):
        """Remove all the peakLines from the scene.
        """
        self._removeLinesFromScene(self.assignmentLines)

    def removeConnectionLinesFromScene(self):
        """Remove all the connection from the scene.
        """
        self._removeLinesFromScene(self.connectingLines)

    def clearAllGuiNmrAtoms(self):
        """clear the displacements in the nmrAtoms.
        """
        # clear the displacement values but keep the dict connections
        for guiAtom in self.guiNmrAtoms.values():
            guiAtom.clearConnectedList()

    def rebuildPeakAssignments(self):
        """Rebuild all the peak assignments in the display after changing the number of spectra.
        """

        # remove all previous assignment lines and reset the dict
        self.removeAssignmentLinesFromScene()
        self.clearAllGuiNmrAtoms()
        for nmrChainId in self.nmrChains.keys():
            self._addAllPeakAssignments(nmrChainId)

        # update the endpoints
        self.updateEndPoints(self.assignmentLines)

    def updateEndPoints(self, lineDict):
        """Update the end points from the dict.
        """
        for lineList in lineDict.values():
            for line in lineList:
                try:
                    line.updateEndPoints()
                except Exception as es:
                    pass

    def updateAssignmentLines(self):
        """Update the endpoints of the assignment lines.
        """
        self.updateEndPoints(self.assignmentLines)

    def updateConnectionLines(self):
        """Update the endpoints of the connection lines.
        """
        self.updateEndPoints(self.connectingLines)

    def updateGuiResiduePositions(self, nmrChainId, updateMainChain=True, updateConnectedChains=True):
        """Update the positions of the residues and connected residues in other chains if required.
        """
        if updateMainChain:
            # update the group positions
            self.updateMainChainPositions(nmrChainId)
        if updateConnectedChains:
            # update crossChainResidue positions
            self.updateConnectedChainPositions(nmrChainId)

        # update the endpoints
        self.updateConnectionLines()
        self.updateAssignmentLines()

    #==========================================================================================

    def getAssignmentLinesFromPeaks(self, peaks):
        """Get the list of assignment lines attached o the given peaks.
        """

    def _addAdjacentResiduesToSet(self, nmrResidue, residueSet):
        """Add the adjacent nmrResidues into the set.
        """
        residueSet.add(nmrResidue)
        nmr = nmrResidue.previousNmrResidue
        if nmr and nmr.mainNmrResidue:
            residueSet.add(nmr.mainNmrResidue)
        nmr = nmrResidue.nextNmrResidue
        if nmr and nmr.mainNmrResidue:
            residueSet.add(nmr.mainNmrResidue)

    def rebuildPeakLines(self, peaks, rebuildPeakLines=False, makeListFromPeak=False):
        """Clear all lines on the display associated with peak.
        """
        # find the current lines associated with the notified peak (or list of peaks)
        peaks = makeIterableList(peaks)

        # make list of peakLines attached to these peaks
        peakLines = [peakLine for peakLineList in self.assignmentLines.values()
                     for peakLine in peakLineList
                     if peakLine._peak in peaks]
        guiNmrAtomSet = set()
        nmrResidueSet = set()

        if not makeListFromPeak:
            # make a list of all the guiNmrAtoms/nmrResidues that are associated with the peak
            for peakLine in peakLines:
                # current associated guiNmrAtoms
                guiNmrAtomSet.add(peakLine.guiAtom1)
                guiNmrAtomSet.add(peakLine.guiAtom2)

                # current associated nmrResidues and their previous/nextNmrResidues
                self._addAdjacentResiduesToSet(peakLine.guiAtom1.nmrAtom.nmrResidue, nmrResidueSet)
                self._addAdjacentResiduesToSet(peakLine.guiAtom2.nmrAtom.nmrResidue, nmrResidueSet)

        else:
            # make a new list for creating a peak; necessary for undo of delete peak as the assignedNmrAtom list exists
            assignmentAtoms = set([nmrAtom for peak in peaks
                                   for assignment in peak.assignments
                                   for nmrAtom in assignment
                                   if nmrAtom in self.guiNmrAtoms])
            for nmrAtom in assignmentAtoms:
                guiNmrAtomSet.add(self.guiNmrAtoms[nmrAtom])
                self._addAdjacentResiduesToSet(nmrAtom.nmrResidue, nmrResidueSet)

        for guiAtom in guiNmrAtomSet:
            for peakLineList in self.assignmentLines.values():
                peakLines = [peakLine for peakLine in peakLineList
                             if peakLine.guiAtom1 is guiAtom or peakLine.guiAtom2 is guiAtom]

                # remove all graphic lines
                for peakLine in peakLines:
                    peakLineList.remove(peakLine)
                    if peakLine in self._scene.items():
                        self._scene.removeItem(peakLine)

            # clear connectivity list of guiNmrAtoms, but don't delete
            guiAtom.clearConnectedList()

        if rebuildPeakLines:
            # now rebuild for the new peak values
            # assumes that the peakAssignments have changed - possibly use different notifier
            if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():

                nmrAtomIncludeList = tuple(guiAtom.nmrAtom for
                                           guiAtom in
                                           guiNmrAtomSet)

                # # create a set of sets ordered by spectra for active lines
                # self.LOCALinterResidueAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
                # self.LOCALinterChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
                # self.LOCALcrossChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())

                for nmrResidue in nmrResidueSet:

                    # only process residues in the current visible chain
                    if nmrResidue.nmrChain is self._module.nmrChain:
                        # add the internally connected Lines
                        internalAssignments, interChainAssignments, crossChainAssignments = \
                            self._getPeakAssignmentsForResidue(nmrResidue,
                                                               nmrAtomIncludeList=nmrAtomIncludeList)

                        self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                        self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                        self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                    self.assignmentLines, self.connectingLines)

                # self._addPeakAssignmentLinesToGroup(self.LOCALinterResidueAtomPairing, self.assignmentLines)
                # self._addPeakAssignmentLinesToGroup(self.LOCALinterChainAtomPairing, self.assignmentLines)
                # self._addPeakAssignmentLinesToAdjacentGroup(self._module.nmrChain, self.LOCALcrossChainAtomPairing,
                #                                             self.assignmentLines, self.connectingLines)

            first = next(iter(nmrResidueSet or []), None)
            thisId = first.nmrChain.pid if first else 'noChainId'
            self.updateGuiResiduePositions(thisId, updateMainChain=True, updateConnectedChains=True)

    #==========================================================================================

    def _addNmrAtomToGuiResidues(self, nmrAtom, atomSpacing=None):
        """Update the guiNmrAtoms for a newly created/undeleted nmrAtom.
        Assumes that the nmrAtom/guiNmrAtom does not exist
        """
        nmrResidue = nmrAtom.nmrResidue

        guiAtoms = {}
        if atomSpacing:
            self.atomSpacing = atomSpacing
        atomNames = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
        residueAtoms = DEFAULT_RESIDUE_ATOMS.copy()

        for k, v in residueAtoms.items():
            if k in atomNames:
                fetchedNmrAtom = nmrResidue.fetchNmrAtom(name=k)
            else:
                fetchedNmrAtom = None
            if fetchedNmrAtom is nmrAtom:
                guiAtoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)

        ii = self.getIndexNmrResidue(nmrResidue)
        if ii is not None:
            res, oldGuiAtoms = self._getNmrResiduePair(ii)
            self._setNmrResiduePair(ii, res, oldGuiAtoms.update(guiAtoms))

        if nmrResidue in self.guiNmrResidues:
            guiResidueGroup = self.guiNmrResidues[nmrResidue]

            # add the guiAtoms to the group and set the reverse link
            for item in guiAtoms.values():
                guiResidueGroup.addToGroup(item)
                item.guiNmrResidueGroup = guiResidueGroup

    def createNmrAtoms(self, nmrAtoms):
        """Create new peak lines associated with the created/undeleted nmrAtoms.
        """
        nmrAtoms = makeIterableList(nmrAtoms)
        peaks = [peak for nmrAtom in nmrAtoms for peak in nmrAtom.assignedPeaks]

        for nmrAtom in nmrAtoms:
            if nmrAtom not in self.guiNmrAtoms:
                # add the new nmrAtoms
                self._addNmrAtomToGuiResidues(nmrAtom)
        self.rebuildPeakLines(peaks, rebuildPeakLines=True, makeListFromPeak=True)

    def _searchPeakLines(self, nmrAtoms, includeDeleted=False):
        """Return a list of the peakLines containing one of the nmrAtoms in the list.
        """
        peakLines = []
        for lineList in self.assignmentLines.values():
            for line in lineList:
                nmrAtom = line.guiAtom1.nmrAtom if line.guiAtom1 else None
                if nmrAtom in nmrAtoms and (includeDeleted or not (nmrAtom.isDeleted or nmrAtom._flaggedForDelete)):
                    peakLines.append(line)
                nmrAtom = line.guiAtom2.nmrAtom if line.guiAtom2 else None
                if nmrAtom in nmrAtoms and (includeDeleted or not (nmrAtom.isDeleted or nmrAtom._flaggedForDelete)):
                    peakLines.append(line)

        return peakLines

    def deleteNmrAtoms(self, nmrAtoms):
        """Delete peak lines associated with the deleted nmrAtoms.
        """
        nmrAtoms = makeIterableList(nmrAtoms)
        peakLines = self._searchPeakLines(nmrAtoms, includeDeleted=True)
        peaks = [peakLine._peak for peakLine in peakLines]

        self.rebuildPeakLines(peaks, rebuildPeakLines=True)

    #==========================================================================================

    def rebuildNmrResidues(self, nmrResidues):
        """Rebuild the peaks of a specified nmrResidue.
        """
        # now rebuild for the new peak values
        # assumes that the peakAssignments have changed - possibly use different notifier
        nmrResidues = makeIterableList(nmrResidues)

        nmrAtomIncludeList = tuple(nmrAtom for nmrResidue in nmrResidues for nmrAtom in nmrResidue.nmrAtoms)
        guiNmrAtomSet = set([self.guiNmrAtoms[nmrAtom] for nmrAtom in nmrAtomIncludeList
                             if nmrAtom in self.guiNmrAtoms])

        for guiAtom in guiNmrAtomSet:
            for peakLineList in self.assignmentLines.values():
                peakLines = [peakLine for peakLine in peakLineList
                             if peakLine.guiAtom1 is guiAtom or peakLine.guiAtom2 is guiAtom]

                # remove all graphic lines
                for peakLine in peakLines:
                    peakLineList.remove(peakLine)
                    if peakLine in self._scene.items():
                        self._scene.removeItem(peakLine)

            # clear connectivity list of guiNmrAtoms
            guiAtom.clearConnectedList()

        if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():

            # # create a set of sets ordered by spectra for active lines
            # self.LOCALinterResidueAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
            # self.LOCALinterChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
            # self.LOCALcrossChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())

            for nmrResidue in nmrResidues:

                # only process residues in the current visible chain
                if nmrResidue is nmrResidue.mainNmrResidue and nmrResidue.nmrChain is self._module.nmrChain:
                    # add the internally connected Lines
                    internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue,
                                                                                                                           nmrAtomIncludeList=nmrAtomIncludeList)

                    self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                self.assignmentLines, self.connectingLines)

                if nmrResidue.nextNmrResidue and nmrResidue.nextNmrResidue.mainNmrResidue:
                    nextNmrResidue = nmrResidue.nextNmrResidue.mainNmrResidue

                    # only process residues in the current visible chain
                    # add the internally connected Lines
                    internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nextNmrResidue,
                                                                                                                           nmrAtomIncludeList=nmrAtomIncludeList)

                    self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToAdjacentGroup(nextNmrResidue, crossChainAssignments,
                                                                self.assignmentLines, self.connectingLines)

            # self._addPeakAssignmentLinesToGroup(self.LOCALinterResidueAtomPairing, self.assignmentLines)
            # self._addPeakAssignmentLinesToGroup(self.LOCALinterChainAtomPairing, self.assignmentLines)
            # self._addPeakAssignmentLinesToAdjacentGroup(self._module.nmrChain, self.LOCALcrossChainAtomPairing,
            #                                             self.assignmentLines, self.connectingLines)

        self.updateGuiResiduePositions(nmrResidues[0].nmrChain.pid, updateMainChain=True, updateConnectedChains=True)

    #==========================================================================================
    #==========================================================================================
    #==========================================================================================


#==========================================================================================
# Sequence Graph Main Module
#==========================================================================================


LINKTOPULLDOWNCLASS = 'linkToPulldownClass'


class SequenceGraphModule(CcpnModule):
    """
    A module for the display of stretches of sequentially linked and assigned stretches of
    NmrResidues.
    """
    className = 'SequenceGraph'

    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'left'

    # consistent with nmrResidueTable - move to generic class later
    activePulldownClass = NmrChain

    def __init__(self, mainWindow=None, name='Sequence Graph', nmrChain=None):

        CcpnModule.__init__(self, mainWindow=mainWindow, name=name)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        if self.mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current
            ###self.setMode('fragment')  # cannot be moved up!
        else:
            self.application = None
            self.project = None
            self.current = None

        self.splitter = Splitter(self.mainWidget, horizontal=False)
        self._sequenceModuleFrame = Frame(None, setLayout=True)
        # self._SequenceGraphFrame = Frame(self.splitter, setLayout=True)
        self.mainWidget.getLayout().addWidget(self.splitter, 1, 0)

        self.thisSequenceModule = SequenceModule(moduleParent=self,
                                                 parent=self._sequenceModuleFrame,
                                                 mainWindow=mainWindow)

        self.colours = getColours()
        self._lineColour = self.colours[SEQUENCEGRAPHMODULE_LINE]
        self._textColour = self.colours[SEQUENCEGRAPHMODULE_TEXT]

        ###frame = Frame(parent=self.mainWidget)
        self._sequenceGraphScrollArea = QtWidgets.QScrollArea()
        # self._sequenceGraphScrollArea.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        self._sequenceGraphScrollArea.setWidgetResizable(True)
        self._sequenceGraphScrollArea.setMinimumHeight(80)

        self.splitter.addWidget(self._sequenceGraphScrollArea)
        self.splitter.addWidget(self._sequenceModuleFrame)
        self.splitter.setStretchFactor(0, 5)
        self.splitter.setChildrenCollapsible(False)

        # add the settings widgets defined from the following orderedDict - test for refactored
        settingsDict = OrderedDict((('peakAssignments', {'label'   : 'Show peak assignments:',
                                                         'tipText' : 'Show peak assignments on display coloured by positiveContourColour.',
                                                         'callBack': self.showNmrChainFromPulldown,
                                                         'enabled' : True,
                                                         'checked' : True,
                                                         '_init'   : None,
                                                         }),
                                    ('treeView', {'label'   : 'Tree view:',
                                                  'tipText' : 'Show peak assignments as a tree below the main backbone.',
                                                  'callBack': None,  #self._updateShowTreeAssignments,
                                                  'enabled' : False,
                                                  'checked' : False,
                                                  '_init'   : None,  #self._updateShowTreeAssignments,
                                                  }),
                                    ('showSideChain', {'label'   : 'Show side chain:',
                                                       'tipText' : 'Show side chain atoms and connections above the main chain.',
                                                       'callBack': None,  #self.showNmrChainFromPulldown,
                                                       'enabled' : False,
                                                       'checked' : False,
                                                       '_init'   : None,
                                                       }),
                                    ('sequentialStrips', {'label'   : 'Show sequential strips:',
                                                          'tipText' : 'Show nmrResidue in all strips.',
                                                          'callBack': None,  #self.showNmrChainFromPulldown,
                                                          'enabled' : True,
                                                          'checked' : False,
                                                          '_init'   : None,
                                                          }),
                                    ('markPositions', {'label'   : 'Mark positions:',
                                                       'tipText' : 'Mark positions in all strips.',
                                                       'callBack': None,  #self.showNmrChainFromPulldown,
                                                       'enabled' : True,
                                                       'checked' : True,
                                                       '_init'   : None,
                                                       }),
                                    ('autoClearMarks', {'label'   : 'Auto clear marks:',
                                                        'tipText' : 'Auto clear all previous marks',
                                                        'callBack': None,
                                                        'enabled' : True,
                                                        'checked' : True,
                                                        '_init'   : None,
                                                        }),
                                    ))
        if self.activePulldownClass:
            settingsDict.update(OrderedDict(((LINKTOPULLDOWNCLASS, {'label'   : 'Link to current %s:' % self.activePulldownClass.className,
                                                       'tipText' : 'Set/update current %s when selecting from pulldown' % self.activePulldownClass.className,
                                                       'callBack': None,
                                                       'enabled' : True,
                                                       'checked' : True,
                                                       '_init'   : None,
                                                       }),
                                )))

        self._SGwidget = SequenceGraphSettings(parent=self.settingsWidget, mainWindow=self.mainWindow,
                                               settingsDict=settingsDict,
                                               grid=(0, 0))

        self.initialiseScene()
        self.residueCount = 0
        self.atomSpacing = 66
        self.nmrResidueList = NmrResidueList(self.mainWindow, self._SGwidget, self._lineColour, self._textColour, self.atomSpacing,
                                             self.scene, self)
        self._deleteStore = {}

        colwidth = 180
        self._MWwidget = Widget(self.mainWidget, setLayout=True,
                                grid=(0, 0), vAlign='top', hAlign='left')

        self.nmrChainPulldown = NmrChainPulldown(self._MWwidget, self.mainWindow, grid=(0, 0), gridSpan=(1, 1),
                                                 showSelectName=True,
                                                 fixedWidths=(colwidth, colwidth, colwidth),
                                                 callback=self.showNmrChainFromPulldown)

        # self.refreshCheckBox = CheckBoxCompoundWidget(self._MWwidget,
        #                                               labelText='Auto refresh NmrChain:',
        #                                               checked=True,
        #                                               fixedWidths=(colwidth, 15),
        #                                               orientation='right', hAlign='left',
        #                                               tipText='Update display when current.nmrChain changes',
        #                                               grid=(0, 1), gridSpan=(1, 1))

        self.sequenceCheckBox = CheckBoxCompoundWidget(self._MWwidget,
                                                       labelText='Show Sequence:',
                                                       checked=True,
                                                       fixedWidths=(colwidth, 15),
                                                       orientation='right', hAlign='left',
                                                       tipText='Show chain sequences',
                                                       callback=self._toggleSequence,
                                                       grid=(0, 2), gridSpan=(1, 1))

        self.nmrResiduesCheckBox = CheckBoxCompoundWidget(self._MWwidget,
                                                          labelText='Show all NmrResidues:',
                                                          checked=True,
                                                          fixedWidths=(colwidth, 15),
                                                          orientation='right', hAlign='left',
                                                          tipText='Show all the NmrResidues in the NmrChain',
                                                          callback=self.showNmrChainFromPulldown,
                                                          grid=(0, 3), gridSpan=(1, 1))

        self._MWwidget.setMinimumWidth(self._MWwidget.sizeHint().width())
        self._MWwidget.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.settingsWidget.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Minimum)

        self.editingToolbar = ToolBar(self._MWwidget, grid=(0, 6), gridSpan=(1, 1), hAlign='right', iconSizes=(24, 24))
        # self.editingToolbar = ToolBar(self._SequenceModuleFrame, grid=(0, 6), gridSpan=(1, 1), hAlign='right', iconSizes=(24,24))

        self.disconnectPreviousAction = self.editingToolbar.addAction("disconnectPrevious", self.disconnectPreviousNmrResidue)
        self.disconnectPreviousIcon = Icon('icons/disconnectPrevious')
        self.disconnectPreviousAction.setIcon(self.disconnectPreviousIcon)
        self.disconnectAction = self.editingToolbar.addAction("disconnect", self.disconnectNmrResidue)
        self.disconnectIcon = Icon('icons/disconnect')
        self.disconnectAction.setIcon(self.disconnectIcon)
        self.disconnectNextAction = self.editingToolbar.addAction("disconnectNext", self.disconnectNextNmrResidue)
        self.disconnectNextIcon = Icon('icons/disconnectNext')
        self.disconnectNextAction.setIcon(self.disconnectNextIcon)

        # add mouse handler for the QGraphicsLineItems
        self._preMouserelease = self.scene.mouseReleaseEvent
        self.scene.mouseReleaseEvent = self._sceneMouseRelease

        # calulate the connections between axes based on experiment types
        self._updateMagnetisationTransfers()

        # stop the mainWidget from squishing during a resize
        self.mainWidget.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)

        # install the event filter to handle maximising from floated dock
        # self.installMaximiseEventHandler(self._maximise, self._closeModule)

        # initialise notifiers
        self._registerNotifiers()

        self.selectSequence(nmrChain)

    def _sceneMouseRelease(self, event):
        """Add a mouse handler to popupa menu from the contained scene
        """
        if event.button() == QtCore.Qt.RightButton:
            object = self.scene.mouseGrabberItem()
            # print('>>>grab', object)
            if object:
                self._raiseContextMenu(object, event)
        self._preMouserelease(event)

    # def _checkLayoutInit(self):
    #     """This is a hack so that the state changes when the layout loads
    #     After the layout initialise, this function is removed
    #     """
    #     self._updateShowTreeAssignments()
    #     self.assignmentsTreeCheckBox.checkBox.stateChanged.disconnect(self._checkLayoutInit)

    def _maximise(self):
        """Maximise the attached table
        """
        pass

    def _updateMagnetisationTransfers(self):
        """Generate the list that defines which couplings there are between the nmrAtoms attached to each peak.
        """
        self.magnetisationTransfers = OrderedDict()
        for spec in self.project.spectra:
            if not spec._flaggedForDelete:
                self.magnetisationTransfers[spec] = {}
                for mt in spec.magnetisationTransfers:
                    self.magnetisationTransfers[spec][mt] = set()

    def _blockEvents(self):
        """Block all updates/signals/notifiers in the scene.
        """
        self.setUpdatesEnabled(False)
        self.scene.blockSignals(True)
        # self.setBlankingAllNotifiers(True)

    def _unblockEvents(self):
        """Unblock all updates/signals/notifiers in the scene.
        """
        # self.setBlankingAllNotifiers(False)
        self.scene.blockSignals(False)
        self.setUpdatesEnabled(True)

    @contextmanager
    def sceneBlocking(self):
        """Context manager to handle blocking, unblocking, resizing of the scene.
        """
        self._blockEvents()
        try:
            # pass control to the calling function
            yield

        finally:
            self._unblockEvents()

            # resize to the new items and spawns a repaint
            self.scene.setSceneRect(self.scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))

    def _updateSpectra(self, data=None):
        """Update list of current spectra and generate new magnetisationTransfer list
        """
        if data:
            trigger = data[Notifier.TRIGGER]

            self._updateMagnetisationTransfers()

            if trigger in [Notifier.CREATE, Notifier.DELETE]:
                nmrChainPid = self.nmrChainPulldown.getText()
                if nmrChainPid:
                    with self.sceneBlocking():
                        self.nmrResidueList.rebuildPeakAssignments()

    def selectSequence(self, nmrChain=None):
        """Manually select a Sequence from the pullDown
        """
        if nmrChain is None:
            # logger.warning('select: No Sequence selected')
            # raise ValueError('select: No Sequence selected')
            self.nmrChainPulldown.selectFirstItem()
        else:
            if not isinstance(nmrChain, NmrChain):
                logger.warning('select: Object is not of type Sequence')
                raise TypeError('select: Object is not of type Sequence')
            else:
                for widgetObj in self.nmrChainPulldown.textList:
                    if nmrChain.pid == widgetObj:
                        self.nmrChainPulldown.select(nmrChain.pid)

    def _registerNotifiers(self):
        """Register the required notifiers
        """
        self._peakNotifier = self.setNotifier(self.project,
                                              [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
                                              Peak.className,
                                              self._updatePeaks,
                                              onceOnly=True)

        # not required
        # self._nmrChainNotifier = self.setNotifier(self.project,
        #                                           [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
        #                                           NmrChain.className,
        #                                           self._updateNmrChains,
        #                                           onceOnly=True)

        self._nmrResidueNotifier = self.setNotifier(self.project,
                                                    [Notifier.CREATE, Notifier.DELETE, Notifier.RENAME],
                                                    NmrResidue.className,
                                                    self._updateNmrResidues,
                                                    onceOnly=True)

        self._nmrResidueChangeNotifier = self.setNotifier(self.project,
                                                          [Notifier.CHANGE],
                                                          NmrResidue.className,
                                                          self._changeNmrResidues,
                                                          onceOnly=True)

        self._nmrAtomNotifier = self.setNotifier(self.project,
                                                 [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
                                                 NmrAtom.className,
                                                 self._updateNmrAtoms,
                                                 onceOnly=True)

        # notifier to change the magnetisationTransfer list when new spectrum added
        self._spectrumListNotifier = self.setNotifier(self.project,
                                                      [Notifier.CREATE, Notifier.DELETE],
                                                      Spectrum.className,
                                                      self._updateSpectra)

        self._currentNmrResidueNotifier = self.setNotifier(self.current,
                                                           [Notifier.CURRENT],
                                                           targetName=NmrResidue._pluralLinkName,
                                                           callback=self._selectCurrentNmrResidues)

        # new notifier to respond to changing current nmrChain
        if self.activePulldownClass:
            self._setCurrentPulldown = Notifier(self.current,
                                                [Notifier.CURRENT],
                                                targetName=self.activePulldownClass._pluralLinkName,
                                                callback=self._selectCurrentPulldownClass)

    def _selectCurrentPulldownClass(self, data):
        """Respond to change in current activePulldownClass
        """
        checkBox = self._SGwidget._getCheckBox(LINKTOPULLDOWNCLASS)
        if self.activePulldownClass and checkBox and checkBox.isChecked() and \
                self.current.nmrChain and self.current.nmrChain != self.nmrChain:
            self.nmrChainPulldown.select(self.current.nmrChain.pid)

    def _repopulateModule(self):
        """CCPN Internal: Repopulate the required widgets in the module
        This is will be attached to GuiNotifiers
        """
        self.showNmrChainFromPulldown()

    # def _updateModule(self, nmrChains=None):
    #     """Update in response to change of current.nmrChains
    #     """
    #     #if nmrChains is None or len(nmrChains)==0: return
    #     nmrChain = self.current.nmrChain
    #     if not nmrChain:
    #         return
    #
    #     if not self.refreshCheckBox.isChecked():
    #         return
    #
    #     # select the chain from the pullDown - should automatically change the display
    #     self.nmrChainPulldown.select(nmrChain.pid)

    # def setMode(self, mode):
    #   if self.project.nmrChains:
    #     self.editingToolbar.hide()
    #     if mode == 'fragment':
    #       self.editingToolbar.show()
    #       #self.nmrChainPulldown.setData([c.pid for c in self.project.nmrChains])
    #       #self.nmrChainLabel.setText('NmrChain: ')
    #     elif mode == 'Assigned - backbone':
    #       pass
    #       #self.nmrChainLabel.setText('Chain: ')
    #       #self.nmrChainPulldown.setData([self.project.getByPid('NC:%s' % chain.shortName).pid for chain in self.project.chains if self.project.getByPid('NC:%s' % chain.shortName)])
    #     self.modePulldown.select(mode)
    #     self.setNmrChainDisplay(self.nmrChainPulldown.getText())
    #   else:
    #     logger.warning('No valid NmrChain is selected.')

    # def _addAdjacentResiduesToSet(self, nmrResidue, residueSet):
    #     residueSet.add(nmrResidue)
    #     nmr = nmrResidue.previousNmrResidue
    #     if nmr:
    #         nmr = nmr.mainNmrResidue
    #         residueSet.add(nmr)
    #     nmr = nmrResidue.nextNmrResidue
    #     if nmr:
    #         nmr = nmr.mainNmrResidue
    #         residueSet.add(nmr)

    def _rebuildNmrResidues(self, nmrChainId, nmrResidues):
        """Rebuild the peaks of a specified nmrResidue.
        """
        # now rebuild for the new peak values
        # assumes that the peakAssignments have changed - possibly use different notifier
        nmrResidues = makeIterableList(nmrResidues)

        nmrAtomIncludeList = tuple(nmrAtom for nmrResidue in nmrResidues for nmrAtom in nmrResidue.nmrAtoms)
        guiNmrAtomSet = set([self.nmrResidueList.guiNmrAtoms[nmrAtom] for nmrAtom in nmrAtomIncludeList])

        for guiAtom in guiNmrAtomSet:
            for peakLineList in self.nmrResidueList.assignmentLines.values():
                peakLines = [peakLine for peakLine in peakLineList
                             if peakLine.guiAtom1 is guiAtom or peakLine.guiAtom2 is guiAtom]

                # remove all graphic lines
                for peakLine in peakLines:
                    peakLineList.remove(peakLine)
                    if peakLine in self.scene.items():
                        self.scene.removeItem(peakLine)

            # clear connectivity list of guiNmrAtoms
            guiAtom.clearConnectedList()

        if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():

            # # create a set of sets ordered by spectra for active lines
            # self.LOCALinterResidueAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
            # self.LOCALinterChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())
            # self.LOCALcrossChainAtomPairing = OrderedDict((spec, set()) for spec in self._module.magnetisationTransfers.keys())

            for nmrResidue in nmrResidues:

                # only process residues in the current visible chain
                if nmrResidue is nmrResidue.mainNmrResidue and nmrResidue.nmrChain is self.nmrChain:
                    # add the internally connected Lines
                    internalAssignments, interChainAssignments, crossChainAssignments = self.nmrResidueList._getPeakAssignmentsForResidue(nmrResidue,
                                                                                                                                          nmrAtomIncludeList=nmrAtomIncludeList)

                    self.nmrResidueList._addPeakAssignmentLinesToGroup(internalAssignments, self.nmrResidueList.assignmentLines)
                    self.nmrResidueList._addPeakAssignmentLinesToGroup(interChainAssignments, self.nmrResidueList.assignmentLines)
                    self.nmrResidueList._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                               self.nmrResidueList.assignmentLines, self.nmrResidueList.connectingLines)

                if nmrResidue.nextNmrResidue and nmrResidue.nextNmrResidue.mainNmrResidue:
                    nextNmrResidue = nmrResidue.nextNmrResidue.mainNmrResidue

                    # only process residues in the current visible chain
                    # add the internally connected Lines
                    internalAssignments, interChainAssignments, crossChainAssignments = self.nmrResidueList._getPeakAssignmentsForResidue(nextNmrResidue,
                                                                                                                                          nmrAtomIncludeList=nmrAtomIncludeList)

                    self.nmrResidueList._addPeakAssignmentLinesToGroup(internalAssignments, self.nmrResidueList.assignmentLines)
                    self.nmrResidueList._addPeakAssignmentLinesToGroup(interChainAssignments, self.nmrResidueList.assignmentLines)
                    self.nmrResidueList._addPeakAssignmentLinesToAdjacentGroup(nextNmrResidue, crossChainAssignments,
                                                                               self.nmrResidueList.assignmentLines, self.nmrResidueList.connectingLines)

            # self._addPeakAssignmentLinesToGroup(self.LOCALinterResidueAtomPairing, self.assignmentLines)
            # self._addPeakAssignmentLinesToGroup(self.LOCALinterChainAtomPairing, self.assignmentLines)
            # self._addPeakAssignmentLinesToAdjacentGroup(self._module.nmrChain, self.LOCALcrossChainAtomPairing,
            #                                             self.assignmentLines, self.connectingLines)

        self.nmrResidueList.updateGuiResiduePositions(nmrChainId, updateMainChain=True, updateConnectedChains=True)

    # def _rebuildPeakLines(self, peaks, rebuildPeakLines=False, makeListFromPeak=False):
    #     """Clear all lines on the display associated with peak.
    #     """
    #     # find the current lines associated with the notified peak (or list of peaks)
    #     peaks = makeIterableList(peaks)
    #
    #     # make list of peakLines attached to these peaks
    #     peakLines = [peakLine for peakLineList in self.assignmentLines.values()
    #                  for peakLine in peakLineList
    #                  if peakLine._peak in peaks]
    #
    #     guiNmrAtomSet = set()
    #     nmrResidueSet = set()
    #
    #     if not makeListFromPeak:
    #         # make a list of all the guiNmrAtoms/nmrResidues that are associated with the peak
    #         for peakLine in peakLines:
    #             # current associated guiNmrAtoms
    #             guiNmrAtomSet.add(peakLine.guiAtom1)
    #             guiNmrAtomSet.add(peakLine.guiAtom2)
    #
    #             # current associated nmrResidues and their previous/nextNmrResidues
    #             self._addAdjacentResiduesToSet(peakLine.guiAtom1.nmrAtom.nmrResidue, nmrResidueSet)
    #             self._addAdjacentResiduesToSet(peakLine.guiAtom2.nmrAtom.nmrResidue, nmrResidueSet)
    #
    #     else:
    #         # make a new list for creating a peak; necessary for undo of delete peak as the assignedNmrAtom list exists
    #         assignmentAtoms = set([nmrAtom for peak in peaks
    #                                for assignment in peak.assignments
    #                                for nmrAtom in assignment
    #                                if nmrAtom in self.guiNmrAtomDict])
    #         for nmrAtom in assignmentAtoms:
    #             # for peak in peaks:
    #             #     for assignment in peak.assignments:
    #             #         for nmrAtom in assignment:
    #             #
    #             #             if nmrAtom in self.guiNmrAtomDict:
    #             guiNmrAtomSet.add(self.guiNmrAtomDict[nmrAtom])
    #             self._addAdjacentResiduesToSet(nmrAtom.nmrResidue, nmrResidueSet)
    #
    #     for guiAtom in guiNmrAtomSet:
    #         for peakLineList in self.assignmentLines.values():
    #             peakLines = [peakLine for peakLine in peakLineList
    #                          if peakLine.guiAtom1 is guiAtom or peakLine.guiAtom2 is guiAtom]
    #
    #             # remove all graphic lines
    #             for peakLine in peakLines:
    #                 peakLineList.remove(peakLine)
    #                 self.scene.removeItem(peakLine)
    #
    #         # clear connectivity list of guiNmrAtoms
    #         guiAtom.clearConnectedList()
    #
    #     if rebuildPeakLines:
    #         # now rebuild for the new peak values
    #         # assumes that the peakAssignments have changed - possibly use different notifier
    #         if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
    #             for nmrResidue in nmrResidueSet:
    #
    #                 # only process residues in the current visible chain
    #                 if nmrResidue is nmrResidue.mainNmrResidue and nmrResidue.nmrChain is self.nmrChain:
    #                     # add the internally connected Lines
    #                     internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue,
    #                                                                                                                            nmrAtomIncludeList=tuple(
    #                                                                                                                                    guiAtom.nmrAtom for
    #                                                                                                                                    guiAtom in
    #                                                                                                                                    guiNmrAtomSet))
    #                     self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
    #                     self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
    #                     self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
    #                                                                 self.assignmentLines, self.connectingLines)
    #
    #         self._updateGuiResiduePositions(updateMainChain=True, updateConnectedChains=True)

    def _updatePeaks(self, data):
        """Update the peaks in the display.
        """
        peak = data[Notifier.OBJECT]
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if trigger == Notifier.DELETE:
                # print('>>>_updatePeaks delete', peak)
                self.nmrResidueList.rebuildPeakLines(peak, rebuildPeakLines=True)

            elif trigger == Notifier.CREATE:
                # print('>>>_updatePeaks create', peak)
                self.nmrResidueList.rebuildPeakLines(peak, rebuildPeakLines=True, makeListFromPeak=True)

            elif trigger == Notifier.CHANGE:
                # print('>>>_updatePeaks change', peak)
                self.nmrResidueList.rebuildPeakLines(peak, rebuildPeakLines=True, makeListFromPeak=True)

    # def _updateNmrChains(self, data):
    #     """Update the nmrChains in the display.
    #     """
    #     nmrChain = data[Notifier.OBJECT]
    #
    #     # print('>>>_updateNmrChains', nmrChain)
    #     trigger = data[Notifier.TRIGGER]
    #
    #     # with self.sceneBlocking():
    #     #     if trigger == Notifier.DELETE:
    #     #         print('>>>delete nmrChain - no action', nmrChain)
    #     #
    #     #     elif trigger == Notifier.CREATE:
    #     #         print('>>>create nmrChain - no action', nmrChain)
    #     #
    #     #     elif trigger == Notifier.CHANGE:
    #     #         print('>>>change nmrChain - no action', nmrChain)

    def _selectCurrentNmrResidues(self, data):
        """
        Notifier Callback for selecting current nmrResidue
        """
        objList = data[CallBack.OBJECT]

        if not self.nmrResiduesCheckBox.isChecked():
            nmrResidue = objList.nmrResidue

            # redraw the nmrResidues if current is in the displayed chain and not already visible
            if nmrResidue.nmrChain == self.nmrChain and nmrResidue not in self.nmrResidueList.guiNmrResidues:
                self.showNmrChainFromPulldown(nmrResidue)

    def _updateNmrResidues(self, data):
        """Update the nmrResidues in the display.
        """
        nmrResidue = data[Notifier.OBJECT]

        # print('>>>_updateNmrResidues', nmrResidue)
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if trigger == Notifier.DELETE:
                # print('>>>delete nmrResidue', nmrResidue)
                self._deleteNmrResidues(nmrResidue)

            elif trigger == Notifier.CREATE:
                # print('>>>create nmrResidue', nmrResidue)
                if not self._createNmrResidues(nmrResidue):
                    # print('>>>error? redraw list')
                    # self.setNmrChainDisplay(self.nmrChain)
                    pass

            elif trigger == Notifier.RENAME:
                oldPid = data[Notifier.OLDPID]
                self._renameNmrResidue(nmrResidue, oldPid)

            # elif trigger == Notifier.CHANGE:
            #     print('>>>change nmrResidue - no action', nmrResidue)
            #
            # elif trigger == Notifier.OBSERVE:
            #     print('>>>observe nmrResidue - no action', nmrResidue)

    def _changeNmrResidues(self, data):
        """Update the nmrResidues in the display.
        """
        nmrResidue = data[Notifier.OBJECT]

        # print('>>>_changeNmrResidues', nmrResidue)
        trigger = data[Notifier.TRIGGER]

        try:
            with self.sceneBlocking():
                if nmrResidue in self.nmrChain.nmrResidues and nmrResidue not in self.nmrResidueList.guiNmrResidues:
                    # print('>>>change nmrResidue - create', nmrResidue)
                    if not self._createNmrResidues(nmrResidue):
                        # print('>>>error? redraw list')
                        # self.setNmrChainDisplay(self.nmrChain)
                        pass

                elif nmrResidue in self.nmrResidueList.guiNmrResidues and nmrResidue not in self.nmrChain.nmrResidues:
                    # not in chain, but in residues as other chain
                    self._deleteGuiNmrResidues(nmrResidue)
                    # self._deleteNmrResidues(nmrResidue)

                else:
                    # print('>>>change2 nmrResidue - create **** rename', nmrResidue)

                    # this is the event that fires on a name change
                    # self._deleteBadNmrResidues(nmrResidue)
                    self._deleteNmrResidues(nmrResidue)
                    if not self._createNmrResidues(nmrResidue):
                        # print('>>>error? redraw list')
                        # self.setNmrChainDisplay(self.nmrChain)
                        pass

        except Exception as es:
            # strange error not traced yet, interesting, but not fatal if trapped - think I've found it
            getLogger().warning(str(es))

    def _updateNmrAtoms(self, data):
        """Update the nmrAtoms in the display.
        """

        # Done
        nmrAtom = data[Notifier.OBJECT]
        trigger = data[Notifier.TRIGGER]

        # only mainNmrResidues are shown on the screen
        with self.sceneBlocking():
            if trigger == Notifier.DELETE:
                # print('>>>delete nmrAtom', nmrAtom)
                self.nmrResidueList.deleteNmrAtoms(nmrAtom)

            elif trigger == Notifier.CREATE:
                # print('>>>create nmrAtom', nmrAtom)
                self.nmrResidueList.createNmrAtoms(nmrAtom)

    def _renameNmrResidue(self, nmrResidue, oldPid: str):
        """Reset pid for NmrResidue and all offset NmrResidues
        """

        with self.sceneBlocking():
            if nmrResidue in self.nmrResidueList.guiNmrResidues:
                self.nmrResidueList.guiNmrResidues[nmrResidue].nmrResidueLabel._update()
            if nmrResidue in self.nmrResidueList.guiGhostNmrResidues:
                self.nmrResidueList.guiGhostNmrResidues[nmrResidue].nmrResidueLabel._update()

            # nmrChainPid = self.nmrChainPulldown.getText()
            # if self.project.getByPid(nmrChainPid):
            #
            #     for nr in [nmrResidue] + list(nmrResidue.offsetNmrResidues):
            #         for guiNmrResidueGroup in self.nmrResidueList.guiNmrResidues.values():
            #             if guiNmrResidueGroup.nmrResidue is nr:
            #                 guiNmrResidueGroup.nmrResidueLabel._update()

    def _createNmrResidues(self, nmrResidues):
        """Create new nmrResidues in the scene
        """
        # print('>>>_createNmrResidues')

        # get the nmrResidue in the current selected chain
        nmrResidues = makeIterableList(nmrResidues)

        for nmrChainId in self.nmrResidueList.nmrChains.keys():
            thisResList = [nmrResidue for nmrResidue in nmrResidues if nmrResidue.nmrChain.pid == nmrChainId]
            if thisResList:
                self._buildNmrResidues(nmrChainId, thisResList)

        return True

    def _deleteNmrResidues(self, nmrResidues):
        """Delete the nmrResidue from the scene
        """
        # print('>>>_deleteNmrResidues')

        # get the nmrResidue in the current selected chain
        nmrResidues = makeIterableList(nmrResidues)

        for nmrChainId in self.nmrResidueList.nmrChains.keys():
            thisResList = [nmrResidue for nmrResidue in nmrResidues if nmrResidue.nmrChain.pid == nmrChainId]
            if thisResList:
                self._removeNmrResidues(nmrChainId, thisResList)

        return True

    def _removeNmrResidues(self, nmrChainId, nmrResidues):
        """Delete the nmrResidue from the scene
        """
        # print('>>>  _removeNmrResidues')

        nmrResidues = makeIterableList(nmrResidues)

        for nmrResidue in nmrResidues:

            # this SHOULD never fail
            if nmrResidue not in self.nmrResidueList.nmrChains[nmrChainId]:
                continue

            ii = self.nmrResidueList.nmrChains[nmrChainId].index(nmrResidue)

            # take the current nmrList, and remove the deleted nmrResidue
            # with the mainNmrResidues - should be correct order or empty
            delResSet = list(OrderedSet(self.nmrResidueList.nmrChains[nmrChainId]) -
                             OrderedSet([nmrResidue]))
            # print('>>>    delResSet', delResSet)
            self.nmrResidueList.nmrChains[nmrChainId] = delResSet

            if ii > 0:
                # rebuild peak lines, etc, for the previous nmrResidue - may need to check for +1 offsets
                previousNmrResidue = self.nmrResidueList.nmrChains[nmrChainId][ii - 1]
                self._rebuildNmrResidues(nmrChainId, previousNmrResidue)

            self._cleanupNmrResidue(nmrResidue)

        # update the prediction in the sequenceModule
        self.predictSequencePosition(self.nmrResidueList.nmrChains[nmrChainId])

        # # ignore if not in the visible chain
        # ii = self.nmrResidueList.getIndexNmrResidue(nmrResidue)
        # if ii is not None:
        #     self.nmrResidueList.deleteNmrResidueIndex(ii)
        #
        #     if ii:
        #         previousNmrResidue, guiAtom = self.nmrResidueList._getNmrResiduePair(ii - 1)
        #         self._rebuildNmrResidues(previousNmrResidue)

        # ii = self.predictedStretch.index(nmrResidue)
        # # nmrResiduesToUpdate = self.predictedStretch[max(0, ii - 1):min(ii + 1, len(self.predictedStretch))]
        #
        # # remove from the visible list
        # self.predictedStretch.pop(ii)
        # self.guiResiduesShown.pop(ii)
        #
        # if ii:
        #     previousNmrResidue = self.predictedStretch[ii - 1]
        #     self._rebuildNmrResidues(previousNmrResidue)

        # if nmrResidue.previousNmrResidue:
        #     previousNmrResidue = nmrResidue.previousNmrResidue.mainNmrResidue
        #     self._rebuildNmrResidues(previousNmrResidue)

        # guiNmrAtomSet = set([self.nmrResidueList.guiNmrAtoms[nmrAtom] for nmrAtom in nmrResidue.nmrAtoms
        #                      if nmrAtom in self.nmrResidueList.guiNmrAtoms])
        #
        # for guiAtom in guiNmrAtomSet:
        #     for guiLineList in self.nmrResidueList.assignmentLines.values():
        #         guiLines = [guiLine for guiLine in guiLineList
        #                     if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]
        #
        #         # remove all graphic lines
        #         for guiLine in guiLines:
        #             guiLineList.remove(guiLine)
        #             if guiLine in self.scene.items():
        #                 self.scene.removeItem(guiLine)
        #
        #     for guiLineList in self.nmrResidueList.connectingLines.values():
        #         guiLines = [guiLine for guiLine in guiLineList
        #                     if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]
        #
        #         # remove all graphic lines
        #         for guiLine in guiLines:
        #             guiLineList.remove(guiLine)
        #             if guiLine in self.scene.items():
        #                 self.scene.removeItem(guiLine)
        #
        #     # clear connectivity list of guiNmrAtoms
        #     guiAtom.clearConnectedList()
        #
        # self.scene.removeItem(self.nmrResidueList.guiNmrResidues[nmrResidue])

    def _cleanupNmrResidue(self, nmrResidue):
        """remove guiNmrAtoms and guiNmrResidue from scene
        Clean up dicts in nmrResidueList
        """
        guiNmrAtomSet = set([self.nmrResidueList.guiNmrAtoms[nmrAtom] for nmrAtom in nmrResidue.nmrAtoms
                             if nmrAtom in self.nmrResidueList.guiNmrAtoms])

        for guiAtom in guiNmrAtomSet:
            for guiLineList in self.nmrResidueList.assignmentLines.values():
                guiLines = [guiLine for guiLine in guiLineList
                            if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]

                # remove all graphic lines
                for guiLine in guiLines:
                    guiLineList.remove(guiLine)
                    if guiLine in self.scene.items():
                        self.scene.removeItem(guiLine)

            for guiLineList in self.nmrResidueList.connectingLines.values():
                guiLines = [guiLine for guiLine in guiLineList
                            if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]

                # remove all graphic lines
                for guiLine in guiLines:
                    guiLineList.remove(guiLine)
                    if guiLine in self.scene.items():
                        self.scene.removeItem(guiLine)

            # clear connectivity list of guiNmrAtoms
            guiAtom.clearConnectedList()

        self.scene.removeItem(self.nmrResidueList.guiNmrResidues[nmrResidue])

        del self.nmrResidueList.guiNmrResidues[nmrResidue]
        for nmrAtom in nmrResidue.nmrAtoms:
            if nmrAtom in self.nmrResidueList.guiNmrAtoms:
                del self.nmrResidueList.guiNmrAtoms[nmrAtom]

    def _deleteGuiNmrResidues(self, nmrResidues):
        """Delete items from an old nmrResidue.
        """
        # print('>>>_deleteGuiNmrResidues')

        nmrResidues = makeIterableList(nmrResidues)

        for nmrResidue in nmrResidues:
            self._cleanupNmrResidue(nmrResidue)

        # nmrChains contain nmrResidue in the wrong place and must be removed before sequence prediction
        for nmrChainId, nmrList in self.nmrResidueList.nmrChains.items():
            thisResList = [nmrResidue for nmrResidue in nmrList if nmrResidue.nmrChain.pid == nmrChainId]
            if thisResList:
                # update the prediction in the sequenceModule
                self.predictSequencePosition(thisResList)

            # put all the guiResidueGroups in the correct positions
            self.nmrResidueList.updateGuiResiduePositions(nmrChainId, updateMainChain=True, updateConnectedChains=True)

            # guiNmrAtomSet = set([self.nmrResidueList.guiNmrAtoms[nmrAtom] for nmrAtom in nmrResidue.nmrAtoms
            #                      if nmrAtom in self.nmrResidueList.guiNmrAtoms])
            #
            # for guiAtom in guiNmrAtomSet:
            #     for guiLineList in self.nmrResidueList.assignmentLines.values():
            #         guiLines = [guiLine for guiLine in guiLineList
            #                     if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]
            #
            #         # remove all graphic lines
            #         for guiLine in guiLines:
            #             guiLineList.remove(guiLine)
            #             if guiLine in self.scene.items():
            #                 self.scene.removeItem(guiLine)
            #
            #     for guiLineList in self.nmrResidueList.connectingLines.values():
            #         guiLines = [guiLine for guiLine in guiLineList
            #                     if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]
            #
            #         # remove all graphic lines
            #         for guiLine in guiLines:
            #             guiLineList.remove(guiLine)
            #             if guiLine in self.scene.items():
            #                 self.scene.removeItem(guiLine)
            #
            #     # clear connectivity list of guiNmrAtoms
            #     guiAtom.clearConnectedList()
            #
            # self.scene.removeItem(self.nmrResidueList.guiNmrResidues[nmrResidue])

            # oldInd = self.nmrResidueList.getIndexNmrResidue(nmrResidue)
            # ii = self.nmrChain.mainNmrResidues
            # self._storeIndex(nmrResidue, oldInd)

            # self.nmrResidueList.deleteNmrResidue(nmrResidue)
            # del self.nmrResidueList.guiNmrResidues[nmrResidue]
            # for nmrAtom in nmrResidue.nmrAtoms:
            #     del self.nmrResidueList.guiNmrAtoms[nmrAtom]

    # def _storeIndex(self, nmrResidue, index):
    #     """Store the index for the undo.
    #     """
    #     if not nmrResidue in self._deleteStore:
    #         self._deleteStore[nmrResidue] = []
    #     self._deleteStore[nmrResidue].append(index)

    # def _restoreIndex(self, nmrResidue):
    #     """Retrieve the index of the deleted item.
    #     """
    #     if nmrResidue in self._deleteStore:
    #         return self._deleteStore[nmrResidue].pop()
    #
    #     return None

    # def _deleteBadNmrResidues(self, nmrResidues):
    #     """Delete the nmrResidue from the scene
    #     """
    #     print('>>>_deleteBadNmrResidues')
    #
    #     nmrResidues = makeIterableList(nmrResidues)
    #
    #     for nmrResidue in nmrResidues:
    #
    #         ii = self.nmrResidueList.getIndexNmrResidue(nmrResidue)
    #         if ii is not None:
    #             self.nmrResidueList.deleteNmrResidueIndex(ii)
    #
    #             # # ignore if not in the visible chain
    #             # if nmrResidue in self.predictedStretch:
    #             #
    #             #     ii = self.predictedStretch.index(nmrResidue)
    #             #     # nmrResiduesToUpdate = self.predictedStretch[max(0, ii - 1):min(ii + 1, len(self.predictedStretch))]
    #             #
    #             #     # remove from the visible list
    #             #     self.predictedStretch.pop(ii)
    #             #     self.guiResiduesShown.pop(ii)
    #
    #             if ii:
    #                 # previousNmrResidue = self.predictedStretch[ii - 1]
    #                 previousNmrResidue, guiAtom = self.nmrResidueList._getNmrResiduePair(ii - 1)
    #                 self._rebuildNmrResidues(previousNmrResidue)
    #
    #             # if nmrResidue.previousNmrResidue:
    #             #     previousNmrResidue = nmrResidue.previousNmrResidue.mainNmrResidue
    #             #     self._rebuildNmrResidues(previousNmrResidue)
    #
    #             guiNmrAtomSet = set([self.nmrResidueList.guiNmrAtoms[nmrAtom] for nmrAtom in nmrResidue.nmrAtoms])
    #
    #             for guiAtom in guiNmrAtomSet:
    #                 for guiLineList in self.nmrResidueList.assignmentLines.values():
    #                     guiLines = [guiLine for guiLine in guiLineList
    #                                 if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]
    #
    #                     # remove all graphic lines
    #                     for guiLine in guiLines:
    #                         guiLineList.remove(guiLine)
    #                         if guiLine in self.scene.items():
    #                             self.scene.removeItem(guiLine)
    #
    #                 for guiLineList in self.nmrResidueList.connectingLines.values():
    #                     guiLines = [guiLine for guiLine in guiLineList
    #                                 if guiLine.guiAtom1 is guiAtom or guiLine.guiAtom2 is guiAtom]
    #
    #                     # remove all graphic lines
    #                     for guiLine in guiLines:
    #                         guiLineList.remove(guiLine)
    #                         if guiLine in self.scene.items():
    #                             self.scene.removeItem(guiLine)
    #
    #                 # clear connectivity list of guiNmrAtoms
    #                 guiAtom.clearConnectedList()
    #
    #             self.scene.removeItem(self.nmrResidueList.guiNmrResidues[nmrResidue])

    # def _createNmrAtoms(self, nmrAtoms):
    #     """Create new peak lines associated with the created/undeleted nmrAtoms.
    #     """
    #     nmrAtoms = makeIterableList(nmrAtoms)
    #     peaks = [peak for nmrAtom in nmrAtoms for peak in nmrAtom.assignedPeaks]
    #
    #     for nmrAtom in nmrAtoms:
    #         if nmrAtom not in self.guiNmrAtomDict:
    #             self._addNmrAtomToGuiResidues(nmrAtom)
    #     self.nmrResidueList.rebuildPeakLines(peaks, rebuildPeakLines=True, makeListFromPeak=True)
    #
    # def _deleteNmrAtoms(self, nmrAtoms):
    #     """Delete peak lines associated with the deleted nmrAtoms.
    #     """
    #     nmrAtoms = makeIterableList(nmrAtoms)
    #     peakLines = self._searchPeakLines(nmrAtoms, includeDeleted=True)
    #     peaks = [peakLine._peak for peakLine in peakLines]
    #
    #     self.nmrResidueList.rebuildPeakLines(peaks, rebuildPeakLines=True)

    # def _removeLinesFromScene(self, lineDist):
    #     """Remove all the peakLines from the scene.
    #     """
    #     for lineList in lineDist.values():
    #         for line in lineList:
    #             self.scene.removeItem(line)
    #     lineDist.clear()

    # def _searchPeakLines(self, nmrAtoms, includeDeleted=False):
    #     """Return a list of the peakLines containing one of the nmrAtoms in the list.
    #     """
    #     peakLines = []
    #     for lineList in self.assignmentLines.values():
    #         for line in lineList:
    #             nmrAtom = line.guiAtom1.nmrAtom if line.guiAtom1 else None
    #             if nmrAtom in nmrAtoms and (includeDeleted or not (nmrAtom.isDeleted or nmrAtom._flaggedForDelete)):
    #                 peakLines.append(line)
    #             nmrAtom = line.guiAtom2.nmrAtom if line.guiAtom2 else None
    #             if nmrAtom in nmrAtoms and (includeDeleted or not (nmrAtom.isDeleted or nmrAtom._flaggedForDelete)):
    #                 peakLines.append(line)
    #
    #     return peakLines

    def _updateEndPoints(self, lineDict):
        """Update the end points from the dict.
        """
        for lineList in lineDict.values():
            for line in lineList:
                try:
                    line.updateEndPoints()
                except Exception as es:
                    pass

    def _buildNmrResidues(self, nmrChainId, nmrResidueList):
        """Build the new residues in the list, inserting into the predicted stretch at the correct index.
        """
        # print('>>>  _buildNmrResidues')

        # this SHOULD be a list of only one item but use like this just to be sure
        for nmrResidue in nmrResidueList:
            if nmrResidue is nmrResidue.mainNmrResidue:

                # take the current nmrList, and the new nmrResidue, and get intersection
                # with the mainNmrResidues - should be correct order or empty
                newResSet = list((OrderedSet(self.nmrResidueList.nmrChains[nmrChainId]) |
                                  OrderedSet([nmrResidue])) & \
                                 OrderedSet(nmrResidue.nmrChain.mainNmrResidues))  # keeps the ordering of the second set
                # print('>>>    newResSet', newResSet)

                # if not showing all nmrResidues then get the connected stretch from the current nmrResidue
                if not self.nmrResiduesCheckBox.isChecked():
                    # get the connected stretch of mainNmrResidues
                    if self.current.nmrResidue:
                        mainNmrRes = self.current.nmrResidue.mainNmrResidue
                        if mainNmrRes in newResSet:
                            indL = indR = newResSet.index(mainNmrRes)
                            while newResSet[indL].previousNmrResidue and indL > 0:
                                indL -= 1
                            while newResSet[indR].nextNmrResidue and indR < len(newResSet):
                                indR += 1
                            newResSet = newResSet[indL:indR + 1]

                if not newResSet:
                    return

                # update to the new list
                self.nmrResidueList.nmrChains[nmrChainId] = newResSet

                # iterate through and and add new residues that don't already exists
                for ii, nmrRes in enumerate(newResSet):
                    # do not _insertNmrRes as the list is built
                    self.nmrResidueList.addNmrResidue(nmrChainId, nmrRes, ii, _insertNmrRes=False)

        # add the connecting lines
        # guiNmrResidues = [self.guiNmrResidues[nmrResidue] for nmrResidue in nmrResidueList if nmrResidue is nmrResidue.mainNmrResidue]

        # add the connecting lines
        # self.nmrResidueList.addConnectionsBetweenGroups()

        # for ii, res in enumerate(guiNmrResidues):
        #     if not self.nmrResiduesCheckBox.isChecked() or ii in connectingLinesNeeded:
        #         self._addConnectingLineToGroup(tuple(self.guiNmrResidues.values())[ii],
        #                                        res['C'], self.guiResiduesShown[ii + 1]['N'],
        #                                        self._lineColour, 1.0, lineList=self.connectingLines, lineId=res)

        # # add connecting lines
        # newList = OrderedSet()
        # connectsNeeded = self._module.nmrResiduesCheckBox.isChecked()
        #
        # for nmr in nmrResidueList:
        #     newList.add(nmr)
        #     if nmr.nextNmrResidue and nmr.nextNmrResidue.mainNmrResidue:
        #         newList.add(nmr.nextNmrResidue.mainNmrResidue)
        #
        #         nextNmr = nmr.nextNmrResidue.mainNmrResidue
        #
        #         indThis = self.nmrResidueList.getIndexNmrResidue(nmr)
        #         indNext = self.nmrResidueList.getIndexNmrResidue(nextNmr)
        #
        #         if indThis is not None and indNext is not None:
        #             # add connecting line
        #
        #             prev = self.nmrResidueList._getNmrResiduePair(indThis)
        #             this = self.nmrResidueList._getNmrResiduePair(indNext)
        #
        #             # add the connection the minus residue and point to the right - may need to change for +1 residues
        #             prevRes, prevGuiAtoms = prev
        #             thisRes, thisGuiAtoms = this
        #
        #             if (prevRes.nextNmrResidue and prevRes.nextNmrResidue is thisRes) and connectsNeeded:
        #                 # connect from this 'N' to the previous 'C'
        #                 self._addConnectingLineToGroup(self.guiNmrResidues[prevRes],
        #                                                prevGuiAtoms['C'], thisGuiAtoms['N'],
        #                                                self._lineColour, 1.0, lineList=self.connectingLines,
        #                                                lineId=thisRes)

        self.nmrResidueList.rebuildNmrResidues(nmrResidueList)

        # update the prediction in the sequenceModule
        self.predictSequencePosition(self.nmrResidueList.nmrChains[nmrChainId])

        return True

        # # add the peakAssignments
        # if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
        #     for nmrResidue in nmrResidueList:
        #         if nmrResidue is nmrResidue.mainNmrResidue:
        #             # add the internally connected Lines
        #             internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue)
        #             self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
        #             self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
        #             self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
        #                                                         self.assignmentLines, self.connectingLines)
        #
        # self._updateGuiResiduePositions(updateMainChain=True, updateConnectedChains=True)

    def removeNmrChainNotifiers(self):
        """Remove notifiers that are set on nmrChains.
        """
        nmrChains = tuple(self.project.nmrChains)
        foundNotifiers = self.searchNotifiers(objects=nmrChains, triggers=[Notifier.OBSERVE], targetName='nmrResidues')
        for notifier in foundNotifiers:
            # print('>>>deleting notifier', notifier)
            self.deleteNotifier(notifier)

    def addNmrChainNotifiers(self):
        """Add new notifiers for all nmrChains in the project.
        """
        for nmrChain in self.project.nmrChains:
            self.setNotifier(nmrChain, triggers=[Notifier.OBSERVE], targetName='nmrResidues',
                             callback=self._changeNmrResidues)

    def resetScene(self):
        """Reset all gui items and data in the scene.
        """
        self.nmrResidueList.reset()
        self.scene.clear()
        self.scene.setSceneRect(self.scene.itemsBoundingRect())
        self.thisSequenceModule._initialiseChainLabels()

    def setNmrChain(self, nmrChain):
        self.nmrResidueList.nmrChain = nmrChain

    def setNmrChainDisplay(self, nmrChainOrPid):

        # print('>>>setNmrChainDisplay')

        if isinstance(nmrChainOrPid, str):
            if not Pid.isValid(nmrChainOrPid):
                self.resetScene()
                return

            nmrChain = self.project.getByPid(nmrChainOrPid)
        else:
            nmrChain = nmrChainOrPid

        # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
        if not nmrChain:
            self.resetScene()
            return

        self.nmrChain = nmrChain
        thisChainId = nmrChain.pid

        with notificationEchoBlocking():

            # currently only handles one visible nmrChain at a time - but changing to a dict
            self.resetScene()
            self.setNmrChain(nmrChain)

            # self.removeNmrChainNotifiers()
            # self.addNmrChainNotifiers()

            nmrList = nmrChain.mainNmrResidues

            if self.nmrResiduesCheckBox.isChecked():
                nmrList = nmrChain.mainNmrResidues

            if not self.nmrResiduesCheckBox.isChecked():

                # get the connected stretch of mainNmrResidues
                if self.current.nmrResidue:
                    mainNmrRes = self.current.nmrResidue.mainNmrResidue
                    if mainNmrRes in nmrList:
                        indL = indR = nmrList.index(mainNmrRes)
                        while nmrList[indL].previousNmrResidue and indL > 0:
                            indL -= 1
                        while nmrList[indR].nextNmrResidue and indR < len(nmrList):
                            indR += 1
                        nmrList = nmrList[indL:indR + 1]

            # add the nmrResidues to the scene
            for ii, nmrRes in enumerate(nmrList):
                self.nmrResidueList.addNmrResidue(thisChainId, nmrRes, index=ii)

            # add the connecting lines
            self.nmrResidueList.addConnectionsBetweenGroups(thisChainId)

            # add the peakAssignment lines
            self.nmrResidueList._addAllPeakAssignments(thisChainId)

            # put all the guiResidueGroups in the correct positions
            self.nmrResidueList.updateGuiResiduePositions(thisChainId, updateMainChain=True, updateConnectedChains=True)

            # update the prediction in the sequenceModule
            if thisChainId in self.nmrResidueList.nmrChains:
                self.predictSequencePosition(self.nmrResidueList.nmrChains[thisChainId])

    def showNmrChainFromPulldown(self, data=None):
        """Clear and redraw the nmrChain selected from the pulldown.
        """
        # print('>>>showNmrChainFromPulldown')

        nmrChainPid = self.nmrChainPulldown.getText()
        if nmrChainPid:
            with self.sceneBlocking():
                self.setNmrChainDisplay(nmrChainPid)

            # check whther to update self.current.nmrChain
            self._setCurrentNmrChain(nmrChainPid)
        else:
            # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
            self.resetScene()

    def _setCurrentNmrChain(self, nmrChainOrPid):

        if isinstance(nmrChainOrPid, str):
            if not Pid.isValid(nmrChainOrPid):
                return

            nmrChain = self.project.getByPid(nmrChainOrPid)
        else:
            nmrChain = nmrChainOrPid

        # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
        if not nmrChain:
            self.resetScene()
            return

        checkBox = self._SGwidget._getCheckBox(LINKTOPULLDOWNCLASS)
        if self.current.nmrChain and self.current.nmrChain != nmrChain and checkBox and checkBox.isChecked():
            self.current.nmrChain = nmrChain

    def resetSequenceGraph(self):
        """Reset the module to the default nmrChain.
        """
        self.nmrChainPulldown.pulldownList.select('NC:@-')

    def _closeModule(self):
        """CCPN-INTERNAL: used to close the module
        """
        # self._unRegisterNotifiers()
        self.thisSequenceModule.close()
        super()._closeModule()

    def close(self):
        """Close the table from the commandline
        """
        self._closeModule()

    def unlinkNearestNmrResidue(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            if self.current.nmrResidue.mainNmrResidue.previousNmrResidue:
                with progressManager(self.mainWindow, 'unlinking Previous NmrResidue to:\n ' + selected):
                    try:
                        self.current.nmrResidue.unlinkPreviousNmrResidue()
                    except Exception as es:
                        showWarning(str(self.windowTitle()), str(es))
                        if self.application._isInDebugMode:
                            raise es
            elif self.current.nmrResidue.mainNmrResidue.previousNmrResidue:
                with progressManager(self.mainWindow, 'unlinking Next NmrResidue to:\n ' + selected):
                    try:
                        self.current.nmrResidue.unlinkNextNmrResidue()
                    except Exception as es:
                        showWarning(str(self.windowTitle()), str(es))
                        if self.application._isInDebugMode:
                            raise es

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()

    def disconnectPreviousNmrResidue(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            with progressManager(self.mainWindow, 'disconnecting Previous NmrResidue to:\n ' + selected):
                try:
                    self.current.nmrResidue.disconnectPrevious()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))
                    if self.application._isInDebugMode:
                        raise es

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()

    def disconnectNmrResidue(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            with progressManager(self.mainWindow, 'disconnecting NmrResidue:\n ' + selected):
                try:
                    self.current.nmrResidue.disconnect()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))
                    if self.application._isInDebugMode:
                        raise es

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()

    def disconnectNextNmrResidue(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            with progressManager(self.mainWindow, 'disconnecting Next NmrResidue to:\n ' + selected):
                try:
                    self.current.nmrResidue.disconnectNext()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))
                    if self.application._isInDebugMode:
                        raise es

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()

    def disconnectAllNmrResidues(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            with progressManager(self.mainWindow, 'disconnecting all NmrResidues connected to:\n ' + selected):
                try:
                    self.current.nmrResidue.disconnectAll()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))
                    if self.application._isInDebugMode:
                        raise es

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()

    def deassignNmrChain(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.nmrChain.pid)
            with progressManager(self.mainWindow, 'deassigning nmrResidues in NmrChain:\n ' + selected):
                try:
                    self.current.nmrResidue.deassignNmrChain()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))
                    if self.application._isInDebugMode:
                        raise es

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()

    def deassignPeak(self, selectedPeak=None, selectedNmrAtom=None):
        """Deassign the peak by removing the assigned nmrAtoms from the list
        """
        selectedPeak = self.project.getByPid(selectedPeak) if isinstance(selectedPeak, str) else selectedPeak
        if not isinstance(selectedPeak, Peak):
            raise TypeError('selectedPeak must be of type Peak')
        selectedNmrAtom = self.project.getByPid(selectedNmrAtom) if isinstance(selectedNmrAtom, str) else selectedNmrAtom
        if not isinstance(selectedNmrAtom, NmrAtom):
            raise TypeError('selectedNmrAtom must be of type NmrAtom')

        if selectedPeak:
            with undoBlock():
                try:
                    newList = []
                    # remove the nmrAtom from the list and replace with None
                    for atomList in selectedPeak.assignedNmrAtoms:
                        atoms = [(atom if atom != selectedNmrAtom else None) for atom in list(atomList)]
                        newList.append(tuple(atoms))

                    selectedPeak.assignedNmrAtoms = tuple(newList)

                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))
                    if self.application._isInDebugMode:
                        raise es

    def initialiseScene(self):
        """Replace the scene with a new one to reset the size of the scrollbars.
        """
        # Only needed to be done the first time, scene is resized at the end of setNmrChainDisplay
        self.scene = QtWidgets.QGraphicsScene(self)
        self.scrollContents = QtWidgets.QGraphicsView(self.scene, self)
        self.scrollContents.setRenderHints(QtGui.QPainter.Antialiasing)
        self.scrollContents.setInteractive(True)
        self.scrollContents.setGeometry(QtCore.QRect(0, 0, 300, 400))
        self.scrollContents.setAlignment(QtCore.Qt.AlignCenter)
        self._sequenceGraphScrollArea.setWidget(self.scrollContents)

    # def deassignNmrAtom(self, selectedNmrAtom=None):
    #     """Remove the selected peaks from the assignedPeaks list
    #     """
    #     selectedNmrAtom = self.project.getByPid(selectedNmrAtom) if isinstance(selectedNmrAtom, str) else selectedNmrAtom
    #     if not isinstance(selectedNmrAtom, NmrAtom):
    #         raise TypeError('selectedNmrAtom must be of type NmrAtom')
    #
    #     if selectedNmrAtom:
    #         with undoBlock():
    #             try:
    #                 atoms = list(selectedNmrAtom.assignedPeaks)
    #                 # selectedPeak.assignedNmrAtoms = ()
    #
    #                 # for peak in selectedNmrAtom.assignedPeaks:
    #                 #
    #                 #   allAtoms = list(peak.dimensionNmrAtoms)
    #                 #   for dim in range(len(peak.dimensionNmrAtoms)):
    #                 #     dimNmrAtoms = list(peak.dimensionNmrAtoms[dim])
    #                 #     if selectedNmrAtom in dimNmrAtoms:
    #                 #       dimNmrAtoms.remove(selectedNmrAtom)
    #                 #
    #                 #       # allAtoms = list(peak.dimensionNmrAtoms)
    #                 #       allAtoms[dim] = dimNmrAtoms
    #                 #
    #                 #   peak.dimensionNmrAtoms = allAtoms
    #
    #             except Exception as es:
    #                 showWarning(str(self.windowTitle()), str(es))
    #                 if self.application._isInDebugMode:
    #                     raise es

    # def _assembleGroupResidue(self, nmrResidue: NmrResidue, atoms: typing.Dict[str, GuiNmrAtom], pos=None, lineList=None):
    #     """Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
    #     creates a graphical representation of a residue in the assigner
    #     """
    #     guiResidueGroup = GuiNmrResidueGroup(self, nmrResidue, atoms['CA'], 0)  #atoms['H'].x())
    #     self.guiNmrResidues[nmrResidue] = guiResidueGroup
    #     self.scene.addItem(guiResidueGroup)
    #
    #     # add the atoms to the group and set the reverse link
    #     for item in atoms.values():
    #         guiResidueGroup.addToGroup(item)
    #         item.guiNmrResidueGroup = guiResidueGroup
    #
    #     # modify the group
    #     nmrAtoms = [atom.name for atom in nmrResidue.nmrAtoms]
    #     if "CB" in list(atoms.keys()):
    #         self._addConnectingLineToGroup(guiResidueGroup, atoms['CA'], atoms['CB'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #
    #     if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
    #         self._addConnectingLineToGroup(guiResidueGroup, atoms['H'], atoms['N'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #
    #     self._addConnectingLineToGroup(guiResidueGroup, atoms['N'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #     self._addConnectingLineToGroup(guiResidueGroup, atoms['C'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #
    #     self._addGroupResiduePredictions(guiResidueGroup, nmrResidue, atoms['CA'])
    #
    #     return guiResidueGroup

    # def _assembleGhostResidue(self, nmrResidue: NmrResidue, atoms: typing.Dict[str, GuiNmrAtom], lineList=None):
    #     """Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
    #     creates a graphical representation of a residue in the assigner
    #     """
    #     guiResidueGroup = GuiNmrResidueGroup(self, nmrResidue, atoms['CA'], 0)  #atoms['H'].x())
    #     self.guiGhostNmrResidues[nmrResidue] = guiResidueGroup
    #     self.scene.addItem(guiResidueGroup)
    #
    #     # add the atoms to the group and set the reverse link
    #     for item in atoms.values():
    #         guiResidueGroup.addToGroup(item)
    #         item.guiNmrResidueGroup = guiResidueGroup
    #
    #     # modify the group
    #     nmrAtoms = [atom.name for atom in nmrResidue.nmrAtoms]
    #     if "CB" in list(atoms.keys()):
    #         self._addConnectingLineToGroup(guiResidueGroup, atoms['CA'], atoms['CB'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #
    #     if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
    #         self._addConnectingLineToGroup(guiResidueGroup, atoms['H'], atoms['N'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #
    #     self._addConnectingLineToGroup(guiResidueGroup, atoms['N'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #     self._addConnectingLineToGroup(guiResidueGroup, atoms['C'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
    #
    #     # self._addGroupResiduePredictions(guiResidueGroup, nmrResidue, atoms['CA'])
    #
    #     return guiResidueGroup

    # def addSideChainAtoms(self, nmrResidue, cbAtom, atoms, colour, lineList):
    #     """Add the sideChain atoms and connecting lines above the backbone line.
    #     """
    #     residue = {}
    #     for k, v in ATOM_POSITION_DICT[nmrResidue.residueType].items():
    #         if k != 'boundAtoms':
    #             position = [cbAtom.x() + v[0], cbAtom.y() + v[1]]
    #             nmrAtom = nmrResidue.fetchNmrAtom(name=k)
    #             newAtom = self._createGuiNmrAtom(k, position, nmrAtom)
    #             self.scene.addItem(newAtom)
    #             residue[k] = newAtom
    #             atoms[k] = newAtom

    # self.guiNmrAtomDict[nmrAtom] = newAtom

    # for boundAtomPair in ATOM_POSITION_DICT[nmrResidue.residueType]['boundAtoms']:
    #     guiAtom1 = residue[boundAtomPair[0]]
    #     guiAtom2 = residue[boundAtomPair[1]]
    #     newLine = AssignmentLine(guiAtom1.x(), guiAtom1.y(), guiAtom2.x(), guiAtom2.y(), colour, 1.0)
    #     self.scene.addItem(newLine)

    # def _addNmrAtomToGuiResidues(self, nmrAtom, atomSpacing=None):
    #     """Update the guiNmrAtoms for a newly created/undeleted nmrAtom.
    #     Assumes that the nmrAtom/guiNmrAtom does not exist
    #     """
    #     nmrResidue = nmrAtom.nmrResidue
    #
    #     atoms = {}
    #     if atomSpacing:
    #         self.atomSpacing = atomSpacing
    #     nmrAtoms = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
    #     residueAtoms = DEFAULT_RESIDUE_ATOMS.copy()
    #
    #     for k, v in residueAtoms.items():
    #         if k in nmrAtoms:
    #             fetchedNmrAtom = nmrResidue.fetchNmrAtom(name=k)
    #         else:
    #             fetchedNmrAtom = None
    #         if fetchedNmrAtom is nmrAtom:
    #             atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)
    #
    #     if nmrResidue in self.predictedStretch:
    #         ii = self.predictedStretch.index(nmrResidue)
    #         self.guiResiduesShown[ii].update(atoms)
    #
    #     if nmrResidue in self.guiNmrResidues:
    #         guiResidueGroup = self.guiNmrResidues[nmrResidue]
    #
    #         # add the atoms to the group and set the reverse link
    #         for item in atoms.values():
    #             guiResidueGroup.addToGroup(item)
    #             item.guiNmrResidueGroup = guiResidueGroup

    # def addResidue(self, nmrResidue: NmrResidue, nmrResidueIndex: int, atomSpacing=None, lineList=None):
    #     """Takes an Nmr Residue and a direction, either '-1 or '+1', and adds a residue to the sequence graph
    #     corresponding to the Nmr Residue.
    #     Nmr Residue name displayed beneath CA of residue drawn and residue type predictions displayed
    #     beneath Nmr Residue name
    #     """
    #     atoms = {}
    #     pos = np.array([0.0, 0.0])
    #     if atomSpacing:
    #         self.atomSpacing = atomSpacing
    #     nmrAtoms = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
    #     residueAtoms = DEFAULT_RESIDUE_ATOMS.copy()
    #
    #     if nmrResidue.residueType == 'GLY':
    #         # GLY doesn't have CB
    #         del residueAtoms['CB']
    #
    #     # add the new nmrResidue to the current list
    #     if not self.guiResiduesShown:
    #         for k, v in residueAtoms.items():
    #             if k in nmrAtoms:
    #                 nmrAtom = nmrResidue.fetchNmrAtom(name=k)
    #             else:
    #                 nmrAtom = None
    #             atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)
    #         self.guiResiduesShown.append(atoms)
    #         self.predictedStretch.append(nmrResidue)
    #
    #         # add the sideChain atoms
    #         if self._SGwidget.checkBoxes['showSideChain']['checkBox'].isChecked():
    #             if 'CB' in residueAtoms and nmrResidue.residueType:
    #                 cbAtom = atoms['CB']
    #                 self.addSideChainAtoms(nmrResidue, cbAtom, atoms, self._lineColour, lineList)
    #
    #     else:
    #         for k, v in residueAtoms.items():
    #             if k in nmrAtoms:
    #                 nmrAtom = nmrResidue.fetchNmrAtom(name=k)
    #             else:
    #                 nmrAtom = None
    #             atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)
    #
    #         # # insert into the list at the correct position
    #         # if direction == '-1':
    #         #     self.guiResiduesShown.insert(0, atoms)
    #         #     self.predictedStretch.insert(0, nmrResidue)
    #         # else:
    #         #     self.guiResiduesShown.append(atoms)
    #         #     self.predictedStretch.append(nmrResidue)
    #
    #         self.guiResiduesShown.insert(nmrResidueIndex, atoms)
    #         self.predictedStretch.insert(nmrResidueIndex, nmrResidue)
    #
    #         # add the sideChain atoms
    #         if self._SGwidget.checkBoxes['showSideChain']['checkBox'].isChecked():
    #             if 'CB' in residueAtoms and nmrResidue.residueType:
    #                 cbAtom = atoms['CB']
    #                 self.addSideChainAtoms(nmrResidue, cbAtom, atoms, self._lineColour, lineList)
    #
    #     newGuiResidueGroup = self._assembleGroupResidue(nmrResidue, atoms, lineList=lineList)  #, pos[0])
    #
    #     return newGuiResidueGroup

    # def _addGroupResiduePredictions(self, group: GuiNmrResidueGroup, nmrResidue: NmrResidue, caAtom: GuiNmrAtom):
    #     """Gets predictions for residue type based on BMRB statistics and determines label positions
    #     based on caAtom position.
    #     """
    #     predictions = list(set(map(tuple, (getNmrResiduePrediction(nmrResidue, self.project.chemicalShiftLists[0])))))
    #     predictions.sort(key=lambda a: float(a[1][:-1]), reverse=True)
    #     for prediction in predictions:
    #         predictionLabel = QtWidgets.QGraphicsTextItem()
    #         predictionLabel.setPlainText(prediction[0] + ' ' + prediction[1])
    #         predictionLabel.setDefaultTextColor(QtGui.QColor(self._textColour))
    #         predictionLabel.setFont(textFontSmallBold)
    #         predictionLabel.setPos(caAtom.x() - caAtom.boundingRect().width() / 2,
    #                                caAtom.y() + (30 * (predictions.index(prediction) + 2)))
    #         group.addToGroup(predictionLabel)

    # def _addResiduePredictions(self, nmrResidue: NmrResidue, caAtom: GuiNmrAtom):
    #     """
    #     Gets predictions for residue type based on BMRB statistics and determines label positions
    #     based on caAtom position.
    #     """
    #     predictions = list(set(map(tuple, (getNmrResiduePrediction(nmrResidue, self.project.chemicalShiftLists[0])))))
    #     predictions.sort(key=lambda a: float(a[1][:-1]), reverse=True)
    #     for prediction in predictions:
    #         predictionLabel = QtWidgets.QGraphicsTextItem()
    #         predictionLabel.setPlainText(prediction[0] + ' ' + prediction[1])
    #         predictionLabel.setDefaultTextColor(QtGui.QColor(self._textColour))
    #         predictionLabel.setFont(textFontSmallBold)
    #         predictionLabel.setPos(caAtom.x() - caAtom.boundingRect().width() / 2,
    #                                caAtom.y() + (30 * (predictions.index(prediction) + 2)))
    #         self.scene.addItem(predictionLabel)

    def predictSequencePosition(self, nmrResidueList: list):
        """
        Predicts sequence position for Nmr residues displayed in the Assigner and highlights appropriate
        positions in the Sequence Module if it is displayed.
        """
        if len(nmrResidueList) < 3:
            self.thisSequenceModule._initialiseChainLabels()
            return

        if self.project.chains and self.project.chemicalShiftLists:

            nmrResidues = nmrResidueList  # [item[0] for item in nmrResidueList]

            matchesDict = {}
            for chainNum, chain in enumerate(self.project.chains):
                matchesDict[chainNum] = []
                for chemList in self.project.chemicalShiftLists:
                    match = getSpinSystemsLocation(self.project, nmrResidues,
                                                   chain, chemList)
                    if match:
                        matchesDict[chainNum].append(match)

            for chainNum in matchesDict.keys():

                # possibleMatches = getSpinSystemsLocation(self.project, nmrResidues,
                #                   self.project.chains[0], self.project.chemicalShiftLists[0])

                self.thisSequenceModule._clearStretches(chainNum)
                possibleMatches = matchesDict[chainNum]

                if possibleMatches:
                    for chemList in possibleMatches:
                        for possibleMatch in chemList:
                            if possibleMatch[0] > 1 and not len(possibleMatch[1]) < len(nmrResidues):
                                # if hasattr(self.application, 'sequenceModule'):
                                # self.application.sequenceModule._highlightPossibleStretches(possibleMatch[1])

                                self.thisSequenceModule._highlightPossibleStretches(chainNum, possibleMatch[1])

    def _toggleSequence(self):
        if not self.sequenceCheckBox.isChecked():
            self._sequenceModuleFrame.hide()
        else:
            self._sequenceModuleFrame.show()

    def _addConnectingLineToGroup(self, group: GuiNmrResidueGroup, guiAtom1: GuiNmrAtom, guiAtom2: GuiNmrAtom,
                                  colour: str, width: float, displacement: float = None, style: str = None,
                                  peak: Peak = None, lineList=None, lineId='default'):
        """Adds a line between two GuiNmrAtoms using the width, colour, displacement and style specified.
        """
        newLine = AssignmentLine(0, 0, 0, 0, colour, width,
                                 parent=self, style=style, peak=peak,
                                 guiAtom1=guiAtom1, guiAtom2=guiAtom2, displacement=displacement)

        # not sure why this is different
        # group.addToGroup(newLine)
        newLine.setParentItem(group)

        itemKey = id(lineId)
        if itemKey not in lineList:
            lineList[itemKey] = []
        lineList[itemKey].append(newLine)
        return newLine

    # def _createGuiNmrAtom(self, atomType: str, position: tuple, nmrAtom: NmrAtom = None) -> GuiNmrAtom:
    #     """Creates a GuiNmrAtom specified by the atomType and graphical position supplied.
    #     GuiNmrAtom can be linked to an NmrAtom by supplying it to the function.
    #     """
    #     atom = GuiNmrAtom(text=atomType, pos=position, nmrAtom=nmrAtom)
    #     self.guiNmrAtomDict[nmrAtom] = atom
    #     return atom
    #
    # def _createGhostGuiNmrAtom(self, atomType: str, position: tuple, nmrAtom: NmrAtom = None) -> GuiNmrAtom:
    #     """Creates a GuiNmrAtom specified by the atomType and graphical position supplied.
    #     GuiNmrAtom can be linked to an NmrAtom by supplying it to the function.
    #     """
    #     atom = GuiNmrAtom(text=atomType, pos=position, nmrAtom=nmrAtom)
    #     self.guiNmrAtomDict[nmrAtom] = atom
    #     return atom

    # def _getPeakAssignmentsForResidue(self, nmrResidue, nmrAtomIncludeList=None):
    #     """Get the list of peak assignments from the nmrAtoms
    #     interResidueAtomPairing is the linking within the same nmrResidue
    #     interChainAtomPairing is the linking within the same chain but to different nmrResidues
    #     crossChainAtomPairing is the linking to different chains
    #     """
    #
    #     # create a set of sets ordered by spectra
    #     interResidueAtomPairing = OrderedDict((spec, set()) for spec in self.magnetisationTransfers.keys())
    #     interChainAtomPairing = OrderedDict((spec, set()) for spec in self.magnetisationTransfers.keys())
    #     crossChainAtomPairing = OrderedDict((spec, set()) for spec in self.magnetisationTransfers.keys())
    #
    #     nmrChain = nmrResidue.nmrChain
    #
    #     for nmrAtom in nmrResidue.nmrAtoms:
    #
    #         if nmrAtom._flaggedForDelete or nmrAtom.isDeleted:
    #             continue
    #
    #         if '.ASP' in str(nmrAtom.nmrResidue.id) and nmrAtom.name == 'CB':
    #             pass
    #
    #         for peak in nmrAtom.assignedPeaks:
    #
    #             # ignore peaks that are due for delete (can probably also use the notifier list)
    #             if peak._flaggedForDelete or peak.isDeleted:
    #                 continue
    #
    #             spec = peak.peakList.spectrum
    #             for assignment in peak.assignments:
    #
    #                 # find the mainNmrResidue for -1 and +1 connections
    #                 newCon = list(assignment)
    #                 for conNum in range(len(assignment)):
    #
    #                     # assignments could be None
    #                     if assignment[conNum] and assignment[conNum].nmrResidue.relativeOffset == -1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:
    #
    #                         # this is a minus residue so find connected, have to traverse to the previousNmrResidue
    #                         # will it always exist?
    #                         conName = assignment[conNum].name
    #                         preN = assignment[conNum].nmrResidue.mainNmrResidue.previousNmrResidue
    #                         if preN:
    #                             newConSwap = [nmrA for nmrA in preN.nmrAtoms if nmrA.name == conName]
    #                             if newConSwap:
    #                                 newCon[conNum] = newConSwap[0]
    #                         else:
    #                             newCon[conNum] = None  # not connected so skip
    #
    #                     elif assignment[conNum] and assignment[conNum].nmrResidue.relativeOffset == +1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:
    #
    #                         # this is a plus residue so find connected, have to traverse to the nextNmrResidue
    #                         # will it always exist?
    #                         conName = assignment[conNum].name
    #                         preN = assignment[conNum].nmrResidue.mainNmrResidue.nextNmrResidue
    #                         if preN:
    #                             newConSwap = [nmrA for nmrA in preN.nmrAtoms if nmrA.name == conName]
    #                             if newConSwap:
    #                                 newCon[conNum] = newConSwap[0]
    #                         else:
    #                             newCon[conNum] = None  # not connected so skip
    #
    #                 assignment = newCon
    #
    #                 # only get the assignments a-b if a and b are defined in the spectrum magnetisationTransfers list
    #                 for mag in self.magnetisationTransfers[spec]:
    #                     nmrAtom0 = assignment[mag[0] - 1]
    #                     nmrAtom1 = assignment[mag[1] - 1]
    #                     nmrAtom0 = nmrAtom0 if nmrAtom0 and not (nmrAtom0.isDeleted or nmrAtom0._flaggedForDelete) else None
    #                     nmrAtom1 = nmrAtom1 if nmrAtom1 and not (nmrAtom1.isDeleted or nmrAtom1._flaggedForDelete) else None
    #
    #                     if not None in (nmrAtom0, nmrAtom1):
    #
    #                         # ignore nmrAtoms that are not in the include list (if specified)
    #                         if nmrAtomIncludeList is not None and not (nmrAtom0 in nmrAtomIncludeList or nmrAtom1 in nmrAtomIncludeList):
    #                             continue
    #
    #                         if (nmrAtom0.nmrResidue is nmrResidue) and (nmrAtom1.nmrResidue is nmrResidue):
    #
    #                             # interResidueAtomPairing
    #                             if (nmrAtom1, nmrAtom0, peak) not in interResidueAtomPairing[spec]:
    #                                 interResidueAtomPairing[spec].add((nmrAtom0, nmrAtom1, peak))
    #
    #                         elif (nmrAtom0.nmrResidue.nmrChain is nmrChain) and (nmrAtom1.nmrResidue.nmrChain is nmrChain):
    #
    #                             # connections within the same chain
    #                             if (nmrAtom1, nmrAtom0, peak) not in interChainAtomPairing[spec]:
    #                                 interChainAtomPairing[spec].add((nmrAtom0, nmrAtom1, peak))
    #
    #                         # elif (nmrAtom0.nmrResidue.nmrChain is nmrChain) and (nmrAtom1.nmrResidue.nmrChain is not nmrChain):
    #                         else:
    #
    #                             # connections to a dif
    #                             if (nmrAtom1, nmrAtom0, peak) not in crossChainAtomPairing[spec]:
    #                                 crossChainAtomPairing[spec].add((nmrAtom0, nmrAtom1, peak))
    #
    #     return interResidueAtomPairing, interChainAtomPairing, crossChainAtomPairing

    # def _addPeakAssignmentLinesToGroup(self, assignments, lineList):
    #     for specAssignments in assignments.values():
    #         for nmrAtomPair in specAssignments:
    #
    #             guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
    #                               self.guiNmrAtomDict.get(nmrAtomPair[1]),
    #                               nmrAtomPair[2]
    #                               )
    #
    #             # skip if not defined
    #             if None in guiNmrAtomPair:
    #                 continue
    #
    #             # get the peak and the spectrum
    #             peak = guiNmrAtomPair[2]
    #             spectrum = peak.peakList.spectrum
    #             displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
    #
    #             # add the internal line to the guiNmrResidueGroup, should now move when group is moved
    #             try:
    #                 group = self.guiNmrResidues[guiNmrAtomPair[0].nmrAtom.nmrResidue]
    #             except Exception as es:
    #                 pass
    #
    #             self._addConnectingLineToGroup(group,
    #                                            guiNmrAtomPair[0],
    #                                            guiNmrAtomPair[1],
    #                                            spectrum.positiveContourColour,
    #                                            2.0, displacement=displacement,
    #                                            peak=peak, lineList=lineList, lineId=peak)
    #
    #             # update displacements for both guiNmrAtoms
    #             guiNmrAtomPair[0].addConnectedList(guiNmrAtomPair[1])
    #             guiNmrAtomPair[1].addConnectedList(guiNmrAtomPair[0])

    # def _addPeakAssignmentLinesToAdjacentGroup(self, nmrResidue, assignments, peaklineList, connectingLineList):
    #     for specAssignments in assignments.values():
    #         for nmrAtomPair in specAssignments:
    #
    #             guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
    #                               self.guiNmrAtomDict.get(nmrAtomPair[1]),
    #                               nmrAtomPair[2]
    #                               )
    #
    #             # get the peak and the spectrum
    #             peak = guiNmrAtomPair[2]
    #             spectrum = peak.peakList.spectrum
    #
    #             if guiNmrAtomPair[0] is None and guiNmrAtomPair[1] is None:
    #                 continue
    #
    #             if guiNmrAtomPair[0] is None:
    #                 if nmrAtomPair[1].nmrResidue.nmrChain is not nmrResidue.nmrChain:
    #                     continue
    #
    #                 newGhostResidue = self.addGhostResidue(nmrAtomPair[0].nmrResidue,
    #                                                        guiNmrAtomPair[1],
    #                                                        nmrAtomPair[1].nmrResidue,
    #                                                        nmrAtomPair[0].name,
    #                                                        nmrAtomPair[1].name,
    #                                                        True,
    #                                                        lineList=connectingLineList)
    #                 guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
    #                                   self.guiNmrAtomDict.get(nmrAtomPair[1]),
    #                                   nmrAtomPair[2]
    #                                   )
    #
    #                 group = self.guiNmrResidues[nmrAtomPair[1].nmrResidue]
    #                 displacement = guiNmrAtomPair[1].getConnectedList(guiNmrAtomPair[0])
    #                 self._addConnectingLineToGroup(group,
    #                                                guiNmrAtomPair[1],
    #                                                guiNmrAtomPair[0],
    #                                                spectrum.positiveContourColour,
    #                                                2.0, displacement=displacement,
    #                                                peak=peak, lineList=peaklineList, lineId=nmrResidue)
    #
    #             elif guiNmrAtomPair[1] is None:
    #                 if nmrAtomPair[0].nmrResidue.nmrChain is not nmrResidue.nmrChain:
    #                     continue
    #
    #                 newGhostResidue = self.addGhostResidue(nmrAtomPair[1].nmrResidue,
    #                                                        guiNmrAtomPair[0],
    #                                                        nmrAtomPair[0].nmrResidue,
    #                                                        nmrAtomPair[1].name,
    #                                                        nmrAtomPair[0].name,
    #                                                        True,
    #                                                        lineList=connectingLineList)
    #                 guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
    #                                   self.guiNmrAtomDict.get(nmrAtomPair[1]),
    #                                   nmrAtomPair[2]
    #                                   )
    #
    #                 group = self.guiNmrResidues[nmrAtomPair[0].nmrResidue]
    #                 displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
    #                 self._addConnectingLineToGroup(group,
    #                                                guiNmrAtomPair[0],
    #                                                guiNmrAtomPair[1],
    #                                                spectrum.positiveContourColour,
    #                                                2.0, displacement=displacement,
    #                                                peak=peak, lineList=peaklineList, lineId=nmrResidue)
    #
    #             else:
    #                 if nmrAtomPair[0].nmrResidue.nmrChain is nmrResidue.nmrChain:
    #                     group = self.guiNmrResidues[nmrAtomPair[0].nmrResidue]
    #                     displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
    #                     self._addConnectingLineToGroup(group,
    #                                                    guiNmrAtomPair[0],
    #                                                    guiNmrAtomPair[1],
    #                                                    spectrum.positiveContourColour,
    #                                                    2.0, displacement=displacement,
    #                                                    peak=peak, lineList=peaklineList, lineId=nmrResidue)
    #
    #                 elif nmrAtomPair[1].nmrResidue.nmrChain is nmrResidue.nmrChain:
    #                     group = self.guiNmrResidues[nmrAtomPair[1].nmrResidue]
    #                     displacement = guiNmrAtomPair[1].getConnectedList(guiNmrAtomPair[0])
    #                     self._addConnectingLineToGroup(group,
    #                                                    guiNmrAtomPair[1],
    #                                                    guiNmrAtomPair[0],
    #                                                    spectrum.positiveContourColour,
    #                                                    2.0, displacement=displacement,
    #                                                    peak=peak, lineList=peaklineList, lineId=nmrResidue)
    #                 else:
    #                     continue
    #
    #             guiNmrAtomPair[0].addConnectedList(guiNmrAtomPair[1])
    #             guiNmrAtomPair[1].addConnectedList(guiNmrAtomPair[0])

    # def addGhostResidue(self, nmrResidueCon1: NmrResidue,
    #                     guiRef: GuiNmrAtom,
    #                     nmrResidueCon0: NmrResidue,
    #                     name1: str, name0: str,
    #                     offsetAdjust,
    #                     atomSpacing=None, lineList=None):
    #     """Takes an Nmr Residue and a direction, either '-1 or '+1', and adds a residue to the sequence graph
    #     corresponding to the Nmr Residue.
    #     Nmr Residue name displayed beneath CA of residue drawn and residue type predictions displayed
    #     beneath Nmr Residue name
    #     """
    #
    #     # need to keep a list of the atoms that have been added so don't repeat
    #     count = 0
    #     if nmrResidueCon0 in self.ghostList:
    #         count = len(self.ghostList[nmrResidueCon0])
    #         if nmrResidueCon1 in self.ghostList[nmrResidueCon0]:
    #             # already exists in the dict so exit
    #             return
    #     else:
    #         self.ghostList[nmrResidueCon0] = ()
    #
    #     nmrResidue = nmrResidueCon1
    #     atoms = {}
    #     if atomSpacing:
    #         self.atomSpacing = atomSpacing
    #     nmrAtoms = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
    #     residueAtoms = DEFAULT_RESIDUE_ATOMS.copy()
    #
    #     if nmrResidue.residueType == 'GLY':
    #         del residueAtoms['CB']
    #
    #     for k, v in residueAtoms.items():
    #         if k in nmrAtoms:
    #             nmrAtom = nmrResidue.fetchNmrAtom(name=k)
    #         else:
    #             nmrAtom = None
    #         atoms[k] = self._createGhostGuiNmrAtom(k, v, nmrAtom)
    #
    #     newGuiResidueGroup = self._assembleGhostResidue(nmrResidue, atoms, lineList=lineList)
    #     newGuiResidueGroup.crossChainCount = count
    #     newGuiResidueGroup.crossChainResidue = nmrResidueCon0
    #
    #     self.ghostList[nmrResidueCon0] += (nmrResidueCon1,)
    #     return atoms

    def _getDisplays(self):
        """Return list of displays to navigate - if needed
        """
        if not self.application:
            return []

        displays = []
        # check for valid displays
        gids = self._SGwidget.displaysWidget.getTexts()
        if len(gids) == 0: return displays
        if ALL in gids:
            displays = self.application.ui.mainWindow.spectrumDisplays
        else:
            displays = [self.application.getByGid(gid) for gid in gids if gid != ALL]
        return displays

    def navigateToNmrResidue(self, selectedNmrResidue=None):
        """Navigate in selected displays to nmrResidue; skip if none defined
        """
        nmrResidue = self.current.nmrResidue
        if not nmrResidue:
            return

        logger.debug('nmrResidue=%s' % (nmrResidue.id))

        displays = self._getDisplays()

        if len(displays) == 0:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        with undoBlock():
            # optionally clear the marks
            if self._SGwidget.checkBoxes['autoClearMarks']['checkBox'].isChecked():
                self.mainWindow.clearMarks()

            # navigate the displays
            for display in displays:
                if display and len(display.strips) > 0 and display.strips[0].spectrumViews:
                    newWidths = None  #_getCurrentZoomRatio(display.strips[0].viewRange())
                    if display.strips[0].spectrumViews[0].spectrum.dimensionCount <= 2:
                        widths = _getCurrentZoomRatio(display.strips[0].viewRange())

                    navigateToNmrResidueInDisplay(nmrResidue, display, stripIndex=0,
                                                  widths=newWidths,  #['full'] * len(display.strips[0].axisCodes),
                                                  showSequentialResidues=(len(display.axisCodes) > 2) and
                                                                         self._SGwidget.checkBoxes['sequentialStrips']['checkBox'].isChecked(),
                                                  markPositions=self._SGwidget.checkBoxes['markPositions']['checkBox'].isChecked()
                                                  )

    def _raiseContextMenu(self, object, event: QtGui.QMouseEvent):
        """Creates and raises a context menu enabling items to be disconnected
        """
        cursor = QtGui.QCursor()
        contextMenu = Menu('', self, isFloatWidget=True)

        pressed = self.scene.mouseGrabberItem()

        if isinstance(object, AssignmentLine):  # self.selectedLine:
            thisLine = pressed  #self.selectedLine

            if thisLine._peak and thisLine._peak.assignedNmrAtoms:
                contextMenu.addAction('deassign nmrAtoms from Peak: %s' % str(thisLine._peak.id))
                contextMenu.addSeparator()

                # add the nmrAtoms to the menu
                for nmrAtomList in thisLine._peak.assignedNmrAtoms:
                    for nmrAtom in nmrAtomList:
                        if nmrAtom:

                            # contextMenu.addAction(nmrAtom.id, partial(self.deassignPeak, thisLine._peak, nmrAtom))

                            if nmrAtom.nmrResidue is None or nmrAtom.nmrResidue.offsetNmrResidues:
                                contextMenu.addAction(nmrAtom.id, partial(self.deassignPeak, thisLine._peak, nmrAtom))
                            else:
                                contextMenu.addAction('(' + nmrAtom.id + ')', partial(self.deassignPeak, thisLine._peak, nmrAtom))

                contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
                contextMenu.exec()
                contextMenu = None

        elif isinstance(pressed, GuiNmrResidue):

            # create the nmrResidue menu
            self._disconnectPreviousActionMenu = contextMenu.addAction(self.disconnectPreviousIcon, 'disconnect Previous nmrResidue',
                                                                       partial(self.disconnectPreviousNmrResidue))
            self._disconnectActionMenu = contextMenu.addAction(self.disconnectIcon, 'disconnect nmrResidue', partial(self.disconnectNmrResidue))
            self._disconnectNextActionMenu = contextMenu.addAction(self.disconnectNextIcon, 'disconnect Next nmrResidue',
                                                                   partial(self.disconnectNextNmrResidue))
            contextMenu.addSeparator()
            self._disconnectAllActionMenu = contextMenu.addAction('disconnect all nmrResidues', partial(self.disconnectAllNmrResidues))
            if object.nmrResidue.residue:
                contextMenu.addSeparator()
                self._deassignNmrChainActionMenu = contextMenu.addAction('deassign nmrChain', partial(self.deassignNmrChain))

                assign = pressed.nmrResidue.residue is not None
                self._deassignNmrChainActionMenu.setEnabled(assign)

            contextMenu.addSeparator()
            self._showActionMenu = contextMenu.addAction('Show nmrResidue', partial(self.showNmrResidue, object))

            prev = pressed.nmrResidue.previousNmrResidue is not None
            nxt = pressed.nmrResidue.nextNmrResidue is not None

            self._disconnectPreviousActionMenu.setEnabled(prev)
            self._disconnectActionMenu.setEnabled(prev or nxt)
            self._disconnectNextActionMenu.setEnabled(nxt)
            self._disconnectAllActionMenu.setEnabled(prev or nxt)

            contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
            contextMenu.exec()
            contextMenu = None

        elif isinstance(pressed, GuiNmrAtom):

            # create the nmrAtom menu
            contextMenu.addAction('deassign nmrAtoms from Peaks')
            contextMenu.addSeparator()

            if pressed.nmrAtom and pressed.nmrAtom.assignedPeaks:

                # add nmrAtoms to the menu
                allPeaks = OrderedSet()

                # add the internal peaks to the list
                for peak in pressed.nmrAtom.assignedPeaks:
                    allPeaks.add(peak)

                self._addPeaksToMenu(allPeaks, contextMenu)

                allAdjacentPeaks = OrderedSet()

                atom = pressed.nmrAtom
                atomName = atom.name
                res = atom.nmrResidue
                nextRes = res.nextNmrResidue.mainNmrResidue if res.nextNmrResidue else None
                prevRes = res.previousNmrResidue.mainNmrResidue if res.previousNmrResidue else None

                # add peaks corresponding to the previous residue
                if prevRes:
                    for offRes in prevRes.offsetNmrResidues:
                        if offRes.relativeOffset == 1:

                            for offAtm in offRes.nmrAtoms:

                                if offAtm.name == atomName:
                                    for peak in offAtm.assignedPeaks:
                                        allAdjacentPeaks.add(peak)

                # add peaks corresponding to the next residue
                if nextRes:
                    for offRes in nextRes.offsetNmrResidues:
                        if offRes.relativeOffset == -1:

                            for offAtm in offRes.nmrAtoms:

                                if offAtm.name == atomName:
                                    for peak in offAtm.assignedPeaks:
                                        allAdjacentPeaks.add(peak)

                if allPeaks and allAdjacentPeaks:
                    contextMenu.addSeparator()

                self._addPeaksToMenu(allAdjacentPeaks, contextMenu)

                contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
                contextMenu.exec()
                contextMenu = None

    def _addPeaksToMenu(self, allPeaks, contextMenu):
        for peak in allPeaks:

            if peak and peak.assignedNmrAtoms:
                subMenu = contextMenu.addMenu(peak.id)

                subMenu.addAction('nmrAtoms')
                subMenu.addSeparator()
                for nmrAtomList in peak.assignedNmrAtoms:
                    for nmrAtom in nmrAtomList:
                        if nmrAtom:

                            if nmrAtom.nmrResidue is None or nmrAtom.nmrResidue.offsetNmrResidues:
                                subMenu.addAction(nmrAtom.id, partial(self.deassignPeak, peak, nmrAtom))
                            else:
                                subMenu.addAction('(' + nmrAtom.id + ')', partial(self.deassignPeak, peak, nmrAtom))

    def showNmrResidue(self, object):
        self.navigateToNmrResidue(selectedNmrResidue=object.nmrResidue)


import math


atomSpacing = 66
cos36 = math.cos(math.pi / 5)
sin36 = math.sin(math.pi / 5)
tan36 = math.tan(math.pi / 5)

cos54 = math.cos(3 * math.pi / 10)
sin54 = math.sin(3 * math.pi / 10)

cos60 = math.cos(math.pi / 3)
sin60 = math.sin(math.pi / 3)
sin72 = math.sin(2 * math.pi / 5)
cos72 = math.cos(2 * math.pi / 5)

DEFAULT_RESIDUE_ATOMS = {'H' : np.array([0, 0]),
                         'N' : np.array([0, -1 * atomSpacing]),
                         'CA': np.array([atomSpacing, -1 * atomSpacing]),
                         'CB': np.array([atomSpacing, -2 * atomSpacing]),
                         'C': np.array([2 * atomSpacing, -1 * atomSpacing])
                         }

# Use loadCompoundPickle from chemBuild to load the structure for these; found in compound.variants.bonds

ATOM_POSITION_DICT = {

    'ALA': {'HB%'       : (0.0, -0.75 * atomSpacing),
            'boundAtoms': ('')},
    'CYS': {'SG': (0.0, -1 * atomSpacing), 'HG': (0, -1.75 * atomSpacing)},
    'ASP': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CG' : (0, -1 * atomSpacing)},
    'ASN': {'HBx' : (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CG'  : (0, -1 * atomSpacing), 'ND2': (0, -2 * atomSpacing),
            'HD2x': (atomSpacing * -0.75, -2 * atomSpacing - (0.75 * atomSpacing * cos60)),
            'HD2y': (atomSpacing * +0.75, -2 * atomSpacing - (0.75 * atomSpacing * cos60)),
            },
    'GLU': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'HGx': (atomSpacing * -0.75, -1 * atomSpacing), 'HGy': (atomSpacing * 0.75, -1 * atomSpacing),
            'CG' : (0, -1 * atomSpacing), 'CD': (0, -2 * atomSpacing)
            },
    'GLN': {'HBx' : (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'HGx' : (atomSpacing * -0.75, -1 * atomSpacing), 'HGy': (atomSpacing * 0.75, -1 * atomSpacing),
            'CG'  : (0, -1 * atomSpacing), 'CD': (0, -2 * atomSpacing), 'NE2': (0, -3 * atomSpacing),
            'HD2x': (atomSpacing * -0.75, -3 * atomSpacing - (0.75 * atomSpacing * cos60)),
            'HD2y': (atomSpacing * +0.75, -3 * atomSpacing - (0.75 * atomSpacing * cos60)),
            },
    'PHE': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CG' : (0, -1 * atomSpacing), 'CD1': (-1 * atomSpacing, (-1 - cos60) * atomSpacing),
            'CD2': (1 * atomSpacing, (-1 - cos60) * atomSpacing),
            'CE1': (-1 * atomSpacing, (-2 - cos60) * atomSpacing),
            'CE2': (1 * atomSpacing, (-2 - cos60) * atomSpacing),
            'HD1': (-1.75 * atomSpacing, (-1 - cos60) * atomSpacing),
            'HD2': (1.75 * atomSpacing, (-1 - cos60) * atomSpacing),
            'HE1': (-1.75 * atomSpacing, (-2 - cos60) * atomSpacing),
            'HE2': (1.75 * atomSpacing, (-2 - cos60) * atomSpacing),
            'CZ' : (0, (-2 - cos60 - sin60) * atomSpacing), 'HZ': (0, (-2 - cos60 - sin60 - 0.75) * atomSpacing)
            },
    'TYR': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CG' : (0, -1 * atomSpacing), 'CD1': (-1 * atomSpacing, (-1 - cos60) * atomSpacing),
            'CD2': (1 * atomSpacing, (-1 - cos60) * atomSpacing),
            'CE1': (-1 * atomSpacing, (-2 - cos60) * atomSpacing),
            'CE2': (1 * atomSpacing, (-2 - cos60) * atomSpacing),
            'HD1': (-1.75 * atomSpacing, (-1 - cos60) * atomSpacing),
            'HD2': (1.75 * atomSpacing, (-1 - cos60) * atomSpacing),
            'HE1': (-1.75 * atomSpacing, (-2 - cos60) * atomSpacing),
            'HE2': (1.75 * atomSpacing, (-2 - cos60) * atomSpacing),
            'CZ' : (0, (-2 - cos60 - sin60) * atomSpacing), 'HH': (0, (-2 - cos60 - sin60 - 0.75) * atomSpacing)
            },
    'SER': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'HG' : (0, -1 * atomSpacing)
            },
    'THR': {'HG1': (atomSpacing * -0.75, 0.0), 'HB': (atomSpacing * 0.75, 0.0),
            'CG2': (0, -1 * atomSpacing), 'HG2%': (0, -1.75 * atomSpacing)
            },
    'MET': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'HGx': (atomSpacing * -0.75, -1 * atomSpacing), 'HGy': (atomSpacing * 0.75, -1 * atomSpacing),
            'CG' : (0, -1 * atomSpacing), 'SD': (0, -2 * atomSpacing), 'CE': (0, -3 * atomSpacing),
            'HE%': (0, -3.75 * atomSpacing)
            },
    'ARG': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'HGx': (atomSpacing * -0.75, -1 * atomSpacing), 'HGy': (atomSpacing * 0.75, -1 * atomSpacing),
            'CG' : (0, -1 * atomSpacing), 'CD': (0, -2 * atomSpacing), 'NE': (0, -3 * atomSpacing),
            'CZ' : (0, -4 * atomSpacing), 'NH1': (atomSpacing * -1, -4 * atomSpacing - (0.75 * atomSpacing * cos60)),
            'NH2': (atomSpacing * +1, -4 * atomSpacing - (0.75 * atomSpacing * cos60)),
            },
    'VAL': {'HBx' : (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CGx' : (-1 * atomSpacing, -1 * (cos60 * atomSpacing)),
            'CGy' : (1 * atomSpacing, -1 * (cos60 * atomSpacing)),
            'HGx%': (atomSpacing * -1, -1 * (cos60 * atomSpacing) - (0.75 * atomSpacing)),
            'HGy%': (atomSpacing * +1, -1 * (cos60 * atomSpacing) - (0.75 * atomSpacing))
            },
    'LEU': {'HBx' : (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'HGx' : (atomSpacing * -0.75, -1 * atomSpacing), 'HGy': (atomSpacing * 0.75, -1 * atomSpacing),
            'CG'  : (0, -1 * atomSpacing),
            'CDx' : (-1 * atomSpacing, (-1 - cos60) * atomSpacing),
            'CDy' : (1 * atomSpacing, (-1 - cos60) * atomSpacing),
            'HDx%': (atomSpacing * -1, ((-1 - cos60) * atomSpacing) - (0.75 * atomSpacing)),
            'HDy%': (atomSpacing * +1, ((-1 - cos60) * atomSpacing) - (0.75 * atomSpacing))
            },
    'ILE': {'HBx' : (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CG1' : (-1 * atomSpacing, -1 * (cos60 * atomSpacing)),
            'CG2%': (1 * atomSpacing, -1 * (cos60 * atomSpacing)),
            'HG1x': (atomSpacing * -1.75, -1 * (cos60 * atomSpacing)),
            'HG1y': (atomSpacing * -0.25, -1 * (cos60 * atomSpacing)),
            'HG2%': (1 * atomSpacing, -1 * (cos60 * atomSpacing) - (0.75 * atomSpacing)),
            'CD1%': (-1 * atomSpacing, -1 * (cos60 * atomSpacing) - atomSpacing),
            'HD1%': (-1 * atomSpacing, -1 * (cos60 * atomSpacing) - (1.75 * atomSpacing)),
            },
    'LYS': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'HGx': (atomSpacing * -0.75, -1 * atomSpacing), 'HGy': (atomSpacing * 0.75, -1 * atomSpacing),
            'CG' : (0, -1 * atomSpacing), 'CD': (0, -2 * atomSpacing),
            'HDx': (atomSpacing * -0.75, -2 * atomSpacing), 'HDy': (atomSpacing * 0.75, -2 * atomSpacing),
            'HEx': (atomSpacing * -0.75, -3 * atomSpacing), 'HEy': (atomSpacing * 0.75, -3 * atomSpacing),
            'CE' : (0, -3 * atomSpacing),
            'NZ' : (0, -4 * atomSpacing), 'HZ%': (0, -4.75 * atomSpacing),
            },
    'HIS': {'HBx': (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CG' : (0, -1 * atomSpacing), 'ND1': (-1 * atomSpacing, -1 * (atomSpacing + (atomSpacing / (2 * tan36)))),
            'CD2': (atomSpacing, -1 * (atomSpacing + (atomSpacing / (2 * tan36)))),
            'NE2': (atomSpacing / 2, -1 * (atomSpacing + (atomSpacing / (2 * sin36)) + (atomSpacing / (2 * tan36)))),
            'CD1': (-0.5 * atomSpacing, -1 * (atomSpacing + (atomSpacing / (2 * sin36)) + (atomSpacing / (2 * tan36)))),
            },

    'TRP': {'HBx'       : (atomSpacing * -0.75, 0.0), 'HBy': (atomSpacing * 0.75, 0.0),
            'CG'        : (0, -1 * atomSpacing), 'CD1': (atomSpacing, -1 * atomSpacing),
            'NE1'       : (atomSpacing + (atomSpacing * cos72), -1 * (atomSpacing + (atomSpacing * sin72))),
            'CE2'       : (atomSpacing + (atomSpacing * cos72) - (atomSpacing * sin54),
                           -1 * (atomSpacing + (atomSpacing * sin72) + (atomSpacing * cos54))),
            'CD2'       : (-1 * (atomSpacing * cos72), -1 * (atomSpacing + (atomSpacing * sin72))),
            'CE3'       : (atomSpacing + (atomSpacing * cos72) - (atomSpacing * sin54) - (2 * (atomSpacing * sin60)),
                           -1 * (atomSpacing + (atomSpacing * sin72) + (atomSpacing * cos54))),
            'CZ2'       : (atomSpacing + (atomSpacing * cos72) - (atomSpacing * sin54),
                           -1 * (2 * atomSpacing + (atomSpacing * sin72) + (atomSpacing * cos54))),
            'CZ3'       : (atomSpacing + (atomSpacing * cos72) - (atomSpacing * sin54) - (2 * (atomSpacing * sin60)),
                           -1 * (2 * atomSpacing + (atomSpacing * sin72) + (atomSpacing * cos54))),
            'CH2'       : (-1 * (atomSpacing * cos72), -1 * (2 * atomSpacing + (atomSpacing * sin72) + (atomSpacing * cos54) + (atomSpacing * cos60))),

            'boundAtoms': (('CG', 'CD1'), ('CG', 'CD2'), ('CD2', 'CE3'), ('CD2', 'CE2'),
                           ('CD1', 'NE1'), ('CE2', 'CZ2'), ('CE3', 'CZ3'), ('CZ3', 'CH2'),
                           ('CZ2', 'CH2'), ('NE1', 'CE2'))
            },

    'PRO': {
        'CB': (atomSpacing * cos72, -1 * (atomSpacing * sin72) + atomSpacing),
        'CG': (-0.5 * atomSpacing, -1 * atomSpacing / (2 * tan36)),
        'CD': (-1 * (atomSpacing + (atomSpacing * cos72)), -1 * (atomSpacing * sin72) + atomSpacing),
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

def greekKey(word):
    greekSort = '0123456789ABGDEZHQIKLMNXOPRSTUFCYWabgdezhqiklmnxoprstufcyw'
    greekLetterCount = len(greekSort)

    key = (0,)
    if word:
        key = (ord(word[0]),)
        key += tuple(greekSort.index(c) if c in greekSort else greekLetterCount for c in word[1:])
    return key


if __name__ == '__main__':
    from ccpn.ui.gui.widgets.Application import TestApplication
    from ccpn.ui.gui.widgets.TextEditor import TextEditor
    from ccpnmodel.ccpncore.lib.assignment.ChemicalShift import PROTEIN_ATOM_NAMES, ALL_ATOMS_SORTED

    isotopes = ['13C', '1H', '15N']
    axisCodes = ['C', 'H', 'N']

    atoms = {'H', 'N', 'CA', 'CB', 'C'}
    for k, v in ALL_ATOMS_SORTED.items():
        atoms = atoms | set(v)

    thisDict = {}
    for axis, isotope in zip(axisCodes, isotopes):
        thisDict[isotope] = sorted([at for at in set(atoms) if at.startswith(axis)], key=greekKey)

    # app = TestApplication()
    #
    # popup = SequenceGraphModule()
    #
    # popup.show()
    # popup.raise_()
    # app.start()
