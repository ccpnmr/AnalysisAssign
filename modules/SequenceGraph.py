"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:21 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b4 $"
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
from itertools import islice
from PyQt5 import QtGui, QtWidgets, QtCore
from collections import OrderedDict
from contextlib import contextmanager

from ccpn.core.lib.Pid import Pid
from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.Peak import Peak
from ccpn.core.Spectrum import Spectrum
from ccpn.core.NmrChain import NmrChain
from ccpn.core.lib.AssignmentLib import getNmrResiduePrediction
from ccpn.core.lib.AssignmentLib import nmrAtomPairsByDimensionTransfer
from ccpn.core.lib.Notifiers import Notifier
from ccpn.ui.gui.lib.Strip import navigateToNmrResidueInDisplay, _getCurrentZoomRatio
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.ui.gui.guiSettings import textFontSmall, textFontSmallBold, textFont
from ccpn.ui.gui.guiSettings import getColours
from ccpn.ui.gui.guiSettings import GUINMRATOM_NOTSELECTED, GUINMRATOM_SELECTED, \
    GUINMRRESIDUE, SEQUENCEGRAPHMODULE_LINE, SEQUENCEGRAPHMODULE_TEXT

from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Icon import Icon
from ccpn.ui.gui.widgets.ToolBar import ToolBar
from ccpn.ui.gui.widgets.Spacer import Spacer
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
from ccpn.core.lib.ContextManagers import logCommandBlock

from ccpn.util.decorators import profile


logger = getLogger()
ALL = '<all>'


class GuiNmrAtom(QtWidgets.QGraphicsTextItem):
    """
    A graphical object specifying the position and name of an atom when created by the Assigner.
    Can be linked to a Nmr Atom.
    """

    def __init__(self, mainWindow, text, pos=None, nmrAtom=None):

        super(GuiNmrAtom, self).__init__()

        self.setPlainText(text)
        self.setPos(QtCore.QPointF((pos[0] - self.boundingRect().x()), (pos[1] - self.boundingRect().y())))

        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current
        self.nmrAtom = nmrAtom

        self.connectedAtoms = 0
        self.connectedList = {}  # ejb - new connection test

        # wb104: not sure why below is needed rather than setFlags() but it is
        self.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        ###self.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable | self.flags())
        #self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable)

        self.colours = getColours()
        if self.isSelected:
            self.setDefaultTextColor(QtGui.QColor(self.colours[GUINMRATOM_SELECTED]))
        else:
            self.setDefaultTextColor(QtGui.QColor(self.colours[GUINMRATOM_NOTSELECTED]))

    def mouseDoubleClickEvent(self, event):
        """
        CCPN INTERNAL - re-implementation of double click event
        """
        pass

    def mousePressEvent(self, event):
        """
        CCPN INTERNAL - re-implementation of mouse press event
        """
        if self.nmrAtom is not None:
            self.current.nmrAtom = self.nmrAtom
            self.current.nmrResidue = self.nmrAtom.nmrResidue
            event.accept()

    def _raiseContextMenu(self, event: QtGui.QMouseEvent):
        """
        Creates and raises a context menu enabling items to be disconnected
        """
        from ccpn.ui.gui.widgets.Menu import Menu
        from functools import partial

        contextMenu = Menu('', event.widget(), isFloatWidget=True)
        contextMenu.addAction('deassign all Peaks', partial(self._deassignAllPeaksFromNmrAtom))
        cursor = QtGui.QCursor()
        contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
        contextMenu.exec()

    def addConnectedList(self, connectedAtom):
        # keyVal = str(connectedAtom.nmrAtom.pid)

        keyVal = connectedAtom
        if keyVal in self.connectedList:
            self.connectedList[keyVal] += 1
        else:
            self.connectedList[keyVal] = 1

    def removeConnectedList(self, connectedAtom):
        # keyVal = str(connectedAtom.nmrAtom.pid)

        keyVal = connectedAtom
        if keyVal in self.connectedList:
            self.connectedList[keyVal] -= 1
        else:
            raise RuntimeError('Connection does not exist')

    def getConnectedList(self, connectedAtom):
        # keyVal = str(connectedAtom.nmrAtom.pid)

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


class GuiNmrResidueGroup(QtWidgets.QGraphicsItemGroup):
    """
    Group item to group all nmrAtoms/connecting lines of nmrResidue
    """

    def __init__(self, parent, nmrResidue, caAtom, pos):
        super(GuiNmrResidueGroup, self).__init__()

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


class GuiNmrResidue(QtWidgets.QGraphicsTextItem):
    """
    Object linking residues displayed in Assigner and Nmr Residues. Contains functionality for drag and
    drop assignment in conjunction with the Sequence Module.
    """

    def __init__(self, parent, nmrResidue, caAtom):

        super(GuiNmrResidue, self).__init__()
        self.setPlainText(nmrResidue.id)

        self.mainWindow = parent.mainWindow
        self.application = self.mainWindow.application
        self.project = self.mainWindow.project
        self.current = self.mainWindow.application.current

        self.setFont(textFontSmall)
        self.colours = getColours()
        self.setDefaultTextColor(QtGui.QColor(self.colours[GUINMRRESIDUE]))

        self.setPos(caAtom.x() - caAtom.boundingRect().width() / 2, caAtom.y() + 30)
        ###self.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable | self.flags())
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable)
        self._parent = parent
        self.nmrResidue = nmrResidue
        # self.mousePressEvent = self._mousePressEvent
        # self.mouseMoveEvent = self._mouseMoveEvent
        # self.mouseReleaseEvent = self._mouseReleaseEvent
        # self.mouseDoubleClickEvent = self._mouseDoubleClickEvent
        self.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)

    def _update(self):
        self.setPlainText(self.nmrResidue.id)

    def _mouseMoveEvent(self, event):
        """create a drag item if left button pressed
        """
        nmrItem = None
        if (event.buttons() == QtCore.Qt.LeftButton):  # and (event.modifiers() & QtCore.Qt.ShiftModifier):

            nmrItem = self
            nmrChain = self.nmrResidue.nmrChain

            if nmrItem:
                drag = QtGui.QDrag(event.widget())
                mimeData = QtCore.QMimeData()
                itemData = json.dumps({'pids': [nmrItem.nmrResidue.pid]})  # nmrChain.pid

                # ejb - added so that itemData works with PyQt5
                tempData = QtCore.QByteArray()
                stream = QtCore.QDataStream(tempData, QtCore.QIODevice.WriteOnly)
                stream.writeQString(itemData)
                mimeData.setData(ccpnmrJsonData, tempData)

                # mimeData.setData(ccpnmrJsonData, itemData)
                mimeData.setText(itemData)
                drag.setMimeData(mimeData)

                dragLabel = QtWidgets.QLabel()
                dragLabel.setText(self.toPlainText())
                dragLabel.setFont(textFont)
                dragLabel.setStyleSheet('color : %s' % (self.colours[GUINMRRESIDUE]))

                # pixmap = QtGui.QPixmap.grabWidget(dragLabel)    # ejb -    this gets the whole window   event.widget())
                pixmap = dragLabel.grab()  # ejb -    this gets the whole window   event.widget())

                painter = QtGui.QPainter(pixmap)
                painter.setCompositionMode(painter.CompositionMode_DestinationIn)
                painter.fillRect(pixmap.rect(), QtGui.QColor(0, 0, 0, 240))
                painter.end()
                drag.setPixmap(pixmap)
                drag.setHotSpot(QtCore.QPoint(dragLabel.width() // 2, dragLabel.height() // 2))

                drag.exec_(QtCore.Qt.MoveAction)  # ejb - same as BackboneAssignment

    def _mousePressEvent(self, event):
        self.current.nmrResidue = self.nmrResidue
        self.setSelected(True)

    def _mouseDoubleClickEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self._parent.showNmrResidue(self)
            event.accept()
        else:
            super(GuiNmrResidue, self).mouseDoubleClickEvent(event)

    def _raiseLineMenu(self, scene, event):
        """
        Creates and raises a context menu enabling items to be disconnected
        """
        from ccpn.ui.gui.widgets.Menu import Menu
        from functools import partial

        cursor = QtGui.QCursor()
        contextMenu = Menu('', event.widget(), isFloatWidget=True)
        pressed = self.scene.mouseGrabberItem()

        if self.selectedLine:
            thisLine = self.selectedLine
            contextMenu.addAction('deassign nmrAtoms from Peak: %s' % str(thisLine._peak.id))
            contextMenu.addSeparator()
            if thisLine._peak:

                # skip if nothing connected
                if not thisLine._peak.assignedNmrAtoms:
                    return

                # add the nmrAtoms to the menu
                for nmrAtomList in thisLine._peak.assignedNmrAtoms:
                    for nmrAtom in nmrAtomList:
                        if nmrAtom:
                            contextMenu.addAction(nmrAtom.id, partial(self._deassignPeak, thisLine._peak, nmrAtom))


class AssignmentLine(QtWidgets.QGraphicsLineItem):
    """
    Object to create lines between GuiNmrAtoms with specific style, width, colour and displacement.
    Displacement allows multiplet peak lines to be shown as adjacent lines between the same nmrAtom.
    """

    def __init__(self, x1, y1, x2, y2, colour, width,
                 parent=None, style=None, peak=None, atom1=None, atom2=None, displacement=None):
        QtWidgets.QGraphicsLineItem.__init__(self)

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
        self.atom1 = atom1
        self.atom2 = atom2
        self.displacement = displacement

        # enable hovering so the current line can be set
        self.setAcceptedMouseButtons(QtCore.Qt.RightButton)
        self.setAcceptHoverEvents(True)

    def updateEndPoints(self):
        """Update the endPoints of the line to point. Co-ordinates are relative to the group to
        which the graphicsItem belongs, in this case the guiNmrResidue group. GuiNmrResidue group is the top level relative to the scene.
        """
        atom1 = self.atom1
        atom2 = self.atom2
        residue1 = atom1.guiNmrResidueGroup
        residue2 = atom2.guiNmrResidueGroup

        rx1 = residue1.x()
        ry1 = residue1.y()
        rx2 = residue2.x()
        ry2 = residue2.y()

        atom1Rect = atom1.boundingRect()
        atom2Rect = atom2.boundingRect()
        w1 = atom1Rect.width()
        h1 = atom1Rect.height()
        w2 = atom2Rect.width()
        h2 = atom2Rect.height()

        x1 = atom1.x()  # + rx1
        y1 = atom1.y()  # + ry1
        x2 = atom2.x() + (rx2 - rx1)
        y2 = atom2.y() + (ry2 - ry1)

        dx = x2 - x1
        dy = y2 - y1
        length = 2.0 * pow(dx * dx + dy * dy, 0.5)
        if self.displacement is not None:
            count = (atom1.connectedList[atom2] - 1) // 2
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

    def paint(self, painter: QtGui.QPainter, option: 'QStyleOptionGraphicsItem', widget: typing.Optional[QtWidgets.QWidget] = ...):
        """Automatically update the end-points of the assignment lines to point to the correct guiNmrAtoms
        """
        self.updateEndPoints()
        super().paint(painter, option, widget)

    def hoverEnterEvent(self, event):
        self._parent.selectedLine = self

    def hoverLeaveEvent(self, event):
        self._parent.selectedLine = None

    def mousePressEvent(self, event):
        pass


class SequenceGraphModule(CcpnModule):
    """
    A module for the display of stretches of sequentially linked and assigned stretches of
    NmrResidues.
    """
    className = 'SequenceGraph'

    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'left'

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

        self._registerNotifiers()

        self.atomSpacing = 66
        self.guiResiduesShown = []
        self.predictedStretch = []
        self.direction = None
        self.selectedStretch = []
        self.guiNmrResidues = OrderedDict()
        self.guiNmrResidueLabels = []
        self.guiNmrAtomDict = {}
        self.guiGhostNmrResidues = OrderedDict()
        self.ghostList = {}
        self.selectedLine = None

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

        self.resetScene()

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
        self._SGwidget = SequenceGraphSettings(parent=self.settingsWidget, mainWindow=self.mainWindow,
                                               settingsDict=settingsDict,
                                               grid=(0, 0))

        self.residueCount = 0

        colwidth = 140
        self._MWwidget = Widget(self.mainWidget, setLayout=True,
                                grid=(0, 0), vAlign='top', hAlign='left')

        self.nmrChainPulldown = NmrChainPulldown(self._MWwidget, self.project, grid=(0, 0), gridSpan=(1, 1),
                                                 showSelectName=True,
                                                 fixedWidths=(colwidth, colwidth, colwidth),
                                                 callback=self.showNmrChainFromPulldown)

        self.refreshCheckBox = CheckBoxCompoundWidget(self._MWwidget,
                                                      labelText='Auto refresh NmrChain:',
                                                      checked=True,
                                                      fixedWidths=(colwidth, 15),
                                                      orientation='right', hAlign='left',
                                                      tipText='Update display when current.nmrChain changes',
                                                      grid=(0, 1), gridSpan=(1, 1))

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

        self._updateMagnetisationTransfers()

        # stop the mainWidget from squishing during a resize
        self.mainWidget.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)

        # install the event filter to handle maximising from floated dock
        # self.installMaximiseEventHandler(self._maximise, self._closeModule)

        if nmrChain is not None:
            self.selectSequence(nmrChain)

    def _sceneMouseRelease(self, event):
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

        except Exception as es:
            raise es
        finally:
            self._unblockEvents()

            # resize to the new items and spawns a repaint
            self.scene.setSceneRect(self.scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))

    def _updateSpectra(self, data=None):
        """Update list of current spectra and generate new magnetisationTransfer list
        """
        print('>>>_updateSpectra')

        if data:
            trigger = data[Notifier.TRIGGER]
            object = data[Notifier.OBJECT]
            print('>>> spectrum:', trigger, object)

            self._updateMagnetisationTransfers()

            if trigger in [Notifier.CREATE, Notifier.DELETE]:
                nmrChainPid = self.nmrChainPulldown.getText()
                if nmrChainPid:
                    with self.sceneBlocking():
                        self._rebuildPeakAssignments(nmrChainPid)

    def selectSequence(self, nmrChain=None):
        """Manually select a Sequence from the pullDown
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

    def _registerNotifiers(self):
        """Register the required notifiers
        """
        self._peakNotifier = self.setNotifier(self.project,
                                              [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
                                              Peak.className,
                                              self._updatePeaks,
                                              onceOnly=True)

        self._nmrChainNotifier = self.setNotifier(self.project,
                                                  [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
                                                  NmrChain.className,
                                                  self._updateNmrChains,
                                                  onceOnly=True)

        self._nmrResidueNotifier = self.setNotifier(self.project,
                                                    [Notifier.CREATE, Notifier.DELETE],
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

    def _repopulateModule(self):
        """CCPN Internal: Repopulate the required widgets in the module
        This is will be attached to GuiNotifiers
        """
        self.showNmrChainFromPulldown()

    def _updateModule(self, nmrChains=None):
        """Update in response to change of current.nmrChains
        """
        #if nmrChains is None or len(nmrChains)==0: return
        nmrChain = self.current.nmrChain
        if not nmrChain:
            return

        if not self.refreshCheckBox.isChecked():
            return

        # select the chain from the pullDown - should automatically change the display
        self.nmrChainPulldown.select(nmrChain.pid)

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

    def _addAdjacentResiduesToSet(self, nmrResidue, residueSet):
        residueSet.add(nmrResidue)
        nmr = nmrResidue.previousNmrResidue
        if nmr:
            nmr = nmr.mainNmrResidue
            residueSet.add(nmr)
        nmr = nmrResidue.nextNmrResidue
        if nmr:
            nmr = nmr.mainNmrResidue
            residueSet.add(nmr)

    def _rebuildNmrResidues(self, nmrResidues):
        """Rebuild the peaks of a specified nmrResidue.
        """
        # now rebuild for the new peak values
        # assumes that the peakAssignments have changed - possibly use different notifier
        nmrResidues = makeIterableList(nmrResidues)

        nmrAtomIncludeList = tuple(nmrAtom for nmrResidue in nmrResidues for nmrAtom in nmrResidue.nmrAtoms)
        guiNmrAtomSet = set([self.guiNmrAtomDict[nmrAtom] for nmrAtom in nmrAtomIncludeList])

        for guiAtom in guiNmrAtomSet:
            for peakLineList in self.assignmentLines.values():
                peakLines = [peakLine for peakLine in peakLineList
                             if peakLine.atom1 is guiAtom or peakLine.atom2 is guiAtom]

                # remove all graphic lines
                for peakLine in peakLines:
                    peakLineList.remove(peakLine)
                    self.scene.removeItem(peakLine)

            # clear connectivity list of guiNmrAtoms
            guiAtom.clearConnectedList()

        if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
            for nmrResidue in nmrResidues:

                # only process residues in the current visible chain
                if nmrResidue is nmrResidue.mainNmrResidue and nmrResidue.nmrChain is self.nmrChain:
                    # add the internally connected Lines
                    internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue,
                                                                                                                           nmrAtomIncludeList=nmrAtomIncludeList)
                    self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                self.assignmentLines, self.connectingLines)

        self._updateGuiResiduePositions(updateMainChain=True, updateConnectedChains=True)

    def _rebuildPeakLines(self, peaks, rebuildPeakLines=False, makeListFromPeak=False):
        """Clear all lines on the display associated with peak.
        """
        # find the current lines associated with the notified peak (or list of peaks)
        peaks = makeIterableList(peaks)

        peakLines = [peakLine for peakLineList in self.assignmentLines.values() for peakLine in peakLineList if peakLine._peak in peaks]
        guiNmrAtomSet = set()
        nmrResidueSet = set()

        if not makeListFromPeak:
            # make a list of all the guiNmrAtoms/nmrResidues that are associated with the peak
            for peakLine in peakLines:
                # current associated guiNmrAtoms
                guiNmrAtomSet.add(peakLine.atom1)
                guiNmrAtomSet.add(peakLine.atom2)

                # current associated nmrResidues and their previous/nextNmrResidues
                self._addAdjacentResiduesToSet(peakLine.atom1.nmrAtom.nmrResidue, nmrResidueSet)
                self._addAdjacentResiduesToSet(peakLine.atom2.nmrAtom.nmrResidue, nmrResidueSet)

        else:
            # make a new list for creating a peak; necessary for undo of delete peak as the assignedNmrAtom list exists
            assignmentAtoms = set([nmrAtom for peak in peaks
                                   for assignment in peak.assignments
                                   for nmrAtom in assignment
                                   if nmrAtom in self.guiNmrAtomDict])
            for nmrAtom in assignmentAtoms:
                # for peak in peaks:
                #     for assignment in peak.assignments:
                #         for nmrAtom in assignment:
                #
                #             if nmrAtom in self.guiNmrAtomDict:
                guiNmrAtomSet.add(self.guiNmrAtomDict[nmrAtom])
                self._addAdjacentResiduesToSet(nmrAtom.nmrResidue, nmrResidueSet)

        for guiAtom in guiNmrAtomSet:
            for peakLineList in self.assignmentLines.values():
                peakLines = [peakLine for peakLine in peakLineList
                             if peakLine.atom1 is guiAtom or peakLine.atom2 is guiAtom]

                # remove all graphic lines
                for peakLine in peakLines:
                    peakLineList.remove(peakLine)
                    self.scene.removeItem(peakLine)

            # clear connectivity list of guiNmrAtoms
            guiAtom.clearConnectedList()

        if rebuildPeakLines:
            # now rebuild for the new peak values
            # assumes that the peakAssignments have changed - possibly use different notifier
            if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
                for nmrResidue in nmrResidueSet:

                    # only process residues in the current visible chain
                    if nmrResidue is nmrResidue.mainNmrResidue and nmrResidue.nmrChain is self.nmrChain:
                        # add the internally connected Lines
                        internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue,
                                                                                                                               nmrAtomIncludeList=tuple(
                                                                                                                                       guiAtom.nmrAtom for
                                                                                                                                       guiAtom in
                                                                                                                                       guiNmrAtomSet))
                        self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                        self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                        self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                    self.assignmentLines, self.connectingLines)

            self._updateGuiResiduePositions(updateMainChain=True, updateConnectedChains=True)

    def _updatePeaks(self, data):
        """Update the peaks in the display.
        """
        peak = data[Notifier.OBJECT]

        print('>>>_updatePeaks', peak)
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if trigger == Notifier.DELETE:
                self._rebuildPeakLines(peak, rebuildPeakLines=True)

            elif trigger == Notifier.CREATE:
                self._rebuildPeakLines(peak, rebuildPeakLines=True, makeListFromPeak=True)

            elif trigger == Notifier.CHANGE:
                self._rebuildPeakLines(peak, rebuildPeakLines=True)

    def _updateNmrChains(self, data):
        """Update the nmrChains in the display.
        """
        nmrChain = data[Notifier.OBJECT]

        print('>>>_updateNmrChains', nmrChain)
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if trigger == Notifier.DELETE:
                print('>>>delete nmrChain - no action', nmrChain)

            elif trigger == Notifier.CREATE:
                print('>>>create nmrChain - no action', nmrChain)

            elif trigger == Notifier.CHANGE:
                print('>>>change nmrChain - no action', nmrChain)

    def _updateNmrResidues(self, data):
        """Update the nmrResidues in the display.
        """
        nmrResidue = data[Notifier.OBJECT]

        print('>>>_updateNmrResidues', nmrResidue)
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if trigger == Notifier.DELETE:
                print('>>>delete nmrResidue', nmrResidue)
                self._deleteNmrResidues(nmrResidue)

            elif trigger == Notifier.CREATE:
                print('>>>create nmrResidue', nmrResidue)
                self._createNmrResidues(nmrResidue)

            elif trigger == Notifier.CHANGE:
                print('>>>change nmrResidue - no action', nmrResidue)

            elif trigger == Notifier.OBSERVE:
                print('>>>observe nmrResidue - no action', nmrResidue)

    def _changeNmrResidues(self, data):
        """Update the nmrResidues in the display.
        """
        nmrResidue = data[Notifier.OBJECT]

        print('>>>_changeNmrResidues', nmrResidue)
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if nmrResidue in self.nmrChain.nmrResidues:
                print('>>>change nmrResidue - no action', nmrResidue)

    def _updateNmrAtoms(self, data):
        """Update the nmrAtoms in the display.
        """
        nmrAtom = data[Notifier.OBJECT]

        print('>>>_updateNmrAtoms', nmrAtom)
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if trigger == Notifier.DELETE:
                print('>>>delete nmrAtom', nmrAtom)
                self._deleteNmrAtoms(nmrAtom)

            elif trigger == Notifier.CREATE:
                print('>>>create nmrAtom', nmrAtom)
                self._createNmrAtoms(nmrAtom)

            elif trigger == Notifier.CHANGE:
                print('>>>change nmrAtom - no action', nmrAtom)

    def _resetNmrResiduePidForAssigner(self, data):  #nmrResidue, oldPid:str):
        """Reset pid for NmrResidue and all offset NmrResidues
        """
        print('>>>_resetNmrResiduePidForAssigner - no action')

        return

        nmrResidue = data['object']
        trigger = data[Notifier.TRIGGER]

        with self.sceneBlocking():
            if trigger == Notifier.RENAME:
                nmrChainPid = self.nmrChainPulldown.getText()
                if self.project.getByPid(nmrChainPid):

                    for nr in [nmrResidue] + list(nmrResidue.offsetNmrResidues):
                        for guiNmrResidueGroup in self.guiNmrResidues.values():
                            if guiNmrResidueGroup.nmrResidue is nr:
                                guiNmrResidueGroup.nmrResidueLabel._update()

    def _createNmrResidues(self, nmrResidues):
        """Create new nmrResidues in the scene
        """
        nmrResidues = makeIterableList(nmrResidues)
        nmrResidues = [nmrResidue for nmrResidue in nmrResidues if nmrResidue.nmrChain is self.nmrChain]
        self._buildNmrResidues(nmrResidues)

    def _deleteNmrResidues(self, nmrResidues):
        """Delete the nmrResidue from the scene
        """
        nmrResidues = makeIterableList(nmrResidues)

        for nmrResidue in nmrResidues:

            # ignore if not in the visible chain
            if nmrResidue in self.predictedStretch:

                ii = self.predictedStretch.index(nmrResidue)
                # nmrResiduesToUpdate = self.predictedStretch[max(0, ii - 1):min(ii + 1, len(self.predictedStretch))]

                # remove from the visible list
                self.predictedStretch.pop(ii)
                self.guiResiduesShown.pop(ii)

                if nmrResidue.previousNmrResidue:
                    previousNmrResidue = nmrResidue.previousNmrResidue.mainNmrResidue
                    self._rebuildNmrResidues(previousNmrResidue)

                self.scene.removeItem(self.guiNmrResidues[nmrResidue])

    def _createNmrAtoms(self, nmrAtoms):
        """Create new peak lines associated with the created/undeleted nmrAtoms.
        """
        nmrAtoms = makeIterableList(nmrAtoms)
        peaks = [peak for nmrAtom in nmrAtoms for peak in nmrAtom.assignedPeaks]

        for nmrAtom in nmrAtoms:
            if nmrAtom not in self.guiNmrAtomDict:
                self._addNmrAtomToGuiResidues(nmrAtom)
        self._rebuildPeakLines(peaks, rebuildPeakLines=True, makeListFromPeak=True)

    def _deleteNmrAtoms(self, nmrAtoms):
        """Delete peak lines associated with the deleted nmrAtoms.
        """
        nmrAtoms = makeIterableList(nmrAtoms)
        peakLines = self._searchPeakLines(nmrAtoms, includeDeleted=True)
        peaks = [peakLine._peak for peakLine in peakLines]

        self._rebuildPeakLines(peaks, rebuildPeakLines=True)

    def _removeLinesFromScene(self, lineDist):
        """Remove all the peakLines from the scene.
        """
        for lineList in lineDist.values():
            for line in lineList:
                self.scene.removeItem(line)
        lineDist.clear()

    def _searchPeakLines(self, nmrAtoms, includeDeleted=False):
        """Return a list of the peakLines containing one of the nmrAtoms in the list.
        """
        peakLines = []
        for lineList in self.assignmentLines.values():
            for line in lineList:
                nmrAtom = line.atom1.nmrAtom if line.atom1 else None
                if nmrAtom in nmrAtoms and (includeDeleted or not (nmrAtom.isDeleted or nmrAtom._flaggedForDelete)):
                    peakLines.append(line)
                nmrAtom = line.atom2.nmrAtom if line.atom2 else None
                if nmrAtom in nmrAtoms and (includeDeleted or not (nmrAtom.isDeleted or nmrAtom._flaggedForDelete)):
                    peakLines.append(line)

        return peakLines

    def _updateEndPoints(self, lineDict):
        """Update the end points from the dict.
        """
        for lineList in lineDict.values():
            for line in lineList:
                try:
                    line.updateEndPoints()
                except Exception as es:
                    pass

    def _updateGuiResiduePositions(self, updateMainChain=True, updateConnectedChains=True):
        """Update the positions of the residues and connected residues in other chains if required.
        """
        if updateMainChain:
            # update the group positions
            for ii, res in enumerate(self.guiNmrResidues.items()):
                # second element of tuple; res is (k, v)
                res[1].setPos(QtCore.QPointF(ii * self.atomSpacing * 3.0, 0.0))

        if updateConnectedChains:
            # update crossChainResidue positions
            for res in self.guiGhostNmrResidues.values():
                if res.crossChainResidue:
                    link = self.guiNmrResidues[res.crossChainResidue]
                    count = res.crossChainCount

                    newPosx = link.x()
                    newPosy = link.y()
                    res.setPos(QtCore.QPointF(newPosx + (count * 0.5 - 1.0) * self.atomSpacing,
                                              newPosy + (count * 2.5 + 5.0) * self.atomSpacing))

        # update the endpoints
        self._updateEndPoints(self.connectingLines)
        self._updateEndPoints(self.assignmentLines)

    def _rebuildPeakAssignments(self, nmrChainOrPid):
        """Rebuild all the peak assignments in the display after changing the number of spectra.
        """
        if isinstance(nmrChainOrPid, str):
            if not Pid.isValid(nmrChainOrPid):
                return
            nmrChain = self.project.getByPid(nmrChainOrPid)
        else:
            nmrChain = nmrChainOrPid

        # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
        if not nmrChain:
            return

        # remove all previous assignment lines and reset the dict
        self._removeLinesFromScene(self.assignmentLines)

        # clear the displacement values but keep the dict connections
        for guiAtom in self.guiNmrAtomDict.values():
            guiAtom.clearConnectedList()

        # add the peakAssignments
        if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
            for nmrResidue in nmrChain.nmrResidues:
                if nmrResidue is nmrResidue.mainNmrResidue:
                    # add the internally connected Lines
                    internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue)
                    self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                self.assignmentLines, self.connectingLines)

        # update the endpoints
        self._updateEndPoints(self.assignmentLines)

    def _buildNmrResidues(self, nmrResidueList):
        """Build the new residues in the list, inserting into the predicted stretch at the correct index.
        """
        connectingLinesNeeded = set()
        if self.nmrResiduesCheckBox.isChecked():
            for nmrResidue in nmrResidueList:
                if nmrResidue is nmrResidue.mainNmrResidue:

                    ii = self.nmrChain.nmrResidues.index(nmrResidue)
                    self.addResidue(nmrResidue, ii, lineList=self.connectingLines)
                    # self.setNotifier(nmrResidue, [Notifier.OBSERVE], 'nmrChain',
                    #                  callback=self._changeNmrResidue)

                    # add a connecting line to the adjacent residue
                    if nmrResidue.nextNmrResidue:
                        connectingLinesNeeded.add(len(self.guiResiduesShown) - 1)
        else:
            # find where to add in the predicted stretch
            raise RuntimeError('error - not implemented yet')

        if len(self.predictedStretch) > 2:
            self.predictSequencePosition(self.predictedStretch)

        # add the connecting lines
        guiNmrResidues = [self.guiNmrResidues[nmrResidue] for nmrResidue in nmrResidueList]

        for ii, res in enumerate(guiNmrResidues):
            if not self.nmrResiduesCheckBox.isChecked() or ii in connectingLinesNeeded:
                self._addConnectingLineToGroup(tuple(self.guiNmrResidues.values())[ii],
                                               res['CO'], self.guiResiduesShown[ii + 1]['N'],
                                               self._lineColour, 1.0, lineList=self.connectingLines, lineId=res)

        # add the peakAssignments
        if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
            for nmrResidue in nmrResidueList:
                if nmrResidue is nmrResidue.mainNmrResidue:
                    # add the internally connected Lines
                    internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue)
                    self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                    self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                self.assignmentLines, self.connectingLines)

        self._updateGuiResiduePositions(updateMainChain=True, updateConnectedChains=True)

    def removeNmrChainNotifiers(self):
        """Remove notifiers that are set on nmrChains.
        """
        nmrChains = tuple(self.project.nmrChains)
        foundNotifiers = self.searchNotifiers(objects=nmrChains, triggers=[Notifier.OBSERVE], targetName='nmrResidues')
        for notifier in foundNotifiers:
            print('>>>deleting notifier', notifier)
            self.deleteNotifier(notifier)

    def addNmrChainNotifiers(self):
        """Add new notifiers for all nmrChains in the project.
        """
        for nmrChain in self.project.nmrChains:
            self.setNotifier(nmrChain, triggers=[Notifier.OBSERVE], targetName='nmrResidues',
                             callback=self._changeNmrResidues)

    def setNmrChainDisplay(self, nmrChainOrPid):

        if isinstance(nmrChainOrPid, str):
            if not Pid.isValid(nmrChainOrPid):
                self.scene.clear()
                self.scene.setSceneRect(self.scene.itemsBoundingRect())
                return

            nmrChain = self.project.getByPid(nmrChainOrPid)
        else:
            nmrChain = nmrChainOrPid

        # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
        if not nmrChain:
            self.scene.clear()
            self.scene.setSceneRect(self.scene.itemsBoundingRect())
            return

        with logCommandBlock(get='self') as log:
            log('setNmrChainDisplay', nmrChainOrPid=repr(nmrChain.pid))

            self.clearAllItems()
            # self.removeNmrChainNotifiers()
            # self.addNmrChainNotifiers()

            ###nmrChain = self.project.getByPid(nmrChainPid)
            ###if self.modePulldown.currentText() == 'fragment':
            if True:

                connectingLinesNeeded = set()
                if self.nmrResiduesCheckBox.isChecked():
                    for nmrResidue in nmrChain.nmrResidues:
                        if nmrResidue is nmrResidue.mainNmrResidue:
                            self.addResidue(nmrResidue, len(self.predictedStretch), lineList=self.connectingLines)
                            # self.setNotifier(nmrResidue, [Notifier.OBSERVE], 'nmrChain',
                            #                  callback=self._changeNmrResidues)

                            # add a connecting line to the adjacent residue
                            if nmrResidue.nextNmrResidue:
                                connectingLinesNeeded.add(len(self.guiResiduesShown) - 1)

                else:
                    nmrResidue = self.current.nmrResidue
                    if nmrResidue in nmrChain.nmrResidues:
                        while nmrResidue.previousNmrResidue:  # go to start of connected stretch
                            nmrResidue = nmrResidue.previousNmrResidue
                    elif nmrChain.isConnected or nmrChain.chain:  # either NC:# or NC:A type nmrChains but not NC:@
                        nmrResidue = nmrChain.mainNmrResidues[0]
                    else:
                        # uses current nmrResidue
                        nmrResidue = nmrChain.nmrResidues[0]

                    while nmrResidue:  # add all of connected stretch
                        self.addResidue(nmrResidue, len(self.predictedStretch), lineList=self.connectingLines)
                        # self.setNotifier(nmrResidue, [Notifier.OBSERVE], 'nmrChain',
                        #                  callback=self._changeNmrResidues)
                        nmrResidue = nmrResidue.nextNmrResidue

                if len(self.predictedStretch) > 2:
                    self.predictSequencePosition(self.predictedStretch)

                # add the connecting lines
                for ii, res in enumerate(self.guiResiduesShown[:-1]):
                    if not self.nmrResiduesCheckBox.isChecked() or ii in connectingLinesNeeded:
                        self._addConnectingLineToGroup(tuple(self.guiNmrResidues.values())[ii],
                                                       res['CO'], self.guiResiduesShown[ii + 1]['N'],
                                                       self._lineColour, 1.0, lineList=self.connectingLines, lineId=res)

                # add the peakAssignments
                if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
                    for nmrResidue in nmrChain.nmrResidues:
                        if nmrResidue is nmrResidue.mainNmrResidue:
                            # add the internally connected Lines
                            internalAssignments, interChainAssignments, crossChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue)
                            self._addPeakAssignmentLinesToGroup(internalAssignments, self.assignmentLines)
                            self._addPeakAssignmentLinesToGroup(interChainAssignments, self.assignmentLines)
                            self._addPeakAssignmentLinesToAdjacentGroup(nmrResidue, crossChainAssignments,
                                                                        self.assignmentLines, self.connectingLines)

                self._updateGuiResiduePositions(updateMainChain=True, updateConnectedChains=True)

        self.nmrChain = nmrChain

    def showNmrChainFromPulldown(self, data=None):
        """Clear and redraw the nmrChain selected from the pulldown.
        """
        print('>>>showNmrChainFromPulldown')

        nmrChainPid = self.nmrChainPulldown.getText()
        if nmrChainPid:
            with self.sceneBlocking():
                self.setNmrChainDisplay(nmrChainPid)

        else:
            # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
            self.scene.clear()
            self.scene.setSceneRect(self.scene.itemsBoundingRect())

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
        if selectedPeak:
            with logCommandBlock(get='self') as log:
                log('deassignPeak', selectedPeak=repr(selectedPeak.pid), selectedNmrAtom=repr(selectedNmrAtom.pid))

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

    def deassignNmrAtom(self, selectedNmrAtom=None):
        """Remove the selected peaks from the assignedPeaks list
        """
        if selectedNmrAtom:
            with logCommandBlock(get='self') as log:
                log('deassignPeak', selectedNmrAtom=repr(selectedNmrAtom.pid))

                try:
                    atoms = list(selectedNmrAtom.assignedPeaks)
                    # selectedPeak.assignedNmrAtoms = ()

                    # for peak in selectedNmrAtom.assignedPeaks:
                    #
                    #   allAtoms = list(peak.dimensionNmrAtoms)
                    #   for dim in range(len(peak.dimensionNmrAtoms)):
                    #     dimNmrAtoms = list(peak.dimensionNmrAtoms[dim])
                    #     if selectedNmrAtom in dimNmrAtoms:
                    #       dimNmrAtoms.remove(selectedNmrAtom)
                    #
                    #       # allAtoms = list(peak.dimensionNmrAtoms)
                    #       allAtoms[dim] = dimNmrAtoms
                    #
                    #   peak.dimensionNmrAtoms = allAtoms

                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))
                    if self.application._isInDebugMode:
                        raise es

    def clearAllItems(self):
        """Removes all displayed residues in the sequence graph and resets items count to zero.
        """
        self.residueCount = 0
        self.predictedStretch = []
        self.guiResiduesShown = []
        self.guiNmrResidues = OrderedDict()
        self.guiNmrResidueLabels = []
        self.guiNmrAtomDict = {}
        self.guiGhostNmrResidues = OrderedDict()
        self.ghostList = {}
        self.connectingLines = {}
        self.assignmentLines = {}
        self.scene.clear()

    def resetScene(self):
        """
        Replace the scene with a new one to reset the size of the scrollbars.
        """
        # ejb - only needed to be done the first time, scene is resized at the end of setNmrChainDisplay
        self.scene = QtWidgets.QGraphicsScene(self)
        self.scrollContents = QtWidgets.QGraphicsView(self.scene, self)
        self.scrollContents.setRenderHints(QtGui.QPainter.Antialiasing)
        self.scrollContents.setInteractive(True)
        self.scrollContents.setGeometry(QtCore.QRect(0, 0, 300, 400))
        self.scrollContents.setAlignment(QtCore.Qt.AlignCenter)
        self._sequenceGraphScrollArea.setWidget(self.scrollContents)

    def _assembleGroupResidue(self, nmrResidue: NmrResidue, atoms: typing.Dict[str, GuiNmrAtom], pos=None, lineList=None):
        """Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
        creates a graphical representation of a residue in the assigner
        """
        guiResidueGroup = GuiNmrResidueGroup(self, nmrResidue, atoms['CA'], 0)  #atoms['H'].x())
        self.guiNmrResidues[nmrResidue] = guiResidueGroup
        self.scene.addItem(guiResidueGroup)

        # add the atoms to the group and set the reverse link
        for item in atoms.values():
            guiResidueGroup.addToGroup(item)
            item.guiNmrResidueGroup = guiResidueGroup

        # modify the group
        nmrAtoms = [atom.name for atom in nmrResidue.nmrAtoms]
        if "CB" in list(atoms.keys()):
            self._addConnectingLineToGroup(guiResidueGroup, atoms['CA'], atoms['CB'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
            self._addConnectingLineToGroup(guiResidueGroup, atoms['H'], atoms['N'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        self._addConnectingLineToGroup(guiResidueGroup, atoms['N'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
        self._addConnectingLineToGroup(guiResidueGroup, atoms['CO'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        self._addGroupResiduePredictions(guiResidueGroup, nmrResidue, atoms['CA'])

        return guiResidueGroup

    def _assembleGhostResidue(self, nmrResidue: NmrResidue, atoms: typing.Dict[str, GuiNmrAtom], lineList=None):
        """Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
        creates a graphical representation of a residue in the assigner
        """
        guiResidueGroup = GuiNmrResidueGroup(self, nmrResidue, atoms['CA'], 0)  #atoms['H'].x())
        self.guiGhostNmrResidues[nmrResidue] = guiResidueGroup
        self.scene.addItem(guiResidueGroup)

        # add the atoms to the group and set the reverse link
        for item in atoms.values():
            guiResidueGroup.addToGroup(item)
            item.guiNmrResidueGroup = guiResidueGroup

        # modify the group
        nmrAtoms = [atom.name for atom in nmrResidue.nmrAtoms]
        if "CB" in list(atoms.keys()):
            self._addConnectingLineToGroup(guiResidueGroup, atoms['CA'], atoms['CB'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
            self._addConnectingLineToGroup(guiResidueGroup, atoms['H'], atoms['N'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        self._addConnectingLineToGroup(guiResidueGroup, atoms['N'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)
        self._addConnectingLineToGroup(guiResidueGroup, atoms['CO'], atoms['CA'], self._lineColour, 1.0, lineList=lineList, lineId=nmrResidue)

        # self._addGroupResiduePredictions(guiResidueGroup, nmrResidue, atoms['CA'])

        return guiResidueGroup

    def addSideChainAtoms(self, nmrResidue, cbAtom, atoms, colour, lineList):
        """Add the sideChain atoms and connecting lines above the backbone line.
        """
        residue = {}
        for k, v in ATOM_POSITION_DICT[nmrResidue.residueType].items():
            if k != 'boundAtoms':
                position = [cbAtom.x() + v[0], cbAtom.y() + v[1]]
                nmrAtom = nmrResidue.fetchNmrAtom(name=k)
                newAtom = self._createGuiNmrAtom(k, position, nmrAtom)
                self.scene.addItem(newAtom)
                residue[k] = newAtom
                atoms[k] = newAtom

                # self.guiNmrAtomDict[nmrAtom] = newAtom

        # for boundAtomPair in ATOM_POSITION_DICT[nmrResidue.residueType]['boundAtoms']:
        #     atom1 = residue[boundAtomPair[0]]
        #     atom2 = residue[boundAtomPair[1]]
        #     newLine = AssignmentLine(atom1.x(), atom1.y(), atom2.x(), atom2.y(), colour, 1.0)
        #     self.scene.addItem(newLine)

    def _addNmrAtomToGuiResidues(self, nmrAtom, atomSpacing=None):
        """Update the guiNmrAtoms for a newly created/undeleted nmrAtom.
        Assumes that the nmrAtom/guiNmrAtom does not exist
        """
        nmrResidue = nmrAtom.nmrResidue

        atoms = {}
        if atomSpacing:
            self.atomSpacing = atomSpacing
        nmrAtoms = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
        residueAtoms = DEFAULT_RESIDUE_ATOMS.copy()

        for k, v in residueAtoms.items():
            if k in nmrAtoms:
                fetchedNmrAtom = nmrResidue.fetchNmrAtom(name=k)
            else:
                fetchedNmrAtom = None
            if fetchedNmrAtom is nmrAtom:
                atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)

        if nmrResidue in self.predictedStretch:
            ii = self.predictedStretch.index(nmrResidue)
            self.guiResiduesShown[ii].update(atoms)

        guiResidueGroup = self.guiNmrResidues[nmrResidue]

        # add the atoms to the group and set the reverse link
        for item in atoms.values():
            guiResidueGroup.addToGroup(item)
            item.guiNmrResidueGroup = guiResidueGroup

    def addResidue(self, nmrResidue: NmrResidue, nmrResidueIndex: int, atomSpacing=None, lineList=None):
        """Takes an Nmr Residue and a direction, either '-1 or '+1', and adds a residue to the sequence graph
        corresponding to the Nmr Residue.
        Nmr Residue name displayed beneath CA of residue drawn and residue type predictions displayed
        beneath Nmr Residue name
        """
        atoms = {}
        pos = np.array([0.0, 0.0])
        if atomSpacing:
            self.atomSpacing = atomSpacing
        nmrAtoms = [nmrAtom.name for nmrAtom in nmrResidue.nmrAtoms]
        residueAtoms = DEFAULT_RESIDUE_ATOMS.copy()

        if nmrResidue.residueType == 'GLY':
            # GLY doesn't have CB
            del residueAtoms['CB']

        # add the new nmrResidue to the current list
        if not self.guiResiduesShown:
            for k, v in residueAtoms.items():
                if k in nmrAtoms:
                    nmrAtom = nmrResidue.fetchNmrAtom(name=k)
                else:
                    nmrAtom = None
                atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)
            self.guiResiduesShown.append(atoms)
            self.predictedStretch.append(nmrResidue)

            # add the sideChain atoms
            if self._SGwidget.checkBoxes['showSideChain']['checkBox'].isChecked():
                if 'CB' in residueAtoms and nmrResidue.residueType:
                    cbAtom = atoms['CB']
                    self.addSideChainAtoms(nmrResidue, cbAtom, atoms, self._lineColour, lineList)

        else:
            for k, v in residueAtoms.items():
                if k in nmrAtoms:
                    nmrAtom = nmrResidue.fetchNmrAtom(name=k)
                else:
                    nmrAtom = None
                atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)

            # # insert into the list at the correct position
            # if direction == '-1':
            #     self.guiResiduesShown.insert(0, atoms)
            #     self.predictedStretch.insert(0, nmrResidue)
            # else:
            #     self.guiResiduesShown.append(atoms)
            #     self.predictedStretch.append(nmrResidue)

            self.guiResiduesShown.insert(nmrResidueIndex, atoms)
            self.predictedStretch.insert(nmrResidueIndex, nmrResidue)

            # add the sideChain atoms
            if self._SGwidget.checkBoxes['showSideChain']['checkBox'].isChecked():
                if 'CB' in residueAtoms and nmrResidue.residueType:
                    cbAtom = atoms['CB']
                    self.addSideChainAtoms(nmrResidue, cbAtom, atoms, self._lineColour, lineList)

        newGuiResidueGroup = self._assembleGroupResidue(nmrResidue, atoms, lineList=lineList)  #, pos[0])

        return newGuiResidueGroup

    def _addGroupResiduePredictions(self, group: GuiNmrResidueGroup, nmrResidue: NmrResidue, caAtom: GuiNmrAtom):
        """Gets predictions for residue type based on BMRB statistics and determines label positions
        based on caAtom position.
        """
        predictions = list(set(map(tuple, (getNmrResiduePrediction(nmrResidue, self.project.chemicalShiftLists[0])))))
        predictions.sort(key=lambda a: float(a[1][:-1]), reverse=True)
        for prediction in predictions:
            predictionLabel = QtWidgets.QGraphicsTextItem()
            predictionLabel.setPlainText(prediction[0] + ' ' + prediction[1])
            predictionLabel.setDefaultTextColor(QtGui.QColor(self._textColour))
            predictionLabel.setFont(textFontSmallBold)
            predictionLabel.setPos(caAtom.x() - caAtom.boundingRect().width() / 2,
                                   caAtom.y() + (30 * (predictions.index(prediction) + 2)))
            group.addToGroup(predictionLabel)

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

    def predictSequencePosition(self, nmrResidues: list):
        """
        Predicts sequence position for Nmr residues displayed in the Assigner and highlights appropriate
        positions in the Sequence Module if it is displayed.
        """
        if self.project.chains and self.project.chemicalShiftLists:

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

    def _addConnectingLineToGroup(self, group: GuiNmrResidueGroup, atom1: GuiNmrAtom, atom2: GuiNmrAtom,
                                  colour: str, width: float, displacement: float = None, style: str = None,
                                  peak: Peak = None, lineList=None, lineId='default'):
        """Adds a line between two GuiNmrAtoms using the width, colour, displacement and style specified.
        """
        newLine = AssignmentLine(0, 0, 0, 0, colour, width,
                                 parent=self, style=style, peak=peak,
                                 atom1=atom1, atom2=atom2, displacement=displacement)

        # not sure why this is different
        # group.addToGroup(newLine)
        newLine.setParentItem(group)

        itemKey = id(lineId)
        if itemKey not in lineList:
            lineList[itemKey] = []
        lineList[itemKey].append(newLine)
        return newLine

    def _createGuiNmrAtom(self, atomType: str, position: tuple, nmrAtom: NmrAtom = None) -> GuiNmrAtom:
        """Creates a GuiNmrAtom specified by the atomType and graphical position supplied.
        GuiNmrAtom can be linked to an NmrAtom by supplying it to the function.
        """
        atom = GuiNmrAtom(self.mainWindow, text=atomType, pos=position, nmrAtom=nmrAtom)
        self.guiNmrAtomDict[nmrAtom] = atom
        return atom

    def _createGhostGuiNmrAtom(self, atomType: str, position: tuple, nmrAtom: NmrAtom = None) -> GuiNmrAtom:
        """Creates a GuiNmrAtom specified by the atomType and graphical position supplied.
        GuiNmrAtom can be linked to an NmrAtom by supplying it to the function.
        """
        atom = GuiNmrAtom(self.mainWindow, text=atomType, pos=position, nmrAtom=nmrAtom)
        self.guiNmrAtomDict[nmrAtom] = atom
        return atom

    def _getPeakAssignmentsForResidue(self, nmrResidue, nmrAtomIncludeList=None):
        """Get the list of peak assignments from the nmrAtoms
        interResidueAtomPairing is the linking within the same nmrResidue
        interChainAtomPairing is the linking within the same chain but to different nmrResidues
        crossChainAtomPairing is the linking to different chains
        """

        # create a set of sets ordered by spectra
        interResidueAtomPairing = OrderedDict((spec, set()) for spec in self.magnetisationTransfers.keys())
        interChainAtomPairing = OrderedDict((spec, set()) for spec in self.magnetisationTransfers.keys())
        crossChainAtomPairing = OrderedDict((spec, set()) for spec in self.magnetisationTransfers.keys())

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
                    for mag in self.magnetisationTransfers[spec]:
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

    def _addPeakAssignmentLinesToGroup(self, assignments, lineList):
        for specAssignments in assignments.values():
            for nmrAtomPair in specAssignments:

                guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
                                  self.guiNmrAtomDict.get(nmrAtomPair[1]),
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
                try:
                    group = self.guiNmrResidues[guiNmrAtomPair[0].nmrAtom.nmrResidue]
                except Exception as es:
                    pass

                self._addConnectingLineToGroup(group,
                                               guiNmrAtomPair[0],
                                               guiNmrAtomPair[1],
                                               spectrum.positiveContourColour,
                                               2.0, displacement=displacement,
                                               peak=peak, lineList=lineList, lineId=peak)

                # update displacements for both guiNmrAtoms
                guiNmrAtomPair[0].addConnectedList(guiNmrAtomPair[1])
                guiNmrAtomPair[1].addConnectedList(guiNmrAtomPair[0])

    def _addPeakAssignmentLinesToAdjacentGroup(self, nmrResidue, assignments, peaklineList, connectingLineList):
        for specAssignments in assignments.values():
            for nmrAtomPair in specAssignments:

                guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
                                  self.guiNmrAtomDict.get(nmrAtomPair[1]),
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

                    newGhostResidue = self.addGhostResidue(nmrAtomPair[0].nmrResidue,
                                                           guiNmrAtomPair[1],
                                                           nmrAtomPair[1].nmrResidue,
                                                           nmrAtomPair[0].name,
                                                           nmrAtomPair[1].name,
                                                           True,
                                                           lineList=connectingLineList)
                    guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
                                      self.guiNmrAtomDict.get(nmrAtomPair[1]),
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

                    newGhostResidue = self.addGhostResidue(nmrAtomPair[1].nmrResidue,
                                                           guiNmrAtomPair[0],
                                                           nmrAtomPair[0].nmrResidue,
                                                           nmrAtomPair[1].name,
                                                           nmrAtomPair[0].name,
                                                           True,
                                                           lineList=connectingLineList)
                    guiNmrAtomPair = (self.guiNmrAtomDict.get(nmrAtomPair[0]),
                                      self.guiNmrAtomDict.get(nmrAtomPair[1]),
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

    def addGhostResidue(self, nmrResidueCon1: NmrResidue,
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

        with logCommandBlock(get='self') as log:
            log('navigateToNmrResidue', selectedNmrResidue=repr(nmrResidue.pid))

            # optionally clear the marks
            if self._SGwidget.checkBoxes['autoClearMarks']['checkBox'].isChecked():
                self.mainWindow.clearMarks()

            # navigate the displays
            for display in displays:
                if len(display.strips) > 0 and display.strips[0].spectrumViews:
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
        from ccpn.ui.gui.widgets.Menu import Menu
        from functools import partial

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
                            contextMenu.addAction(nmrAtom.id, partial(self.deassignPeak, thisLine._peak, nmrAtom))

                contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
                contextMenu.exec()
                contextMenu = None

        elif isinstance(pressed, GuiNmrResidue):

            # create the nmrResidue menu
            contextMenu.addAction(self.disconnectPreviousIcon, 'disconnect Previous nmrResidue', partial(self.disconnectPreviousNmrResidue))
            contextMenu.addAction(self.disconnectIcon, 'disconnect nmrResidue', partial(self.disconnectNmrResidue))
            contextMenu.addAction(self.disconnectNextIcon, 'disconnect Next nmrResidue', partial(self.disconnectNextNmrResidue))
            contextMenu.addSeparator()
            contextMenu.addAction('disconnect all nmrResidues', partial(self.disconnectAllNmrResidues))
            if object.nmrResidue.residue:
                contextMenu.addSeparator()
                contextMenu.addAction('deassign nmrChain', partial(self.deassignNmrChain))

            contextMenu.addSeparator()
            contextMenu.addAction('Show nmrResidue', partial(self.showNmrResidue, object))
            contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
            contextMenu.exec()
            contextMenu = None

        elif isinstance(pressed, GuiNmrAtom):

            # create the nmrAtom menu
            contextMenu.addAction('deassign nmrAtoms from Peaks')
            contextMenu.addSeparator()
            if pressed.nmrAtom and pressed.nmrAtom.assignedPeaks:

                # add nmrAtoms to the menu
                for peak in pressed.nmrAtom.assignedPeaks:
                    if peak and peak.assignedNmrAtoms:
                        subMenu = contextMenu.addMenu(peak.id)

                        subMenu.addAction('nmrAtoms')
                        subMenu.addSeparator()
                        for nmrAtomList in peak.assignedNmrAtoms:
                            for nmrAtom in nmrAtomList:
                                if nmrAtom:
                                    subMenu.addAction(nmrAtom.id, partial(self.deassignPeak, peak, nmrAtom))

                contextMenu.move(cursor.pos().x(), cursor.pos().y() + 10)
                contextMenu.exec()
                contextMenu = None

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
                         'CO': np.array([2 * atomSpacing, -1 * atomSpacing])
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

if __name__ == '__main__':
    from ccpn.ui.gui.widgets.Application import TestApplication
    from ccpn.ui.gui.widgets.TextEditor import TextEditor


    app = TestApplication()

    popup = SequenceGraphModule()

    popup.show()
    popup.raise_()
    app.start()
