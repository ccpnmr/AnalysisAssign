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

    def getConnectedList(self, connectedAtom):
        # keyVal = str(connectedAtom.nmrAtom.pid)

        keyVal = connectedAtom
        if keyVal in self.connectedList:
            return self.connectedList[keyVal]
        else:
            return 0

    def clearConnectedList(self):
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

        if atom2.x() < atom1.x():
            x1 = atom1.x()  # + rx1
            y1 = atom1.y()  # + ry1
            x2 = atom2.x() + (rx2 - rx1)
            y2 = atom2.y() + (ry2 - ry1)
        else:
            x1 = atom2.x() + (rx2 - rx1)
            y1 = atom2.y() + (ry2 - ry1)
            x2 = atom1.x()  # + rx1
            y2 = atom1.y()  # + ry1

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

        # # add a callback that fires when the layout changes the state of the checkbox
        # # move this into the settings widget as _init callback
        # self._SGwidget.assignmentsTreeCheckBox.checkBox.stateChanged.connect(self._checkLayoutInit)

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

        if nmrChain is not None:
            self.selectSequence(nmrChain)

        # stop the mainWidget from squishing during a resize
        self.mainWidget.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)

        # install the event filter to handle maximising from floated dock
        # self.installMaximiseEventHandler(self._maximise, self._closeModule)

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
                    self.setUpdatesEnabled(False)
                    self.scene.blockSignals(True)

                    self._rebuildPeakAssignments(nmrChainPid)

                    self.scene.blockSignals(False)
                    self.setUpdatesEnabled(True)

                    # resize to the new items and spawns a repaint
                    self.scene.setSceneRect(self.scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))

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
                        self.setNmrChainDisplay(nmrChain)

    def _registerNotifiers(self):
        """Register the required notifiers
        """
        # self.current.registerNotify(self._updateModule, 'nmrChains')      # doesn't work

        # self._nmrResidueNotifier = self.setNotifier(self.project,
        #                                             [Notifier.RENAME, Notifier.CHANGE],
        #                                             NmrResidue.__name__,
        #                                             self._resetNmrResiduePidForAssigner,
        #                                             onceOnly=True)
        #
        # self._peakNotifier = self.setNotifier(self.project,
        #                                       [Notifier.CHANGE],
        #                                       Peak.__name__,
        #                                       self._updateShownAssignments,
        #                                       onceOnly=True)

        self._peakNotifier = self.setNotifier(self.project,
                                      [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
                                      Peak.className,
                                      self._updatePeaks,
                                      onceOnly=True)

        self._nmrAtomNotifier = self.setNotifier(self.project,
                                      [Notifier.CHANGE, Notifier.DELETE],
                                      NmrAtom.className,
                                      self._updateNmrAtoms,
                                      onceOnly=True)

        # self._spectrumNotifier = self.setNotifier(self.project,
        #                                           [Notifier.CHANGE],
        #                                           Spectrum.__name__,
        #                                           self._updateShownAssignments)

        # notifier to change the magnetisationTransfer list when new spectrum added
        self._spectrumListNotifier = self.setNotifier(self.project,
                                                      [Notifier.CREATE, Notifier.DELETE],
                                                      Spectrum.className,
                                                      self._updateSpectra)

        # # notifier for changing the selected chain - draw new display
        # self._nmrChainNotifier = self.setNotifier(self.project,
        #                                           [Notifier.CHANGE, Notifier.DELETE],
        #                                           NmrChain.__name__,
        #                                           self._updateChain,
        #                                           onceOnly=True)

    # def _unRegisterNotifiers(self):
    #     # use the new notifier class
    #     if self._nmrResidueNotifier:
    #         self._nmrResidueNotifier.unRegister()
    #     if self._peakNotifier:
    #         self._peakNotifier.unRegister()
    #     if self._spectrumNotifier:
    #         self._spectrumNotifier.unRegister()
    #     if self._nmrChainNotifier:
    #         self._nmrChainNotifier.unRegister()

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
        # self.showNmrChainFromPulldown(nmrChain)

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

    def _deleteNmrResidue(self, data):
        """
        delete the nmrResidue form the scene
        """
        nmrResidue = data[Notifier.OBJECT]
        print('>>>_deleteNmrResidue', nmrResidue)

        if nmrResidue in self.predictedStretch:
            # print('  >>>delete items')

            sceneItems = self.scene.items()
            for item in sceneItems:
                if isinstance(item, GuiNmrResidueGroup):
                    if item.nmrResidue == nmrResidue:
                        # print('    >>>delete groupItem:', item)
                        self.scene.removeItem(item)
                # elif isinstance(item, GuiNmrAtom):
                #   if item.nmrAtom and item.nmrAtom.nmrResidue == nmrResidue:
                #     print('    >>>delete atomItem:', item)
                #     self.scene.removeItem(item)

            self.predictedStretch.remove(nmrResidue)
        # print('>>>end')

    def _updatePeaks(self, data):
        """
        update the peaks in the display
        :param data:
        """
        peak = data[Notifier.OBJECT]

        print('>>>_updatePeaks', peak)
        trigger = data[Notifier.TRIGGER]

        if trigger == Notifier.DELETE:

            print ('>>>delete peak', peak)

            items = [item for item in self.assignmentLines if item._peak is peak]
            for item in items:
                print ('  removing', item)
                self.scene.removeItem(item)

        elif trigger == Notifier.CREATE:

            print ('>>>create peak - no action', peak)

        elif trigger == Notifier.CHANGE:

            print ('>>>change peak', peak)

    def _updateNmrAtoms(self, data):
        """
        update the nmrAtoms in the display
        :param data:
        """
        nmrAtom = data[Notifier.OBJECT]

        print('>>>_updateNmrAtoms', nmrAtom)
        trigger = data[Notifier.TRIGGER]

        if trigger == Notifier.DELETE:

            print ('>>>delete nmrAtom', nmrAtom)

            items = [item for item in self.assignmentLines if item.atom1.nmrAtom is nmrAtom or item.atom2.nmrAtom is nmrAtom]
            for item in items:
                print ('  removing', item)
                self.scene.removeItem(item)

        elif trigger == Notifier.CREATE:

            print ('>>>create nmrAtom - no action', nmrAtom)

        elif trigger == Notifier.CHANGE:

            print ('>>>change nmrAtom', nmrAtom)

    def _rebuildPeakAssignments(self, nmrChainOrPid):

        if isinstance(nmrChainOrPid, str):
            if not Pid.isValid(nmrChainOrPid):
                return
            nmrChain = self.project.getByPid(nmrChainOrPid)
        else:
            nmrChain = nmrChainOrPid

        # nmrChainOrPid could be '<Select>' in which case nmrChain would be None
        if not nmrChain:
            return

        # update the endpoints
        for peakLine in self.assignmentLines:
            self.scene.removeItem(peakLine)

        self.assignmentLines = []

        # reset the displacement values
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
        for peakLine in self.assignmentLines:
            peakLine.updateEndPoints()

    # @profile
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

            ###nmrChain = self.project.getByPid(nmrChainPid)
            ###if self.modePulldown.currentText() == 'fragment':
            if True:

                connectingLinesNeeded = set()
                if self.nmrResiduesCheckBox.isChecked():
                    for nmrResidue in nmrChain.nmrResidues:
                        if nmrResidue is nmrResidue.mainNmrResidue:
                            self.addResidue(nmrResidue, '+1', lineList=self.connectingLines)

                            # add a connecting line to the adjacent residue
                            if nmrResidue.nextNmrResidue:
                                connectingLinesNeeded.add(len(self.guiResiduesShown) - 1)

                            # if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
                            #     # add the internally connected Lines
                            #     internalAssignments, interChainAssignments = self._getPeakAssignmentsForResidue(nmrResidue)
                            #     self._addPeakAssignmentLinesToGroup(internalAssignments)
                            #     self._addPeakAssignmentLinesToAdjacentGroup(interChainAssignments)

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
                        self.addResidue(nmrResidue, '+1', lineList=self.connectingLines)
                        nmrResidue = nmrResidue.nextNmrResidue

                if len(self.predictedStretch) > 2:
                    # TODO:ED causes a crash from here GuiNmrResidue has been deleted
                    self.predictSequencePosition(self.predictedStretch)

                ###elif self.modePulldown.currentText() == 'Assigned - backbone':
                ###  self._showBackboneAssignments(nmrChain)

                # add the connecting lines
                for ii, res in enumerate(self.guiResiduesShown[:-1]):
                    if not self.nmrResiduesCheckBox.isChecked() or ii in connectingLinesNeeded:
                        self._addConnectingLineToGroup(tuple(self.guiNmrResidues.values())[ii],
                                                       res['CO'], self.guiResiduesShown[ii + 1]['N'],
                                                       self._lineColour, 1.0, lineList=self.connectingLines)

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

                # update the group positions
                for ii, res in enumerate(self.guiNmrResidues.items()):
                    # second element of tuple; res is (k, v)
                    res[1].setPos(QtCore.QPointF(ii * self.atomSpacing * 3.0, 0.0))

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
                for connectingLine in self.connectingLines:
                    connectingLine.updateEndPoints()
                for peakLine in self.assignmentLines:
                    peakLine.updateEndPoints()

            # if self._SGwidget.assignmentsCheckBox.isChecked():
            # if self._SGwidget.checkBoxes['peakAssignments']['checkBox'].isChecked():
            #     self._getAssignmentsFromSpectra()

        # self.scene.setSceneRect(self.scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))  # resize to the new items
        self.nmrChain = nmrChain

    def resetSequenceGraph(self):

        self.nmrChainPulldown.pulldownList.select('NC:@-')

    def _closeModule(self):
        """
        CCPN-INTERNAL: used to close the module
        """
        # self._unRegisterNotifiers()
        self.thisSequenceModule.close()
        super()._closeModule()

    def close(self):
        """
        Close the table from the commandline
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
            elif self.current.nmrResidue.mainNmrResidue.previousNmrResidue:
                with progressManager(self.mainWindow, 'unlinking Next NmrResidue to:\n ' + selected):
                    try:
                        self.current.nmrResidue.unlinkNextNmrResidue()
                    except Exception as es:
                        showWarning(str(self.windowTitle()), str(es))

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

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()
                # self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
            ###self.updateNmrResidueTable()

    def disconnectNmrResidue(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            with progressManager(self.mainWindow, 'disconnecting NmrResidue:\n ' + selected):
                try:
                    self.current.nmrResidue.disconnect()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()
                # self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
            #self.updateNmrResidueTable()

    def disconnectNextNmrResidue(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            with progressManager(self.mainWindow, 'disconnecting Next NmrResidue to:\n ' + selected):
                try:
                    self.current.nmrResidue.disconnectNext()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()
                # self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
            #self.updateNmrResidueTable()

    def disconnectAllNmrResidues(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.pid)
            with progressManager(self.mainWindow, 'disconnecting all NmrResidues connected to:\n ' + selected):
                try:
                    self.current.nmrResidue.disconnectAll()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()
                # self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
            #self.updateNmrResidueTable()

    def deassignNmrChain(self, selectedNmrResidue=None):
        if self.current.nmrResidue:
            selected = str(self.current.nmrResidue.nmrChain.pid)
            with progressManager(self.mainWindow, 'deassigning nmrResidues in NmrChain:\n ' + selected):
                try:
                    self.current.nmrResidue.deassignNmrChain()
                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

            if self.current.nmrResidue:
                self.showNmrChainFromPulldown()
                # self.setNmrChainDisplay(self.current.nmrResidue.nmrChain.pid)
            #self.updateNmrResidueTable()

    def _resetNmrResiduePidForAssigner(self, data):  #nmrResidue, oldPid:str):
        """Reset pid for NmrResidue and all offset NmrResidues
        """
        print('>>>_resetNmrResiduePidForAssigner')

        nmrResidue = data['object']

        nmrChainPid = self.nmrChainPulldown.getText()
        if self.project.getByPid(nmrChainPid):

            for nr in [nmrResidue] + list(nmrResidue.offsetNmrResidues):
                for guiNmrResidueGroup in self.guiNmrResidues.values():
                    if guiNmrResidueGroup.nmrResidue is nr:
                        guiNmrResidueGroup.nmrResidueLabel._update()

    def deassignPeak(self, selectedPeak=None, selectedNmrAtom=None):
        """Deassign the peak by removing the assigned nmrAtoms from the list
        """
        if selectedPeak:
            with logCommandBlock(get='self') as log:
                log('deassignPeak', selectedPeak=repr(selectedPeak.pid), selectedNmrAtom=repr(selectedNmrAtom.pid))

                try:
                    newList = []
                    for atomList in selectedPeak.assignedNmrAtoms:
                        atoms = [atom for atom in list(atomList) if atom != selectedNmrAtom]
                        newList.append(tuple(atoms))

                    selectedPeak.assignedNmrAtoms = tuple(newList)

                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

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

    def clearAllItems(self):
        """
        Removes all displayed residues in the sequence graph and resets items count to zero.
        """
        # qDeleteAll(self.sene.items)
        # if self.scene:
        #   for item in self.scene.items():
        #     self.scene.removeItem(item)
        self.residueCount = 0
        self.predictedStretch = []
        self.guiResiduesShown = []
        self.guiNmrResidues = OrderedDict()
        self.guiNmrResidueLabels = []
        self.guiNmrAtomDict = {}
        self.guiGhostNmrResidues = OrderedDict()
        self.ghostList = {}
        self.connectingLines = []
        self.assignmentLines = []
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

        # self.scrollContents.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)

    # def _assembleResidue(self, nmrResidue: NmrResidue, atoms: typing.Dict[str, GuiNmrAtom]):
    #     """
    #     Takes an Nmr Residue and a dictionary of atom names and GuiNmrAtoms and
    #     creates a graphical representation of a residue in the assigner
    #     """
    #
    #     for item in atoms.values():
    #         self.scene.addItem(item)
    #
    #     nmrAtoms = [atom.name for atom in nmrResidue.nmrAtoms]
    #     if "CB" in list(atoms.keys()):
    #         self._addConnectingLine(atoms['CA'], atoms['CB'], self._lineColour, 1.0, 0)
    #
    #     if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
    #         self._addConnectingLine(atoms['H'], atoms['N'], self._lineColour, 1.0, 0)
    #
    #     # if nmrResidue.residueType != 'PRO':
    #     #     self._addConnectingLine(atoms['H'], atoms['N'], self._lineColour, 1.0, 0)
    #     # else:
    #     #     self.scene.removeItem(atoms['H'])
    #
    #     # if not 'CB' in nmrAtoms:
    #     #   self.scene.removeItem(atoms['CB'])
    #     #   self.scene.removeItem(cbLine)
    #
    #     self._addConnectingLine(atoms['N'], atoms['CA'], self._lineColour, 1.0, 0)
    #     self._addConnectingLine(atoms['CO'], atoms['CA'], self._lineColour, 1.0, 0)
    #     self.nmrResidueLabel = GuiNmrResidue(self, nmrResidue, atoms['CA'])
    #     self.guiNmrResidues.append(self.nmrResidueLabel)
    #     self.scene.addItem(self.nmrResidueLabel)
    #     self._addResiduePredictions(nmrResidue, atoms['CA'])

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
            self._addConnectingLineToGroup(guiResidueGroup, atoms['CA'], atoms['CB'], self._lineColour, 1.0, lineList=lineList)

        if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
            self._addConnectingLineToGroup(guiResidueGroup, atoms['H'], atoms['N'], self._lineColour, 1.0, lineList=lineList)

        self._addConnectingLineToGroup(guiResidueGroup, atoms['N'], atoms['CA'], self._lineColour, 1.0, lineList=lineList)
        self._addConnectingLineToGroup(guiResidueGroup, atoms['CO'], atoms['CA'], self._lineColour, 1.0, lineList=lineList)

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
            self._addConnectingLineToGroup(guiResidueGroup, atoms['CA'], atoms['CB'], self._lineColour, 1.0, lineList=lineList)

        if "H" in list(atoms.keys()) and nmrResidue.residueType != 'PRO':
            self._addConnectingLineToGroup(guiResidueGroup, atoms['H'], atoms['N'], self._lineColour, 1.0, lineList=lineList)

        self._addConnectingLineToGroup(guiResidueGroup, atoms['N'], atoms['CA'], self._lineColour, 1.0, lineList=lineList)
        self._addConnectingLineToGroup(guiResidueGroup, atoms['CO'], atoms['CA'], self._lineColour, 1.0, lineList=lineList)

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

    def addResidue(self, nmrResidue: NmrResidue, direction: str, atomSpacing=None, lineList=None):
        """
        Takes an Nmr Residue and a direction, either '-1 or '+1', and adds a residue to the sequence graph
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

                # if direction == '-1':
                #
                #     # # use the 'H' as the reference
                #     # pos = np.array([self.guiResiduesShown[0]['H'].x() - 3 * self.atomSpacing, self.guiResiduesShown[0]['H'].y()])
                #     # atoms[k] = self._createGuiNmrAtom(k, v + pos, nmrAtom)
                #     atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)
                #
                # else:
                #     # pos = np.array([self.guiResiduesShown[-1]['H'].x() + 3 * self.atomSpacing, self.guiResiduesShown[-1]['H'].y()])
                #     # atoms[k] = self._createGuiNmrAtom(k, v + pos, nmrAtom)
                #     atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)

                atoms[k] = self._createGuiNmrAtom(k, v, nmrAtom)

            # insert into the list at the correct position
            if direction == '-1':
                self.guiResiduesShown.insert(0, atoms)
                self.predictedStretch.insert(0, nmrResidue)
            else:
                self.guiResiduesShown.append(atoms)
                self.predictedStretch.append(nmrResidue)

            # add the sideChain atoms
            if self._SGwidget.checkBoxes['showSideChain']['checkBox'].isChecked():
                if 'CB' in residueAtoms and nmrResidue.residueType:
                    cbAtom = atoms['CB']
                    self.addSideChainAtoms(nmrResidue, cbAtom, atoms, self._lineColour, lineList)

        newGuiResidueGroup = self._assembleGroupResidue(nmrResidue, atoms, lineList=lineList)  #, pos[0])

        # # move to the required position
        # newGuiResidueGroup.setPos(QtCore.QPointF(pos[0], pos[1]))

        return newGuiResidueGroup

    def _addGroupResiduePredictions(self, group: GuiNmrResidueGroup, nmrResidue: NmrResidue, caAtom: GuiNmrAtom):
        """
        Gets predictions for residue type based on BMRB statistics and determines label positions
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

    # def _updateShowTreeAssignments(self, peak=None):
    #     nmrChainPid = self.nmrChainPulldown.getText()
    #     if nmrChainPid:
    #         self.setUpdatesEnabled(False)
    #         self.scene.blockSignals(True)
    #
    #         self.setNmrChainDisplay(nmrChainPid)
    #
    #         # test moving the groups
    #         for ii, res in enumerate(self.guiNmrResidues.items()):
    #             print('>>>', res[0])
    #             res[1].setPos(QtCore.QPointF(ii * self.atomSpacing * 3.5, 0))
    #
    #         self.scene.blockSignals(False)
    #         self.setUpdatesEnabled(True)

    def showNmrChainFromPulldown(self, data=None):
        ###if self.current.nmrChain is not None:
        ###  self.setNmrChainDisplay(self.current.nmrChain.pid)
        print('>>>showNmrChainFromPulldown')

        nmrChainPid = self.nmrChainPulldown.getText()
        if nmrChainPid:
            self.setUpdatesEnabled(False)
            self.scene.blockSignals(True)

            self.setNmrChainDisplay(nmrChainPid)

            self.scene.blockSignals(False)
            self.setUpdatesEnabled(True)

            # resize to the new items and spawns a repaint
            self.scene.setSceneRect(self.scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))

    # def _updateChain(self, data):
    #     ###if self.current.nmrChain is not None:
    #     ###  self.setNmrChainDisplay(self.current.nmrChain.pid)
    #     print('>>>_updateChain')
    #
    #     nmrChainPid = self.nmrChainPulldown.getText()
    #     if nmrChainPid:
    #         self.setUpdatesEnabled(False)
    #         self.scene.blockSignals(True)
    #
    #         self.setNmrChainDisplay(nmrChainPid)
    #
    #         # test moving the groups
    #         for ii, res in enumerate(self.guiNmrResidues.items()):
    #             print('>>>', res[0])
    #             res[1].setPos(QtCore.QPointF(ii * self.atomSpacing * 3.5, 0))
    #
    #         self.scene.blockSignals(False)
    #         self.setUpdatesEnabled(True)

    # def _addConnectingLine(self, atom1: GuiNmrAtom, atom2: GuiNmrAtom,
    #                        colour: str, width: float, displacement: float = None, style: str = None,
    #                        peak: Peak = None, lineList=None):
    #     """Adds a line between two GuiNmrAtoms using the width, colour, displacement and style specified.
    #     """
    #     newLine = AssignmentLine(0, 0, 0, 0, colour, width,
    #                              parent=self, style=style, peak=peak,
    #                              atom1=atom1, atom2=atom2, displacement=displacement)
    #
    #     self.scene.addItem(newLine)
    #     lineList.append(newLine)
    #     return newLine

    def _addConnectingLineToGroup(self, group: GuiNmrResidueGroup, atom1: GuiNmrAtom, atom2: GuiNmrAtom,
                                  colour: str, width: float, displacement: float = None, style: str = None,
                                  peak: Peak = None, lineList=None):
        """Adds a line between two GuiNmrAtoms using the width, colour, displacement and style specified.
        """
        newLine = AssignmentLine(0, 0, 0, 0, colour, width,
                                 parent=self, style=style, peak=peak,
                                 atom1=atom1, atom2=atom2, displacement=displacement)

        # not sure why this is different
        # group.addToGroup(newLine)
        newLine.setParentItem(group)

        lineList.append(newLine)
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

    def _getPeakAssignmentsForResidue(self, nmrResidue):
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

            for peak in nmrAtom.assignedPeaks:
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

                        # if assignment[conNum].nmrResidue.relativeOffset == +1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:
                        #
                        #     # this is a plus residue so find connected, have to traverse to the nextNmrResidue
                        #     # will it always exist?
                        #     conName = assignment[conNum].name
                        #     preN = assignment[conNum].nmrResidue.mainNmrResidue.nextNmrResidue
                        #     if preN:
                        #         newConSwap = [nmrA for nmrA in preN.nmrAtoms if nmrA.name == conName]
                        #         if newConSwap:
                        #             newCon[conNum] = newConSwap[0]
                        #     else:
                        #         newCon[conNum] = None                   # not connected so skip

                    assignment = newCon

                    # only get the assignments a-b if a and b are defined in the spectrum magnetisationTransfers list
                    for mag in self.magnetisationTransfers[spec]:
                        nmrAtom0 = assignment[mag[0] - 1]
                        nmrAtom1 = assignment[mag[1] - 1]
                        nmrAtom0 = nmrAtom0 if nmrAtom0 else None
                        nmrAtom1 = nmrAtom1 if nmrAtom1 else None

                        if not None in (nmrAtom0, nmrAtom1):

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

                            # if (nmrAtom1, nmrAtom0, peak) not in atomPairList:
                            #     atomPairList.add((nmrAtom0, nmrAtom1, peak))

                        # need to sort by spectra

                    # nmrAtoms = [nA for nA in assignment if nA and nA.nmrResidue is nmrResidue and nA is not nmrAtom]
                    #
                    # for nA in nmrAtoms:
                    #     if (nA, nmrAtom, peak) not in atomPairList:
                    #         atomPairList.add((nmrAtom, nA, peak))

        return interResidueAtomPairing, interChainAtomPairing, crossChainAtomPairing

    # def _getPeakAssignmentsForResidueExternal(self, nmrResidue, mTList=None):
    #     """Get the list of peak assignments from the nmrAtoms, only peaks exterior to the nmrResidue
    #     """
    #
    #     # generate the list that defines which couplings there are between the nmrAtoms attached to each peak - move outside
    #     mTList = OrderedDict()
    #     for spec in self.project.spectra:
    #         mTList[spec] = {}
    #         for mt in spec.magnetisationTransfers:
    #             mTList[spec][mt] = set()
    #
    #     # create a set of sets ordered by spectra
    #     atomPairing = OrderedDict((spec, set()) for spec in mTList.keys())
    #     atomPairingCount = 0
    #
    #     for nmrAtom in nmrResidue.nmrAtoms:
    #
    #         for peak in nmrAtom.assignedPeaks:
    #             spec = peak.peakList.spectrum
    #             for assignment in peak.assignments:
    #
    #                 # find the mainNmrResidue for -1 and +1 connections
    #                 newCon = list(assignment)
    #                 for conNum in range(len(assignment)):
    #                     if assignment[conNum].nmrResidue.relativeOffset == -1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:
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
    #                     # if assignment[conNum].nmrResidue.relativeOffset == +1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:
    #                     #
    #                     #     # this is a plus residue so find connected, have to traverse to the nextNmrResidue
    #                     #     # will it always exist?
    #                     #     conName = assignment[conNum].name
    #                     #     preN = assignment[conNum].nmrResidue.mainNmrResidue.nextNmrResidue
    #                     #     if preN:
    #                     #         newConSwap = [nmrA for nmrA in preN.nmrAtoms if nmrA.name == conName]
    #                     #         if newConSwap:
    #                     #             newCon[conNum] = newConSwap[0]
    #                     #     else:
    #                     #         newCon[conNum] = None                   # not connected so skip
    #
    #                 assignment = newCon
    #
    #                 # only get the assignments a-b if a and b are defined in the spectrum magnetisationTransfers list
    #                 for mag in mTList[spec]:
    #                     nmrAtom0 = assignment[mag[0] - 1]
    #                     nmrAtom1 = assignment[mag[1] - 1]
    #                     nmrAtom0 = nmrAtom0 if nmrAtom0 else None
    #                     nmrAtom1 = nmrAtom1 if nmrAtom1 else None
    #
    #                     if not None in (nmrAtom0, nmrAtom1):
    #                         if (nmrAtom0.nmrResidue is not nmrResidue) or (nmrAtom1.nmrResidue is not nmrResidue):
    #                             if (nmrAtom1, nmrAtom0, peak) not in atomPairing[spec]:
    #                                 atomPairing[spec].add((nmrAtom0, nmrAtom1, peak))
    #                                 atomPairingCount += 1
    #
    #                         # if (nmrAtom1, nmrAtom0, peak) not in atomPairList:
    #                         #     atomPairList.add((nmrAtom0, nmrAtom1, peak))
    #
    #                     # need to sort by spectra
    #
    #                 # nmrAtoms = [nA for nA in assignment if nA and nA.nmrResidue is nmrResidue and nA is not nmrAtom]
    #                 #
    #                 # for nA in nmrAtoms:
    #                 #     if (nA, nmrAtom, peak) not in atomPairList:
    #                 #         atomPairList.add((nmrAtom, nA, peak))
    #
    #     return atomPairing if atomPairingCount else None

    # def _getAssignmentsFromSpectra(self):
    #     """Get the peak assignments attached to the nmrResidues and display as lines, coloured by spectrum, to
    #     the nmrResidues.
    #     """
    #     for spectrum in self.project.spectra:
    #         connections = [x for y in list(nmrAtomPairsByDimensionTransfer(spectrum.peakLists).values())
    #                        for x in y]
    #
    #         # print ('  >>>connections')
    #         # for ii in connections:
    #         #   print ('        >>>', ii)
    #
    #         # find the minus links and update the links to the previousNmrResidue
    #         minusResList = []
    #         for inCon in connections:
    #             newCon = list(inCon)
    #             for conNum in range(0, 2):
    #                 if inCon[conNum].nmrResidue.relativeOffset == -1:  # and inCon[conNum].nmrResidue.nmrChain.isConnected:
    #
    #                     # this is a minus residue so find connected, have to traverse to the previousNmrResidue
    #                     # will it always exist?
    #                     conName = inCon[conNum].name
    #                     preN = inCon[conNum].nmrResidue.mainNmrResidue.previousNmrResidue
    #                     if preN:
    #                         newConSwap = [nmrA for nmrA in preN.nmrAtoms if nmrA.name == conName]
    #                         if newConSwap:
    #                             newCon[conNum] = newConSwap[0]
    #
    #                 # if newCon:
    #                 #   cc = (newCC[0], cc[1])    # replace the minus residue
    #
    #             minusResList.append(newCon)
    #
    #             if ('A.31' in inCon[0].pid):
    #                 pass
    #
    #         # the original routine to add the links to adjacent atoms
    #         # sometimes the link maybe a huge distance away on the scene
    #         for ii, connection in enumerate(minusResList):  # ejb - was connections
    #
    #             # nmrAtomPair = [self.project._data2Obj.get(connection[0]).nmrAtom,
    #             #                self.project._data2Obj.get(connection[1]).nmrAtom]
    #             # sorting makes sure drawing is done properly
    #             guiNmrAtomPair = [self.guiNmrAtomDict.get(a) for a in sorted(connection[0:2], reverse=True)]
    #             guiNmrResiduePair = [a for a in sorted(connection[0:2], reverse=True)]
    #
    #             if None not in guiNmrAtomPair:
    #                 # displacement = min(guiNmrAtomPair[0].connectedAtoms, guiNmrAtomPair[1].connectedAtoms)
    #
    #                 displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
    #                 # peak0 = guiNmrAtomPair[0].nmrAtom.assignedPeaks[0]  # must be true to draw the line
    #                 # peak1 = guiNmrAtomPair[1].nmrAtom.assignedPeaks[0]  # must be true to draw the line
    #
    #                 peak0 = connection[2]
    #                 peak1 = connection[2]
    #
    #                 # if (abs(guiNmrAtomPair[0].x() - guiNmrAtomPair[1].x()) < 6 * self.atomSpacing) or self.assignmentsTreeCheckBox.isChecked() is False:
    #                 if (abs(guiNmrAtomPair[0].x() - guiNmrAtomPair[1].x()) < 6 * self.atomSpacing) or \
    #                         self._SGwidget.checkBoxes['treeView']['checkBox'].isChecked() is False:
    #
    #                     if guiNmrAtomPair[0].nmrAtom.nmrResidue is guiNmrAtomPair[1].nmrAtom.nmrResidue:
    #
    #                         # add the internal line to the guiNmrResidueGroup, should now move when group is moved
    #                         group = self.guiNmrResidues[guiNmrAtomPair[0].nmrAtom.nmrResidue]
    #                         self._addConnectingLineToGroup(group,
    #                                                        guiNmrAtomPair[0],
    #                                                        guiNmrAtomPair[1],
    #                                                        spectrum.positiveContourColour,
    #                                                        2.0, displacement=displacement,
    #                                                        peak=peak0)
    #                     else:
    #
    #                         # do not connect to a group yet, possibly add to another group that requires updating on create/delete nmrResidues
    #                         self._addConnectingLine(guiNmrAtomPair[0],
    #                                                 guiNmrAtomPair[1],
    #                                                 spectrum.positiveContourColour,
    #                                                 2.0, displacement=displacement,
    #                                                 peak=peak0)
    #
    #                     guiNmrAtomPair[0].addConnectedList(guiNmrAtomPair[1])
    #                     guiNmrAtomPair[1].addConnectedList(guiNmrAtomPair[0])
    #
    #                 elif (guiNmrAtomPair[0].x() - guiNmrAtomPair[1].x()) > 0:
    #                     # add a new 'ghost' atom below the line and link to it instead
    #                     # only goes to the right so far...
    #
    #                     # print ('>>>right ', guiNmrResiduePair[0].nmrResidue.pid, guiNmrResiduePair[1].nmrResidue.pid)
    #                     tempAtoms = self.addGhostResidue(guiNmrResiduePair[1].nmrResidue,
    #                                                      guiNmrAtomPair[0],
    #                                                      guiNmrResiduePair[0].nmrResidue,
    #                                                      guiNmrResiduePair[1].name,
    #                                                      guiNmrResiduePair[0].name,
    #                                                      True)
    #
    #                     displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
    #
    #                     self._addConnectingLine(guiNmrAtomPair[0],
    #                                             tempAtoms[guiNmrResiduePair[1].name],
    #                                             spectrum.positiveContourColour,
    #                                             2.0, displacement,
    #                                             peak=peak0)
    #
    #                     guiNmrAtomPair[0].addConnectedList(guiNmrAtomPair[1])
    #                     guiNmrAtomPair[1].addConnectedList(guiNmrAtomPair[0])
    #
    #                     # make a duplicate going the other way
    #                     # print ('>>>left  ', guiNmrResiduePair[0].nmrResidue.pid, guiNmrResiduePair[1].nmrResidue.pid)
    #                     tempAtoms = self.addGhostResidue(guiNmrResiduePair[0].nmrResidue,
    #                                                      guiNmrAtomPair[1],
    #                                                      guiNmrResiduePair[1].nmrResidue,
    #                                                      guiNmrResiduePair[0].name,
    #                                                      guiNmrResiduePair[1].name,
    #                                                      False)
    #
    #                     displacement = guiNmrAtomPair[1].getConnectedList(guiNmrAtomPair[0])
    #
    #                     self._addConnectingLine(guiNmrAtomPair[1],
    #                                             tempAtoms[guiNmrResiduePair[0].name],
    #                                             spectrum.positiveContourColour,
    #                                             2.0, displacement,
    #                                             peak=peak1)

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
                group = self.guiNmrResidues[guiNmrAtomPair[0].nmrAtom.nmrResidue]
                self._addConnectingLineToGroup(group,
                                               guiNmrAtomPair[0],
                                               guiNmrAtomPair[1],
                                               spectrum.positiveContourColour,
                                               2.0, displacement=displacement,
                                               peak=peak, lineList=lineList)

                # print('>>>_group', group)
                # for child in group.childItems():
                #     print('   >>>_child', child)

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
                                                   peak=peak, lineList=peaklineList)

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
                                                   peak=peak, lineList=peaklineList)

                else:
                    if nmrAtomPair[0].nmrResidue.nmrChain is nmrResidue.nmrChain:
                        group = self.guiNmrResidues[nmrAtomPair[0].nmrResidue]
                        displacement = guiNmrAtomPair[0].getConnectedList(guiNmrAtomPair[1])
                        self._addConnectingLineToGroup(group,
                                                       guiNmrAtomPair[0],
                                                       guiNmrAtomPair[1],
                                                       spectrum.positiveContourColour,
                                                       2.0, displacement=displacement,
                                                       peak=peak, lineList=peaklineList)

                    elif nmrAtomPair[1].nmrResidue.nmrChain is nmrResidue.nmrChain:
                        group = self.guiNmrResidues[nmrAtomPair[1].nmrResidue]
                        displacement = guiNmrAtomPair[1].getConnectedList(guiNmrAtomPair[0])
                        self._addConnectingLineToGroup(group,
                                                       guiNmrAtomPair[1],
                                                       guiNmrAtomPair[0],
                                                       spectrum.positiveContourColour,
                                                       2.0, displacement=displacement,
                                                       peak=peak, lineList=peaklineList)
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
        newGuiResidueGroup.crossChainResidue = nmrResidueCon0       # self.guiNmrResidues[nmrResidueCon0]

        # self.ghostList.append((nmrResidueCon1.pid, nmrResidueCon0.pid, atoms))
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
        """
        Creates and raises a context menu enabling items to be disconnected
        """
        from ccpn.ui.gui.widgets.Menu import Menu
        from functools import partial

        cursor = QtGui.QCursor()
        contextMenu = Menu('', self, isFloatWidget=True)
        # self.scene.update()
        pressed = self.scene.mouseGrabberItem()

        if isinstance(object, AssignmentLine):  # self.selectedLine:
            # print('>>>', object._peak, pressed._peak)
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
