"""
AnalysisAssign Program
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
__dateModified__ = "$dateModified: 2017-07-07 16:32:20 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b5 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

from ccpn.framework.Framework import Framework
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets import MessageDialog
from ccpn.util.Logging import getLogger


class Assign(Framework):
    """Root class for Assign application"""

    def __init__(self, applicationName, applicationVersion, commandLineArguments):
        Framework.__init__(self, applicationName, applicationVersion, commandLineArguments)
        # self.components.add('Assignment')

    def setupMenus(self):
        super().setupMenus()
        menuSpec = ('Assign', [("Set up NmrResidues", self.showSetupNmrResiduesPopup, [('shortcut', 'sn')]),
                               ("Pick and Assign", self.showPickAndAssignModule, [('shortcut', 'pa')]),
                               (),
                               ("Backbone Assignment", self.showBackboneAssignmentModule, [('shortcut', 'bb')]),
                               ("Sidechain Assignment", self.showSidechainAssignmentModule, [('shortcut', 'sc'), ('enabled', False)]),
                               (),
                               ("Peak Assigner", self.showPeakAssigner, [('shortcut', 'aa')]),
                               ("Assignment Inspector", self.showAssignmentInspectorModule, [('shortcut', 'ai')]),
                               # ("Residue Information", self.showResidueInformation, [('shortcut', 'ri')]),
                               ])
        self.addApplicationMenuSpec(menuSpec)

        viewMenuItems = [("Sequence Graph", self.showSequenceGraph, [('shortcut', 'sg')]),
                         ("NmrAtom Assigner", self.showAtomSelector, [('shortcut', 'as')]),
                         ()
                         ]
        self.addApplicationMenuItems('View', viewMenuItems, position=8)

    # overrides superclass
    def _closeExtraWindows(self):

        # remove links to modules when closing them
        for attr in ('sequenceGraph', 'backboneModule', 'sidechainAssignmentModule'):
            if hasattr(self, attr):
                delattr(self, attr)

        Framework._closeExtraWindows(self)

    def showSetupNmrResiduesPopup(self):
        if not self.project.peakLists:
            getLogger().warning('No peaklists in project. Cannot assign peaklists.')
            MessageDialog.showWarning('No peaklists in project.', 'Cannot assign peaklists.')
        else:
            from ccpn.ui.gui.popups.SetupNmrResiduesPopup import SetupNmrResiduesPopup

            popup = SetupNmrResiduesPopup(parent=self.ui.mainWindow, mainWindow=self.ui.mainWindow)
            popup.exec_()

    def showPickAndAssignModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """
        Displays Pick and Assign module.
        """
        from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea  # ejb
        self.pickAndAssignModule = PickAndAssignModule(mainWindow=mainWindow)
        mainWindow.moduleArea.addModule(self.pickAndAssignModule, position=position, relativeTo=relativeTo)
        mainWindow.pythonConsole.writeConsoleCommand("application.showPickAndAssignModule()")
        getLogger().info("application.showPickAndAssignModule()")
        return self.pickAndAssignModule

    def showBackboneAssignmentModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """
        Displays Backbone Assignment module.
        """
        from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea  # ejb
        self.backboneModule = BackboneAssignmentModule(mainWindow=mainWindow)
        mainWindow.moduleArea.addModule(self.backboneModule, position=position, relativeTo=relativeTo)
        mainWindow.pythonConsole.writeConsoleCommand("application.showBackboneAssignmentModule()")
        getLogger().info("application.showBackboneAssignmentModule()")
        return self.backboneModule

    def showSidechainAssignmentModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """
        Displays Backbone Assignment module.
        """
        MessageDialog.showWarning('Not implemented',
                                  'Sidechain Assignment Module\n'
                                  'is not implemented yet')
        return

        # from ccpn.AnalysisAssign.modules.SideChainAssignmentModule import SideChainAssignmentModule
        #
        # if hasattr(self, 'sidechainAssignmentModule'):
        #     return
        #
        # mainWindow = self.ui.mainWindow
        #
        # if not relativeTo:
        #     relativeTo = mainWindow.moduleArea  # ejb
        # self.sidechainAssignmentModule = SideChainAssignmentModule(mainWindow=mainWindow)  # ejb self, self.project)
        # mainWindow.moduleArea.addModule(self.sidechainAssignmentModule, position=position, relativeTo=relativeTo)
        # mainWindow.pythonConsole.writeConsoleCommand("application.showSidechainAssignmentModule()")
        # getLogger().info("application.showSidechainAssignmentModule()")
        #
        # return self.sidechainAssignmentModule

    def showPeakAssigner(self, position='bottom', relativeTo=None):
        """Displays peak assignment module."""
        from ccpn.AnalysisAssign.modules.PeakAssigner import PeakAssigner

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea  # ejb
        self.assignmentModule = PeakAssigner(mainWindow=mainWindow)
        mainWindow.moduleArea.addModule(self.assignmentModule, position=position, relativeTo=relativeTo)
        mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentModule()")
        getLogger().info("application.showAssignmentModule()")

    def showAssignmentInspectorModule(self, nmrAtom=None, position: str = 'bottom', relativeTo: CcpnModule = None):
        from ccpn.AnalysisAssign.modules.AssignmentInspectorModule import AssignmentInspectorModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea  # ejb
        self.assignmentInspectorModule = AssignmentInspectorModule(mainWindow=mainWindow)
        mainWindow.moduleArea.addModule(self.assignmentInspectorModule, position=position, relativeTo=relativeTo)
        mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentInspectorModule()")
        getLogger().info("application.showAssignmentInspectorModule()")

    def showSequenceGraph(self, position: str = 'bottom', relativeTo: CcpnModule = None, nmrChain=None):
        """
        Displays sequence graph at the bottom of the screen, relative to another module if nextTo is specified.
        """
        from ccpn.AnalysisAssign.modules.SequenceGraph import SequenceGraphModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea  # ejb
        self.sequenceGraphModule = SequenceGraphModule(mainWindow=mainWindow, nmrChain=nmrChain)
        mainWindow.moduleArea.addModule(self.sequenceGraphModule, position=position, relativeTo=relativeTo)
        mainWindow.pythonConsole.writeConsoleCommand("application.showSequenceGraph()")
        getLogger().info("application.showSequenceGraph()")
        return self.sequenceGraphModule

    def showAtomSelector(self, position: str = 'bottom', relativeTo: CcpnModule = None, nmrAtom=None):
        """Displays Atom Selector."""
        from ccpn.AnalysisAssign.modules.NmrAtomAssigner import NmrAtomAssignerModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea  # ejb
        self.nmrAtomAssigner = NmrAtomAssignerModule(mainWindow=mainWindow, nmrAtom=nmrAtom)
        mainWindow.moduleArea.addModule(self.nmrAtomAssigner, position=position, relativeTo=relativeTo)
        mainWindow.pythonConsole.writeConsoleCommand("application.showAtomSelector()")
        getLogger().info("application.showAtomSelector()")
        return self.nmrAtomAssigner
