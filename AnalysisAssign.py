"""
AnalysisAssign Program
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
__dateModified__ = "$dateModified: 2017-07-07 16:32:20 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b3 $"
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
    menuSpec = ('Assign', [("Setup NmrResidues", self.showSetupNmrResiduesPopup, [('shortcut', 'sn')]),
                           ("Pick and Assign", self.showPickAndAssignModule, [('shortcut', 'pa')]),
                           (),
                           ("Backbone Assignment", self.showBackboneAssignmentModule, [('shortcut', 'bb')]),
                           ("Sidechain Assignment", self.showSidechainAssignmentModule, [('shortcut', 'sc'), ('enabled', False)]),
                           (),
                           ("Peak Assigner", self.showPeakAssigner, [('shortcut', 'aa')]),
                           ("Assignment Inspector", self.showAssignmentInspectorModule, [('shortcut', 'ai')]),
                           ("Residue Information", self.showResidueInformation, [('shortcut', 'ri')]),
                          ])
    self.addApplicationMenuSpec(menuSpec)

    viewMenuItems = [ ("Sequence Graph", self.showSequenceGraph, [('shortcut', 'sg')]),
                      ("Atom Selector", self.showAtomSelector, [('shortcut', 'as')]),
                      ()
                    ]
    self.addApplicationMenuItems('View', viewMenuItems, position=9)

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

  def showPickAndAssignModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Pick and Assign module.
    """
    from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule

    mainWindow = self.ui.mainWindow
    #FIXME:ED - crashes sometimes opening a module
    if not relativeTo:
      relativeTo = mainWindow.moduleArea    # ejb
    self.pickAndAssignModule = PickAndAssignModule(mainWindow=mainWindow)
    mainWindow.moduleArea.addModule(self.pickAndAssignModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showPickAndAssignModule()")
    getLogger().info("application.showPickAndAssignModule()")
    return self.pickAndAssignModule


  def showBackboneAssignmentModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Backbone Assignment module.
    """
    from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule

    mainWindow = self.ui.mainWindow
    #FIXME:ED - crashes sometimes opening a module
    if not relativeTo:
      relativeTo = mainWindow.moduleArea    # ejb
    self.backboneModule = BackboneAssignmentModule(mainWindow=mainWindow)
    mainWindow.moduleArea.addModule(self.backboneModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showBackboneAssignmentModule()")
    getLogger().info("application.showBackboneAssignmentModule()")
    return self.backboneModule


  def showSidechainAssignmentModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Backbone Assignment module.
    """
    MessageDialog.showWarning('Not implemented',
                              'Sidechain Assignment Module\n'
                              'is not implemented yet')
    return

    from ccpn.AnalysisAssign.modules.SideChainAssignmentModule import SideChainAssignmentModule

    if hasattr(self, 'sidechainAssignmentModule'):
      return

    mainWindow = self.ui.mainWindow
    #FIXME:ED - crashes sometimes opening a module
    if not relativeTo:
      relativeTo = mainWindow.moduleArea    # ejb
    self.sidechainAssignmentModule = SideChainAssignmentModule(mainWindow=mainWindow)   # ejb self, self.project)
    mainWindow.moduleArea.addModule(self.sidechainAssignmentModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showSidechainAssignmentModule()")
    getLogger().info("application.showSidechainAssignmentModule()")

    return self.sidechainAssignmentModule


  def showPeakAssigner(self, position='bottom', relativeTo=None):
    """Displays peak assignment module."""
    from ccpn.AnalysisAssign.modules.PeakAssigner import PeakAssigner

    mainWindow = self.ui.mainWindow
    #FIXME:ED - crashes sometimes opening a module
    if not relativeTo:
      relativeTo = mainWindow.moduleArea    # ejb
    self.assignmentModule = PeakAssigner(mainWindow=mainWindow)
    mainWindow.moduleArea.addModule(self.assignmentModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentModule()")
    getLogger().info("application.showAssignmentModule()")


  def showResidueInformation(self, position: str='bottom', relativeTo:CcpnModule=None):
    """Displays Residue Information module."""
    from ccpn.ui.gui.modules.ResidueInformation import ResidueInformation
    if not self.project.residues:
      getLogger().warning('No Residues in project. Residue Information Module requires Residues in the project to launch.')
      MessageDialog.showWarning('No Residues in project.',
                                'Residue Information Module requires Residues in the project to launch.')
      return

    mainWindow = self.ui.mainWindow
    #FIXME:ED - crashes sometimes opening a module
    if not relativeTo:
      relativeTo = mainWindow.moduleArea    # ejb
    self.residueModule = ResidueInformation(mainWindow=mainWindow)
    mainWindow.moduleArea.addModule(self.residueModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showResidueInformation()")
    getLogger().info("application.showResidueInformation()")


  def showAssignmentInspectorModule(self, nmrAtom=None, position: str='bottom', relativeTo:CcpnModule=None):
    from ccpn.AnalysisAssign.modules.AssignmentInspectorModule import AssignmentInspectorModule

    mainWindow = self.ui.mainWindow
    #FIXME:ED - crashes sometimes opening a module
    if not relativeTo:
      relativeTo = mainWindow.moduleArea    # ejb
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
    # FIXME:ED - sometimes crashes
    if not relativeTo:
      relativeTo = mainWindow.moduleArea  # ejb
    self.sequenceGraphModule = SequenceGraphModule(mainWindow=mainWindow, nmrChain=nmrChain)
    mainWindow.moduleArea.addModule(self.sequenceGraphModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showSequenceGraph()")
    getLogger().info("application.showSequenceGraph()")
    return self.sequenceGraphModule


  def showAtomSelector(self, position: str = 'bottom', relativeTo: CcpnModule = None, nmrAtom=None):
    """Displays Atom Selector."""
    from ccpn.AnalysisAssign.modules.AtomSelector import AtomSelectorModule

    mainWindow = self.ui.mainWindow
    # FIXME:ED - sometimes crashes
    if not relativeTo:
      relativeTo = mainWindow.moduleArea  # ejb
    self.atomSelectorModule = AtomSelectorModule(mainWindow=mainWindow, nmrAtom=nmrAtom)
    mainWindow.moduleArea.addModule(self.atomSelectorModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showAtomSelector()")
    getLogger().info("application.showAtomSelector()")
    return self.atomSelectorModule

