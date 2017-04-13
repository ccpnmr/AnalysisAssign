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
__author__ = "$Author: Geerten Vuister $"
__date__ = "$Date: 2017-04-13 12:24:47 +0100 (Thu, April 13, 2017) $"

#=========================================================================================
# Start of code
#=========================================================================================

# from ccpn.ui.gui.AppBase import AppBase, defineProgramArguments
# from ccpn.ui.gui.lib.Window import MODULE_DICT
# from ccpn.ui.gui.modules import GuiStrip

# from ccpn.framework.lib.SvnRevision import applicationVersion

from ccpn.framework.Framework import Framework
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets import MessageDialog


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
                           ("Sidechain Assignment", self.showSidechainAssignmentModule, [('shortcut', 'sc')]),
                           (),
                           ("Peak Assigner", self.showPeakAssigner, [('shortcut', 'aa')]),
                           ("Assignment Inspector", self.showAssignmentInspectorModule, [('shortcut', 'ai')]),
                           ("Residue Information", self.showResidueInformation, [('shortcut', 'ri')]),
                          ])
    self.addApplicationMenuSpec(menuSpec)

  # overrides superclass
  def _closeExtraWindows(self):

    # remove links to modules when closing them
    for attr in ('sequenceGraph', 'backboneModule', 'sidechainAssignmentModule'):
      if hasattr(self, attr):
        delattr(self, attr)

    Framework._closeExtraWindows(self)

  def showSetupNmrResiduesPopup(self):
    from ccpn.ui.gui.popups.SetupNmrResiduesPopup import SetupNmrResiduesPopup
    popup = SetupNmrResiduesPopup(self.ui.mainWindow, self.project)
    popup.exec_()


  def showSequenceGraph(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays sequence graph at the bottom of the screen, relative to another module if nextTo is specified.
    """
    from ccpn.AnalysisAssign.modules.SequenceGraph import SequenceGraph

    if hasattr(self, 'sequenceGraph'):
      return
    self.sequenceGraph = SequenceGraph(self, project=self.project)
    #if hasattr(self, 'backboneModule'):
    #  self.backboneModule._connectSequenceGraph(self.sequenceGraph)

    if relativeTo is not None:
      self.ui.mainWindow.moduleArea.addModule(self.sequenceGraph, position=position, relativeTo=relativeTo)
    else:
      self.ui.mainWindow.moduleArea.addModule(self.sequenceGraph, position=position)
    self.ui.mainWindow.pythonConsole.writeConsoleCommand("application.showSequenceGraph()")
    self.project._logger.info("application.showSequenceGraph()")
    return self.sequenceGraph


  def showPickAndAssignModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule

    """Displays Pick and Assign module."""
    mainWindow = self.ui.mainWindow
    self.pickAndAssignModule = PickAndAssignModule(mainWindow.moduleArea, self.project)
    mainWindow.moduleArea.addModule(self.pickAndAssignModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showPickAndAssignModule()")
    self.project._logger.info("application.showPickAndAssignModule()")
    return self.pickAndAssignModule


  def showBackboneAssignmentModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Backbone Assignment module.
    """
    from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule

    #if hasattr(self, 'backboneModule'):
    #  return

    self.backboneModule = BackboneAssignmentModule(self)

    mainWindow = self.ui.mainWindow
    mainWindow.moduleArea.addModule(self.backboneModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showBackboneAssignmentModule()")
    self.project._logger.info("application.showBackboneAssignmentModule()")
    #if hasattr(self, 'sequenceGraph'):
    #  self.backboneModule._connectSequenceGraph(self.sequenceGraph)

    return self.backboneModule


  def showSidechainAssignmentModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Backbone Assignment module.
    """
    from ccpn.AnalysisAssign.modules.SideChainAssignmentModule import SideChainAssignmentModule

    if hasattr(self, 'sidechainAssignmentModule'):
      return

    self.sidechainAssignmentModule = SideChainAssignmentModule(self, self.project)

    mainWindow = self.ui.mainWindow
    mainWindow.moduleArea.addModule(self.sidechainAssignmentModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showSidechainAssignmentModule()")
    self.project._logger.info("application.showSidechainAssignmentModule()")

    return self.sidechainAssignmentModule


  def showPeakAssigner(self, position='bottom', relativeTo=None):
    """Displays peak assignment module."""
    from ccpn.ui.gui.modules.PeakAssigner import PeakAssigner

    mainWindow = self.ui.mainWindow
    self.assignmentModule = PeakAssigner(mainWindow)
    mainWindow.moduleArea.addModule(self.assignmentModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentModule()")
    self.project._logger.info("application.showAssignmentModule()")


  def showResidueInformation(self, position: str='bottom', relativeTo:CcpnModule=None):
    """Displays Residue Information module."""
    from ccpn.ui.gui.modules.ResidueInformation import ResidueInformation
    if not self.project.residues:
      self.project._logger.warn('No Residues in  project. Residue Information Module requires Residues in the project to launch.')
      MessageDialog.showWarning('No Residues in  project.',
                                'Residue Information Module requires Residues in the project to launch.')
      return

    mainWindow = self.ui.mainWindow
    mainWindow.moduleArea.addModule(ResidueInformation(self, self.project), position=position,
                              relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showResidueInformation()")
    self.project._logger.info("application.showResidueInformation()")


  def showAssignmentInspectorModule(self, nmrAtom=None, position: str='bottom', relativeTo:CcpnModule=None):
    from ccpn.AnalysisAssign.modules.AssignmentInspectorModule import AssignmentInspectorModule

    if not nmrAtom and len(self.project.nmrAtoms) == 0:
      self.project._logger.warn('No NmrAtom selected or defined. The Modify Assignments Module requires an NmrAtom to be functional')
      #MessageDialog.showWarning('No NmrAtom selected or defined.',
      #                          'The Modify Assignments Module requires an NmrAtom to be functional')
      #return

    mainWindow = self.ui.mainWindow
    self.assignmentInspectorModule = AssignmentInspectorModule(mainWindow.moduleArea)
    mainWindow.moduleArea.addModule(self.assignmentInspectorModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentInspectorModule()")
    self.project._logger.info("application.showAssignmentInspectorModule()")
