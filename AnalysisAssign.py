"""
AnalysisAssign Program

Please cite:  Skinner et al, J Biomol NMR (2016) 66:111–124; DOI 10.1007/s10858-016-0060-y 

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2017"
__credits__ = ("Wayne Boucher, Ed Brooksbank, Rasmus H Fogh, Luca Mureddu, Timothy J Ragan"
               "Simon P Skinner & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for licence text")
__reference__ = ("For publications, please use reference from http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2017-04-07 11:40:21 +0100 (Fri, April 07, 2017) $"
__version__ = "$Revision: 3.0.b1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"

__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
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
    if not self.project.peakLists:
      self.project._logger.warn('No peaklists in project. Cannot assign peaklists.')
      MessageDialog.showWarning('No peaklists in project.', 'Cannot assign peaklists.')
    else:
      from ccpn.ui.gui.popups.SetupNmrResiduesPopup import SetupNmrResiduesPopup
      popup = SetupNmrResiduesPopup(self.ui.mainWindow, self.project)
      popup.exec_()

  def showPickAndAssignModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Pick and Assign module.
    """
    from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule

    mainWindow = self.ui.mainWindow
    self.pickAndAssignModule = PickAndAssignModule(mainWindow=mainWindow, name='Pick and Assign')
    mainWindow.moduleArea.addModule(self.pickAndAssignModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showPickAndAssignModule()")
    self.project._logger.info("application.showPickAndAssignModule()")
    return self.pickAndAssignModule


  def showBackboneAssignmentModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Backbone Assignment module.
    """
    from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule

    mainWindow = self.ui.mainWindow
    self.backboneModule = BackboneAssignmentModule(mainWindow=mainWindow)
    mainWindow.moduleArea.addModule(self.backboneModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showBackboneAssignmentModule()")
    self.project._logger.info("application.showBackboneAssignmentModule()")
    return self.backboneModule


  def showSidechainAssignmentModule(self, position:str='bottom', relativeTo:CcpnModule=None):
    """
    Displays Backbone Assignment module.
    """
    MessageDialog.showWarning('PickandAssignModule',
                              'SideChainAssignmentModule.py (16-20)\n'
                              'Not implemented yet')
    return

    from ccpn.AnalysisAssign.modules.SideChainAssignmentModule import SideChainAssignmentModule

    if hasattr(self, 'sidechainAssignmentModule'):
      return

    mainWindow = self.ui.mainWindow
    self.sidechainAssignmentModule = SideChainAssignmentModule(mainWindow=mainWindow)   # ejb self, self.project)
    mainWindow.moduleArea.addModule(self.sidechainAssignmentModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showSidechainAssignmentModule()")
    self.project._logger.info("application.showSidechainAssignmentModule()")

    return self.sidechainAssignmentModule


  def showPeakAssigner(self, position='bottom', relativeTo=None):
    """Displays peak assignment module."""
    from ccpn.ui.gui.modules.PeakAssigner import PeakAssigner

    mainWindow = self.ui.mainWindow
    self.assignmentModule = PeakAssigner(mainWindow=mainWindow)
    mainWindow.moduleArea.addModule(self.assignmentModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentModule()")
    self.project._logger.info("application.showAssignmentModule()")


  def showResidueInformation(self, position: str='bottom', relativeTo:CcpnModule=None):
    """Displays Residue Information module."""
    from ccpn.ui.gui.modules.ResidueInformation import ResidueInformation
    if not self.project.residues:
      self.project._logger.warn('No Residues in project. Residue Information Module requires Residues in the project to launch.')
      MessageDialog.showWarning('No Residues in project.',
                                'Residue Information Module requires Residues in the project to launch.')
      return

    mainWindow = self.ui.mainWindow
    mainWindow.moduleArea.addModule(ResidueInformation(self, self.project), position=position,
                              relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showResidueInformation()")
    self.project._logger.info("application.showResidueInformation()")


  def showAssignmentInspectorModule(self, nmrAtom=None, position: str='bottom', relativeTo:CcpnModule=None):
    from ccpn.AnalysisAssign.modules.AssignmentInspectorModule import AssignmentInspectorModule

    mainWindow = self.ui.mainWindow
    self.assignmentInspectorModule = AssignmentInspectorModule(mainWindow=mainWindow)
    mainWindow.moduleArea.addModule(self.assignmentInspectorModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentInspectorModule()")
    self.project._logger.info("application.showAssignmentInspectorModule()")
