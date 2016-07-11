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
__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Revision$"

#=========================================================================================
# Start of code
#=========================================================================================

# from ccpn.ui.gui.AppBase import AppBase, defineProgramArguments
# from ccpn.ui.gui.lib.Window import MODULE_DICT
# from ccpn.ui.gui.modules import GuiStrip

# from ccpn.framework.lib.SvnRevision import applicationVersion

from ccpn.framework.Framework import Framework
# from ccpn.ui.gui.modules import GuiStripNd
# from ccpn.ui.gui.modules import GuiSpectrumDisplay
# from ccpn.ui.gui.modules import GuiStripDisplayNd
from ccpn.ui.gui.widgets import MessageDialog
from ccpn.ui.gui.widgets.Module import CcpnModule


# applicationName = 'AnalysisAssign'

# class Assign(AppBase):
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
                           # ("Sidechain Assignment", self.showSetupNmrResiduesPopup, 'sc'),
                           (),
                           ("Peak Assigner", self.showPeakAssigner, [('shortcut', 'aa')]),
                           ("Modify Assignments", self.showModifyAssignmentModule, [('shortcut', 'ma')]),
                           ("Residue Information", self.showResidueInformation, [('shortcut', 'ri')]),
                          ])
    self.addApplicationMenuSpec(menuSpec)


  # def initGraphics(self):
  #   """Set up graphics system after loading"""
  #
  #   # Initialise strips
  #   project = self.project
  #   for strip in project.strips:
  #     GuiStrip._setupGuiStrip(project, strip._wrappedData)
  #
  #     # if isinstance(strip, GuiStripNd) and not strip.haveSetupZWidgets:
  #     #   strip.setZWidgets()
  #
  #   # Initialise Rulers
  #   for task in project.tasks:
  #     for apiMark in task._wrappedData.sortedMarks():
  #       for apiRuler in apiMark.sortedRulers():
  #         GuiStrip._rulerCreated(project, apiRuler)
  #
  #   # Initialise SpectrumViews
  #   for spectrumDisplay in project.spectrumDisplays:
  #     for strip in spectrumDisplay.strips:
  #       for spectrumView in strip.spectrumViews:
  #         spectrumView._createdSpectrumView()
  #         for peakList in spectrumView.spectrum.peakLists:
  #           strip.showPeaks(peakList)
  #
  #   self.initLayout()
  #
  # def initLayout(self):
  #   """
  #   Restore layout of modules from previous save after graphics have been set up.
  #   """
  #   import yaml, os
  #   if os.path.exists(os.path.join(self.project.path, 'layouts', 'layout.yaml')):
  #     with open(os.path.join(self.project.path, 'layouts', 'layout.yaml')) as f:
  #       layout = yaml.load(f)
  #       typ, contents, state = layout['main']
  #
  #       # TODO: When UI has a main window, change the call below (then move the whole function!)
  #       containers, modules = self.ui.mainWindow.moduleArea.findAll()
  #       flatten = lambda *n: (e for a in n
  #       for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))
  #       flatContents = list(flatten(contents))
  #       for item in flatContents:
  #         if item in list(MODULE_DICT.keys()):
  #           obj = modules.get(item)
  #           if not obj:
  #            func = getattr(self, MODULE_DICT[item])
  #            func()
  #       for s in layout['float']:
  #         typ, contents, state = s[0]['main']
  #         containers, modules = self.ui.mainWindow.moduleArea.findAll()
  #         for item in contents:
  #           if item[0] == 'dock':
  #             obj = modules.get(item[1])
  #             if not obj:
  #               func = getattr(self, MODULE_DICT[item[1]])
  #               func()
  #       self.ui.mainWindow.moduleArea.restoreState(layout)


  def showSetupNmrResiduesPopup(self):
    from ccpn.ui.gui.popups.SetupNmrResiduesPopup import SetupNmrResiduesPopup
    popup = SetupNmrResiduesPopup(self.ui.mainWindow, self.project)
    popup.exec_()


  def showPickAndAssignModule(self, position:str= 'bottom', relativeTo:CcpnModule=None):
    from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule

    """Displays Pick and Assign module."""
    mainWindow = self.ui.mainWindow
    self.pickAndAssignModule = PickAndAssignModule(mainWindow.moduleArea, self.project)
    mainWindow.moduleArea.addModule(self.pickAndAssignModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showPickAndAssignModule()")
    self.project._logger.info("application.showPickAndAssignModule()")
    return self.pickAndAssignModule


  def showBackboneAssignmentModule(self, position:str= 'bottom', relativeTo:CcpnModule=None):
    """
    Displays Backbone Assignment module.
    """
    from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule

    if hasattr(self, 'backboneModule'):
      return

    self.backboneModule = BackboneAssignmentModule(self, self.project)

    mainWindow = self.ui.mainWindow
    mainWindow.moduleArea.addModule(self.backboneModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showBackboneAssignmentModule()")
    self.project._logger.info("application.showBackboneAssignmentModule()")
    if hasattr(self, 'assigner'):
      self.backboneModule._connectSequenceGraph(self.assigner)

    return self.backboneModule


  def showPeakAssigner(self, position='bottom', relativeTo=None):
    """Displays assignment module."""
    from ccpn.ui.gui.modules.PeakAssigner import PeakAssigner

    mainWindow = self.ui.mainWindow
    self.assignmentModule = PeakAssigner(self, self.project, self.current.peaks)
    mainWindow.moduleArea.addModule(self.assignmentModule, position=position, relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showAssignmentModule()")
    self.project._logger.info("application.showAssignmentModule()")


  def showResidueInformation(self, position: str='bottom', relativeTo:CcpnModule=None):
    """Displays Residue Information module."""
    from ccpn.ui.gui.modules.ResidueInformation import ResidueInformation

    mainWindow = self.ui.mainWindow
    mainWindow.moduleArea.addModule(ResidueInformation(self, self.project), position=position,
                              relativeTo=relativeTo)
    mainWindow.pythonConsole.writeConsoleCommand("application.showResidueInformation()")
    self.project._logger.info("application.showResidueInformation()")


  def showModifyAssignmentModule(self, nmrAtom=None, position: str='bottom', relativeTo:CcpnModule=None):
    from ccpn.AnalysisAssign.modules.ModifyAssignmentModule import ModifyAssignmentModule
    if not nmrAtom:
      self.project._logger.warn('No NmrAtom selected. The Modify Assignments Module requires an NmrAtom to launch')
      MessageDialog.showWarning('No NmrAtom selected.', 'The Modify Assignments Module requires an NmrAtom to launch',
                                colourScheme=self.preferences.general.colourScheme)
      return
    mainWindow = self.ui.mainWindow
    self.modifyAssignmentsModule = ModifyAssignmentModule(mainWindow.moduleArea, self.project, nmrAtom=nmrAtom)
    mainWindow.moduleArea.addModule(self.modifyAssignmentsModule, position=position,
                              relativeTo=relativeTo)
