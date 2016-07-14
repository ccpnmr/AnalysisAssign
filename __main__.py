__author__ = 'TJ'

from ccpn.core.lib.Version import applicationVersion
from ccpn.framework import Framework
from ccpn.AnalysisAssign.AnalysisAssign import Assign as Application


if __name__ == '__main__':

  # argument parser
  parser = Framework.defineProgramArguments()

  # add any additional commandline argument here
  commandLineArguments = parser.parse_args()

  application = Application('AnalysisAssign', applicationVersion, commandLineArguments)
  Framework._getApplication = lambda: application
  application.start()
