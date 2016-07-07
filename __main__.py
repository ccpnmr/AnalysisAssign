__author__ = 'TJ'

from ccpn.core.lib.Version import applicationVersion
from ccpn.framework.Framework import defineProgramArguments
from ccpn.AnalysisAssign.AnalysisAssign import Assign as Application


if __name__ == '__main__':

  # argument parser
  parser = defineProgramArguments()

  # add any additional commandline argument here
  commandLineArguments = parser.parse_args()

  application = Application('Assign', applicationVersion, commandLineArguments)
  application.start()
