__author__ = 'TJ'

from ccpn.framework import Framework
from ccpn.AnalysisAssign.AnalysisAssign import Assign as Application
from ccpn.framework.Version import applicationVersion

if __name__ == '__main__':
  from ccpn.util.GitTools import getAllRepositoriesGitCommit
  applicationVersion = 'development: {AnalysisAssign:.8s}'.format(**getAllRepositoriesGitCommit())

  # argument parser
  parser = Framework.defineProgramArguments()

  # add any additional commandline argument here
  commandLineArguments = parser.parse_args()

  application = Application('AnalysisAssign', applicationVersion, commandLineArguments)
  Framework._getApplication = lambda: application
  application.start()
