import glob
import sys

import LineRegress


def getFiles(pattern):
    files =glob.glob(pattern)
    return files


if __name__ == "__main__":
    neFile = sys.argv[1]
    configFile = None
    if len(sys.argv) ==3:
        configFile = sys.argv[2]
    LineRegress.neGrapher(neFile,configFile)
