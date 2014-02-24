# Python initialization
import sys, os

# Name of the executable we're running
executable = sys.argv[0]

# Dereference links
while os.path.islink(executable):
    executable = os.readlink(executable)

# Pyraf-specific startup
if os.path.split(executable)[1] == "pyraf":
    from pyraf import iraf
    sys.path.append("/Users/scenko/python/")

# Clean up the namespace
del executable
