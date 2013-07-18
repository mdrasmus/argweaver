
import sys, os

def load_deps(dirname="deps"):
    sys.path.append(os.path.realpath(
        os.path.join(os.path.dirname(__file__), dirname)))

