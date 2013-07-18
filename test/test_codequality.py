
import os
import sys
import subprocess


def run(cmd):
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out = pipe.stdout.read()
    retcode = pipe.wait()
    if retcode != 0:
        print out
    return retcode



def run_pyflakes(filenames, key=lambda line: True):

    cmd = " ".join(["pyflakes"] + filenames)
    print cmd
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lines = [line for line in pipe.stdout if key(line)]
    retcode = pipe.wait()
    return lines


def pyflakes_filter(line):

    # Allow arghmm to be imported without being used since it
    # imports arghmm.deps
    if "'arghmm' imported but unused" in line:
        return False

    # allow summon.core to use 'import *'
    if 'from summon.core import *' in line:
        return False

    return True


def test_import_arghmm():
    """
    Ensure arghmm library can be imported.
    """

    assert os.system("PYTHONPATH= python -c 'import arghmm'") == 0
    assert os.system("PYTHONPATH= python -c 'import arghmm.popsize'") == 0


def get_python_scripts(path):
    """
    Return the python scripts in a directory
    """
    filenames = sorted(os.listdir("bin"))
    for filename in filenames:
        filename = os.path.join("bin", filename)
        if filename.endswith(".py"):
            yield filename
            continue

        with open(filename) as infile:
            if "python" in infile.readline():
                yield filename


def test_bin():
    """
    Ensure all scripts can run without external PYTHONPATH.
    """

    filenames = os.listdir("bin")
    errors = []
    for filename in filenames:
        filename = os.path.join("bin", filename)
        if not os.path.isfile(filename):
            continue

        cmd = (("export PYTHONPATH=\n(head -n1 %s | grep python -q -v) || "
                "%s --help < /dev/null 2>&1 ") % (filename, filename))
        if run(cmd) != 0:
            print "ERROR>", filename
            errors.append(filename)

    if len(errors) > 0:
        print "scripts with erroneous imports:"
        print "\n".join(errors)
        raise Exception()


def test_pyflakes():
    """
    Run pyflakes on python code base.
    """

    filenames = list(get_python_scripts("bin"))
    lines = run_pyflakes(filenames, key=pyflakes_filter)

    if len(lines) > 0:
        print "pyflakes errors:"
        print "".join(lines)
        raise Exception()
