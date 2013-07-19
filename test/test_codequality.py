
import os
import subprocess


def run(cmd, shell=True):
    """
    Run a command and check the return code.
    """
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=shell)
    out = pipe.stdout.read()
    retcode = pipe.wait()
    if retcode != 0:
        print out
    return retcode


def run_pyflakes(filenames, key=lambda line: True):
    """
    Run pyflakes and return all errors.
    """
    cmd = " ".join(["pyflakes"] + filenames)
    print cmd
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lines = [line for line in pipe.stdout if key(line)]
    pipe.wait()
    return lines


def run_pep8(filenames):
    """
    Run pep8 and return all errors.
    """
    cmd = " ".join(["pep8"] + filenames)
    print cmd
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lines = [line for line in pipe.stdout]
    pipe.wait()
    return lines


def pyflakes_filter(line):
    """
    Standard filter for pyflakes.
    """

    # Allow arghmm to be imported without being used since it
    # imports arghmm.deps
    if "'arghmm' imported but unused" in line:
        return False

    # allow summon.core to use 'import *'
    if 'from summon.core import *' in line:
        return False

    if 'arghmm/debug.py' in line:
        return False

    if 'arghmm/bottle.py' in line:
        return False

    if 'arghmm/arghmmc.py' in line:
        return False

    return True


def test_import_arghmm():
    """
    Ensure arghmm library can be imported.
    """

    assert os.system("PYTHONPATH= python -c 'import arghmm'") == 0
    assert os.system("PYTHONPATH= python -c 'import arghmm.popsize'") == 0


def get_python_scripts(*paths):
    """
    Return the python scripts in a directory
    """
    filenames = []
    for path in paths:
        files = sorted(os.listdir(path))
        filenames.extend(os.path.join(path, filename) for filename in files)
    for filename in filenames:
        # Skip directories
        if not os.path.isfile(filename):
            continue

        # Return filenames ending in *.py
        if filename.endswith(".py"):
            yield filename
            continue

        # Return filenames containing 'python' in the first line
        with open(filename) as infile:
            line = infile.readline()
            if "python" in line and "python-i" not in line:
                yield filename


def test_bin():
    """
    Ensure all scripts can run without external PYTHONPATH.
    """
    filenames = get_python_scripts("bin")
    errors = []
    for filename in filenames:
        cmd = "export PYTHONPATH=\n %s --help < /dev/null" % filename
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
    filenames = list(get_python_scripts("bin", "arghmm", "test"))
    lines = run_pyflakes(filenames, key=pyflakes_filter)

    if len(lines) > 0:
        print "pyflakes errors:"
        print "".join(lines)
        raise Exception()


def test_pep8():
    """
    Ensure pep8 compliance on python code base.
    """
    filenames = list(get_python_scripts("bin", "test"))
    lines = run_pep8(filenames)

    if len(lines) > 0:
        print "pep8 errors:"
        print "".join(lines)
        raise Exception()
