
import os
import sys
import subprocess

def test_import_arghmm():

    assert os.system("PYTHONPATH= python -c 'import arghmm'") == 0
    assert os.system("PYTHONPATH= python -c 'import arghmm.popsize'") == 0

def test_bin():
    filenames = os.listdir("bin")
    errors = []
    for filename in filenames:
        filename = os.path.join("bin", filename)
        if not os.path.isfile(filename):
            continue

        cmd = (("export PYTHONPATH=\n(head -n1 %s | grep python -q -v) || "
                "%s --help < /dev/null 2>&1 ") % (filename, filename))

        pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        out = pipe.stdout.read()
        if pipe.wait() != 0:
            print "ERROR>", filename
            print out
            errors.append(filename)

    if len(errors) > 0:
        print "scripts with erroneous imports:"
        print "\n".join(errors)
        raise Exception()
