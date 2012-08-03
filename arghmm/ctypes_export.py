"""

   Functions for automating the export of c functions from a library using
   ctypes.

"""

import sys
import os
from ctypes import *


def c_list(c_type, lst):
    """Make a C array from a list"""
    list_type = c_type * len(lst)
    return list_type(* lst)

def c_matrix(c_type, mat):
    """Make a C matrix from a list of lists (mat)"""

    row_type = c_type * len(mat[0])
    mat_type = POINTER(c_type) * len(mat)
    mat = mat_type(* [row_type(* row) for row in mat])
    return cast(mat, POINTER(POINTER(c_type)))

class c_out (object):
    """
    This wrapper object specifies that an argument should be used for output
    """
    def __init__(self, c_type):
        self.c_type = c_type


# additional primary C types
c_bool_p = POINTER(c_bool)
c_bool_p_p = POINTER(POINTER(c_bool))
c_int_p = POINTER(c_int)
c_int_p_p = POINTER(POINTER(c_int))
c_float_p = POINTER(c_float)
c_float_p_p = POINTER(POINTER(c_float))
c_double_p = POINTER(c_double)
c_double_p_p = POINTER(POINTER(c_double))
c_char_p_p = POINTER(c_char_p)

# basic auto conversion list and matrix types
c_bool_list = (c_int_p, lambda x: c_list(c_bool, x))
c_bool_matrix = (c_int_p_p, lambda x: c_matrix(c_bool, x))
c_char_list = (c_char_p, lambda x: c_list(c_char, x))
c_char_matrix = (c_char_p_p, lambda x: c_matrix(c_char, x))
c_int_list = (c_int_p, lambda x: c_list(c_int, x))
c_int_matrix = (c_int_p_p, lambda x: c_matrix(c_int, x))
c_float_list = (c_float_p, lambda x: c_list(c_float, x))
c_float_matrix = (c_float_p_p, lambda x: c_matrix(c_float, x))
c_double_list = (c_double_p, lambda x: c_list(c_double, x))
c_double_matrix = (c_double_p_p, lambda x: c_matrix(c_double, x))
c_char_p_list = (c_char_p_p, lambda x: c_list(c_char_p, x))


class Exporter (object):

    def __init__(self, env):
        self._env = env


    def export(self, lib, funcname, return_type, prototypes, newname=None):
        """
        Exports a C function with documentation

        lib         -- module to export from (returned from ctypes)
        funcname    -- function to export
        return_type -- return type of function
        prototypes  -- a list defining the function prototype
                       e.g. [type1, name1, type2, name2, type3, name3, ...]
        newname     -- if given, the name of the function after export
                       (default: funcname)
        """

        converts = []

        if newname is None:
            newname = funcname

        # get c arguments
        arg_types = prototypes[0::2]
        for i, arg in enumerate(arg_types):
            if isinstance(arg, c_out):
                # this a output argument
                arg = arg.c_type
            
            if isinstance(arg, tuple):
                # use special conversion function
                arg_types[i] = arg[0]
                converts.append(arg[1])
            else:
                # use type as conversion function
                converts.append(lambda x: x)

        # tell ctypes the prototype of the function
        cfunc = lib.__getattr__(funcname)
        cfunc.restype = return_type
        cfunc.argtypes = arg_types


        # wrapper function for exported function
        def wrapper(*args):

            # record array sizes for outputs
            sizes = {}
            for i, argtype in enumerate(prototypes[0::2]):
                if isinstance(argtype, c_out):
                    sizes[i] = len(args[i])

            # convert arguments to c types
            cargs = [f(a) for f, a in zip(converts, args)]

            # call c function
            ret = cfunc(*cargs)        

            # pass back arguments
            # used for pointers to arrays
            for i, size in sizes.iteritems():
                args[i][:] = cargs[i][:sizes[i]]

            return ret

        # install wrapper function into environment
        self._env[newname] = wrapper

        # set documentation
        args_doc = prototypes[1::2]
        self._env[newname].__doc__ = "%s(%s)" % (funcname, ",".join(args_doc))


def load_library(path, lib):
    try:
        # use library from source path
        libdir = os.path.join(os.path.dirname(__file__), *path)
        return cdll.LoadLibrary(os.path.join(libdir, lib))
    except Exception, e:
        # search for lib in library path
        try:
            return cdll.LoadLibrary(lib)
        except Exception, e2:
            print >>sys.stderr, e
            print >>sys.stderr, e2
            return None

