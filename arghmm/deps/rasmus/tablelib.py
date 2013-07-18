"""
tablelib.py

Portable Tabular Format (PTF)

Implements and standardizes Manolis style tab-delimited table file format.


--Example----------------------------------------------------
##types:string int
name num
mike 23
alex 12
matt 7
-------------------------------------------------------------

File is tab delimited.

Directives are on a single line and begin with two hashes '##'
No space after colon is allowed.


Table can also handle custom types.  Custom types must do the following

 1. default value: 
   default = mytype()  
   returns default value

 2. convert from string
   val = mytype(string)  
   converts from string to custom type

 3. convert to string
   string = str(val)
   converts val of type 'mytype' to a string
   TODO: I could change this interface...
   I could just use  mytype.__str__(val)

 4. type inference (optional)
   type(val)
   returns instance of 'mytype'
   TODO: I could not require this (only map() really needs it and __init__())


"""

# python libs
import copy
import StringIO
import sys
import os
import itertools

# rasmus libs
from rasmus import util


# table directives
DIR_TYPES    = 1

# a special unique null type (more 'null' than None)
NULL = object()


class TableException (Exception):
    """Exception class for Table"""
    def __init__(self, errmsg, filename=None, lineno=None):
        msg = ""
        add_space = False
        add_semicolon = False
        
        if filename:
            msg += "%s" % filename
            add_space = True
            add_semicolon = True
            
        if lineno:
            add_semicolon = True        
            if add_space:
                msg += " "
            msg += "line %d" % lineno
        
        if add_semicolon:
            msg += ": "
        
        msg = msg + errmsg
        
        Exception.__init__(self, msg)


#===========================================================================
# Types handling
#

def guess_type(text):
    """Guesses the type of a value encoded in a string"""
    
    if text.isdigit():
        return int

    try:
        float(text)
        return float
    except ValueError:
        pass
    
    try:
        str2bool(text)
        return bool
    except ValueError:
        pass
    
    return str


def str2bool(text=None):
    """Will parse every way manolis stores a boolean as a string"""
    
    if text is None:
        # default value
        return False
    
    text2 = text.lower()
    
    if text2 == "false":
        return False
    elif text2 == "true":
        return True
    else:
        raise ValueError("unknown string for bool '%s'" % text)
        




#=============================================================================

_type_definitions = [["string", str],
                     ["unknown", str], # backwards compatiable name
                     ["str",    str],  # backwards compatiable name
                     ["int",    int],
                     ["float",  float],
                     ["bool",   bool]]

# NOTE: ordering of name-type pairs is important
#   the first occurrence of a type gives the perferred name for writing
    

def parse_type(text):
    for name, t in _type_definitions:
        if text == name:
            return t
    raise Exception("unknown type '%s'" % text)


def format_type(t):
    for name, t2 in _type_definitions:
        if t == t2:
            return name
    raise Exception("unknown type '%s'" % t)



#===========================================================================
# Table class
#

class Table (list):
    """Class implementing the Portable Table Format"""

    def __init__(self, rows=None, 
                       headers=None,
                       types={},
                       filename=None):
        
        # set table info
        self.headers = copy.copy(headers)
        self.types = copy.copy(types)
        self.filename = filename
        self.comments = []
        self.delim = "\t"
        self.nheaders = 1
                
        
        # set data
        if rows is not None:
            it = iter(rows)
            try:
                first_row = it.next()

                # data is a list of dicts
                if isinstance(first_row, dict):
                    self.append(first_row)
                    for row in it:
                        self.append(dict(row))

                    if self.headers is None:
                        self.headers = sorted(self[0].keys())

                # data is a list of lists
                elif isinstance(first_row, (list, tuple)):
                    if self.headers is None:
                        self.headers = range(len(first_row))
                        self.nheaders = 0
                    for row in itertools.chain([first_row], it):
                        self.append(dict(zip(self.headers, row)))


                # set table info
                for key in self.headers:
                    # guess any types not specified
                    if key not in self.types:
                        self.types[key] = type(self[0][key])
                
            except StopIteration:
                pass
            
        
            
    def clear(self, headers=None, delim="\t", nheaders=1, types=None):
        """Clears the contents of the table"""
        
        self[:] = []
        self.headers = copy.copy(headers)
        if types is None:
            self.types = {}
        else:
            self.types = copy.copy(types)
        self.comments = []
        self.delim = delim
        self.nheaders = nheaders
    
    
    def new(self, headers=None):
        """
        return a new table with the same info but no data
            
        headers - if specified, only a subset of the headers will be copied
        """
        
        if headers is None:
            headers = self.headers
        
        tab = type(self)(headers=headers)
        
        tab.types = util.subdict(self.types, headers)
        tab.comments = copy.copy(self.comments)
        tab.delim = self.delim
        tab.nheaders = self.nheaders
        
        return tab
    
    
    #===================================================================
    # Input/Output
    #
    
    def read(self, filename, delim="\t", nheaders=1,
                   headers=None, types=None, guess_types=True):
        for row in self.read_iter(filename, delim=delim, nheaders=nheaders,
                                  headers=headers, types=types,
                                  guess_types=guess_types):
            self.append(row)
            
    
    def read_iter(self, filename, delim="\t", nheaders=1,
                 headers=None, types=None, guess_types=True):
        """Reads a character delimited file and returns a list of dictionaries
            
           notes:
           Lines that start with '#' are treated as comments and are skiped
           Blank lines are skipped.

           If the first comment starts with '#Types:' the following tokens
           are interpreted as the data type of the column and values in that
           column are automatically converted.

           supported datatypes:
           - string
           - int
           - float
           - bool
           - unknown (no conversion is done, left as a string)

        """
        
        infile = util.open_stream(filename)
        
        # remember filename for later saving
        if isinstance(filename, str):
            self.filename = filename
        

        # clear table
        self.clear(headers, delim, nheaders, types)

                
        # temps for reading only
        self.tmptypes = None

        
        # line number for error reporting
        lineno = 0
        
        
        try:
            for line in infile:
                line = line.rstrip()        
                lineno += 1

                # skip blank lines
                if len(line) == 0:
                    continue

                # handle comments
                if line[0] == "#":
                    if not self._read_directive(line):
                        self.comments.append(line)
                    continue

                # split row into tokens
                tokens = line.split(delim)
                
                # if no headers read yet, use this line as a header
                if not self.headers:
                    # parse headers
                    if self.nheaders > 0:
                        self._parse_header(tokens)
                        continue
                    else:
                        # default headers are numbers
                        self.headers = range(len(tokens))

                assert len(tokens) == len(self.headers), tokens

                # populate types
                if not self.types:
                    if self.tmptypes:
                        assert len(self.tmptypes) == len(self.headers)
                        self.types = dict(zip(self.headers, self.tmptypes))
                    else:
                        # default types
                        if guess_types:
                            for token, header in zip(tokens, self.headers):
                                self.types.setdefault(header,
                                                      guess_type(token))
                        else:
                            for header in self.headers:
                                self.types.setdefault(header, str)
                
                # parse data
                row = {}
                for i in xrange(len(tokens)):
                    key = self.headers[i]
                    t = self.types[key]
                    if t is bool:
                        row[key] = str2bool(tokens[i])
                    else:
                        row[key] = t(tokens[i])
                
                # return completed row
                yield row
                
                
        except Exception, e:
            # report error in parsing input file
            raise TableException(str(e), self.filename, lineno)
            #raise
        
        
        # clear temps
        del self.tmptypes
        
        raise StopIteration

    
    
    def _parse_header(self, tokens):
        """Parse the tokens as headers"""
        
        self.headers = tokens

        # check that headers are unique
        check = set()
        for header in self.headers:
            if header in check:
                raise TableException("Duplicate header '%s'" % header)
            check.add(header)

    
    
    def write(self, filename=sys.stdout, delim="\t"):
        """Write a table to a file or stream.
           
           If 'filename' is a string it will be opened as a file.
           If 'filename' is a stream it will be written to directly.
        """
        
        # remember filename for later saving
        if isinstance(filename, str):
            self.filename = filename
    
        out = util.open_stream(filename, "w")
        
        self.write_header(out, delim=delim)
        
        # tmp variable
        types = self.types
        
        # write data
        for row in self:
            # code is inlined here for speed
            rowstr = []
            for header in self.headers:
                if header in row:
                    rowstr.append(types[header].__str__(row[header]))
                else:
                    rowstr.append('')
            print >>out, delim.join(rowstr)
    
    
    def write_header(self, out=sys.stdout, delim="\t"):
        # ensure all info is complete
        for key in self.headers:
            if key not in self.types:
                if len(self) > 0:
                    self.types[key] = type(self[0][key])
                else:
                    self.types[key] = str
        
                    
        # ensure types are in directives
        if DIR_TYPES not in self.comments:
            self.comments = [DIR_TYPES] + self.comments
        

        # write comments
        for line in self.comments:
            if isinstance(line, str):
                print >>out, line
            else:
                self._write_directive(line, out, delim)
        
        
        # write header
        if self.nheaders > 0:
            print >>out, delim.join(self.headers)

    
    
    def write_row(self, out, row, delim="\t"):
        rowstr = []
        types = self.types
        for header in self.headers:
            if header in row:
                rowstr.append(types[header].__str__(row[header]))
            else:
                rowstr.append('')
        out.write(delim.join(rowstr))
        out.write("\n")

    # NOTE: back-compat
    writeRow = write_row
    
    
    def save(self):
        """Writes the table to the last used filename for the read() or write()
           function"""
        
        if self.filename is not None:
            self.write(self.filename)
        else:
            raise Exception("Table has no filename")
    
    
    #===================================================================
    # Input/Output: Directives
    #
    
    def _determine_directive(self, line):
        if line.startswith("#Types:") or \
           line.startswith("#types:") or \
           line.startswith("##types:"):
            # backwards compatible
            return DIR_TYPES
            
        else:
            return None
    
    
    
    def _read_directive(self, line):
        """Attempt to read a line with a directive"""
        
        directive = self._determine_directive(line)
        
        if directive is None:
            return False
        
        rest = line[line.index(":")+1:]         
        self.comments.append(directive)
        
        if directive == DIR_TYPES:
            self.tmptypes = map(parse_type, rest.rstrip().split(self.delim))
            return True
        
        else:
            return False
    
    
    def _write_directive(self, line, out, delim):
        """Write a directive"""
        
        if line == DIR_TYPES:
            out.write("##types:" + delim.join(format_type(self.types[h])
                                              for h in self.headers) + "\n")
                
        else:
            raise "unknown directive:", line
    

    #===================================================================
    # Table manipulation
    #
    
    def add(self, **kargs):
        """Add a row to the table
           
           tab.add(col1=val1, col2=val2, col3=val3)
        """
        self.append(kargs)
    
    
    def add_col(self, header, coltype=None, default=NULL, pos=None, data=None):
        """Add a column to the table.  You must populate column data yourself.
        
           header  - name of the column
           coltype - type of the values in that column
           default - default value of the column
           pos     - position to insert column (default: right-end)
        """
        # ensure header is unique
        if header in self.headers:
            raise Exception("header '%s' is already in table" % header)
        
        # default column position is last column
        if pos is None:
            pos = len(self.headers)
        
        # default coltype is guessed from data
        if coltype is None:
            if data is None:
                raise Exception("must specify data or coltype")
            else:
                coltype = type(data[0])
        
        # default value is inferred from column type
        if default is NULL:
            default = coltype()
        
        # update table info
        self.headers.insert(pos, header)
        self.types[header] = coltype
        
        # add data
        if data is not None:
            for i in xrange(len(self)):
                self[i][header] = data[i]

    
    def remove_col(self, *cols):
        """Removes a column from the table"""
        
        for col in cols:
            self.headers.remove(col)
            del self.types[col]
            
            for row in self:
                del row[col]
    
    
    def rename_col(self, oldname, newname):
        """Renames a column"""
        
        # change header
        col = self.headers.index(oldname)
        
        if col == -1:
            raise Exception("column '%s' is not in table" % oldname)
        
        self.headers[col] = newname
        
        # change info
        self.types[newname] = self.types[oldname]
        del self.types[oldname]
        
        # change data
        for row in self:
            row[newname] = row[oldname]
            del row[oldname]

       
    def get_matrix(self, rowheader="rlabels"):
        """Returns mat, rlabels, clabels
        
           where mat is a copy of the table as a 2D list
                 rlabels are the row labels
                 clabels are the column labels
        """
        
        # get labels
        if rowheader is not None and rowheader in self.headers:
            rlabels = self.cget(rowheader)
            clabels = copy.copy(self.headers)
            clabels.remove(rowheader)
        else:
            rlabels = range(len(self))
            clabels = copy.copy(self.headers)

        # get data
        mat = []
        for row in self:
            mat.append(util.mget(row, clabels))
        
        return mat, rlabels, clabels
    
    
    def filter(self, cond):
        """Returns a table with a subset of rows such that cond(row) == True"""
        tab = self.new()
        
        for row in self:
            if cond(row):
                tab.append(row)
        
        return tab


    def map(self, func, headers=None):
        """Returns a new table with each row mapped by function 'func'"""

        if len(self) == 0:
            # handle case of zero length table
            return self.new()

        # determine what table will look like from first row
        first_row = func(self[0])

        # determine headers of new table
        if headers is None:
            # try order new headers the same way as old headers
            headers = first_row.keys()
            lookup = util.list2lookup(self.headers)
            top = len(headers)            
            headers.sort(key=lambda x: (lookup.get(x, top), x))
        
        tab = type(self)(
            itertools.chain([first_row], (func(x) for x in self[1:])),
            headers=headers)
        tab.delim = self.delim
        tab.nheaders = self.nheaders
        
        return tab


    def uniq(self, key=None, col=None):
        """
        Returns a copy of this table with consecutive repeated rows removed
        """

        tab = self.new()

        if len(self) == 0:
            return tab

        if col is not None:
            key = lambda x: x[col]

        if key is None:
            last_row = self[0]
            for row in self[1:]:
                if row != last_row:
                    tab.append(row)
                last_row = row
        else:
            last_row = key(self[0])
            for row in self[1:]:
                key_row = key(row)
                if key_row != last_row:
                    tab.append(row)
                last_row = key_row
            

        return tab
        
    
    def groupby(self, key=None):
        """Groups the row of the table into separate tables based on the 
           function key(row).  Returns a dict where the keys are the values
           retruned from key(row) and the values are tables.
           
           Ex:
           tab = Table([{'name': 'matt', 'major': 'CS'},
                        {'name': 'mike', 'major': 'CS'},
                        {'name': 'alex', 'major': 'bio'}])
           lookup = tab.groupby(lambda x: x['major'])
           
           lookup ==> {'CS': Table([{'name': 'matt', 'major': 'CS'},
                                    {'name': 'mike', 'major': 'CS'}]),
                       'bio': Table([{'name': 'alex', 'major': 'bio'}])}
            
           Can also use a column name such as:
           tab.groupby('major')
            
        """
           
           
        groups = {}
        
        if isinstance(key, str):
            keystr = key
            key = lambda x: x[keystr]
        
        if key is None:
            raise Exception("must specify keyfunc")
        
        
        for row in self:
            key2 = key(row)
            
            # add new table if necessary
            if key2 not in groups:
                groups[key2] = self.new()
            
            groups[key2].append(row)
        
        return groups
    
    
    def lookup(self, *keys, **options):
        """Returns a lookup dict based on a column 'key'
           or multiple keys
           
           extra options:
           default=None
           uselast=False    # allow multiple rows, just use last
        """
        
        options.setdefault("default", None)
        options.setdefault("uselast", False)
        lookup = util.Dict(dim=len(keys), default=options["default"])
        uselast = options["uselast"]
        
        for row in self:
            keys2 = util.mget(row, keys)
            ptr = lookup
            for i in xrange(len(keys2) - 1):
                ptr = lookup[keys2[i]]
            if not uselast and keys2[-1] in ptr:
                raise Exception("duplicate key '%s'" % str(keys2[-1]))
            ptr[keys2[-1]] = row
        
        lookup.insert = False
        return lookup
    
    
    def get(self, rows=None, cols=None):
        """Returns a table with a subset of the rows and columns"""
        
        # determine rows and cols
        if rows is None:
            rows = range(len(self))
        
        if cols is None:
            cols = self.headers
            
        tab = self.new(cols)
        
        # copy data        
        for i in rows:
            row = {}
            for j in cols:
                row[j] = self[i][j]
            tab.append(row)
        
        return tab
    
    
    def cget(self, *cols):
        """Returns columns of the table as separate lists"""
        
        ret = []
        
        for col in cols:
            newcol = []
            ret.append(newcol)
            
            for row in self:
                newcol.append(row[col])
        
        if len(ret) == 1:
            return ret[0]
        else:    
            return ret


    def get_row(self, *rows):
        """Returns row(s) as list(s)"""

        if len(rows) == 1:
            # return one row
            row = self[rows[0]]
            return [row[j] for j in self.headers]

        else:
            # return multiple rows (or zero)
            return [[self[i][j] for j in self.headers]
                    for i in rows]
                
            
        
    
    
    def sort(self, cmp=None, key=None, reverse=False, col=None):
        """Sorts the table inplace"""
        
        if col is not None:
            key = lambda row: row[col]
        elif cmp is None and key is None:
            # sort by first column
            key = lambda row: row[self.headers[0]]
        
        list.sort(self, cmp=cmp, key=key, reverse=reverse)
    
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            # return another table if key is a slice
            tab = self.new()
            tab[:] = list.__getitem__(self, key)
            return tab
        else:
            return list.__getitem__(self, key)
    
    
    def __getslice__(self, a, b):
        # for python version compatibility
        return self.__getitem__(slice(a, b))
    

    def __repr__(self):
        s = StringIO.StringIO("w")
        self.write_pretty(s)
        return s.getvalue()
    
    
    def write_pretty(self, out=sys.stdout, spacing=2):
        mat2, rlabels, clabels = self.get_matrix(rowheader=None)

        mat = []
        
        # get headers
        mat.append(clabels)
        
        # get data
        mat.extend(mat2)
        
        util.printcols(mat, spacing=spacing, out=out)

    
    def __str__(self):
        s = StringIO.StringIO("w")
        self.write(s)
        return s.getvalue()
    
    



#===========================================================================
# convenience functions
#

def read_table(filename, delim="\t", headers=None,
               nheaders=1, types=None,
               guess_types=True):
    """Read a Table from a file written in PTF"""
    
    table = Table()
    table.read(filename, delim=delim, headers=headers,
               nheaders=nheaders, types=types,
               guess_types=guess_types)
    return table

# NOTE: back-compat
readTable = read_table


def iter_table(filename, delim="\t", nheaders=1):
    """Iterate through the rows of a Table from a file written in PTF"""
    
    table = Table()
    return table.read_iter(filename, delim=delim, nheaders=nheaders)

# NOTE: back-compat
iterTable = iter_table


def histtab(items, headers=["item", "count", "percent"]):
    h = util.hist_dict(items)
    tab = Table(headers=headers)
    tot = float(sum(h.itervalues()))

    if len(headers) == 2:    
        for key, val in h.items():
            tab.append({headers[0]: key,
                        headers[1]: val})
    
    elif len(headers) == 3:
        for key, val in h.items():
            tab.append({headers[0]: key,
                        headers[1]: val,
                        headers[2]: val / tot})
    
    else:
        raise Exception("Wrong number of headers (2 or 3 only)")
    
    tab.sort(col=headers[1], reverse=True)
    
    return tab


def join_tables(* args, **kwargs):
    """Join together tables into one table.
       Each argument is a tuple (table_i, key_i, cols_i)
       
       key_i is either a column name or a function that maps a 
       table row to a unique key
    """
    
    if len(args) == 0:
        return Table()
    
    # determine common keys
    tab, key, cols = args[0]
    if isinstance(key, str):
        keys = tab.cget(key)
        lookups = [tab.lookup(key)]        
    else:
        keys = map(key, tab)
        lookup = {}
        for row in tab:
            lookup[key(row)] = row
        lookups = [lookup]
        
    keyset = set(keys)
    

    for tab, key, cols in args[1:]:
        if isinstance(key, str):
            keyset = keyset & set(tab.cget(key))
            lookups.append(tab.lookup(key))            
        else:
            keyset = keyset & set(map(key, tab))
            lookup = {}
            for row in tab:
                lookup[key(row)] = row
            
            lookups.append(lookup)
    
    keys = filter(lambda x: x in keyset, keys)
    
    
    # build new table
    if "headers" not in kwargs:
        headers = util.concat(*util.cget(args, 2))
    else:
        headers = kwargs["headers"]
    tab = Table(headers=headers)
    
    for key in keys:
        row = {}
        for (tab2, key2, cols), lookup in zip(args, lookups):
            row.update(util.subdict(lookup[key], cols))
        tab.append(row)
    
    return tab


def showtab(tab, name='table'):
    """Show a table in a new xterm"""
    
    name = name.replace("'", "")
    tmp = util.tempfile(".", "tmp", ".tab")
    tab.write_pretty(file(tmp, "w"))
    os.system("(xterm -T '%s' -n '%s' -e less -S %s; rm %s) &" %
              (name, name, tmp, tmp))


def sqlget(dbfile, query, maxrows=None, headers=None, headernum=False):
    """Get a table from a sqlite file"""
    try:
        from pysqlite2 import dbapi2 as sqlite    
    except ImportError:
        try:
            from sqlite3 import dbapi2 as sqlite
        except ImportError:
            import sqlite

    # open database
    if hasattr(dbfile, "cursor"):
        con = dbfile
        cur = con.cursor()
        auto_close = False
    else:
        con = sqlite.connect(dbfile, isolation_level="DEFERRED")
        cur = con.cursor()
        auto_close = True
    
    cur.execute(query)

    # infer header names
    if headers is None and not headernum:
        headers = [x[0] for x in cur.description]
    
    if maxrows is not None:
        lst = []
        try:
            for i in xrange(maxrows):
                lst.append(cur.next())
        except StopIteration:
            pass
        tab = Table(lst, headers=headers)
    else:
        tab = Table(list(cur), headers=headers)

    if auto_close:
        con.close()
    return tab


def sqlexe(dbfile, sql):

    try:
        from pysqlite2 import dbapi2 as sqlite    
    except ImportError:
        try:
            from sqlite3 import dbapi2 as sqlite
        except ImportError:
            import sqlite

    # open database
    if hasattr(dbfile, "cursor"):
        con = dbfile
        cur = con.cursor()
        auto_close = False
    else:
        con = sqlite.connect(dbfile, isolation_level="DEFERRED")
        cur = con.cursor()
        auto_close = True

    cur.execute(sql)

    if auto_close:
        con.close()


def sql_create_table(cur, table_name, tab, overwrite=True):
    """Create an SQL based on a tab"""

    def issubclass2(t1, t2):
        if type(t1) != type:
            return False
        return issubclass(t1, t2)    

    # drop old table if needed
    if overwrite:
        cur.execute("DROP TABLE IF EXISTS %s;" % table_name)

    # build columns
    cols = []
    for header in tab.headers:

        t = tab.types[header]


        if issubclass2(t, basestring):
            cols.append("%s TEXT" % header)
        elif issubclass2(t, int):
            cols.append("%s INTEGER" % header)
        elif issubclass2(t, float):
            cols.append("%s FLOAT" % header)
        elif issubclass2(t, bool):
            cols.append("%s BOOLEAN" % header)
        else:
            # default is text
            cols.append("%s TEXT" % header)

    cols = ",".join(cols)

    # create table
    cur.execute("""CREATE TABLE %s (%s);""" %
                    (table_name, cols))

    

#def sql_insert_rows(cur, headers, types, rows

def sqlput(dbfile, table_name, tab, overwrite=True, create=True):
    """Insert a table into a sqlite file"""

    try:
        from pysqlite2 import dbapi2 as sqlite    
    except ImportError:
        try:
            from sqlite3 import dbapi2 as sqlite
        except ImportError:
            import sqlite

    # open database
    if hasattr(dbfile, "cursor"):
        con = dbfile
        cur = con.cursor()
        auto_close = False
    else:
        con = sqlite.connect(dbfile, isolation_level="DEFERRED")
        cur = con.cursor()
        auto_close = True

    # read table from file
    if not isinstance(tab, Table):
        filename = tab
        tab = Table()
        it = tab.read_iter(filename)
        
        try:
            # force a reading of the headers
            row = it.next()
            rows = itertools.chain([row], it)
        except StopIteration:
            rows = []
            pass
    else:
        rows = tab


    if create:
        sql_create_table(cur, table_name, tab, overwrite=overwrite)

    # determine text columns
    def issubclass2(t1, t2):
        if type(t1) != type:
            return False
        return issubclass(t1, t2)    


    text = set()
    for header in tab.headers:
        t = tab.types[header]
        
        if issubclass2(t, basestring) or not (
           issubclass2(t, int) or
           issubclass2(t, float) or
           issubclass2(t, bool)):            
            text.add(header)


    # insert rows
    for row in rows:
        vals = []
        for header in tab.headers:
            if header in text:
                vals.append('"%s"' % row[header])
            else:
                vals.append(tab.types[header].__str__(row[header]))
        vals = ",".join(vals)
        cur.execute("INSERT INTO %s VALUES (%s);" % (table_name, vals))
    
    con.commit()
    if auto_close:
        con.close()
    

#===========================================================================
# Matrix functions
#

def matrix2table(mat, rlabels=None, clabels=None, rowheader="rlabels"):
    """
    convert a matrix into a table
    
    use table.get_matrix()  to convert back to a matrix
    
    """
    
    if clabels is None:
        clabels = range(len(mat[0]))
        nheaders = 0
    else:
        nheaders = 1
    
    if rlabels is None:
        tab = Table(headers=clabels)
    else:
        tab = Table(headers=[rowheader] + clabels)
    tab.nheaders = nheaders
    
   
    for i, row in enumerate(mat):
        if rlabels is not None:
            row2 = {rowheader: rlabels[i]}
        else:
            row2 = {}
            
        for j in xrange(len(mat[i])):
            row2[clabels[j]] = mat[i][j]
        
        tab.append(row2)
    
    return tab


def write_matrix(filename, mat, rlabels=None, clabels=None, rowheader="rlabels"):
    tab = matrix2table(mat,
                       rlabels=rlabels,
                       clabels=clabels,
                       rowheader=rowheader)
    tab.write(filename)



def read_matrix(filename, rowheader="rlabels"):
    tab = read_table(filename)    
    mat, rlabels, clabels = tab.get_matrix(rowheader=rowheader)
    return mat, rlabels, clabels


#===========================================================================
# testing
#


if __name__ == "__main__":
    import StringIO
    

    
    #################################################
    text="""\
##types:str	int	int
#
# hello
#
name	0	1
matt	123	3
alex	456	2
mike	789	1
"""

    tab = read_table(StringIO.StringIO(text), nheaders=0)    
    
    print tab
    print tab[0][1]
    
    
    tab.add_col('extra', bool, False)
    for row in tab:
        row['extra'] = True
    
    

    #################################################
    text="""\
##types:str	int	int
name	num	num2
matt	123	3
alex	456	2
mike	789	1
"""

    tab = read_table(StringIO.StringIO(text))
    tab.sort()
    
    print repr(tab)
    print tab
    print tab.cget('name', 'num')


    #################################################
    # guess types
    text="""\
name	num	num2	status
matt	11123	3.0	false
alex	456	2.0	true
mike	789	1.0	false
"""

    tab = read_table(StringIO.StringIO(text))
    tab.sort()
    
    print repr(tab)


    
'''
    #################################################
    # catch parse error
    if 0:
        text="""\
##types:str	int	int
name	num	num
matt	123	0
alex	456	2
mike	789	1
"""

        tab = readTable(StringIO.StringIO(text))
        tab.sort()

        print repr(tab)
        print tab
        print tab.cget('name', 'num')


    #################################################
    # timing
    if 0:
        from rasmus import util
        
        text=["##types:" + "int\t" * 99 + "int",
              "\t".join(map(str, range(100))) ]

        for i in range(10000):
            text.append("1\t" * 99 + "1")
        text = "\n".join(text)
        
        stream = StringIO.StringIO(text)
        
        util.tic("read table")
        tab = readTable(stream)
        util.toc()
    
    
    #################################################
    # specialized types
    if 1:
        text="""\
##types:str	int	strand_type
name	num	strand
matt	123	+
alex	456	-
mike	789	+
john	0	+
"""
        
        
       
        
        class strand_type:
            def __init__(self, text=None):
                if text is None:
                    self.val = True
                else:
                    if text == "+":
                        self.val = True
                    elif text == "-":
                        self.val = False
                    else:
                        raise Exception("cannot parse '%s' as strand_type" % 
                                        str(text))
                
            
            def __str__(self):
                if self.val:
                    return "+"
                else:
                    return "-"
        

        def strand_parser(text=None):
            if text is None:
                return True
            else:
                if text == "+":
                    return True
                elif text == "-":
                    return False
                else:
                    raise Exception("cannot parse '%s' as strand_type" % 
                                    str(text))
        
        def strand_formatter(val):
            if val:
                return "+"
            else:
                return "-"
        
        strand_type = TableType(strand_parser, strand_formatter)
                    

        stream = StringIO.StringIO(text)
        tab = readTable(stream, type_lookup=[["strand_type", strand_type]])
        print tab.types
        print tab
    
    #################################################
    # quoted strings
    if 1:
        text=\
r"""##types:str	bool	quoted_string
name	num	blah
matt	True	hello\tthere
alex	False	hello\nthere
mike	True	hello\\there
john	False	hello\n\\\nthere
"""                    

        stream = StringIO.StringIO(text)
        tab = readTable(stream)
        print tab.types
        print tab
    
    
    #################################################
    # python data structures/code
    if 1:
        def eval2(text=None):
            if text is None:
                return None
            else:
                return eval(text)
        
        python_type = TableType(eval2, str)
    
    
    
        tab = Table(headers=["name", "list"],
                    types={"list": python_type},
                    type_lookup=[["python", python_type]])
                    
        
        tab.append({"name": "matt", "list": [1,2,3]})
        tab.append({"name": "mike", "list": [4,5,6]})
        tab.append({"name": "alex", "list": [7,8,9]})
        
        tab.write()
    
    ##################################################
    # join tables
    if 1:
        tab1 = Table([[0, 1, 2],
                      [1, 3, 4],
                      [2, 5, 6],
                      [3, 7, 8]],
                     headers=['a', 'b', 'c'])
        tab2 = Table([[0, 6, 6],
                      [1, 7, 7],
                      [3, 8, 8]],
                     headers=['a2', 'b2', 'c2'])
        
        tab3 = joinTables((tab1, lambda x: x['a']+1, ['c', 'b']), (tab2, 'a2', ['b2']))
        
        print tab3
    
'''
