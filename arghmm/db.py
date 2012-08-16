
import sqlite3  as sqlite

import arghmm


class SitesDB (object):

    def __init__(self):
        self.con = None
        self._names = None

    def connect(self, filename):
        if isinstance(filename, basestring):
            self.con = sqlite.connect(filename, isolation_level="DEFERRED")
        else:
            self.con = filename
        self.init_tables()
        return self

    def init_tables(self):
        self.con.execute(u"""CREATE TABLE IF NOT EXISTS Sites (
                             chrom TEXT,
                             pos INTEGER, 
                             col TEXT,
                             UNIQUE(chrom, pos) ON CONFLICT REPLACE);
                             """)
        self.con.execute(u"""CREATE INDEX IF NOT EXISTS IdxSites 
                             ON Sites (chrom, pos);""")
        
        self.con.execute(u"""CREATE TABLE IF NOT EXISTS SitesSequences (
                             name TEXT,
                             seq_order INTEGER)""")

    def add_seq_names(self, names):
        
        names2 = self._query_names()

        if len(names2) == 0:
            # insert names
            for i, name in enumerate(names):
                self.con.execute("INSERT INTO SitesSequences VALUES (?, ?)",
                                 (name, i))
        else:
            assert names == names2

    def _query_names(self):
        res = self.con.execute("SELECT name FROM SitesSequences ORDER BY seq_order")
        return [row[0] for row in res]
        

    def get_names(self):
        if self._names is None:
            self._names = self._query_names()
        return self._names

    def add_sites_file(self, filename):
        it = arghmm.iter_sites(filename)
        header = it.next()
        self.add_seq_names(header["names"])
        
        self.add_sites(header["chrom"], it)

    def add_sites(self, chrom, sites_iter):
        for pos, col in sites_iter:
            self.con.execute("INSERT INTO Sites VALUES (?, ?, ?);",
                             (chrom, pos, col))

    def commit(self):
        self.con.commit()

    def close(self):
        self.con.close()

    def execute(self, *args, **kwargs):
        return self.con.execute(*args, **kwargs)
        
    def get_sites(self, chrom, start, end, names=None):
        res = self.con.execute("""SELECT pos, col FROM Sites
                                  WHERE chrom = ? and ? <= pos and pos < ?
                                  ORDER BY pos;""",
                               (chrom, start, end))
        
        def subset_by_names(res, ind):
            for pos, col in res:
                yield pos, "".join(col[i] for i in ind)

        if names is None:
            return res
        else:
            lookup = dict((n, i) for i, n in enumerate(self.get_names()))
            ind = [lookup[n] for n in names]
            return subset_by_names(res, ind)
    
