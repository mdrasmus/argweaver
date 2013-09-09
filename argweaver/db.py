import sqlite3 as sqlite

import arghmm
import arghmm.vis


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

    def close(self):
        self.con.close()

    def commit(self):
        self.con.commit()

    def execute(self, *args, **kwargs):
        return self.con.execute(*args, **kwargs)

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
            self.con.commit()
        else:
            assert names == names2

    def _query_names(self):
        res = self.con.execute(
            "SELECT name FROM SitesSequences ORDER BY seq_order")
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

    def get_sites(self, chrom, start, end, names=None):
        res = self.con.execute("""SELECT pos, col FROM Sites
                                  WHERE chrom = ? and ? <= pos and pos <= ?
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


class ArgDB (object):

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

    def close(self):
        self.con.close()

    def commit(self):
        self.con.commit()

    def execute(self, *args, **kwargs):
        return self.con.execute(*args, **kwargs)

    def init_tables(self):
        self.con.execute(u"""CREATE TABLE IF NOT EXISTS ArgTrees (
                             sample INTEGER,
                             chrom TEXT,
                             start INTEGER,
                             end INTEGER,
                             tree TEXT,
                             UNIQUE(sample, chrom, start, end)
                                 ON CONFLICT REPLACE);
                             """)
        self.con.execute(u"""CREATE INDEX IF NOT EXISTS IdxArgTrees
                             ON ArgTrees (sample, chrom, start);""")

        self.con.execute(u"""CREATE TABLE IF NOT EXISTS ArgSprs (
                             sample INTEGER,
                             chrom TEXT,
                             pos INTEGER,
                             recomb_node INTEGER,
                             recomb_time FLOAT,
                             coal_node INTEGER,
                             coal_time FLOAT,
                             UNIQUE(sample, chrom, pos) ON CONFLICT REPLACE);
                             """)
        self.con.execute(u"""CREATE INDEX IF NOT EXISTS IdxArgSprs
                             ON ArgSprs (sample, chrom, pos);""")

        self.con.execute(u"""CREATE TABLE IF NOT EXISTS ArgSequences (
                             name TEXT,
                             seq_order INTEGER)""")

    def add_seq_names(self, names):

        names2 = self._query_names()

        if len(names2) == 0:
            # insert names
            for i, name in enumerate(names):
                self.con.execute("INSERT INTO ArgSequences VALUES (?, ?)",
                                 (name, i))
            self.con.commit()
        else:
            assert names == names2

    def _query_names(self):
        res = self.con.execute(
            "SELECT name FROM ArgSequences ORDER BY seq_order")
        return [row[0] for row in res]

    def get_names(self):
        if self._names is None:
            self._names = self._query_names()
        return self._names

    def add_smc_file(self, filename, sample=0):
        it = arghmm.iter_smc_file(filename)
        header = it.next()
        self.add_seq_names(header["names"])
        region = it.next()
        self.add_smc(region["chrom"], it, sample=sample)

    def add_smc(self, chrom, items, sample=0):
        for item in items:
            if item["tag"] == "TREE":
                self.con.execute(
                    "INSERT INTO ArgTrees VALUES (?, ?, ?, ?, ?);",
                    (sample, chrom, item["start"], item["end"], item["tree"]))
            elif item["tag"] == "SPR":
                self.con.execute(
                    "INSERT INTO ArgSprs VALUES (?, ?, ?, ?, ?, ?, ?);",
                    (sample, chrom, item["pos"],
                     item["recomb_node"], item["recomb_time"],
                     item["coal_node"], item["coal_time"]))

    def get_trees(self, chrom, start, end, sample=0):
        for row in self.con.execute(
            """SELECT start, end, tree FROM ArgTrees
               WHERE sample = ? and chrom = ? and ? < end and start < ?
               ORDER BY start;""",
                (sample, chrom, start, end)):
            yield {"tag": "TREE",
                   "start": row[0], "end": row[1], "tree": row[2]}

    def get_sprs(self, chrom, start, end, sample=0):
        for row in self.con.execute(
                """SELECT pos, recomb_node, recomb_time, coal_node, coal_time
                   FROM ArgSprs
                   WHERE sample = ? and chrom = ? and ? <= pos and pos < ?
                   ORDER BY pos;""",
                (sample, chrom, start, end)):
            yield {"tag": "SPR",
                   "pos": row[0],
                   "recomb_node": row[1],
                   "recomb_time": row[2],
                   "coal_node": row[3],
                   "coal_time": row[4]}


class ArgLayoutDB(object):

    def __init__(self):
        self.files = []

    def add_layout_file(self, filename):
        self.files.append(filename)

    def get_layout(self, chrom, start, end):

        for filename in self.files:
            for block_layout in arghmm.vis.query_arg_layout(
                    filename, chrom, start, end):
                yield block_layout
