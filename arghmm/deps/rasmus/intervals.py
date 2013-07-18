#!/usr/bin/env python
#
# common functions for dealing with intervals
#

import heapq
from itertools import chain

from rasmus import util

"""
A region is a list or tuple with the format
 [start, end, ...]

where start <= end and ... can be any additional data

"""


def overlap(region1, region2, inc=True):
    """
    Returns True if range region1=[a,b] overlaps region2=[x,y]
    
    inc -- if True, treat [a,b] and [x,y] as inclusive
    """
    if inc:
        return (region2[1] >= region1[0]) and (region2[0] <= region1[1])
    else:
        return (region2[1] > region1[0]) and (region2[0] < region1[1])


def iter_groups(items, key):
    """
    Iterates through groups of consequentive items x that have the same key(x)

    items -- iterable of values
    key   -- function of one argument
    """
    
    NULL = object()
    last_key = NULL
    group = []

    for item in items:
        k = key(item)
        if k != last_key:
            if group:
                yield group
            
            # start new group
            group = []
            last_key = k
        group.append(item)

    if group:
        yield group


def iter_union_ids(regions):
    """
    Iterate over union groups

    Yields (groupnum, region) for each region in regions

    NOTE: regions must be sorted by start
    """

    # TODO: add inclusive option
    
    start = -util.INF
    end = -util.INF
    group = None
    groupnum = -1

    for reg in regions:        
        if reg[0] > end:
            # start new group
            start = reg[0]
            end = reg[1]
            groupnum += 1

        else:
            # append to current group
            if reg[1] > end:
                end = reg[1]

        yield (groupnum, reg)


def groupby_unions(regions):
    """
    Iterate over union groups
    
    NOTE: regions must be sorted by start
    """

    # TODO: add inclusive option
    
    for group in iter_groups(iter_union_ids(regions), lambda x: x[0]):
        # remove group index from each region
        yield [x[1] for x in group]


def iter_unions(regions):
     """
     Iterate over union groups

     Yields (start, end, [region1, region2, ...]) for each union group

     NOTE: regions must be sorted by start
     """

     # TODO: add inclusive option

     for group in groupby_unions(regions):
         start = min(r[0] for r in group)
         end = max(r[1] for r in group)
         yield (start, end, group)    


def iter_intersections(regions):
    """
    Iterate over intersection groups

    Yields (start, end, [region1, region2]) for each intersection group
    
    NOTE: regions must be sorted by start
    """

    # TODO: think about whether this is inclusive or not
    # useful for DNA coordinates

    # endpoints queue
    group = []

    start = None
    end = None

    for reg in regions:

        while len(group) > 0 and group[0][0] <= reg[0]:
            # process end points

            # yield group upto this endpoint
            end = group[0][0]
            if end != start:
                yield (start, end, [x[1] for x in group])
            heapq.heappop(group)
            start = end
        else:
            # process new start point

            # yield group before new region
            end = reg[0]
            if start != end and len(group) > 0:
                yield (start, end, [x[1] for x in group])

            # add region end to group
            heapq.heappush(group, (reg[1], reg))
            start = end

    # yield remaining groups
    while len(group) > 0:
        end = group[0][0]
        if end != start:
            yield (start, end, [x[1] for x in group])
        heapq.heappop(group)
        start = end


def iter_substract(regions1, regions2):
    """
    Substract regions2 from regions1

    Yields (start, end, original_region_from_regions1)

    NOTE: regions must be sorted by start and
    regions1 and regions2 should each be non-overlaping.
    """

    # add a tag to regions
    regions1 = [(reg[0], reg[1], 1, reg) for reg in regions1]
    regions2 = [(reg[0], reg[1], 2, reg) for reg in regions2]
                
    # combine all regions into one list
    regions = sorted(chain(regions1, regions2), key=lambda x: x[0])

    for a, b, group in iter_intersections(regions):
        if len(group) == 1 and group[0][2] == 1:
            yield (a, b, group[0][3])
            
'''
def iter_combine_regions(*regionsets):
    """
    Combine two or more region sets into one sorted region set

    NOTE: region sets must be sorted by start
    """

    regionsets = map(iter, regionsets)
    next = set()

    for regions in regionsets:
        pass
'''        
    


def query_point_regions(point, regions, inc=True):

    ind = util.sortindex(regions, key=lambda r: r[1])
    rind = util.mget(range(len(regions)), ind)
    regions_by_end = util.mget(regions, ind)

    end = util.binsearch([r[0] for r in regions], x)[1]
    start = util.binsearch([r[1] for r in regions_by_end], x)[0]

    if start is None:
        start = 0
    if end is None:
        end = len(regions)

    if inc:
        for i in xrange(start, end):
            if regions[i][0] <= x <= regions[i][1]:
                yield regions[i]
    else:
        for i in xrange(start, end):
            if regions[i][0] < x < regions[i][1]:
                yield regions[i]


def query_regions_regions(query_regions, regions, inc=True):

    pass



    


if __name__ == "__main__":
    
    print "union"
    print list(iter_union_ids([[1, 10], [2, 4], [2, 5],
                               [12, 20], [13, 22]]))

    print list(iter_unions([[1, 10], [2, 4], [2, 5],
                            [12, 20], [13, 22]]))


    print list(groupby_unions([[1, 10], [2, 4], [2, 5],
                               [12, 20], [13, 22]]))

    print "intersect"
    print list(iter_intersections( 
            [[1, 10], [2, 4], [2, 5],
             [12, 20], [13, 22]]))

    print "union"
    print list(query_point_regions(3, [[1, 10], [2, 4], [2, 5],
                                       [12, 20], [13, 22]]))
