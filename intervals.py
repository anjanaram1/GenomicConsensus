from sys import maxint

def intersect(i1, i2):
    s1, e1 = i1
    s2, e2 = i2
    return not((s1 >= e2) or (s2 >= e1))

def hull(i1, i2):
    s1, e1 = i1
    s2, e2 = i2
    return (min(s1, s2), max(e1, e2))

def hullMany(intervals):
    h = (maxint,-maxint)
    for v in intervals:
        h = hull(h, v)
    return h

def within(needle, haystack):
    ns, ne = needle
    hs, he = haystack
    return hs <= ns <= ne <= he

def ilen(interval):
    return interval[1]-interval[0]
