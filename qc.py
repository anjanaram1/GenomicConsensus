from intervals import *

class QCFailure(Exception):
    pass

def mappingsAreConcordant(templateIntervals):
    # Are the mappings concordant?
    FUDGE=30
    s, e = max(templateIntervals, key=ilen)
    s -= FUDGE
    e += FUDGE
    return all(within(i, (s, e)) for i in templateIntervals)


def checkMapping(alnSubreads):
    # Require unique mapping for now.
    if not all(sr.MapQV == 254 for sr in alnSubreads): return False
    if not len(set(sr.referenceId for sr in alnSubreads)) == 1: return False
    templateIntervals = [ (a.tStart, a.tEnd) for a in alnSubreads ]
    if not mappingsAreConcordant(templateIntervals): return False
    return True
