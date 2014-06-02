import sys, numpy as np, ipdb
from pbcore.io import *
import GenomicConsensus.quiver as Q
import ConsensusCore as cc

from pbcore.io.BasH5IO import (HQ_REGION,
                               ADAPTER_REGION,
                               INSERT_REGION)

def doCCS(zmw):
    pass

def doCCS(subreads):
    pass

def areAdjacent(rl, rr):
    return rl.regionEnd == rr.regionStart

def isAdapter(r):
    return r.regionType == ADAPTER_REGION

def isInsert(r):
    return r.regionType == INSERT_REGION

def regionLength(r):
    return r.regionEnd - r.regionStart

def alternatingTypes(rl, rr):
    return (isAdapter(rl) and isInsert(rr) or \
            isInsert(rl)  and isAdapter(rr))

def scoreRun(regs):
    if len(regs) < 3: return 0
    if (not isAdapter(regs[0])) or (not isAdapter(regs[-1])): return 0
    if not all((alternatingTypes(regs[i], regs[i+1]) and
                areAdjacent(regs[i], regs[i+1]))
               for i in xrange(0, len(regs)-1)): return 0
    if not all(areAdjacent(regs[i], regs[i+1])
               for i in xrange(0, len(regs)-1)): return 0
    return sum([regionLength(reg) for reg in regs
                if isInsert(reg)])

def bracketedSubreads(zmw):
    """
    Get maximal list of subreads from the maximal alternating run
    adapter-insert-adapter-insert-....-insert-adapter (beginning and
    ending with adapters), where regions in the run are all adjacent.
    """
    rt = zmw.regionTable
    rt = rt[(rt.regionType == ADAPTER_REGION) |
            (rt.regionType == INSERT_REGION)]
    rt = np.sort(rt, order="regionStart")

    L = range(len(rt)+1)
    bestScoreRun = max([(i,j) for i in L for j in L if i < j],
                       key=lambda (i,j): scoreRun(rt[i:j]))
    regions = rt[bestScoreRun[0]:bestScoreRun[1]]
    if scoreRun(regions) == 0:
        return []
    else:
        insertRegions = filter(isInsert, regions)
        return [ zmw.read(region.regionStart, region.regionEnd)
                 for region in insertRegions ]


# def main():
#     basFname = sys.argv[1]
#     bas = BasH5Reader(basFname)

#     for zmw in bas.sequencingZmws:
#         print zmw.regions


#if __name__ == '__main__':


# bas = BasH5Reader("~/Data/singleZmwBas/zmw139.bax.h5")
# z = bas[139]
# rt = z.regionTable
# rt = rt[(rt.regionType == ADAPTER_REGION) |
#         (rt.regionType == INSERT_REGION)]
# rt = np.sort(rt, order="regionStart")


bas = BasH5Reader("~/Data/m140518_125000_42140_c110040372550000001823110806241437_s1_p0.1.bax.h5")
for z in bas:
    srs = bracketedSubreads(z)
    srLens = map(len, srs)
    if (len(srLens) >= 2 and max(srLens)-min(srLens) > 0.2*max(srLens)):
        print "Odd zmw:"
        print z.holeNumber, bracketedSubreads(z)
        ipdb.set_trace()
