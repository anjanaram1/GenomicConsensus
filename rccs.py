import sys, numpy as np, ipdb
from pbcore.io import *
import GenomicConsensus.quiver as Q
import ConsensusCore as cc

from pbcore.io.BasH5IO import (HQ_REGION,
                               ADAPTER_REGION,
                               INSERT_REGION)

from intervals import *
from qc import *

class ConsensusSequence(object):
    def __init__(self, seq, qual, coverage, estimatedAccuracy):
        self.seq = seq
        self.qual = qual
        self.coverage = coverage
        self.estimatedAccuracy = estimatedAccuracy

def checkUserPreFilters(alnSubreads):
    return len(alnSubreads) >= 3

def checkUserPostFilters(alnSubreads):
    return True


def rccsOn(quiverConfig, alnSubreads):
    # Raise exceptions based on different QC failure type
    if not checkMapping(alnSubreads):
        raise QCFailure("Subread mappings are ambiguous or discordant")
    refId = alnSubreads[0].referenceId
    v = hullMany([(a.tStart, a.tEnd) for a in alnSubreads])
    window = (refId, v[0], v[1])
    windowLen = v[1] - v[0]
    refSeq = "N" * windowLen
    css = Q.utils.consensusForAlignments(window,
                                         refSeq,
                                         alnSubreads,
                                         quiverConfig)
    return ConsensusSequence("foo", "bar", "baz", "quux")


c = CmpH5Reader("~/Data/lambdaJob/aligned_subreads.cmp.h5")
c.attach("~/Data/lambdaJob/input.1.fofn")

for movieID, holeNumber in sorted(set(zip(c.MovieID,
                                          c.HoleNumber))):

    quiverConfig = Q.model.loadQuiverConfig("P4-C2.AllQVsMergingByChannelModel")
    quiverConfig.minPoaCoverage = 1


    molecule = (movieID, holeNumber)
    alns = c[((c.MovieID    == movieID) &
              (c.HoleNumber == holeNumber))]
    try:
        if not checkUserPreFilters(alns): continue
        css = rccsOn(quiverConfig, alns)
        if not checkUserPostFilters(css): continue
        print "OK: %r" % (molecule,)
    except QCFailure as e:
        print "Skipping %r: %s" % (molecule, str(e))


    #print molecule, alns
    #raw_input()
