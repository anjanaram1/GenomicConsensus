import sys, numpy as np, ipdb, logging
from pbcore.io import *
import GenomicConsensus.quiver as Q
import ConsensusCore as cc

from pbcore.io.BasH5IO import (HQ_REGION,
                               ADAPTER_REGION,
                               INSERT_REGION)

from intervals import *
from qc import *



c = CmpH5Reader("~/Data/lambdaJob/aligned_subreads.cmp.h5")
c.attach("~/Data/lambdaJob/input.1.fofn")
referenceGenome = FastaTable("~/Data/lambdaNEB.fa")[0].sequence[:]


logging.basicConfig(level=logging.WARN)

class ConsensusSequence(object):
    def __init__(self, name, seq, qual, coverage, predictedAccuracy):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.coverage = coverage
        self.predictedAccuracy = predictedAccuracy

def checkUserPreFilters(alnSubreads):
    return len(alnSubreads) >= 3

def checkUserPostFilters(cssObj):
    return (cssObj.coverage >= 3 and
            cssObj.predictedAccuracy >= 0.9)

def estimateAccuracy(confidence):
    def unphred(q):
        return 10**(-q/10)
    def phred(p):
        return -10*np.log10(p)
    errorProbs = unphred(np.array(confidence, dtype=float))
    errorRate = np.mean(errorProbs)
    errorRate = (errorRate*len(confidence) + 1.) / (len(confidence) + 1)
    return 1.-errorRate

def rccsOn(quiverConfig, alnSubreads):
    # Raise exceptions based on different QC failure type
    if not checkMapping(alnSubreads):
        raise QCFailure("Subread mappings are ambiguous or discordant")
    refId = alnSubreads[0].referenceId
    v = hullMany([(a.tStart, a.tEnd) for a in alnSubreads])
    window = (refId, v[0], v[1])
    windowLen = v[1] - v[0]
    refSeq = referenceGenome[v[0]:v[1]]
    css = Q.utils.consensusForAlignments(window,
                                         refSeq,
                                         alnSubreads,
                                         quiverConfig)
    cssName = "/".join(alnSubreads[0].readName.split("/")[:-1]) + "/ccs"
    coverage = sum((a.referenceSpan > 0.8*windowLen) for a in alnSubreads)

    cssObj = ConsensusSequence(cssName,
                               css.sequence,
                               css.confidence,
                               len(alnSubreads),
                               estimateAccuracy(css.confidence))
    return cssObj



quiverConfig = Q.model.loadQuiverConfig("P4-C2.AllQVsMergingByChannelModel")
fastaOut = FastaWriter("css.fasta")
csvOut = open("css.csv", "w")
csvOut.write("Name,NumPasses,PredictedAccuracy\n")

for movieID, holeNumber in sorted(set(zip(c.MovieID,
                                          c.HoleNumber))):

    molecule = (movieID, holeNumber)
    alns = c[((c.MovieID    == movieID) &
              (c.HoleNumber == holeNumber))]
    try:
        if not checkUserPreFilters(alns): continue
        css = rccsOn(quiverConfig, alns)
        print "%r OK: " % (molecule,), css.coverage, css.predictedAccuracy
        if not checkUserPostFilters(css): continue
        fastaOut.writeRecord(css.name, css.seq)
        csvOut.write("%s,%d,%0.3f\n" % (css.name, css.coverage, css.predictedAccuracy))

    except QCFailure as e:
        print "%r Skipping: %s" % (molecule, str(e))


    #print molecule, alns
    #raw_input()
