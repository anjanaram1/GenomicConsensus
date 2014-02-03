from GenomicConsensus.io import BamReader
from pbcore.io import CmpH5Reader
import pbcore.data as D
import string

from nose.tools import assert_equal as EQ
from numpy.testing import assert_array_equal as AEQ

def compareReads(aln1, aln2):
    for o in ["genomic", "native"]:
        for a in [True, False]:
            EQ(aln1.read(orientation=o, aligned=a),
               aln2.read(orientation=o, aligned=a))

TT = string.maketrans("ATGC", "NNNN")

def compareReferences(cmpAln, bamAln):
    for o in ["genomic", "native"]:
        for a in [True, False]:
            cmpRef = cmpAln.reference(orientation=o, aligned=a)
            bamRef = bamAln.reference(orientation=o, aligned=a)
            EQ(cmpRef.translate(TT), bamRef)

class test_BamReader(object):

    def __init__(self):
        self.c = CmpH5Reader(D.getCmpH5())
        self.b = BamReader  (D.getBam())

        bAlns = list(self.b)
        # FWD alns:
        self.baf = bAlns[1]
        self.caf = self.c[1]
        # REV alns:
        self.bar = bAlns[0]
        self.car = self.c[0]

    def test_read(self):
        compareReads(self.car, self.bar)
        compareReads(self.caf, self.baf)

    def test_reference(self):
        compareReferences(self.car, self.bar)
        compareReferences(self.caf, self.baf)

    def test_readCoords(self):
        EQ(self.caf.readStart, self.baf.readStart)
        EQ(self.car.readStart, self.bar.readStart)
        EQ(self.caf.readEnd, self.baf.readEnd)
        EQ(self.car.readEnd, self.bar.readEnd)

    def test_readName(self):
        EQ(self.caf.readName, self.baf.readName)
        EQ(self.car.readName, self.bar.readName)

    def test_clippedTo(self):
        clippedBaf1 = self.baf.clippedTo(100, 120)
        clippedCaf1 = self.caf.clippedTo(100, 120)
        compareReads(clippedCaf1, clippedBaf1)

        clippedBar1 = self.bar.clippedTo(100, 120)
        clippedCar1 = self.car.clippedTo(100, 120)
        compareReads(clippedCar1, clippedBar1)

        # Clipping off the end
        clippedBaf2 = self.baf.clippedTo(360, 380)
        clippedCaf2 = self.caf.clippedTo(360, 380)
        compareReads(clippedCaf2, clippedBaf2)

        clippedBar2 = self.bar.clippedTo(280, 320)
        clippedCar2 = self.car.clippedTo(280, 320)
        compareReads(clippedCar2, clippedBar2)

    def test_referencePositions(self):
        for o in ["native", "genomic"]:
            AEQ(self.caf.referencePositions(o), self.baf.referencePositions(o))
            AEQ(self.car.referencePositions(o), self.bar.referencePositions(o))

    def test_reads_in_range(self):
        for wStart in xrange(0, 50000, 500):
            bReads = list(self.b.readsInRange(0, wStart, wStart+500))
            cReads = list(self.c.readsInRange(1, wStart, wStart+500))
            EQ(len(cReads), len(bReads))
            EQ([c.readName for c in cReads],
               [b.readName for b in bReads])


if __name__ == '__main__':
    c = CmpH5Reader(D.getCmpH5())
    b = BamReader  (D.getBam())

    bAlns = list(b.dual)
    print bAlns[0].seq
    print c[0].read(orientation="genomic", aligned=False)


    print c[0]
    print c[1]
