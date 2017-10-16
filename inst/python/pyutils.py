import sys, os, pysam, statistics
from collections import defaultdict
#import pdb

# bamfile = "/home/gerhard/Projects/graphmap/output/common1/lbc4/shotgun/bwamem.illumina.A.default.100pct.7MAPQ.sorted.bam"
# inpos = [704, 5525]
# ppos = [478, 705, 706, 707, 708, 709, 710, 718, 742, 811]
# nucA = ["A", "T", "A", "T", "G", "C", "T", "G", "T", "T"]
# nucB = ["G", "-", "-", "-", "-", "-", "-", "-", "A", "C"]

def bam_open(fname, mode='r', *args, **kwargs):
    if fname.lower()[-4:] == '.bam':
        return pysam.AlignmentFile(fname, '%sb' % mode, *args, **kwargs)
    if fname.lower()[-5:] == '.cram':
        return pysam.AlignmentFile(fname, '%sc' % mode, *args, **kwargs)
    return pysam.AlignmentFile(fname, '%s' % mode, *args, **kwargs)


def bam_iter(bam):
    for read in bam:
        yield read


def py_get_reads_at_ppos(bamfile, ppos, nucA = None, nucB = None):
    if isinstance(ppos, list):
        ppos = [int(x) for x in ppos]
    else:
        ppos = [int(ppos)]

    if nucA and nucB:
        if isinstance(nucA, list):
            nucA = [str(x) for x in nucA]
        else:
            nucA = [str(nucA)]
        if isinstance(nucB, list):
            nucB = [str(x) for x in nucB]
        else:
            nucB = [str(nucB)]
        assert len(nucA) == len(nucB)
        assert len(ppos) == len(nucA)
        posdict = {pos: {A if B != "+" else "-": [], B if A != "+" else "-": []} for pos, A, B in zip(ppos, nucA, nucB)}
    else:
        posdict = {pos: {"A": [], "C": [], "G": [], "T": [], "-": [], "+": []} for pos in ppos}
    bam = pysam.AlignmentFile(bamfile, "rb")
    ref = bam.references[0]
    pileup = bam.pileup(ref)
    for pileupCol in pileup:
        pos = int(pileupCol.reference_pos)
        if pos not in posdict.keys():
            pass
        else:
            sys.stdout.write("Extracting read IDs for " + " and ".join(posdict[pos].keys()) + " at " + str(pos + 1) + "\n")
            for read in pileupCol.pileups:
                refpos = read.alignment.get_reference_positions(full_length = True)
                try:
                    i = refpos.index(pos)
                    base = read.alignment.query_sequence[i]
                    if "+" not in posdict[pos].keys(): # We don't expect an insertion here
                        try:
                            posdict[pos][base].append(read.alignment.query_name)
                        except KeyError:
                            pass
                    else: # We do expect an insertion here
                        j = refpos.index(pos + 1)
                        ins = read.alignment.query_sequence[i:j][1:]
                        if ins:
                            posdict[pos]["+"].append(read.alignment.query_name)
                        else:
                            posdict[pos]["-"].append(read.alignment.query_name)
                except ValueError:
                    try:
                        posdict[pos]["-"].append(read.alignment.query_name)
                    except KeyError:
                        pass
    return posdict


def py_get_insertions(bamfile, inpos):
    if isinstance(inpos, list):
        inpos = sorted([int(x) for x in inpos])
    else:
        inpos = [int(inpos)]
    sys.stdout.write("Extracting insertions at positions " + ", ".join([str(x) for x in inpos]) + "\n")
    bam = pysam.AlignmentFile(bamfile, "rb")
    ref = bam.references[0]
    pileup = bam.pileup(ref)
    out = defaultdict(list)
    for pileupCol in pileup:
        pos = pileupCol.reference_pos
        nins = 0
        if pos not in inpos:
            pass
        else:
            for read in pileupCol.pileups:
                refpos = read.alignment.get_reference_positions(full_length = True)
                try:
                    i = refpos.index(pos - 1)
                    j = refpos.index(pos)
                    ins = read.alignment.query_sequence[i:j][1:]
                    if ins:
                        nins += 1
                        out[pos].append(ins)
                    else:
                        out[pos].append("-")
                except ValueError:
                    pass
            sys.stdout.write("Pos: " + str(pos) + " Read depth: " + str(pileupCol.nsegments) + " Ins: " + str(nins) + "\n")

    bam.close()
    return out


def prune_read(read):
    if read.is_unmapped:
        return 1

    newcigar = []
    clip_5 = 0
    clip_3 = 0

    changed = False
    inseq = False
    for op, length in read.cigar:
        if op == 5:  # H
            changed = True
        elif op == 4:  # S
            changed = True
            if not inseq:
                clip_5 = length
            else:
                clip_3 = length
        else:
            inseq = True
            newcigar.append((op, length))

    if not changed:
        return 0

    read.cigar = newcigar
    orig_length = len(read.seq)

    s = read.seq
    q = read.qual

    if clip_3:
        read.seq = s[clip_5:-clip_3]
        if q:
            read.qual = q[clip_5:-clip_3]
    else:
        read.seq = s[clip_5:]
        if q:
            read.qual = q[clip_5:]

    newtags = []
    if clip_5:
        newtags.append(('ZA', clip_5))
    if clip_3:
        newtags.append(('ZB', clip_3))

    newtags.append(('ZC', float(clip_5 + clip_3) / orig_length))

    read.tags = read.tags + newtags

    return 2

def py_bam_prune_clipping(infile, outfile):
    bam = pysam.Samfile(infile, "rb")
    out = pysam.Samfile(outfile, "wb", template=bam)
    total = 0
    count = 0
    unmapped = 0
    for read in bam_iter(bam):
        code = prune_read(read)

        if code == 1:
            unmapped += 1
        elif code == 2:
            count += 1

        total += 1
        out.write(read)

    bam.close()
    out.close()
    sys.stderr.write('Wrote %s reads\nAltered: %s\nUnmapped: %s\n' % (total, count, unmapped))

# if __name__=='__main__':
