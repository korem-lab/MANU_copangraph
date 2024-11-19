FASTA_HEADER_PREF = '>'
FIELD_DELIMITER = ';'
FASTG_NODE_DELIMITER = ':'
import sys

class GFATypes:
    S = 'S'
    L = 'L'

class GFA:
    def __init__(self, t):
        self.type = t

class GFASegment(GFA):
    def __init__(self, line):
        t, nid, seq, info = line.strip().split()
        super().__init__(t)
        self.nid = nid
        self.seq = seq
        self.info = info
        self.sample_name = int(self.info.split(':')[3])
        self.contig_name = self.info.split(':')[2]
        self.lb = int(self.info.split(':')[5]) 
        self.rb = int(self.info.split(':')[6]) 

    def __repr__(self):
        return f'GFASegment(nid={self.nid}, seq={self.seq[:10]})'


class GFALine(GFA):
    def __init__(self, line):
        t, l_nid, l_ori, r_nid, r_ori, cigar = line.strip().split()
        super().__init__(t)
        self.r_nid = r_nid
        self.l_nid = l_nid
        self.r_ori = 1 if r_ori == '+' else 0
        self.l_ori = 1 if l_ori == '+' else 0

    def __repr__(self):
        return f'GFALine(l_nid={self.l_nid}, l_ori={" +" if self.l_ori==1 else " -"}, ' + \
               f'r_nid={self.r_nid}, r_ori={" +" if self.r_ori==1 else " -"})'


class Fasta:
    def __init__(self, hdr, seq, fold=80):
        self.hdr = hdr
        self.seq = seq
        self.fold = fold

    def write(self, fs):
        fs.write(
            f'>{self.hdr}\n' +
            '\n'.join([self.seq[i:i+self.fold] for i in range(0, len(self.seq), self.fold)]) +
            '\n'
        )

    def getlen(self):
        return self.seq.__len__()

    def __repr__(self):
        return f'Fasta(hdr={self.hdr}, seq={self.seq[:10]})'


class Fastg(Fasta):
    def __init__(self, hdr, seq, fold=80):
        super().__init__(hdr, seq, fold)
        if FASTG_NODE_DELIMITER in hdr:
            self.out_node, self.in_node = hdr.split(FASTG_NODE_DELIMITER)
            self.in_node = self.in_node[:-1]
        else:
            self.out_node = hdr[:-1]
            self.in_node = None


class Genome(Fasta):
    def __init__(self, hdr, seq, fold=80):
        super().__init__(hdr, seq, fold)


def get_gfa_element(line):
    if line[0] == GFATypes.S:
        return GFASegment(line)
    elif line[0] == GFATypes.L:
        return GFALine(line)


def parse_gfa(fs):
    for line in fs:
        yield get_gfa_element(line)
    fs.close()


def parse(fs, fa_cls, *args, **kwargs):
    '''Input: list of contig fasta file names.'''
    line = next(fs)
    while line[0] == FASTA_HEADER_PREF:
        # strip header, get rid of '>', and split by ':'
        hdr = line.strip()[1:]
        line = next(fs)  # move to seq
        seq = list()
        while line[0] != FASTA_HEADER_PREF:
            seq.append(line.strip())
            line = next(fs, False)
            if not line:
                break
        # reached header for next elem
        yield fa_cls(hdr, ''.join(seq), *args, **kwargs)
        if not line:
            break
    