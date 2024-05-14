FASTA_HEADER_PREF = '>'
FIELD_DELIMITER = ';'
FASTG_NODE_DELIMITER = ':'


class GFATypes:
    S = 'S'
    L = 'L'


class GFA:
    def __init__(self, t):
        self.type = t


class GFASegment(GFA):
    def __init__(self, line):
        t, nid, seq, dp, kc = line.strip().split()
        super().__init__(t)
        self.nid = int(nid)
        self.seq = seq
        self.kc = int(kc.split(':')[-1])
        self.dp = float(dp.split(':')[-1])


class GFALine(GFA):
    def __init__(self, line):
        t, o_nid, o_ori, i_nid, i_ori, cigar = line.strip().split()
        super().__init__(t)
        self.i_nid = int(i_nid)
        self.o_nid = int(o_nid)
        self.i_ori = 1 if i_ori == '+' else -1
        self.o_ori = 1 if o_ori == '+' else -1


def get_gfa_element(line):
    if line[0] == GFATypes.S:
        return GFASegment(line)
    elif line[0] == GFATypes.L:
        return GFALine(line)


def parse_gfa(fs):
    for line in fs:
        yield get_gfa_element(line)


def parse(fs, fa_cls):
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
        yield fa_cls(hdr, ''.join(seq))
        if not line:
            break


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


class RGNodeSeq(Fasta):
    def __init__(self, hdr, seq, fold=80):
        super().__init__(hdr, seq, fold)
        self.node_id, self.sample_id, self.contig_id, self.start_pos, self.end_pos = hdr.split(FIELD_DELIMITER)
        self.start_pos = int(float(self.start_pos))
        self.end_pos = int(float(self.end_pos))


class PEExtContig(Fasta):
    def __init__(self, hdr, seq, fold=80):
        super().__init__(hdr, seq, fold)
        self.sample_id, self.contig_id, self.bwd_ext, self.fwd_ext = hdr.split(FIELD_DELIMITER)
        self.bwd_ext = int(self.bwd_ext)
        self.fwd_ext = int(self.fwd_ext)
