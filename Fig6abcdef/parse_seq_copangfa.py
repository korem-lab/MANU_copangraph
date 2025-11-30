FASTA_HEADER_PREF = '>'
FIELD_DELIMITER = ';'
FASTG_NODE_DELIMITER = ':'
import sys


class GFATypes:
    S = 'S'
    L = 'L'
    H = 'H'


class GFA:
    def __init__(self, t):
        self.type = t

    def write(self, *args, **kwargs):
        NotImplementedError
        return f'GFASegment(nid={self.nid}, seq={self.seq[:10]})'

class GFASegment(GFA):
    def __init__(self, t, nid, seq, *rest):
        super().__init__(t)
        self.nid = nid
        self.seq = seq
        if type(rest) == tuple and len(rest) == 1:
            self.rest = rest[0]
        else:
            self.rest = rest
        
        self.sample_name = self.rest.split(':')[3]
        self.contig_name = self.rest.split(':')[2]
        self.lb = int(self.rest.split(':')[5])
        self.rb = int(self.rest.split(':')[6])

    def __repr__(self):
        return f'GFASegment(nid={self.nid}, seq={self.seq[:10]})'

    def write(self, f):
        f.write(
            '\t'.join([self.type, self.nid, self.seq,self.rest]) + '\n'
        )


class GFALine(GFA):
    def __init__(self, t, l_nid, l_ori, r_nid, r_ori, cigar):
        super().__init__(t)
        self.r_nid = r_nid
        self.l_nid = l_nid
        self.r_ori = 1 if r_ori == '+' else 0
        self.l_ori = 1 if l_ori == '+' else 0
        self.cigar = cigar

    def __repr__(self):
        return f'GFALine(l_nid={self.l_nid}, l_ori={" +" if self.l_ori==1 else " -"}, ' + \
               f'r_nid={self.r_nid}, r_ori={" +" if self.r_ori==1 else " -"})'

    def write(self, f):
        f.write(
            '\t'.join([
                self.type,
                self.l_nid, '+' if self.l_ori == 1 else '-',
                self.r_nid, '+' if self.r_ori == 1 else '-',
                self.cigar
            ]) + '\n'
        )


def get_gfa_element(line):
    if line[0] == GFATypes.S:
        return GFASegment(*line.strip().split())
    elif line[0] == GFATypes.L:
        return GFALine(*line.strip().split())

def parse_gfa(fs, as_fasta=False):
    for line in fs:
        if line[0] == GFATypes.H:
            continue
        if not as_fasta:
            yield get_gfa_element(line)
        else:
            e = get_gfa_element(line)
            if e.type == GFATypes.S:
                yield Fasta(e.nid, e.seq)



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
        if fa_cls == CopanNodeSeq:
            yield fa_cls(*hdr.split(":"), ''.join(seq), *args, **kwargs)
        else:
            yield fa_cls(hdr, ''.join(seq), *args, **kwargs)
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

    def __repr__(self):
        return f'Fasta(hdr={self.hdr}, seq={self.seq[:10]})'

    def __len__(self):
        return len(self.seq)


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
        self.node_id, self.sample_id, self.contig_id, self.start_pos, self.end_pos = hdr.split(':')
        self.start_pos = int(float(self.start_pos))
        self.end_pos = int(float(self.end_pos))

class CopanNodeSeq(Fasta):
    def __init__(self, nid, node_name, sid, cid, contig_name, lp, rp, ori, is_tag, seq, fold=80):
        self.nid = nid
        self.nn = node_name
        self.sid = sid
        self.cid = cid
        self.lp = lp
        self.rp = rp
        self.ori = ori
        self.cn = contig_name
        self.is_tag = is_tag
        hdr = f'{nid}:{nn}:{sid}:{cid}:{contig_name}:{lp}:{rp}:{ori}:{is_tag}'
        super().__init__(hdr, seq, fold)
        print(hdr)
        sys.exit()

class PEExtContig(Fasta):
    def __init__(self, hdr, seq, fold=80):
        super().__init__(hdr, seq, fold)
        self.sample_id, self.contig_id, self.bwd_ext, self.fwd_ext = hdr.split(':')
        self.bwd_ext = int(self.bwd_ext)
        self.fwd_ext = int(self.fwd_ext)

    def __key(self):
        return self.seq

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, PEExtContig):
            return self.__key() == other.__key()
        return NotImplemented

# random access fastq
class RandomAccessFastq:
    def __init__(self, fq_file):
        self.file = fq_file
        self.fs = open(self.file, 'rb')
        self.fq_starts = self.get_fq_starts()

    def at(self, idx):
        self.fs.seek(0)
        assert idx < len(self.fq_starts)-1, Exception('There are only {len(self.fq_starts)-1} elements in the file.')
        # move it to the start if the idx'th element
        self.fs.seek(self.fq_starts[idx]) 
        # write the idx'th element
        return self.fs.read(self.fq_starts[idx+1] - self.fq_starts[idx]-1).decode('utf-8')
        

    def get_fq_starts(self):
        starts = [self.fs.tell()]
        self.fs.readline()
        cur = self.fs.tell()
        while starts[-1] != cur:
            starts.append(cur)
            self.fs.readline()
            cur = self.fs.tell()
        return starts[::4]


if __name__ == '__main__':
    raf = RandomAccessFastq('dump')
    print(raf.at(0))
    print(raf.at(4))
    print(raf.at(2))
    print(raf.at(7))



