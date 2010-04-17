from itertools import groupby

def grouper(qstr):
    """
    parse the mismatch portion of the samline followign the MD tag
    """
    g = ["".join(li) for p, li in groupby(qstr, str.isdigit)]
    for item in g:
        if item.isdigit(): yield int(item)
        else: yield item

class SamLine(object):
    __slots__ = ("qname", "read_id", "flag", "ref", "ref_pos", "seq",
                 "fields", "mismatches", "mismatch_locs", "line")
    def __init__(self, fh):
        line = fh.readline().strip()
        if not line:
            self.read_id = None
            return None
        self.line = line
        line = line.split("\t")
        self.qname = self.read_id = int(line[0])
        self.flag = int(line[1])
        self.ref = line[2]
        # 1-based
        self.ref_pos = int(line[3])
        """ unused for now.
        self.mapq = int(line[4])
        self.cigar = line[5]
        self.mrnm = line[6]
        self.mpos = line[7]
        self.isize = line[8]
        self.qual = line[10]
        """
        self.seq = line[9]
        self.fields = line[11:]
        self.parse_fields()

    def parse_fields(self):
        for field in self.fields:
            name, type, val = field.split(":") 
            if name == "NM": 
                self.mismatches = int(val) 
                continue
            if name != "MD": continue
            li = list(grouper(val))
            self.mismatch_locs = {}
            if li == []: return
            # this is 0-based
            if len(li) % 2: li = li[:-1]
            for i in range(0, len(li), 2):
                self.mismatch_locs[li[i]] = li[i + 1]
            

if __name__ == "__main__":
    data = """\
5	0	3	18765064	255	76M	*	0	0	TGATTTTGAAGTGATTTTTTAAATGGTTGTTTTTTTTTTTGTTTTTGGTTAAATATATTAGAGTTGTTTTAGTGAG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:76	NM:i:0
6	0	c	56152	255	76M	*	0	0	ATATTATAATTTGGTGGAGGAATTTTAGGTTATTTTTGGGGAAATGTATTGGGTGTTGTAGTTAATTGAGTAGTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:0G75	NM:i:1
9	0	m	251933	255	76M	*	0	0	TGTTTATGAATTTGTTTTATTTGGTGGTTAATAAGTGTTTTTTTGGTGTGGTTTGATTTGATGTAGAGTTTTTATT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:2	MD:Z:63A11G0	NM:i:2"""
    for line in data.split("\n"):
        s = SamLine(line)
        print s.mismatch_locs
