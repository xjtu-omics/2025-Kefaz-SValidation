import re


class CIGAR():
    def __init__(self, op, op_len, ref_bases, alt_bases, chrom, ref_start, ref_end, supp_reads, source):
        self.op = op
        self.op_len = op_len
        self.ref_bases = ref_bases
        self.alt_bases = alt_bases

        self.ref_start = ref_start
        self.ref_end = ref_end

        self.chrom = chrom
        self.supp_reads = supp_reads

        self.source = source

    def to_string(self,):
        return "{}-{}-{}->{}-{}-{}-{}-{}".format(self.op, self.op_len, self.ref_bases, self.alt_bases, self.chrom, self.ref_start, self.ref_start, ",".join(self.supp_reads))


def cigar_to_list(cigar, rm_clip=True):
    """
    convert cigar string to list
    :param cigar:
    :return:
    """

    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    if rm_clip == True:
        if ops[0] == "S" or ops[0] == "H":
            ops = ops[1: ]
            lengths = lengths[1: ]

        if ops[-1] == "S" or ops[-1] == "H":
            ops = ops[: -1]
            lengths = lengths[: -1]

    return ops, lengths


def reverse_cigar_str(cigar_str):
    """
    reverse cigar string
    :return:
    """

    rev_cigar_str = ""

    cigar_ops, cigar_op_lengths = cigar_to_list(cigar_str, rm_clip=True)


    for i in range(len(cigar_ops) - 1, -1, -1):
        rev_cigar_str += "{}{}".format(cigar_ops[i], cigar_op_lengths[i])

    return rev_cigar_str


def calculate_ref_and_read_end_by_cigar(ref_start, read_start, cigar_ops, cigar_ops_length):
    """
    given ref and read start, calculate ref and read end by cigar ops
    :param ref_start:
    :param read_start:
    :param ops:
    :param ops_len:
    :return:
    """

    ref_end = ref_start
    read_end = read_start
    for i in range(len(cigar_ops)):
        op = cigar_ops[i]
        op_len = cigar_ops_length[i]

        if op == "N":
            read_end += op_len

        elif op == "I":
            read_end += op_len

        elif op == "D":
            ref_end += op_len

        elif op in ["M", '=', "X", "E"]:
            ref_end += op_len
            read_end += op_len
        else:
            continue

    return ref_end - 1, read_end - 1


def calculate_read_stop_by_ref_stop(expected_ref_stop, origin_ref_start, origin_read_start, cigar_str):
    """
    calculate read end position by given cigar and expected ref cord
    :param expected_ref_stop:
    :param origin_ref_start:
    :param origin_read_start:
    :return:
    """
    cigar_ops, cigar_ops_len = cigar_to_list(cigar_str, rm_clip=True)

    ref_pointer = origin_ref_start
    read_pointer = origin_read_start

    for i in range(len(cigar_ops)):
        op = cigar_ops[i]
        op_len = cigar_ops_len[i]

        if op == "N":
            read_pointer += op_len

        elif op == "I":
            read_pointer += op_len

        elif op == "D":
            ref_pointer += op_len

        elif op in ["M", '=', "X", "E"]:

            if ref_pointer + op_len > expected_ref_stop:
                diff = expected_ref_stop - ref_pointer
                return read_pointer + diff

            ref_pointer += op_len
            read_pointer += op_len

        else:
            continue


