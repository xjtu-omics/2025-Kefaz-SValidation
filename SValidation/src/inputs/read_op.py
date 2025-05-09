from src.inputs.cigar_op import cigar_to_list, calculate_ref_and_read_end_by_cigar, calculate_read_stop_by_ref_stop


class Read:
    def __init__(self, read_name, included_aligns, read_seq):

        self.read_name = read_name
        self.included_aligns = sorted(included_aligns, key=lambda x: x.read_start)

        self.included_aligns_num = len(included_aligns)

        self.read_seq = read_seq

        self.read_seq_length = len(read_seq)

        self.ref_chrom = self.included_aligns[0].ref_chrom
        self.ref_start = min([align.ref_start for align in self.included_aligns])
        self.ref_end = max([align.ref_end for align in self.included_aligns])


    def cur_read_by_ref_cords(self, expected_ref_start, expected_ref_end):
        """
        cur read sequence by given ref start and end cords
        :param expected_ref_start:
        :param expected_ref_end:
        :return:
        """

        # # STEP: initialization
        self.cutted_read_start = 0
        self.cutted_read_end = len(self.read_seq) + 1

        self.cutted_seq = self.read_seq
        self.cutted_ref_start = expected_ref_start
        self.cutted_ref_end = expected_ref_end

        # # STEP: process ref end
        # # current ref_end is shorter than expect (all aligns locate inside the exprect end), then we update the final cords
        if self.ref_end < expected_ref_end:
            self.cutted_ref_end = self.ref_end
        else:
            # # traverse starting from the last inputs
            for index in range(self.included_aligns_num - 1, -1, -1):
                cur_align = self.included_aligns[index]

                # # the whole inputs lies at the right of expected end cord
                if cur_align.ref_start > expected_ref_end:
                    continue

                # # the inputs span the expected end cord
                else:
                    self.cutted_read_end = calculate_read_stop_by_ref_stop(expected_ref_end, cur_align.ref_start, cur_align.read_start, cur_align.cigar_str)

                    # print("cutted_read_end", cutted_read_end)
                    self.cutted_seq = self.cutted_seq[: self.cutted_read_end]

                    break

        # # STEP: process ref start
        # # current ref start is shorted than expect (all aligns locate inside the exprect start), then we update the final cords
        if self.ref_start > expected_ref_start:
            self.cutted_ref_start = self.ref_start
        else:
            # # traverse starting from the first inputs
            for index in range(self.included_aligns_num):
                cur_align = self.included_aligns[index]

                # # the whole inputs lies at the left of expected end cord
                if cur_align.ref_end < expected_ref_start:
                    continue
                # # the inputs span the expected end cord
                else:
                    self.cutted_read_start = calculate_read_stop_by_ref_stop(expected_ref_start, cur_align.ref_start, cur_align.read_start, cur_align.cigar_str)

                    self.cutted_seq = self.cutted_seq[self.cutted_read_start: ]

                    break


    def split_and_cut_aligns(self):
        # # STEP: update cutted inputs's info
        aligns_copy = self.included_aligns.copy()

        self.split_aligns = []

        for align in aligns_copy:
            self.split_aligns.extend(split_align_by_cigar(align))

        self.cutted_split_aligns = []

        for align in self.split_aligns:

            # print(inputs.to_string(), self.cutted_read_start, self.cutted_read_end)
            if align.read_end < self.cutted_read_start:
                continue

            if align.read_start > self.cutted_read_end:
                continue

            if self.cutted_read_start < align.read_start < align.read_end < self.cutted_read_end:
                self.cutted_split_aligns.append(align)

            if align.read_start < self.cutted_read_start < align.read_end:
                align.read_start = self.cutted_read_start
                align.ref_start = self.cutted_ref_start
                self.cutted_split_aligns.append(align)

            if align.read_start < self.cutted_read_end < align.read_end:
                align.read_end = self.cutted_read_end
                align.ref_end = self.cutted_ref_end
                self.cutted_split_aligns.append(align)

    def set_origin_plot_info(self):
        pass

    def set_clear_plot_info(self):
        pass

    def to_string(self):
        pass


class Align:
    def __init__(self, ref_chrom, ref_start, ref_end, read_start, read_end, strand, type, align_seq, cigar_str):
        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_start = read_start
        self.read_end = read_end
        self.strand = strand

        self.type = type

        self.align_seq = align_seq

        self.cigar_str = cigar_str

    def to_string(self):
        return "{}-{}-{}-{}-{}-{}-{}".format(self.ref_chrom, self.ref_start, self.ref_end, self.read_start, self.read_end, self.strand, self.type)


def generate_align_obj_from_tag(tag, type="No-type", seq="No-seq"):
    """
    generate inputs object from tags
    :param tag:
    :return:
    """

    # # STEP: extract info from tag
    ref_chrom = tag[0]
    ref_start = int(tag[1])
    strand = tag[2]  #比对方向
    cigar_str = tag[3]

    # # STEP: calculate end cords by cigar
    cigar_ops, cigar_ops_len = cigar_to_list(cigar_str, rm_clip=False)
    left_soft_clip_length = cigar_ops_len[0] if cigar_ops[0] == "S" else 0

    read_start = left_soft_clip_length + 1

    ref_end, read_end = calculate_ref_and_read_end_by_cigar(ref_start, read_start, cigar_ops, cigar_ops_len)

    return Align(ref_chrom, ref_start, ref_end, read_start, read_end, strand, type, seq, cigar_str)


def generate_align_obj_from_bam(align):
    """
    generate inputs object from aligns in BAM file
    :param align:
    :return:
    """

    return Align(align.reference_name, align.reference_start + 1, align.reference_end,
                 align.query_alignment_start + 1, align.query_alignment_end, "-" if align.is_reverse else "+",
                 "SA" if align.is_supplementary else "PM", align.query, align.cigarstring)


def generate_read_obj_from_aligns(read_name, cur_align, cur_read_seq, other_aligns):
    """
    generate read object from tags
    :param cur_tag:
    :param cur_tag_seq:
    :return:
    """

    # # STEP: current inputs is primary
    if cur_align.type == "PM":

        pm_align = cur_align

        read_seq = cur_read_seq

        # # STEP: update and adjust supp aligns's info
        for align in other_aligns:
            align.type = "SA"

            # # re-set read start cord by strand
            ops, ops_len = cigar_to_list(align.cigar_str, rm_clip=False)
            right_soft_clip_length = ops_len[-1] if ops[-1] == "S" else 0

            # # different strand from pm_strand, then we adjust the cord
            if align.strand != pm_align.strand:
                align.read_start = right_soft_clip_length + 1

                # # re-set end cords by cigar
                align.ref_end, align.read_end = calculate_ref_and_read_end_by_cigar(align.ref_start, align.read_start, ops, ops_len)
                # # re-set strand
                align.strand = "-"

            else:
                align.strand = "+"

        pm_align.strand = "+"

    # # STEP: current inputs is supplementary, them pm_align must in other aligns
    else:
        # # STEP: use the longest inputs as primary
        included_aligns_lengths = [align.ref_end - align.ref_start for align in other_aligns]
        pm_index = included_aligns_lengths.index(max(included_aligns_lengths))

        # # STEP: pm_align found
        pm_align = other_aligns[pm_index]
        pm_align.type = "PM"

        # # STEP: re-calculate seq
        if pm_align.strand != cur_align.strand:
            read_seq = reverse_seq(cur_read_seq)
        else:
            read_seq = cur_read_seq

        # # STEP: update and adjust supp aligns's info
        other_aligns[pm_index] = cur_align
        for align in other_aligns:
            align.type = "SA"

            # # re-set read start cord by strand
            ops, ops_len = cigar_to_list(align.cigar_str, rm_clip=False)
            right_soft_clip_length = ops_len[-1] if ops[-1] == "S" else 0

            # # different strand from pm_strand, then we adjust the cord
            if align.strand != pm_align.strand:
                align.read_start = right_soft_clip_length + 1

                # # re-set end cords by cigar
                align.ref_end, align.read_end = calculate_ref_and_read_end_by_cigar(align.ref_start, align.read_start, ops, ops_len)
                # # re-set strand
                align.strand = "-"
            else:
                align.strand = "+"

        pm_align.strand = "+"

    return Read(read_name, [pm_align] + other_aligns, read_seq)


def reverse_seq(seq):
    """
    reverse sequence (reversed complemented seq)
    :return:
    """

    rev_seq = ""

    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == "A":
            rev_seq += "T"
        elif seq[i] == "T":
            rev_seq += "A"
        elif seq[i] == "C":
            rev_seq += "G"
        elif seq[i] == "G":
            rev_seq += "C"
        else:
            rev_seq += seq[i]

    return rev_seq


def split_align_by_cigar(align):
    """
    split inputs by its cigar, cut when meeting I, D
    :return:
    """

    split_aligns = []
    cigar_ops, cigar_ops_length = cigar_to_list(align.cigar_str, rm_clip=True)

    if align.strand == "+":
        ref_pointer = align.ref_start
        read_pointer = align.read_start
    else:
        ref_pointer = align.ref_start
        read_pointer = align.read_end

    for i in range(len(cigar_ops)):
        op = cigar_ops[i]
        op_len = cigar_ops_length[i]

        # if op_len < 20:
        #     if align.strand == "+":
        #         ref_pointer += op_len
        #         read_pointer += op_len
        #     else:
        #         ref_pointer += op_len
        #         read_pointer -= op_len
        #     continue

        if op == "N" or op == "S":
            continue

        elif op == "I":
            if align.strand == "+":
                read_pointer += op_len
            else:
                read_pointer -= op_len

        elif op == "D":
            if align.strand == "+":
                ref_pointer += op_len
            else:
                ref_pointer += op_len

        elif op in ["M", '=']:
            if align.strand == "+":
                split_aligns.append(Align(align.ref_chrom, ref_pointer, ref_pointer + op_len - 1, read_pointer, read_pointer + op_len - 1, align.strand, "{}M".format(op_len), "seq", align.type))
            else:
                split_aligns.append(Align(align.ref_chrom, ref_pointer, ref_pointer + op_len - 1, read_pointer - op_len + 1, read_pointer, align.strand, "{}M".format(op_len), "seq", align.type))

            if align.strand == "+":
                ref_pointer += op_len
                read_pointer += op_len
            else:
                ref_pointer += op_len
                read_pointer -= op_len

        elif op in ["X", "E"]:
            if align.strand == "+":
                ref_pointer += op_len
                read_pointer += op_len
            else:
                ref_pointer += op_len
                read_pointer -= op_len
        else:
            continue
    return split_aligns
