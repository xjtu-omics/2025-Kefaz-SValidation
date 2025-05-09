from src.inputs.read_op import reverse_seq
import random
import pysam

def fetch_ref_seq(ref_file, chrom, start, end):
    """
    fetch ref seq
    :param ref_file:
    :param chrom:
    :param start:
    :param end:
    :return:
    """
    # -1: convert 1-based to 0-based; +1: convert [ , ) to [ , ]
    return pysam.FastaFile(ref_file).fetch(chrom, start - 1, end - 1 + 1)


def alter_ref_seq_by_sv(ref_seq, ref_start, ref_end, sv):
    """
    alter ref seq by sv, add sv into ref, designed for non-SVision tools, such as Sniffles and cuteSV
    :param ref_seq:
    :param ref_start:
    :param ref_end:
    :param sv:
    :return:
    """

    # # STEP:
    sv_type = sv.origin_type

    altered_refs = []

    # # STEP: normalize cords, make them starting at 0
    sv_start = sv.origin_start - ref_start
    sv_end = sv.origin_end - ref_start
    ref_start -= ref_start
    ref_end -= ref_start

    # # STEP: alter ref seq by different SV types
    if sv_type == "DEL":
        altered_refs.append(ref_seq[: sv_start] + ref_seq[sv_end + 1: ])

    elif sv_type == "INS":
        inserted_seq = sv.origin_read_bases
        altered_refs.append(ref_seq[: sv_start] + inserted_seq + ref_seq[sv_start: ])

    elif sv_type == "INV":
        altered_refs.append(ref_seq[: sv_start] + reverse_seq(ref_seq[sv_start: sv_end + 1]) + ref_seq[sv_end + 1: ])

    elif sv_type == "DUP":
        # # for DUP, there are many conditions since there are no detialed breakpoints, such as tDUP, dDUP
        # # so, we need generate more than 1 altered ref

        # # generate tDUP ref
        duplicated_seq = ref_seq[sv_start: sv_end + 1]
        altered_refs.append(ref_seq[: sv_end + 1] + duplicated_seq + ref_seq[sv_end + 1: ])

        # # generate dDUP ref (from left breakpoint and inserted into right breakpoint)
        # # we use various dup size to simulated the dup length
        dup_size_gradients = [100, 200, 500, 1000]
        for dup_size in dup_size_gradients:
            duplicated_seq = ref_seq[sv_end - dup_size: sv_end + 1]
            altered_refs.append(ref_seq[: sv_start] + duplicated_seq + ref_seq[sv_start: ])

        # # generate dDUP ref (from right breakpoint and inserted into left breakpoint)
        for dup_size in dup_size_gradients:
            duplicated_seq = ref_seq[sv_start: sv_start + dup_size]
            altered_refs.append(ref_seq[: sv_end] + duplicated_seq + ref_seq[sv_end: ])

    elif "/" in sv_type:
        pass

    return altered_refs


def alter_ref_seq_by_simulated_sv(ref_seq, ref_start, ref_end, sv):
    """

    :param ref_seq:
    :param ref_start:
    :param ref_end:
    :param sv:
    :return:
    """

    valid_bases = 'ACGT'

    # # STEP: parse detailed bkps and types
    sv_start = sv.origin_start - ref_start
    sv_end = sv.origin_end - ref_start

    # # STEP: sort bkps, then we can start from the last bkp
    detailed_bkps = sorted(sv.origin_detailed_bkps, key=lambda x: x.insert_pos, reverse=True)
    # print([bkp.to_string() for bkp in detailed_bkps])

    # # STEP: generate altered ref seq by various types
    altered_ref = ref_seq
    for bkp in detailed_bkps:

        bkp_type, bkp_len, bkp_start, bkp_end, insert_pos = bkp.type, bkp.length, bkp.start - ref_start, bkp.end - ref_start, bkp.insert_pos - ref_start

        if bkp_type == "DEL":
            altered_ref = altered_ref[: bkp_start] + altered_ref[bkp_end + 1:]

        elif bkp_type == "INS":
            inserted_seq = "".join(random.choices(valid_bases, k=bkp_len))

            altered_ref = altered_ref[: bkp_start] + inserted_seq + altered_ref[bkp_start: ]

        elif bkp_type == "INV":
            reversed_seq = reverse_seq(altered_ref[bkp_start: bkp_end + 1])

            altered_ref = altered_ref[: bkp_start] + reversed_seq + altered_ref[bkp_end + 1:]

        elif bkp_type == "tDUP" or bkp_type == "dDUP":
            duplicated_seq = ref_seq[bkp_start: bkp_end + 1]

            altered_ref = altered_ref[: insert_pos] + duplicated_seq + altered_ref[insert_pos: ]

        elif bkp_type == 'itDUP' or bkp_type == "idDUP":
            duplicated_seq = reverse_seq(ref_seq[bkp_start: bkp_end + 1])
            altered_ref = altered_ref[: insert_pos] + duplicated_seq + altered_ref[insert_pos: ]

    altered_refs = [altered_ref]

    return altered_refs


def alter_ref_seq_by_detailed_sv(ref_seq, ref_start, ref_end, sv):
    """
    alter ref seq by sv, add sv into ref, this is design for SVision, whose output includes detailed breakpoints
    :param ref_seq:
    :param ref_start:
    :param ref_end:
    :param sv:
    :return:
    """

    valid_bases = 'ACGT'

    # # STEP: parse detailed bkps and types
    sv_start = sv.origin_start - ref_start
    sv_end = sv.origin_end - ref_start

    # # STEP: sort bkps, then we can start from the last bkp
    detailed_bkps = sorted(sv.origin_detailed_bkps, key=lambda x: x.start, reverse=True)

    print([bkp.to_string() for bkp in detailed_bkps])

    # # STEP: generate altered ref seq by various types
    altered_ref = ref_seq

    for bkp in detailed_bkps:
        bkp_type, bkp_len, bkp_start, bkp_end, insert_pos = bkp.type, bkp.length, bkp.start - ref_start, bkp.end - ref_start, bkp.insert_pos - ref_start

        if bkp_type == "DEL":
            altered_ref = altered_ref[: bkp_start] + altered_ref[bkp_end + 1:]

        elif bkp_type == "INS":
            inserted_seq = "".join(random.choices(valid_bases, k=bkp_len))

            altered_ref = altered_ref[: bkp_start] + inserted_seq + altered_ref[bkp_start: ]

        elif bkp_type == "INV":
            reversed_seq = reverse_seq(altered_ref[bkp_start: bkp_end + 1])

            altered_ref = altered_ref[: bkp_start] + reversed_seq + altered_ref[bkp_end + 1:]

        elif bkp_type == "tDUP" or bkp_type == "DUP":
            duplicated_seq = altered_ref[bkp_start: bkp_end + 1]

            altered_ref = altered_ref[: insert_pos] + duplicated_seq + altered_ref[insert_pos: ]

            # duplicated_seq = altered_ref[bkp_start: bkp_end + 1]
            #
            # altered_ref = altered_ref[: bkp_end + 1] + duplicated_seq + altered_ref[bkp_end + 1:]
        # elif bkp_type == "DUP":
        #     duplicated_seq = altered_ref[bkp_start: bkp_end + 1]
        #
        #     altered_ref = altered_ref[: insert_pos] + duplicated_seq + altered_ref[insert_pos: ]

        elif bkp_type == "invDUP":
            duplicated_seq = reverse_seq((altered_ref[bkp_start: bkp_end + 1]))

            altered_ref = altered_ref[: insert_pos] + duplicated_seq + altered_ref[insert_pos: ]

        else:
            print("SV subtype must in [DEL, INS, INV, tDUP, DUP, invDUP]")
            exit()

    altered_refs = [altered_ref]

    return altered_refs