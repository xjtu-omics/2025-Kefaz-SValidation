import os
from src.plots.line_detect import Line
import pysam
import subprocess

class Segment:
    def __init__(self, seg_id, seg_ref_start, seg_ref_end, seg_read_start, seg_read_end, seg_strand):
        self.id = seg_id
        self.ref_start = seg_ref_start
        self.ref_end = seg_ref_end
        self.read_start = seg_read_start
        self.read_end = seg_read_end
        self.strand = seg_strand

        self.flag = None

    def set_seg_flag(self, flag):
        self.flag = flag

    def update_seg_id(self, new_id):
        self.id = new_id

    def to_string(self):
        return "Segment: {}, ref: {}-{}, read: {}-{}, {}, {}".format(self.id, self.ref_start, self.ref_end, self.read_start, self.read_end, self.strand, self.flag)


def convert_lines_to_segments(lines, ref_start, ref_end, ref_seq, read_seq, segment_seq_out_prefix, options):
    """
    convert lines to segments, and perform single base level adjust

    """

    shrink_ratio = max(len(ref_seq), len(read_seq)) / options.pix2pix_img_size

    segments = []

    for line in lines:
        # seg_ref_start = int(line.ref_start * shrink_ratio + ref_start)
        # seg_ref_end = int(line.ref_end * shrink_ratio + ref_start)
        seg_ref_start = int(line.ref_start * shrink_ratio)
        seg_ref_end = int(line.ref_end * shrink_ratio)

        seg_read_start = int(line.read_start * shrink_ratio)
        seg_read_end = int(line.read_end * shrink_ratio)

        seg_strand = "+"

        segments.append(Segment(line.id, seg_ref_start, seg_ref_end, seg_read_start, seg_read_end, seg_strand))

    print(ref_start, ref_end, shrink_ratio, len(ref_seq), len(read_seq), shrink_ratio)

    for seg in segments:
        print(seg.to_string())

    # single_base_segments_adjust(segments, ref_seq, read_seq, segment_seq_out_prefix)

    return segments


def single_base_segments_adjust(segments, ref_seq, read_seq, segment_seq_out_prefix):
    """
    perform single base level bkp adjustment
    """

    for seg_index in range(len(segments)):
        seg = segments[seg_index]

        # # STEP: for first segment, we only need to deal with its right bkp
        if seg_index == 0:
            single_base_bkp_adjust(seg.ref_end, ref_seq, read_seq, "right", seg.id, segment_seq_out_prefix)

        # # STEP: for last segment, we only need to deal with its left bkp
        elif seg_index == len(segments) - 1:
            single_base_bkp_adjust(seg.ref_start, ref_seq, read_seq, "left", seg.id, segment_seq_out_prefix)

        # # STEP: for other segments, we need to deal with both its left and right bkps
        else:
            single_base_bkp_adjust(seg.ref_start, ref_seq, read_seq, "left", seg.id, segment_seq_out_prefix)
            single_base_bkp_adjust(seg.ref_end, ref_seq, read_seq, "right", seg.id, segment_seq_out_prefix)


def single_base_bkp_adjust(bkp, ref_seq, read_seq, bkp_type, seg_id, segment_seq_out_prefix, extension=200):
    """
    perform single base adjustment for left and right bkps of a segment

    """

    # # STEP: fetch ref and read seq around left bkps
    ref_seq_around_bkp = ref_seq[max(0, bkp - extension): bkp + extension]
    read_seq_around_bkp = read_seq[max(0, bkp - extension): bkp + extension]

    ref_seq_around_bkp_out_file = segment_seq_out_prefix + ".{}.{}.fa".format(seg_id, bkp_type)

    with open(ref_seq_around_bkp_out_file, "w") as fout:
        fout.write(">ref\n")
        fout.write(ref_seq_around_bkp + "\n")
        fout.write(">read\n")
        fout.write(read_seq_around_bkp + "\n")

    match_pairs = run_and_parse_mafft(ref_seq_around_bkp_out_file)

    # print(bkp_type, [m[4] for m in match_pairs])
    # print(bkp_type, match_pairs)

def run_and_parse_mafft(input_file):
    """
    run mafft for input file
    """

    # # STEP: run mafft
    output_file = input_file + ".mafft"

    cmd = "mafft --auto {} > {}".format(input_file, output_file)
    cmd_run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=-1)
    cmd_run.wait()

    # # STEP: fetch mafft results
    mafft_res = pysam.FastaFile(output_file)

    ref_match_res_seq = mafft_res["ref"]
    read_match_res_seq = mafft_res["read"]

    # # STEP: search match pairs
    match_pairs = []

    ref_match_bases = ""
    read_match_bases = ""
    match_start_index = -1
    for index_on_seq in range(len(ref_match_res_seq)):

        cur_ref_base = ref_match_res_seq[index_on_seq]
        cur_read_base = read_match_res_seq[index_on_seq]

        if cur_ref_base == "-" or cur_read_base == "-" or index_on_seq == len(ref_match_res_seq) - 1:
            if match_start_index != -1:
                # # save
                match_end_index = index_on_seq
                match_pairs.append([match_start_index, match_end_index, ref_match_bases, read_match_bases, calculate_match_dice(ref_match_bases, read_match_bases)])

                # # init
                match_start_index = -1
                ref_match_bases = ""
                read_match_bases = ""
        else:
            if match_start_index == -1:
                match_start_index = index_on_seq

            # # update
            ref_match_bases += cur_ref_base
            read_match_bases += cur_read_base

    return match_pairs


def calculate_match_dice(seq1, seq2):
    """
    calculate the match dice (same / total) between two seqs

    """

    match_num = 0
    tot_num = len(seq1)

    for i in range(tot_num):
        # # skip "-"
        if seq1[i] == "-" or seq2[i] == "-" :
            tot_num -= 1
            continue

        if seq1[i] == seq2[i]:
            match_num += 1

    return round(match_num / tot_num, 2)

