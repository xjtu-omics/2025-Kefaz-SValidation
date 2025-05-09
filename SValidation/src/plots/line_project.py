import os
import cv2
import numpy as np
from line_detect import Line
from line_to_segment import Segment
ref_span_flags = ["A", "B", "C", "D", "E", 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']


class PJPoint:
    def __init__(self, object, pos, content):
        self.object = object
        self.pos = pos
        self.content = content

    def to_string(self):
        return 'PJPoint: {0} {1} {2}'.format(self.object, self.pos, ';'.join(self.content))


class Span:
    def __init__(self, ref_span_id, ref_span_start, ref_span_end, included_seg_ids):

        self.span_id = ref_span_id
        self.span_start = ref_span_start
        self.span_end = ref_span_end
        self.included_seg_ids = included_seg_ids

        self.span_flag = ref_span_flags[self.span_id]

    def updata_span_flag(self, flag):
        self.span_flag = flag

    def to_string(self):
        return "id={},{},{}-{} {}".format(self.span_id, self.span_flag, self.span_start, self.span_end, self.included_seg_ids)


class Pattern:
    def __init__(self, ref_start, ref_end, pattern_flag):
        self.start = ref_start
        self.end = ref_end
        self.flag = pattern_flag

    def set_pattern_flag(self, pattern_flag):
        self.flag = pattern_flag


def collect_PJPoints(segments, type):
    # if pj already in pj_sets, then add content to this one
    pj_points = []
    for seg in segments:
        if type == 'ref':
            pj_start = seg.ref_start
            pj_end = seg.ref_end
        elif type == 'read':
            pj_start = seg.read_start
            pj_end = seg.read_end
        else:
            pj_start = 0
            pj_end = 0
            print('[ERROR]: type must be ref or read')
            exit()

        seg_id = seg.id

        # both
        start_pj_point = PJPoint(type, pj_start, ['{0}-{1}'.format(seg_id, 'h')])
        end_pj_point = PJPoint(type, pj_end, ['{0}-{1}'.format(seg_id, 't')])

        # add to ref_projection
        start_found_flag = 0
        end_found_flag = 0

        for i in range(len(pj_points)):
            target_pj_point = pj_points[i]

            # check start cutoff, already in:
            if start_pj_point.pos == target_pj_point.pos:
                target_pj_point.content.append(start_pj_point.content[0])
                start_found_flag = 1

            # check end cutoff, already in:
            if end_pj_point.pos == target_pj_point.pos:
                target_pj_point.content.append(end_pj_point.content[0])
                end_found_flag = 1
        # not found, then:
        if start_found_flag == 0:
            pj_points.append(start_pj_point)
        if end_found_flag == 0:
            pj_points.append(end_pj_point)

    return pj_points


def split_read_by_pjpoints(read_pj, segments):
    """
    split read by pjpoints

    """
    all_read_spans = []
    sorted_read_pj = sorted(read_pj, key=lambda aln: aln.pos)

    read_span_id = 0
    for i in range(len(sorted_read_pj) - 1):
        read_span_start = sorted_read_pj[i].pos
        read_span_end = sorted_read_pj[i + 1].pos

        first_pfj_contents = sorted_read_pj[i].content
        second_pfj_contents = sorted_read_pj[i + 1].content

        # if len(first_pfj_contents) != 1 or len(second_pfj_contents) != 1:
        #     print('[ERROR]: Segments overlap on read')
        #     exit()

        if first_pfj_contents[0][0] == second_pfj_contents[0][0]:
            seg_id = int(first_pfj_contents[0][0])
            included_seg_ids = [seg_id]
            all_read_spans.append(Span(read_span_id, read_span_start, read_span_end, included_seg_ids))
        else:
            insert_span = Span(read_span_id, read_span_start, read_span_end, [])
            insert_span.updata_span_flag("I")
            all_read_spans.append(insert_span)

        read_span_id += 1

    return all_read_spans


def split_ref_by_pjpoints(ref_pj):
    """
    split ref by pj points to ref_spans
    ref span format: [ref_span_start, ref_span_end, included_line_ids, ref_span_flags[ref_span_id]]

    """
    sorted_ref_pj = sorted(ref_pj, key=lambda aln: aln.pos)

    # if sorted_ref_pj[0].pos != 0:
    #     print('ERROR: the first pj point not start at 0')
    #     exit(-1)

    # # STEP: split ref by pj points. ref segments are between each two pj point
    all_ref_spans = []
    ref_span_id = 0
    no_end_seg = []

    for i in range(len(sorted_ref_pj) - 1):
        ref_span_start = sorted_ref_pj[i].pos
        ref_span_end = sorted_ref_pj[i + 1].pos

        first_pfj_contents = sorted_ref_pj[i].content
        second_pfj_contents = sorted_ref_pj[i + 1].content

        included_seg_ids = []  # store line's ids that covered by this ref_span

        for c in first_pfj_contents:
            line_id = int(c[0])

            if 'h' in c:
                target_c = c[0: -1] + 't'
                if target_c in second_pfj_contents:
                    included_seg_ids.append(line_id)
                else:   # if only head found but tail not found, then this seg is a no end seg
                    no_end_seg.append(line_id)

            elif 't' in c:
                # found a tail in first pj point, then remove it from a no end set
                if line_id in no_end_seg:
                    no_end_seg.remove(line_id)
            else:
                pass

        included_seg_ids.extend(no_end_seg)

        all_ref_spans.append(Span(ref_span_id, ref_span_start, ref_span_end, included_seg_ids))
        ref_span_id += 1

    return all_ref_spans


def split_seg_by_ref_spans(ref_spans, segments):
    """
    by the given ref spans, split the origin segments, mainly for duplicated segments
    """

    new_segments = []
    new_seg_id = 0

    # # STEP: for each ref span, split its included segments
    for ref_span in ref_spans:
        # # extract ref span info
        ref_span_start, ref_span_end = ref_span.span_start, ref_span.span_end
        included_seg_ids = ref_span.included_seg_ids
        ref_span_flag = ref_span.span_flag

        ref_span_length = ref_span_end - ref_span_start

        for seg_id in included_seg_ids:
            seg = segments[seg_id]

            seg_read_start = seg.read_start
            seg_ref_start = seg.ref_start

            # # split segments
            diff_on_ref = ref_span_start - seg_ref_start

            if seg.strand == "+":
                new_seg_read_start = seg_read_start + diff_on_ref
            else:
                new_seg_read_start = seg_read_start - diff_on_ref

            # # generate new segment
            new_seg = Segment(seg.id, ref_span_start, ref_span_end, new_seg_read_start, new_seg_read_start + ref_span_length, seg.strand)
            if seg.strand == "+":
                new_seg.set_seg_flag(ref_span_flag)
            else:
                new_seg.set_seg_flag(ref_span_flag + '*')
            new_seg_id += 1

            new_segments.append(new_seg)

    return new_segments


def get_ref_pattern(ref_spans, pass_size):
    I_num = 0
    ref_patterns = []
    for span in ref_spans:
        ref_start = span[0]
        ref_end = span[1]
        flag = span[3]
        # too short, then abandon
        if not ref_end - ref_start > pass_size:
            continue

        if flag == 'I':
            flag = 'I{0}'.format(I_num)
            I_num += 1
        ref_patterns.append(Pattern(ref_start, ref_end, flag))

    return ref_patterns


def get_read_pattern(read_spans, pass_size):
    I_num = 0

    read_patterns = []
    for span in read_spans:
        read_start = span[0]
        read_end = span[1]
        flag = span[3]
        # too short, then abandon
        if not read_end - read_start > pass_size:
            continue
        if flag == 'I':
            flag = 'I{0}'.format(I_num)
            I_num += 1
        read_patterns.append(Pattern(read_start, read_end, flag))

    return read_patterns


def project_to_ref_and_read(segments):
    """
    project lines (segments) to both ref and read axis

    """
    # # ------------do projection on ref-------------
    # # STEP: collect pj points for ref
    ref_pjpoints = collect_PJPoints(segments, 'ref')

    for pj in ref_pjpoints:
        print(pj.to_string())

    # # STEP: split ref by pj points
    all_ref_spans = split_ref_by_pjpoints(ref_pjpoints)
    print('Ref spans: ', [ref_span.to_string() for ref_span in all_ref_spans])

    # # STEP: split origin lines by ref spans
    splited_segments = split_seg_by_ref_spans(all_ref_spans, segments)
    splited_segments = sorted(splited_segments, key=lambda aln: (aln.read_start, aln.read_end))

    # assign new seg id
    for i in range(len(splited_segments)):
        splited_segments[i].update_seg_id(i)

    for segment in splited_segments:
        print(segment.to_string())
    #
    # # # ------------do projection on read-------------
    # # # STEP: split read by pj points
    # read_pjpoints = collect_PJPoints(splited_segments, 'read')
    #
    # for pj in read_pjpoints:
    #     print(pj.to_string())
    #
    # all_read_spans = split_read_by_pjpoints(read_pjpoints, splited_segments)
    # print('Read spans: ', [ref_span.to_string() for ref_span in all_read_spans])
    #
    # exit()
    #
    # ref_pattern = get_ref_pattern(all_ref_spans, 20)
    # read_pattern = get_read_pattern(all_read_spans, 20)
    #
    # ref_pattern.insert(0, Pattern(-1, -1, 'anchor1'))
    # ref_pattern.append(Pattern(-1, -1, 'anchor2'))
    # read_pattern.insert(0, Pattern(-1, -1, 'anchor1'))
    # read_pattern.append(Pattern(-1, -1, 'anchor2'))
    #
    # return ref_pattern, read_pattern, splited_segments


if __name__ == '__main__':

    segment1 = Segment(0, 0, 94, 33, 123, "+")

    segment2 = Segment(1, 75, 94, 128, 147, "+")

    segment3 = Segment(2, 179, 254, 148, 219, "+")

    segments = [segment1, segment2, segment3]

    project_to_ref_and_read(segments)

