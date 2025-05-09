import os
from src.plots.line_detect import lsd_line_detection
from src.plots.line_to_segment import convert_lines_to_segments


def valid_dotplots_for_each_sv(sv, origin_dotpot_out_path, clear_dotpot_out_path, segment_seq_out_path, options):
    """
    valid clear dotplots by detecting and projecting line from them
    """
    ref2read_dotplots = sv.ref2read_dotplots
    altref2read_dotplots = sv.altref2read_dotplots

    # # STEP: deal with alter_ref2read, mainly for validation
    valid_right_dotplot_index = []
    for dotplot_index in range(len(altref2read_dotplots)):

        dotplot = altref2read_dotplots[dotplot_index]

        dotplot.set_clear_dotplot_file(dotplot.origin_out_prefix.replace(origin_dotpot_out_path, clear_dotpot_out_path + "/images") + ".dotplot-outputs.png")

        if not os.path.exists(dotplot.clear_dotplot_file):
            print("[WARNING]: no corresponding cleared dotplot")
            continue

        # # detect lines
        lines = sorted(lsd_line_detection(dotplot.clear_dotplot_file), key=lambda x: x.read_start)

        # # search gaps between lines
        if len(lines) == 0:
            print("[WARNING]: no line detected")
        elif len(lines) == 1:
            valid_right_dotplot_index.append(dotplot_index)
        else:
            # # detect non-linear aps between lines
            valid_wrong_flag = False
            for line_index in range(len(lines) - 1):
                cur_line = lines[line_index]
                next_line = lines[line_index + 1]

                gap_on_ref = next_line.ref_start - cur_line.ref_end
                gap_on_read = next_line.read_start - cur_line.read_end

                # # find non-linear gaps
                if abs(gap_on_read - gap_on_ref) > 10:
                    valid_wrong_flag = True

            if valid_wrong_flag is False:
                valid_right_dotplot_index.append(dotplot_index)

    print(valid_right_dotplot_index)

    # # STEP: deal with ref2read, mainly for refinement
    for dotplot in ref2read_dotplots:
        dotplot.set_clear_dotplot_file(dotplot.origin_out_prefix.replace(origin_dotpot_out_path, clear_dotpot_out_path + "/images") + ".dotplot-outputs.png")

        print(dotplot.clear_dotplot_file)
        if not os.path.exists(dotplot.clear_dotplot_file):
            print("[WARNING]: no corresponding cleared dotplot")
            continue

        # # STEP: perform line detection
        lines = sorted(lsd_line_detection(dotplot.clear_dotplot_file), key=lambda x: x.read_start)

        for line in lines:
            print(line.to_string())

        segment_seq_out_prefix = dotplot.origin_out_prefix.replace(origin_dotpot_out_path, segment_seq_out_path)
        # # STEP: convert lines to origin segments and do single base level adjust
        segments = convert_lines_to_segments(lines, dotplot.ref_seq_start, dotplot.ref_seq_end, dotplot.ref_seq, dotplot.read_seq, segment_seq_out_prefix, options)




