import numpy as np
import pandas as pd
import subprocess
from op_dotplot_gepard import Dotplot_Gepard
import pysam
#sys.path.insert(0, r"D:\SV_validation\v6.3_new_train\src")
#sys.path.append('/mnt/d/SV_validation/v6.3_new_train/src')
def get_parent_path(path, level):
    for _ in range(level):
        path = os.path.dirname(path)
    return path



from src.inputs.read_op import generate_align_obj_from_tag, generate_read_obj_from_aligns
from src.inputs.ref_op import alter_ref_seq_by_simulated_sv, fetch_ref_seq
import os

# class Dotplot:
#     def __init__(self, type, read_name, ref_seq, read_seq, zoom, origin_out_prefix, options):
#
#         self.read_name = read_name
#
#         self.type = type
#         self.ref_seq = ref_seq
#         self.read_seq = read_seq
#
#         self.word_size = options.word_size
#
#         self.zoom = zoom
#
#         self.origin_out_prefix = origin_out_prefix
#
#         self.origin_dotplot_file = self.origin_out_prefix + "_matrix.h5"
#         self.clear_dotplot_file = None
#
#         # self.ref_seq_start = int(origin_out_prefix.split("/")[-1].split(".")[5].split("_")[1])
#         # self.ref_seq_end = int(origin_out_prefix.split("/")[-1].split(".")[5].split("_")[2])
#         self.ref_seq_start = int(origin_out_prefix.split("/")[-1].split(".")[2].split("_")[1])
#         self.ref_seq_end = int(origin_out_prefix.split("/")[-1].split(".")[2].split("_")[2])
#         self.create_matrix()
#
#     def set_clear_dotplot_file(self, file_path):
#         self.clear_dotplot_file = file_path
#
#     def create_matrix(self):
#         """
#         create dot matrix by gepard
#         :param work_dir:
#         :return:
#         """
#         # # STEP: output seq to file
#         read_seq_out_file = self.origin_out_prefix + ".read_seq.fa"
#         ref_seq_out_file = self.origin_out_prefix + ".ref_seq.fa"
#
#         with open(ref_seq_out_file, 'w') as fout:
#             fout.write(">ref\n")
#             fout.write("{}\n".format(self.ref_seq))
#
#         with open(read_seq_out_file, 'w') as fout:
#             fout.write(">read\n")
#             fout.write("{}\n".format(self.read_seq))
#
#         # # STEP: use gepard to create dot matrix
#         ref2read_out_file = self.origin_out_prefix + ".matrix"
#
#         # # STEP: run gepard
#         cmd = "java -cp ./utils/gepard/Gepard_altered.jar org.gepard.client.cmdline.CommandLine -seq {} {} -matrix ./utils/gepard/edna.mat -outfile {} -silent -word {} -zoom {}".format(ref_seq_out_file, read_seq_out_file, ref2read_out_file,
#                                                                                                                                                                                          self.word_size, self.zoom)
#         cmd_run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=-1)
#         cmd_run.wait()
#         # os.system("java -cp ./utils/gepard/Gepard_altered.jar org.gepard.client.cmdline.CommandLine -seq {} {} -matrix ./utils/gepard/edna.mat -outfile {} -silent -word {} -zoom {}".format(ref_seq_out_file, read_seq_out_file, ref2read_out_file, self.word_size, self.zoom))
#
#         # # STEP: loading gepard's output as  matrix
#         self.ref2read_dot_matrix = pd.read_csv(ref2read_out_file, sep='\t', header=None).fillna(0)
#
#         # # STEP: normalize the origin matrix
#         self.ref2read_dot_matrix_norm = 255 * abs(self.ref2read_dot_matrix - self.ref2read_dot_matrix.values.max()) / (self.ref2read_dot_matrix.values.max() - self.ref2read_dot_matrix.values.min())
#
#         self.ref2read_dot_matrix_norm = self.ref2read_dot_matrix_norm.round()
#
#         # # STEP: output three channel img
#         self.img_height, self.img_width = np.shape(self.ref2read_dot_matrix_norm)
#
#         channel_values = self.ref2read_dot_matrix_norm.values
#
#         self.img = np.ones((self.img_height, self.img_width, 3))
#
#         self.img[:, :, 0] = channel_values
#         self.img[:, :, 1] = channel_values
#         self.img[:, :, 2] = channel_values
#
#         os.remove(ref_seq_out_file)
#         os.remove(read_seq_out_file)
#         os.remove(ref2read_out_file)
#
#     def output(self):
#         cv.imwrite(self.origin_out_prefix + ".origin_dotplots.png", self.img)



def generate_dotplots(sv, dotplot_out_path, options, mode="test", visual=False, exten=1, reads_limit=True):
    """
    parse partial bam from origin bam file by the SV's start and end
    model choose from [train, test]
    """
    #向两边扩展SV长度的1.5-2倍
    sv_length = sv.origin_length
    if sv_length < 100:
        if exten == 1:
            print(f"extension is mode '{exten}' 1.5 length")
            extension = min(int(sv_length * 1.5), 200)
            # extension = 2500
        elif exten == 2:
            print(f"extension is mode '{exten}' 1.2 length")
            extension = min(int(sv_length * 1.2), 200)
        elif exten == 3:
            extension = 200
        elif exten == 4:
            extension = 500
        elif exten == 5:
            extension = 1000

    elif 100 <= sv_length < 200:
        if exten == 1:
            extension = min(int(sv_length * 1.5), 300)
        elif exten == 2:
            extension = min(int(sv_length * 2.0), 400)
        elif exten == 3:
            extension = 500
        elif exten == 4:
            extension = 1000
        elif exten == 5:
            extension = 2000

    elif  200 <= sv_length < 500:
        if exten == 1:
            extension = min(int(sv_length * 1.5), 1000)
        elif exten == 2:
            extension = min(int(sv_length * 2), 1000)
        elif exten == 3:
            extension = 1500
        elif exten == 4:
            extension = 3000
        elif exten == 5:
            extension = 5000

    else:
        if exten == 1:
            extension = min(int(sv_length * 2), 1000)
        elif exten == 2:
            extension = min(int(sv_length * 3), 2000)
        elif exten == 3:
            extension = 5000
        elif exten == 4:
            extension = 8000
        elif exten == 5:
            extension = 10000

    # extension = int(sv_length * 1.5)
    # extension = 500
    print(f"extension {extension}")


    expected_start = sv.origin_start - extension
    expected_end = sv.origin_end + extension

    # # STEP: fetch partial bam
    bam_file = pysam.AlignmentFile(options.bam_path)

    processed_reads = []

    altref2read_dotplots = []
    ref2read_dotplots = []

    partial_bam = bam_file.fetch(sv.origin_chrom, expected_start, expected_end)
    # bam_list = list(partial_bam)
    # print(f"align number {len(bam_list)}")
    # return
    for i, align in enumerate(partial_bam):
        #在这里添加截断，对于长度在1000bp以上的，限制最大的align 数量为50，超过50的reads数量不再进行后续的dotplot创建操作
        # if sv.origin_length > 1000 and i > 50:
        if reads_limit and i > 50:
            print(f"sv {sv.id} has analysed more than 50 reads, killed dotplot creating")
            break
        # # STEP: only keep supporting reads
        read_name = align.qname

        if read_name in processed_reads:
            continue

        processed_reads.append(read_name)

        # # STEP: refine read name, convert m54329U_190701_222759/11403722/ccs to m54329U_190701_222759_11403722_ccs.
        read_name = read_name.replace("/", "_")
        current_read_seq = align.query_sequence

        # # STEP: collect all inputs tags
        current_align_tags = [align.reference_name, align.reference_start, "-" if align.is_reverse else "+", align.cigarstring, align.mapq, "NM"]
        current_align = generate_align_obj_from_tag(current_align_tags, type="SA" if align.is_supplementary else "PM")

        try:
            other_align_tags = [tag.split(",") for tag in align.get_tag("SA").split(";") if tag != ""]
            other_aligns = [generate_align_obj_from_tag(tag) for tag in other_align_tags]
        except KeyError:
            other_aligns = []

        # # STEP: generate read obj
        read = generate_read_obj_from_aligns(read_name, current_align, current_read_seq, other_aligns)

        # # STEP TMP: apply filters
        include_chroms = []
        for tmp_align in read.included_aligns:
            include_chroms.append(tmp_align.ref_chrom)

        if len(set(include_chroms)) != 1:
            continue

        # # STEP: cut read by expected ref cords
        read.cur_read_by_ref_cords(expected_start, expected_end)

        try:
            cutted_ref_seq = fetch_ref_seq(options.ref_path, read.ref_chrom, read.cutted_ref_start, read.cutted_ref_end)
            cutted_read_seq = read.cutted_seq
            ref_start = read.cutted_ref_start
            ref_end = read.cutted_ref_end
            # print("ref length:{}".format(len(cutted_ref_seq)))
            # print("read length:{}".format(len(cutted_read_seq)))
            # zoom = int((read.cutted_ref_end - read.cutted_ref_start) / options.pix2pix_img_size)
            # print(f"zoom {zoom}")
            zoom = 1
            # print(f"zoom {zoom}")
            name = "{}.{}.{}.{}.{}_{}_{}.ref2read".format(sv.origin_chrom, sv.id, read_name, exten, read.ref_chrom, read.cutted_ref_start, read.cutted_ref_end)

            # out_prefix = os.path.join(dotplot_out_path, "{}.{}.{}.{}_{}_{}.ref2read".format(sv.origin_chrom, sv.id, read_name, read.ref_chrom, read.cutted_ref_start, read.cutted_ref_end))

            if mode == "train":
                pass

            else:
                # # STEP: generate ref2read_dotplot
                dotplot = Dotplot_Gepard(name, cutted_ref_seq, cutted_read_seq, zoom, 10, dotplot_out_path, ref_start, ref_end, visual)
                # ref2read_dotplot = Dotplot_Gepard("ref2read", read_name, cutted_ref_seq, cutted_read_seq, zoom, out_prefix, options)
                # ref2read_dotplot.output()

                # ref2read_dotplot_resize = custom_resize(ref2read_dotplot.img, re_size=options.pix2pix_img_size)

                # # output and save to list
                ref2read_dotplots.append(dotplot)

        except:
            out_prefix = os.path.join(dotplot_out_path, "{}.{}.{}.{}.{}_{}_{}.ref2read".format(sv.origin_chrom, sv.id, read_name, exten, read.ref_chrom, read.cutted_ref_start, read.cutted_ref_end))
            # out_prefix = os.path.join(dotplot_out_path, "{}.{}.{}_{}_{}".format(sv.id, read_name, read.ref_chrom, read.cutted_ref_start, read.cutted_ref_end))
            os.system("rm {}*".format(out_prefix))


    # print("ref2read_dotplots_num:{}".format(len(ref2read_dotplots)))
    return ref2read_dotplots
