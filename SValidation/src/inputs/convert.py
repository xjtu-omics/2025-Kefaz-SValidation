import pysam
from src.inputs.cigar_op import calculate_ref_and_read_end_by_cigar, cigar_to_list


def convert_svision_sim_bed_to_vcf(bed_path):

    out_vcf = open(bed_path.replace(".bed", ".vcf"), "w")

    for line in open(bed_path):
        line_split = line.strip().split("\t")

        sv_chrom, sv_start, sv_end = line_split[0], int(line_split[1]), int(line_split[2])

        sv_type = line_split[3]
        sv_details = line_split[4]

        sv_type_reformat = []

        for type in sv_type.split(";"):
            if type == 'deletion':
                sv_type_reformat.append("DEL")
            elif type == 'insertion':
                sv_type_reformat.append("INS")
            elif type == "inversion":
                sv_type_reformat.append("INV")
            elif type == "tandem duplication":
                sv_type_reformat.append("tDUP")
            elif type == "dispersed duplication":
                sv_type_reformat.append("dDUP")
            elif type == "dispersed inverted duplication":
                sv_type_reformat.append("idDUP")
            elif type == "inverted tandem duplication":
                sv_type_reformat.append("itDUP")
            else:
                print("error type", type)
                exit()

        ref_base, alt_base = "N", "N"
        if ";" not in sv_type:
            if sv_type in ["dispersed duplication", "dispersed inverted duplication"]:
                insert_pos = int(sv_details.split(":")[2])

                bkp_start = sv_start
                bkp_end = sv_end

                sv_start = min([sv_start, sv_end, insert_pos])
                sv_end = max([sv_start, sv_end, insert_pos])

            else:
                if sv_type == "insertion":
                    alt_base = sv_details

                if sv_type in ["tandem duplication", "inverted tandem duplication"]:
                    insert_pos = sv_end
                else:
                    insert_pos = sv_start
                bkp_start = sv_start
                bkp_end = sv_end

            bkp = "{}:{}_{}_{}".format("+".join(sv_type_reformat), bkp_start, bkp_end, insert_pos)
            bkps = [bkp]
            out_vcf.write("{}\t{}\t0\t{}\t{}\t.\t.\tEND={};SVLEN={};SVTYPE={};BKPS={}\n".format(sv_chrom, sv_start, ref_base, alt_base, sv_end, sv_end - sv_start, "+".join(sv_type_reformat), ",".join(bkps)))

        else:
            sv_details = sv_details.split("><")

            if sv_type == "tandem duplication;insertion":
                bkps = []
                first_tdup = sv_details[0].replace("<", "")
                first_tdup = first_tdup.replace(">", "")
                first_tdup = first_tdup.split(",")

                bkp_type = "tDUP"
                bkp_start = int(first_tdup[1])
                bkp_end = int(first_tdup[2])

                tdup_len = bkp_end - bkp_start

                bkp = "{}:{}_{}_{}".format(bkp_type, bkp_start, bkp_end, bkp_end)
                bkps.append(bkp)

                bkp = "{}:{}_{}_{}".format(bkp_type, bkp_start, bkp_end, bkp_end + 1)
                bkps.append(bkp)

                sv_type_reformat = ["tDUP", "tDUP"]

            elif sv_type in ["dispersed inverted duplication;deletion", "dispersed duplication;deletion"]:
                bkps = []

                del_detail = sv_details[1]
                del_detail = del_detail.replace("<", "")
                del_detail = del_detail.replace(">", "")
                del_detail = del_detail.split(",")

                del_type = "DEL"
                del_start = int(del_detail[1])
                del_end = int(del_detail[2])


                ddup_detail = sv_details[0]
                ddup_detail = ddup_detail.replace("<", "")
                ddup_detail = ddup_detail.replace(">", "")
                ddup_detail = ddup_detail.split(",")
                ddup_type = "dDUP" if ddup_detail[0] == "dispersed duplication" else "idDUP"
                ddup_start = int(ddup_detail[1])
                ddup_end = int(ddup_detail[2])
                ddup_insert_pos = del_end + 1

                bkp = "{}:{}_{}_{}".format(ddup_type, ddup_start, ddup_end, ddup_insert_pos)
                bkps.append(bkp)

                bkp = "{}:{}_{}_{}".format(del_type, del_start, del_end, del_start)
                bkps.append(bkp)

            else:
                bkps = []


                for detail in sv_details:
                    detail = detail.replace("<", "")
                    detail = detail.replace(">", "")
                    detail = detail.split(",")

                    bkp_type = detail[0]
                    bkp_start = int(detail[1])
                    bkp_end = int(detail[2])

                    if bkp_type == 'deletion':
                        bkp_type = "DEL"
                    elif bkp_type == 'insertion':
                        bkp_type = "INS"
                    elif bkp_type == "inversion":
                        bkp_type = "INV"
                    elif bkp_type == "tandem duplication":
                        bkp_type = "tDUP"
                    elif bkp_type == "dispersed duplication":
                        bkp_type = "dDUP"
                    elif bkp_type == "dispersed inverted duplication":
                        bkp_type = "idDUP"
                    elif bkp_type == "inverted tandem duplication":
                        bkp_type = "itDUP"
                    else:
                        print("error bkp type", bkp_type)
                        exit()

                    if bkp_type in ["tDUP", "itDUP"]:
                        bkp_insert_pos = bkp_end
                    elif bkp_type in ["dDUP", "idDUP"]:
                        if bkp_start == sv_start:
                            bkp_insert_pos = sv_end
                        elif bkp_end == sv_end:
                            bkp_insert_pos = sv_start - 1
                        else:
                            bkp_insert_pos = -1
                            print("[Error], no insert pos")
                    else:
                        bkp_insert_pos = bkp_start

                    bkp = "{}:{}_{}_{}".format(bkp_type, bkp_start, bkp_end, bkp_insert_pos)
                    bkps.append(bkp)

            out_vcf.write("{}\t{}\t0\t{}\t{}\t.\t.\tEND={};SVLEN={};SVTYPE={};BKPS={}\n".format(sv_chrom, sv_start, ref_base, alt_base, sv_end, sv_end - sv_start,  "+".join(sv_type_reformat), ",".join(bkps)))



def find_supp_reads_for_svision_sim(vcf_path, bam_path):
    """
    bam must be homo
    Parameters
    ----------
    vcf_path
    bam_path

    Returns
    -------

    """

    for record in pysam.VariantFile(vcf_path):
        sv_chrom = record.contig
        sv_start = record.start
        sv_end = record.stop

        supp_reads = []

        # # fetch reads from bam
        for pm_align in pysam.AlignmentFile(bam_path).fetch(sv_chrom, sv_start - 2000, sv_end + 2000):
            if pm_align.is_supplementary:
                continue

            else:

                all_cords = [pm_align.reference_start, pm_align.reference_end]
                try:
                    other_align_tags = [tag.split(",") for tag in pm_align.get_tag("SA").split(";") if tag != ""]
                except KeyError:
                    other_align_tags = []

                for align_tag in other_align_tags:
                    align_chrom = align_tag[0]
                    if align_chrom == sv_chrom:
                        align_start = int(align_tag[1])
                        align_cigar = align_tag[3]
                        cigar_ops, cigar_ops_length = cigar_to_list(align_cigar, rm_clip=True)
                        align_end, _ = calculate_ref_and_read_end_by_cigar(align_start, 0, cigar_ops, cigar_ops_length)

                        all_cords.append(align_start)
                        all_cords.append(align_end)

                if min(all_cords) <= sv_start <= sv_end <= max(all_cords):

                    supp_reads.append(pm_align.qname)


        print(str(record).strip() + ";SUPPORT={};READS={}".format(len(supp_reads), ",".join(supp_reads)))


def split_svision_sim_bed(bed_path, ref_path):
    """
    split the nested bed, since the original file is missed
    Returns
    -------

    """

    ref_file = pysam.FastaFile(ref_path)

    out_bed = open(bed_path.replace(".bed", ".split.bed"), "w")

    for line in open(bed_path):
        line_split = line.strip().split("\t")

        sv_chrom, sv_start, sv_end = line_split[0], int(line_split[1]), int(line_split[2])

        sv_type = line_split[3]
        sv_details = line_split[4]

        if ";" not in sv_type:
            out_bed.write(line.strip() + "\t2\n")

        else:
            sv_details = sv_details.split("><")

            if sv_type == "tandem duplication;insertion":

                first_tdup = sv_details[0].replace("<", "")
                first_tdup = first_tdup.replace(">", "")
                first_tdup = first_tdup.split(",")

                bkp_type = first_tdup[0]
                bkp_start = int(first_tdup[1])
                bkp_end = int(first_tdup[2])

                dup_seq = ref_file.fetch(sv_chrom, bkp_start - 1, bkp_end - 1 + 1)

                if "NN" in dup_seq:
                    continue

                out_bed.write("{}\t{}\t{}\t{}\t2\t2\n".format(sv_chrom, bkp_start, bkp_end, bkp_type))

                out_bed.write("{}\t{}\t{}\t{}\t{}\t2\n".format(sv_chrom, bkp_end + 1, bkp_end + 2, "insertion", dup_seq))

            elif sv_type in ["dispersed inverted duplication;deletion", "dispersed duplication;deletion"]:
                del_detail = sv_details[1]
                del_detail = del_detail.replace("<", "")
                del_detail = del_detail.replace(">", "")
                del_detail = del_detail.split(",")

                del_start = int(del_detail[1])
                del_end = int(del_detail[2])

                ddup_detail = sv_details[0]
                ddup_detail = ddup_detail.replace("<", "")
                ddup_detail = ddup_detail.replace(">", "")
                ddup_detail = ddup_detail.split(",")
                ddup_type = ddup_detail[0]
                ddup_start = int(ddup_detail[1])
                ddup_end = int(ddup_detail[2])
                strand = "forward" if ddup_type == "dispersed duplication" else "reverse"

                bkp_insert_pos = del_end + 2
                out_bed.write("{}\t{}\t{}\t{}\th1:{}:{}:{}\t2\n".format(sv_chrom, ddup_start, ddup_end, ddup_type, sv_chrom, bkp_insert_pos, strand))

                out_bed.write("{}\t{}\t{}\t{}\tNone\t2\n".format(sv_chrom, del_start, del_end, "deletion"))

            else:

                for detail in sv_details:
                    detail = detail.replace("<", "")
                    detail = detail.replace(">", "")
                    detail = detail.split(",")

                    bkp_type = detail[0]
                    bkp_start = int(detail[1])
                    bkp_end = int(detail[2])

                    if bkp_type in ["deletion", "inversion"]:
                        out_bed.write("{}\t{}\t{}\t{}\tNone\t2\n".format(sv_chrom, bkp_start, bkp_end, bkp_type))
                    elif bkp_type in ["tandem duplication", "inverted tandem duplication"]:
                        out_bed.write("{}\t{}\t{}\t{}\t2\t2\n".format(sv_chrom, bkp_start, bkp_end, bkp_type))

                    else:
                        print("ERROR: 1", bkp_type, sv_type)


if __name__ == '__main__':
    # convert_svision_sim_bed_to_vcf("/mnt/d/data/ccs/sim_svision/benchmark_origin_clone0_csv_300.bed")
    # find_supp_reads_for_svision_sim("/mnt/d/data/ccs/sim_svision/benchmark_origin_clone0_csv_300.vcf", "/mnt/d/data/ccs/sim_svision/Sim-homo.ngmlr.srt.bam")

    split_svision_sim_bed("/mnt/d/data/ccs/sim_svision/benchmark_origin_clone0.bed", "/mnt/d/data/ref/grch38/GRCh38.d1.vd1.fa")