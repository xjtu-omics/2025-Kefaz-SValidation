import os


class SV():

    def __init__(self, origin_id, origin_chrom, origin_start, origin_end, origin_type, sv_reptype, origin_ref_bases, origin_read_bases, origin_supp_reads, origin_detailed_bkps):

        self.id = origin_id

        # self.origin_chrom = "chr" + origin_chrom
        self.origin_chrom = origin_chrom
        self.origin_start = origin_start
        self.origin_end = origin_end
        self.origin_type = origin_type
        self.reptype = sv_reptype
        self.origin_ref_bases = origin_ref_bases
        self.origin_read_bases = origin_read_bases

        self.origin_detailed_bkps = origin_detailed_bkps

        self.origin_supp_reads = origin_supp_reads

        self.origin_length = self.origin_end - self.origin_start

        self.refined_chrom = None
        self.refined_start = None
        self.refined_end = None
        self.refined_type = None
        self.clear_dotplots = None
        self.ref2read_dotplots = None
        self.altref2read_dotplots = None

    # def set_dotplots(self, ref2read_dotplots, altref2read_dotplots):
    #     self.ref2read_dotplots = ref2read_dotplots
    #     self.altref2read_dotplots = altref2read_dotplots
    def set_dotplots(self, ref2read_dotplots):
        self.ref2read_dotplots = ref2read_dotplots
    def set_clear_dotplots(self, clear_dotplots):
        self.clear_dotplots = clear_dotplots

    def set_id(self, id):
        self.id_num = id

    def update_id(self, id):
        self.id_num = id

    def to_string(self):
        return "{}-{}-{}-{}-{}-{}-{}".format(self.origin_chrom, self.id, self.origin_length, self.origin_start, self.origin_end, self.origin_type, self.reptype)


class BKP:

    def __init__(self, type, length, start, end):
        self.type = type
        self.length = length
        self.start = start
        self.end = end


        self.insert_pos = -1

    def to_string(self):
        return "{}_{}_{}_{}_{}".format(self.type, self.length, self.start, self.end, self.insert_pos)

def refine_detailed_bkps_for_svision(sv_type, sv_start, sv_end, detailed_bkps):
    """
    do some refinement for SVision, for better adaptation and analysis
    :param sv_type:
    :param detailed_bkps:
    :return:
    """

    if "INV" in sv_type and "DUP" in sv_type:
        inversion_start = -1
        inversion_end = -1
        inversion_bkp = -1

        duplication_start = -1
        duplication_end = -1
        duplication_bkp = -1

        for i in range(len(detailed_bkps)):
            bkp = detailed_bkps[i]
            bkp_type = bkp.type
            bkp_length = bkp.length
            bkp_start = bkp.start
            bkp_end = bkp.end

            if bkp_type == "INV":
                inversion_start = bkp_start
                inversion_end = bkp_end

                inversion_bkp = bkp

            if bkp_type == "DUP" or bkp_type == "tDUP":
                duplication_start = bkp_start
                duplication_end = bkp_end

                duplication_bkp = bkp

        if inversion_start == duplication_start and inversion_end == duplication_end:
            detailed_bkps.remove(inversion_bkp)
            detailed_bkps.remove(duplication_bkp)
            detailed_bkps.append(BKP("invDUP", inversion_end - inversion_start, inversion_start, inversion_end))


    if sv_type == "INS+INV":

        insertion_index = -1
        inversion_start = -1
        inversion_end = -1

        inversion_site = 'left'
        inversion_length = -1
        for i in range(len(detailed_bkps)):
            bkp = detailed_bkps[i]

            bkp_type = bkp.type
            bkp_length = bkp.length
            bkp_start = bkp.start
            bkp_end = bkp.end

            if bkp_type == "INS":

                insertion_index = i

            if bkp_type == "INV":
                inversion_start = bkp_start
                inversion_end = bkp_end
                inversion_length = bkp_length
                inversion_site = 'left' if abs(sv_start - inversion_start) < abs(sv_end - inversion_end) else 'right'

        if inversion_site == 'left':
            detailed_bkps[insertion_index].start = inversion_end
        else:
            detailed_bkps[insertion_index].end = inversion_start

        detailed_bkps[insertion_index].length -= inversion_length

    if sv_type == "DEL+INV":
        deletion_start = -1
        deletion_end = -1
        deletion_index = -1

        inversion_start = -1
        inversion_end = -1

        inversion_site = 'left'

        for i in range(len(detailed_bkps)):
            bkp = detailed_bkps[i]

            bkp_type = bkp.type
            bkp_length = bkp.length
            bkp_start = bkp.start
            bkp_end = bkp.end

            if bkp_type == "DEL":
                deletion_start = bkp_start
                deletion_end = bkp_end

                deletion_index = i

            if bkp_type == "INV":
                inversion_start = bkp_start
                inversion_end = bkp_end

                inversion_site = 'left' if abs(sv_start - inversion_start) < abs(sv_end - inversion_end) else 'right'

        if inversion_site == 'left':
            detailed_bkps[deletion_index].start = inversion_end
        else:
            detailed_bkps[deletion_index].end = inversion_start


    ins_index = -1
    dup_indexes = []

    for index in range(len(detailed_bkps)):
        bkp = detailed_bkps[index]

        bkp_type = bkp.type
        bkp_length = bkp.length
        bkp_start = bkp.start
        bkp_end = bkp.end

        if "DUP" in bkp_type:
            dup_indexes.append(index)
            insert_pos = sv_start if abs(sv_start - bkp_start) > abs(sv_end - bkp_end) else sv_end

            bkp.insert_pos = insert_pos

        else:
            if bkp_type == "INS":
                ins_index = index


    # # STEP: set inserted pos for DUP
    if ins_index != -1:
        insert_pos = detailed_bkps[ins_index].start

        if dup_indexes != []:
            for index in dup_indexes:
                detailed_bkps[index].insert_pos = insert_pos


def parse_vcf_record_sim(record):
    """
    parse simulated vcf
    :param record:
    :return:
    """

    # # parse basic info
    sv_chrom = record.contig
    sv_start = record.start
    sv_end = record.stop

    sv_type = record.info["SVTYPE"]

    # # parse info from record_split
    record_split = str(record).split("\t")

    id = record_split[2]

    ref_base = record_split[3]
    alt_base = record_split[4]

    if sv_type == "INS":
        sv_end = sv_start + len(alt_base)
        alt_base = record.info["BKPS"][0].replace(":", "")


    # # parse support reads
    supp_reads = list(record.info["READS"])
    # supp_reads = [i for i in range(10)]

    detailed_bkps_objs = []


    bkp_cords = record.info["BKPS"]

    for bkp_cord in bkp_cords:
        bkp_type = bkp_cord.split(":")[0]

        insert_pos = int(bkp_cord.split(":")[1].split("_")[2])

        bkp_start = int(bkp_cord.split(":")[1].split("_")[0])
        bkp_end = int(bkp_cord.split(":")[1].split("_")[1])

        bkp = BKP(bkp_type, bkp_end - bkp_start, bkp_start, bkp_end)
        bkp.insert_pos = insert_pos

        detailed_bkps_objs.append(bkp)


    return SV(id, sv_chrom, sv_start, sv_end, sv_type, ref_base, alt_base, supp_reads, detailed_bkps_objs)



def parse_vcf_record(record):
    """
    parse vcf record to generate SV object

    :return:
    """

    # # parse basic info
    sv_chrom = record.contig
    sv_start = record.start
    sv_end = record.stop

    sv_type = record.info["SVTYPE"]

    # # parse info from record_split
    record_split = str(record).split("\t")

    id = record_split[2]

    ref_base = record_split[3]
    alt_base = record_split[4]

    if sv_type == "INS":
        sv_end = sv_start + len(alt_base)

    # # parse support reads
    info_fields = record.info.keys()
    if "REPTYPE" in info_fields:
        sv_reptype = record.info["REPTYPE"]
    else:
        sv_reptype = None

    # if "RNAMES" in info_fields:
    #     supp_reads = list(record.info["RNAMES"])
    # elif "READS" in info_fields:
    #     supp_reads = list(record.info["READS"])
    # else:
    #     supp_reads = []
    if "SUPPORT" in info_fields:
        supp_reads = record.info["SUPPORT"]
    elif "SUPPORT_INLINE" in info_fields:
        supp_reads = record.info["SUPPORT_INLINE"]
    elif "SUPPORT_LONG" in info_fields:
        supp_reads = record.info["SUPPORT_LONG"]
    else:
        supp_reads = 0


    if "BKPS" in info_fields:

        detailed_bkps = record.info["BKPS"]
        #print(detailed_bkps)
        try:
            # # string, then convert to list
            detailed_bkps = detailed_bkps.split(",")

        except:
            # # tuple, then no convert
            detailed_bkps = detailed_bkps
    else:
        detailed_bkps = ["{}:{}-{}-{}".format(sv_type, sv_end - sv_start, sv_start, sv_end)]

    # # convert to object
    detailed_bkps_objs = []

    for bkp in detailed_bkps:
        try:
            bkp_type = bkp.split(":")[0]
            bkp_length = int(bkp.split(":")[1].split("-")[0])
            bkp_start = int(bkp.split(":")[1].split("-")[1])
            bkp_end = int(bkp.split(":")[1].split("-")[2])
        except:
            bkp_type = sv_type
            bkp_length = sv_end - sv_start
            bkp_start = sv_start
            bkp_end = sv_end
        finally:
            detailed_bkps_objs.append(BKP(bkp_type, bkp_length, bkp_start, bkp_end))

    refine_detailed_bkps_for_svision(sv_type, sv_start, sv_end, detailed_bkps_objs)

    return SV(id, sv_chrom, sv_start, sv_end, sv_type, sv_reptype, ref_base, alt_base, supp_reads, detailed_bkps_objs)


def parse_vcf_record_for_train(record):
    """
    parse vcf record to generate SV object

    :return:
    """

    # # parse basic info
    sv_chrom = record.contig
    sv_start = record.start
    sv_end = record.stop

    sv_type = record.info["SVTYPE"]
    # sv_af = float(record.info["AF"][0])
    sv_af = 1

    # # parse info from record_split
    record_split = str(record).split("\t")

    id = record_split[2]

    ref_base = record_split[3]
    alt_base = record_split[4]

    if sv_type == "INS":
        sv_end = sv_start + len(alt_base)

    # # parse support reads
    info_fields = record.info.keys()

    if "RNAMES" in info_fields:
        supp_reads = list(record.info["RNAMES"])
    elif "READS" in info_fields:
        supp_reads = list(record.info["READS"])
    else:
        supp_reads = []

    return SV(id, sv_chrom, sv_start, sv_end, sv_type, ref_base, alt_base, supp_reads, sv_af)