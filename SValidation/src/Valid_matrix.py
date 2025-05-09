import pysam
from src.inputs.sv_op import parse_vcf_record

import datetime
import traceback
import os

import time
from multiprocessing import Process

from src.plots.dotplot_create import generate_dotplots

from src.clear.clearance import clear_dotplots
from src.valid.line_segment import line2Segments


def get_analysed_sv_id(analysed_sv_result):
    analysed_sv_id = []
    with open(analysed_sv_result, 'r') as file:
        lines = file.readlines()
    for line in lines:
        if line is not None:
            analysed_sv_id.append(line.strip())
    return  analysed_sv_id

class Options:
    def __init__(self):

        self.out_path = "/mnt/d/SV_Validation/Valid_Simu/data/hg002_cutesv_429/"
        os.makedirs(self.out_path, exist_ok=True)
        self.bam_path = "/mnt/e/Data/GIAB.HG002.GRCh37.HiFi.minimap2.bam"
        # self.vcf_path = "/mnt/d/SV_Validation/Valid_Simu/data/HG002_SVs_Tier1_v0.6.vcf"
        self.vcf_path = "/mnt/d/SV_Validation/Valid_Simu/data/hg002_sample.svision_pro_v2.3.s5.vcf"
        self.ref_path = "/mnt/e/Data/GRCh37.fasta"
        self.visual = True
        self.word_size = 10
        self.min_mapq = 20
        self.pix2pix_img_size = 256
        self.threads = 1


def process_single_sv(sv, origin_matrix_out_path, clear_matrix_out_path, result_path, options, sv_id, exten=1, retry_cont=0, max_retries=4):
    """
    处理单个 SV 的逻辑，包括矩阵生成、清除重复、生成线段等操作。
    """
    try:
        # 生成矩阵
        print()
        print(f"[Generating matrix]: {sv.to_string()}")
        sv_ref2read_matrix = generate_dotplots(sv, origin_matrix_out_path, options, mode="test", visual=True, exten=exten)

        # 更新 SV 信息
        sv.set_dotplots(sv_ref2read_matrix)
        sv.update_id(sv_id)

        # 清除重复
        clear_dotplots(sv, origin_matrix_out_path, clear_matrix_out_path, visual=True, options=options)

        # 生成线段, 如果识别成功，则返回false， 不在继续试探，如果识别失败，则返回true，然后以新的参数重新执行
        flag = line2Segments(sv, result_path, options)
        print(flag)
        print(retry_cont)
        print(max_retries)
        if flag and retry_cont < max_retries:
            print(f"Error occurred in line2Segments for SV {sv.to_string()}, retrying with different parameters (retry count: {retry_cont + 1})")
            retry_cont += 1
            sv.set_dotplots(None)
            process_single_sv(sv, origin_matrix_out_path, clear_matrix_out_path, result_path, options, sv_id, exten=retry_cont + 1, retry_cont=retry_cont)


    except Exception as e:

        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        error_info = f"Time: {current_time}, Error: {e}"
        error_message = traceback.format_exc()

        error_log_path = './error_log/' + 'error_log_mian_' + result_path.split('/')[-1]

        # 使用'a'模式打开文件以追加内容
        with open(error_log_path, 'a') as file:
            file.write(error_info + "\n")
            file.write("Traceback:\n")
            sv_information = sv.to_string()
            file.write(sv_information + "\n")
            file.write(error_message + "\n")

        print(f"Error processing SV {sv.id}: {e}")
        print(sv.to_string())



# 记录超时信息
def log_timeout(sv_id, result_path, timeout):
    """
    记录超时信息到日志文件
    """
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    timeout_log_path = './error_log/' + 'timeout_log_' + result_path.split('/')[-1]

    with open(timeout_log_path, 'a') as file:
        file.write(f"Time: {current_time}, SV {sv_id} exceeded timeout of {timeout} seconds.\n")
        file.write(f"  \n")


# 主函数
def main():
    result_path = './results/sv_result_hg002_cuteSV_429.txt'

    with open(result_path, 'w') as f:
        pass

    options = Options()
    origin_matrix_out_path = os.path.join(options.out_path, "matrix_origin")
    os.makedirs(origin_matrix_out_path, exist_ok=True)

    clear_matrix_out_path = os.path.join(options.out_path, "matrix_clear")
    os.makedirs(clear_matrix_out_path, exist_ok=True)

    # 读取 VCF 文件
    vcf_file = pysam.VariantFile(options.vcf_path)
    sv_list = []
    sv_id = 0
    inv_cnt = 0

    for record in vcf_file:
        sv = parse_vcf_record(record)

        # if len(sv.origin_read_bases) <= 50 and len(sv.origin_ref_bases) <= 50:
        #     continue

        # if sv.id != "6":
        #     continue
        if sv.origin_type == "INV":
            inv_cnt += 1
            continue
        if sv.origin_type not in ["INS", "DEL", "DUP", "INV"]:
            continue
        if sv.origin_length > 100000 or sv.origin_length <= 50:
            continue
        sv_list.append((sv, origin_matrix_out_path, clear_matrix_out_path, result_path, options, sv_id))
        sv_id += 1
        if sv_id >= 1:
            break

    timeout = 3600  # 设置超时时间 1200s  20min
    max_processes = options.threads  # 设置最大并行进程数

    # 用来保存当前正在运行的进程
    active_processes = []

    def terminate_process(proc, sv_id, result_path, timeout):
        print(f"Process {sv_id} exceeded timeout and will be terminated.")
        proc.terminate()
        proc.join()
        log_timeout(sv_id, result_path, timeout)

    # 按照最大进程数进行调度
    while sv_list or active_processes:
        # 启动新的进程（如果未达到最大进程数限制）
        while len(active_processes) < max_processes and sv_list:
            sv_data = sv_list.pop(0)
            sv, origin_matrix_out_path, clear_matrix_out_path, result_path, options, sv_id = sv_data
            proc = Process(target=process_single_sv,
                           args=(sv, origin_matrix_out_path, clear_matrix_out_path, result_path, options, sv_id))
            proc.start()
            active_processes.append((proc, sv.id, time.time()))  # 记录进程和开始时间

        # 检查是否有超时进程
        for proc, sv_id, start_time in active_processes:
            if not proc.is_alive():
                active_processes.remove((proc, sv_id, start_time))  # 移除已完成的进程
            elif time.time() - start_time > timeout:
                terminate_process(proc, sv_id, result_path, timeout)
                active_processes.remove((proc, sv_id, start_time))  # 移除超时的进程

        time.sleep(3)  # 休眠以减少 CPU 使用率

    print("All processes completed or terminated.")


if __name__ == '__main__':
    main()


# if __name__ == '__main__':
#     """
#     run validation process
#     :return:
#     """
#     # build sv_result.txt
#     result_path = './results/sv_result_829.txt'
#     # os.makedirs(result_path, exist_ok=True)
#     with open(result_path, 'w') as f:
#         pass
#
#     options = Options()
#
#     # # # STEP: check output path
#     origin_matrix_out_path = os.path.join(options.out_path, "matrix_origin")
#     os.makedirs(origin_matrix_out_path, exist_ok=True)
#
#     clear_matrix_out_path = os.path.join(options.out_path, "matrix_clear")
#     os.makedirs(clear_matrix_out_path, exist_ok=True)
#
#     vcf_file = pysam.VariantFile(options.vcf_path)
#     sv_id = 0
#     cnt = 0
#     for record in vcf_file:
#         sv = parse_vcf_record(record)
#         # print(sv.origin_chrom)
#         # if sv.id != "HG2_PB_SVrefine2PBcRplusDovetail_7787":
#             # print(sv.id)
#             # continue
#
#         if sv.origin_chrom != "3":
#             # print(sv.origin_chrom)
#             # print(sv.id)
#             continue
#
#         if sv.origin_type == "INV":
#             cnt += 1
#             print(cnt)
#             continue
#
#         if sv.origin_type not in ["INS", "DEL", "DUP", "INV"]:
#
#             continue
#
#         # if sv.origin_supp_reads < 10:
#         #     print("{} support reads less than 10".format(sv.to_string()))
#         #     continue
#         # if sv.origin_supp_reads > 500:
#         #     print("{} support reads are more than 500".format(sv.to_string()))
#         #     continue
#
#         if sv.origin_length > 100000 or sv.origin_length < 50:
#             print("{} sv length more than 100000bp".format(sv.to_string()))
#             continue
#
#         # # STEP: generate matrix for each SV
#         print("[Generating matrix]: {}".format((sv.to_string())))
#         sv_ref2read_matrix = generate_dotplots(sv, origin_matrix_out_path, options, mode="test", visual=True)
#
#         # # update sv info
#         sv.set_dotplots(sv_ref2read_matrix)
#         sv.update_id(sv_id)
#         sv_id += 1
#
#         clear_dotplots(sv, origin_matrix_out_path, clear_matrix_out_path, options)
#
#         line2Segments(sv, origin_matrix_out_path, clear_matrix_out_path, result_path)
#
#         break
#     print("results analysis")
#     # result_analysis(sv_result)
#     print("finished results analysis")

    # main()