import gc
import subprocess
import pandas as pd
import numpy as np
import os
import cv2 as cv
import cv2
import h5py
import time
class Dotplot_Gepard:
    def __init__(self, name, seq_x, seq_y, zoom, word_size, origin_out_prefix_img, ref_start, ref_end, visual=False):
        self.name = name
        self.seq_x = seq_x
        self.seq_y = seq_y
        self.shrink_ratio = None
        self.word_size = word_size
        self.zoom = zoom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.origin_out_prefix = origin_out_prefix_img
        self.origin_matrix_file = self.origin_out_prefix + '/' + self.name + "_matrix.h5"
        self.origin_dotplot_file = None
        self.visual = visual
        self.matrix_norm = None
        self.img = None
        self.clear_matrix_file = None
        self.clear_dotplot_file = None
        self.create_matrix_gepard()

    def set_clear_matrix_file(self, file_path):
        self.clear_matrix_file = file_path

    def set_original_dotplot_file(self, original_path):
        self.origin_dotplot_file = original_path

    def set_clear_dotplot_file(self, clear_path):
        self.clear_dotplot_file = clear_path

    def set_shrink_ratio(self, shrink_ratio):
        self.shrink_ratio = shrink_ratio

    def create_matrix_gepard(self):
        """
        create dot matrix by gepard
        :param work_dir:
        :return:
        """
        # # STEP: output seq to file
        seq_x_out_file =  f"{self.origin_out_prefix}_{self.name}.seq_x.fa"
        seq_y_out_file = f"{self.origin_out_prefix}_{self.name}.seq_y.fa"

        with open(seq_x_out_file, 'w') as fout:
            fout.write(">seq_x\n")
            fout.write("{}\n".format(self.seq_x))

        with open(seq_y_out_file, 'w') as fout:
            fout.write(">seq_y\n")
            fout.write("{}\n".format(self.seq_y))

        # # STEP: use gepard to create dot matrix
        matrix_out_file = f"{self.origin_out_prefix}_{self.name}.matrix"

        # # STEP: run gepard
        # cmd = "java -cp ./src/pack_dotplot/Gepard_altered.jar org.gepard.client.cmdline.CommandLine -seq {} {} -matrix ./src/pack_dotplot/edna.mat -outfile {} -silent -word {} -zoom {}".format(seq_x_out_file, seq_y_out_file, matrix_out_file, self.word_size, self.zoom)
        cmd = "java -cp ./Gepard_altered.jar org.gepard.client.cmdline.CommandLine -seq {} {} -matrix ./edna.mat -outfile {} -silent -word {} -zoom {}".format(seq_x_out_file, seq_y_out_file, matrix_out_file, self.word_size, self.zoom)
        # cmd = "java -cp ./src/pack_dotplot/Gepard_altered.jar org.gepard.client.cmdline.CommandLine -seq {} {} -matrix ./src/pack_dotplot/edna.mat -outfile {} -silent -word {} ".format(ref_seq_out_file, alt_seq_out_file, ref2alt_out_file, self.word_size, self.zoom)
        try:
            start_time = time.time()
            cmd_run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=-1)

            try:
                cmd_run.wait(timeout=600)  # 设置超时时间为 10 min, 超出时间限制则进行结束子进程

            except subprocess.TimeoutExpired:
                print(f"Command timed out after 10 minutes: {cmd}")
                cmd_run.kill()  # 终止子进程
                raise Exception("Timeout while running Gepard")

            out, err = cmd_run.communicate()

            if err:
                raise Exception(err.decode())

            # print("Command finished in", time.time() - start_time, "seconds.")

            # # STEP: loading gepard's output as  matrix
            self.matrix = pd.read_csv(matrix_out_file, sep='\t', header=None).fillna(0)

            self.matrix_project_x = self.matrix.apply(lambda x: sum(x), axis=0)
            self.matrix_project_y = self.matrix.apply(lambda x: sum(x), axis=1)

            # print(self.ref2alt_ref_project)

            # # STEP: normalize the origin matrix
            self.matrix_norm = abs(self.matrix - self.matrix.values.max()) / (self.matrix.values.max() - self.matrix.values.min())

            matrix_value = self.matrix_norm.values

            with h5py.File(self.origin_matrix_file, 'w') as hf:
                hf.create_dataset('matrix', data=matrix_value, compression='gzip')

            if self.visual:
                visual_path = self.origin_out_prefix.replace("matrix_origin", "dotplots_visual") + '/original'
                os.makedirs(visual_path, exist_ok=True)
                visual_file = visual_path + '/' + self.name + "_dotplots.png"
                self.set_original_dotplot_file(visual_file)
                # STEP: output three channel img
                self.img_height, self.img_width = np.shape(self.matrix_norm)

                channel_values = (self.matrix_norm.values * 255).astype(np.uint8)
                self.img = np.ones((self.img_height, self.img_width))
                self.img = channel_values
                cv.imwrite(visual_file, self.img)
                self.img = None


        except Exception as e:
            print("op dotplot gepard error:", e)

        finally:
            self.matrix_project_x = None
            self.matrix_project_y = None
            self.matrix = None
            os.remove(matrix_out_file)
            os.remove(seq_x_out_file)
            os.remove(seq_y_out_file)
            gc.collect()



    def to_png(self):
        pass
        # cv.imwrite(self.origin_dotplot_file, self.img)

        # cv.imwrite(self.origin_dotplot_file, self.img_resize)

        # #############
        # print(self.origin_dotplot_file)
        # image = cv.imread(self.origin_dotplot_file)
        # gray = cv.cvtColor(image, cv.COLOR_BGR2GRAY)
        # gray[gray > 0] = 1
        # gray[gray == 0] = 255
        # gray[gray == 1] = 0
        # min_line_length = max_gap_length = 50 / self.zoom
        #
        # lines = cv.HoughLinesP(gray, 1, np.pi / 180, 50, minLineLength=min_line_length, maxLineGap=max_gap_length)
        #
        # for line in lines:
        #     # print(line[0])
        #     x1, y1, x2, y2 = line[0]
        #     cv.line(image, (x1, y1), (x2, y2), (0, 255, 0), 2)
        #
        # cv.imwrite(self.origin_dotplot_file, image)

def custom_resize(src, re_size):

    target_img = np.ones((re_size, re_size, 3)) * 255

    target_img_size = np.shape(target_img)
    base_img_size = np.shape(src)

    target_img[0: min(target_img_size[0], base_img_size[0]), 0: min(target_img_size[1], base_img_size[1]), :] = src[0: min(target_img_size[0], base_img_size[0]), 0: min(target_img_size[1], base_img_size[1]), :]

    return target_img


def custom_resize_2(src, re_size):
    # 计算缩放比例。取宽度和高度缩放比例中较小的一个来确保等比例缩放
    scale = min(re_size / src.shape[0], re_size / src.shape[1])

    # 计算新尺寸
    new_size = (int(src.shape[1] * scale), int(src.shape[0] * scale))

    # 使用cv2.resize进行等比例缩放
    resized_img = cv2.resize(src, new_size, interpolation=cv2.INTER_AREA)

    # 创建一个新的目标图像矩阵，初始化为白色
    target_img = np.ones((re_size, re_size, 3), dtype=np.uint8) * 255

    # 计算画布上用于放置已缩放图像的起始点（若图像未填满新尺寸）
    x_offset = (re_size - new_size[0]) // 2
    y_offset = (re_size - new_size[1]) // 2

    # 将已缩放图像置于画布的相应位置
    target_img[y_offset:y_offset + new_size[1], x_offset:x_offset + new_size[0]] = resized_img

    return target_img