
import os
import cv2 as cv
import numpy as np


def calculate_abundance_of_A_against_B(a_img_path, b_img_path, search_around=None):

    a_img = cv.imread(a_img_path, cv.IMREAD_GRAYSCALE)
    b_img = cv.imread(b_img_path, cv.IMREAD_GRAYSCALE)

    # # STEP: collect pixel in target, if the pix val is less than 255, means there is a match point (black)
    b_pixels = np.where(b_img < 255)

    if search_around is None:
        # # STEP: set match points from target img to 255 (white)
        a_img[b_pixels] = 255

    else:
        for i in range(np.shape(b_pixels)[1]):
            pos_x, pos_y = b_pixels[0][i], b_pixels[1][i]

            # for tmp_pox_x in range(pos_x - search_around, pos_x + search_around):
            for tmp_pox_x in range(pos_x, pos_x + 1):

                for tmp_pox_y in range(pos_y - search_around, pos_y + search_around):

                    if tmp_pox_x < 0 or tmp_pox_y < 0:
                        continue

                    a_img_size = np.shape(a_img)[1]
                    if tmp_pox_x >= a_img_size or tmp_pox_y >= a_img_size:
                        continue

                    a_img[tmp_pox_x, tmp_pox_y] = 255

    # # STEP: collect left pixs in output, if the pix val is less than 127, means there is still a match point
    a_abundant_pixels = np.where(a_img < 127)
    a_abundant_pixels_num = np.shape(a_abundant_pixels)[1]

    return a_abundant_pixels_num

if __name__ == '__main__':

    out_path = "/mnt/c/workspace/test/hifi/val_res_batch8_epoch4/images"
    out_file = "/mnt/c/workspace/test/hifi/val_res_batch8_epoch4.txt"

    out_file = open(out_file, "w")

    cnt = 0
    for input_img in os.listdir(out_path):


        if "inputs" not in input_img:
            continue

        cnt += 1

        if cnt % 1000 == 0:
            print(cnt)

        # if "821959-828880" not in input_img:
        #     continue

        output_img = os.path.join(out_path, input_img.replace("inputs", "outputs"))
        target_img = os.path.join(out_path, input_img.replace("inputs", "targets"))

        # # no corresponding output and target imgs
        if not os.path.exists(output_img) or not os.path.exists(target_img):
            continue


        output_againt_target_abundance = calculate_abundance_of_A_against_B(output_img, target_img, search_around=50)

        target_againt_output_abundance = calculate_abundance_of_A_against_B(target_img, output_img, search_around=50)


        target_img = cv.imread(target_img, cv.IMREAD_GRAYSCALE)
        output_img = cv.imread(output_img, cv.IMREAD_GRAYSCALE)

        target_pixels = np.where(target_img < 255)
        output_pixels = np.where(output_img < 255)


        # print("img: {}, output_againt_target: {}, target_againt_output: {}, output_total: {}, target_total: {}".format(input_img, output_againt_target_abundance, target_againt_output_abundance, np.shape(output_pixels), np.shape(target_pixels)))

        if np.shape(output_pixels)[1] != 0 and np.shape(target_pixels)[1] != 0:
            ratio_1 = round(output_againt_target_abundance / np.shape(output_pixels)[1], 5)
            ratio_2 = round(target_againt_output_abundance / np.shape(target_pixels)[1], 5)

            out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(input_img, output_againt_target_abundance, target_againt_output_abundance, np.shape(output_pixels)[1], np.shape(target_pixels)[1], ratio_1, ratio_2))

    out_file.close()


