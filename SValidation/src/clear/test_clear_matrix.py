import os
import torch
import h5py
from Unet import UNet
import torch.nn.functional as F
from math import ceil
import numpy as np
import cv2

def pad_to_nearest_256_multiple(matrix):

    max_dim = max(matrix.shape[0], matrix.shape[1])
    square_dim = 256 * ceil(max_dim / 256)

    # 计算距离新的方阵边长需要填充的行数和列数
    height_pad = square_dim - matrix.shape[0]
    width_pad = square_dim - matrix.shape[1]

    # 对矩阵进行填充
    padded_matrix = np.pad(matrix, ((0, height_pad), (0, width_pad)), 'constant', constant_values=1)
    return padded_matrix

def downsample_to_256(matrix):
    # 分别为每个维度计算下采样步长
    rows, cols = matrix.shape
    height_step_size = rows // 256
    width_step_size = cols // 256

    downsampled_matrix = np.ones((256, 256))
    for i in range(256):
        for j in range(256):
            cols_start = j * width_step_size
            cols_end = (j + 1) * width_step_size
            if np.any(matrix[i * height_step_size, cols_start:cols_end] == 0):
                downsampled_matrix[i, j] = 0
            else:
                downsampled_matrix[i, j] = matrix[i * height_step_size, cols_start]

    return downsampled_matrix
def test_model(original_path, clear_path, original_dotplots):
    device = torch.device("cpu")
    print("using {} device".format(device))

    # 加载模型
    model = UNet().to(device)
    model.load_state_dict(torch.load("../../src/clear/trained_model_87_trf_Unet_pd_1_Conbined_L1pd0_l2edge_epo3_stop_early.pth", map_location="cpu"))

    original_files = os.listdir(original_path)
    for original_file in original_files:

        if original_file.endswith(".h5"):
            print(original_file)
            file_path = os.path.join(original_path, original_file)
            with h5py.File(file_path, 'r') as f:
                data = f['matrix'][:]
            data2save = (data * 255).astype(np.uint8)
            data_pad = pad_to_nearest_256_multiple(data)
            down_sample = downsample_to_256(data_pad)
            data_down2save = (down_sample * 255).astype(np.uint8)
            data_tensor = torch.from_numpy(down_sample).float().to(device)
            output = model(data_tensor)

            output_np = output.cpu().detach().squeeze().numpy()
            output_np_trans = (output_np * 255).astype(np.uint8)


            # output_np_filtered = (output_np * 255).astype(np.uint8)
            # output_np_filtered = cv2.adaptiveThreshold(output_np_filtered, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
            #                                            cv2.THRESH_BINARY, 11, 2)
            #
            # save_original_name = original_file.replace("_matrix.h5", "_filter.png")
            # save_original_path = os.path.join(original_dotplots, save_original_name)
            # cv2.imwrite(save_original_path, output_np_filtered)

            save_original_name = original_file.replace("_matrix.h5", "_original.png")
            save_original_path = os.path.join(original_dotplots, save_original_name)
            cv2.imwrite(save_original_path, data2save)

            save_original_name = original_file.replace("_matrix.h5", "_downsample.png")
            save_original_path = os.path.join(original_dotplots, save_original_name)
            cv2.imwrite(save_original_path, data_down2save)

            save_clear_name = original_file.replace("_matrix.h5", "_clear.png")
            save_clear_path = os.path.join(clear_path, save_clear_name)
            cv2.imwrite(save_clear_path, output_np_trans)
            break


if __name__ == "__main__":
    original_path = '../../data/matrix_origin/'
    clear_path = '../../data/test_dotplots_down_filter/'
    original_dotplots = '../../data/test_dotplots_down_filter/'

    os.makedirs(clear_path,exist_ok=True)
    os.makedirs(original_dotplots, exist_ok=True)

    test_model(original_path, clear_path,original_dotplots)

