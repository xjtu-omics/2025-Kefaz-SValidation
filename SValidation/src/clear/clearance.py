import gc
import os
import torch
import h5py
import numpy as np
import cv2
# from UNet_PLUS_gpu1 import UNet_PP
# from src.clear.Unet import UNet
from src.clear.generate import GeneratorUNet
import torch.nn.functional as F
from math import ceil
import time

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
        rows_start = i * height_step_size
        # rows_end = (i + 1) * height_step_size
        for j in range(256):
            cols_start = j * width_step_size
            cols_end = (j + 1) * width_step_size

            if np.any(matrix[rows_start, cols_start:cols_end] == 0):
                downsampled_matrix[i, j] = 0
            else:
                downsampled_matrix[i, j] = matrix[rows_start, cols_start]

    return downsampled_matrix, (height_step_size, width_step_size)



def safe_h5py_write_with_retries(clear_path, output_np, max_retries=5, delay=2):
    """
    带有重试机制的h5py写入函数。
    :param dotplot: 当前处理的dotplot对象
    :param output_np: 需要写入的数据
    :param max_retries: 最大重试次数
    :param delay: 每次重试之间的延迟（秒）
    """
    for attempt in range(max_retries):
        try:
            with h5py.File(clear_path, 'w') as f:
                f.create_dataset('matrix', data=output_np, compression='gzip')
            return  # 成功写入后直接返回
        except OSError as e:
            print(f"写入文件失败: {e}, 尝试第 {attempt+1} 次重试...")
            time.sleep(delay)  # 等待一段时间后重试

    # 如果超过重试次数仍然失败，则抛出异常
    raise RuntimeError(f"多次重试失败，放弃写入文件: {clear_path}")




def clear_dotplots(sv, origin_dotpot_out_path, clear_dotpot_out_path, visual=False, options=None):
    """
    run denosing model to clear repeats
    Parameters
    ----------
    img_path

    Returns
    -------
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # device = torch.device("cpu")
    print("using {} device clearing {}".format(device, sv.to_string()))

    # 加载模型
    model = GeneratorUNet().to(device)
    # 加载模型检查点
    checkpoint = torch.load("./src/clear/best_model.pth", map_location=device)
    model.load_state_dict(checkpoint['generator_state_dict'])
    model.eval()
    imgs = sv.ref2read_dotplots

    # batch_tensors = []
    # dotplots = []
    # shrink_ratios = []

    # 遍历 SV 中的所有 dotplot，准备批量处理
    for img_index, dotplot in enumerate(imgs):
        if dotplot.matrix_norm is not None:

            dotplot_origin = dotplot.matrix_norm
            dotplot_origin_pad = pad_to_nearest_256_multiple(dotplot_origin)
            dotplot_origin_down, shrink_ratio = downsample_to_256(dotplot_origin_pad)

            # 转换为 Tensor 并添加到批量处理列表中
            dotplot_origin_down_tensor = torch.from_numpy(dotplot_origin_down).unsqueeze(0).unsqueeze(0).float().to(device)
            # batch_tensors.append(dotplot_origin_down_tensor)
            # dotplots.append(dotplot)
            # shrink_ratios.append(shrink_ratio)

            # 如果有图像需要处理，进行批量推理
            # if len(batch_tensors) > 0:
            #     batch_tensor = torch.cat(batch_tensors, dim=0)  # 堆叠成一个批次

            with torch.no_grad():
                output = model(dotplot_origin_down_tensor)

            dotplot.set_clear_matrix_file(dotplot.origin_matrix_file.replace(origin_dotpot_out_path, clear_dotpot_out_path).replace("_matrix.h5", "_matrix_clear.h5"))
            dotplot.set_shrink_ratio(shrink_ratio)

            # 将处理后的结果转换回numpy数组
            output_np = output.cpu().detach().squeeze().numpy()

            # 保存处理后的结果为新的 h5 文件
            # print(dotplot.clear_dotplot_file)
            safe_h5py_write_with_retries(dotplot.clear_matrix_file, output_np)
            # with h5py.File(dotplot.clear_dotplot_file, 'w') as f:
            #     f.create_dataset('matrix', data=output_np, compression='gzip')
            dotplot.matrix_norm = None
            dotplot.matrix = None

            if visual:
                output_np_trans = (output_np * 255).astype(np.uint8)
                visual_path = dotplot.origin_out_prefix.replace("matrix_origin", "dotplots_visual") + "/clear"
                os.makedirs(visual_path, exist_ok=True)
                visual_file = visual_path + '/' + dotplot.name + "_dotplots_clear.png"
                dotplot.set_clear_dotplot_file(visual_file)
                cv2.imwrite(visual_file, output_np_trans)

    ##函数执行完毕后，释放空间
    gc.collect()