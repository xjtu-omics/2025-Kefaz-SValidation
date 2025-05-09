import gc
import os
import numpy as np
import traceback
from src.plots.line_detect import lsd_line_detection, lsd_line_detection_for_inv
# from src.valid.parse_segments import parse_segments
from sklearn.cluster import KMeans, DBSCAN
from src.valid.search_segments import search_segments
import datetime
from collections import Counter
from sklearn.metrics import silhouette_score

class Segment:
    def __init__(self, seg_id, seg_ref_start, seg_ref_end, seg_read_start, seg_read_end, seg_strand, shrink_ratio):
        self.id = seg_id
        self.ref_start = seg_ref_start
        self.ref_end = seg_ref_end
        self.read_start = seg_read_start
        self.read_end = seg_read_end
        self.strand = seg_strand
        self.shrink_ratio = shrink_ratio
        self.flag = None

    def set_seg_flag(self, flag):
        self.flag = flag

    def update_seg_id(self, new_id):
        self.id = new_id

    def to_string(self):
        return "Segment: {}, ref: {}-{}, read: {}-{}, {}, {}, {}".format(self.id, self.ref_start, self.ref_end,
                                                                     self.read_start, self.read_end, self.strand,
                                                                     self.flag, self.shrink_ratio)


def convert_lines_to_segments(lines, shrink_ratio, pix2pix_img_size=256, segment_seq_out_prefix=None,
                              options=None):
    """
    convert lines to segments
    """

    # shrink_ratio = max(len(ref_seq), len(read_seq)) / pix2pix_img_size

    segments = []

    for line in lines:
        # seg_ref_start = int(line.ref_start * shrink_ratio + ref_start)
        # seg_ref_end = int(line.ref_end * shrink_ratio + ref_start)
        seg_ref_start = int(line.ref_start * shrink_ratio[1])
        seg_ref_end = int(line.ref_end * shrink_ratio[1])

        seg_read_start = int(line.read_start * shrink_ratio[0])
        seg_read_end = int(line.read_end * shrink_ratio[0])

        seg_strand = "+"

        segments.append(Segment(line.id, seg_ref_start, seg_ref_end, seg_read_start, seg_read_end, seg_strand, shrink_ratio))

    # for seg in segments:
    #     print(seg.to_string())
    #     print("--------------------------------------------")

    return segments

def filter_points(points, dist_threshold = 1.0, pos_threshold = 1.0):
    '位置和距离双重过滤，去除异常点对于聚类结果的影响'
    distances = np.sqrt(np.sum((points[0::2] - points[1::2]) ** 2, axis=1))
    centers = (points[0::2] + points[1::2]) / 2
    # 计算所有距离的平均值和标准差
    mean_dist = np.mean(distances)
    std_dist = np.std(distances)
    # 计算中心点坐标的平均值和标准差
    mean_center = np.mean(centers, axis=0)
    std_center = np.std(centers, axis=0)
    # 找出异常的点对
    # 过滤异常点对
    outlier_indices = np.where((np.abs(distances - mean_dist) > dist_threshold * std_dist) |
                               (np.sqrt(np.sum((centers - mean_center) ** 2, axis=1)) > pos_threshold * np.max(
                                   std_center)))[0]

    filtered_points = np.delete(points, np.concatenate([2 * outlier_indices, 2 * outlier_indices + 1]), axis=0)
    #如果过滤的比例超过一半，则放弃过去，认为此时的过滤方法不合适
    if len(filtered_points) < 0.5 * len(points):
        return points
    else:
        return  filtered_points


def filter_points_with_clustering(points, eps=10.0, min_samples=5):
    """
    使用DBSCAN聚类算法过滤异常点，DBSCAN能够自动识别离群点。

    Parameters:
    - points: np.array, 点的坐标集
    - eps: float, 两个点之间的最大距离，DBSCAN参数, 小于该距离被视为邻居
    - min_samples: int, 定义一个点至少需要多少邻居才能被视为核心点

    Returns:
    - filtered_points: 过滤后的点
    """
    # 将点对中的每个点组合为一个2D坐标
    centers = (points[0::2] + points[1::2]) / 2

    # 使用DBSCAN聚类
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(centers)

    # 过滤噪声点（label=-1表示噪声点）
    labels = clustering.labels_
    non_noise_indices = np.where(labels != -1)[0]

    # 获取非噪声点的坐标
    filtered_points = np.concatenate([points[2 * i: 2 * i + 2] for i in non_noise_indices])

    # 如果过滤的点太多，返回原始点集
    if len(filtered_points) < 0.4 * len(points):
        return points
    else:
        return filtered_points

def segments_cluster(Segments, clusters, sv):
    # 提取到的所有的坐标值,对于每一个segment的两个坐标分别聚类，不能直接对多个坐标聚类，这样会改变sv的对应信息
    # clusters记录了sv对应直线的平均数

    coordinates = [[] for _ in range(clusters)]
    # print(len(coordinates))
    # 外循环：Segment保存了每一张图像中的所有segment
    # 内循环：segment保存了一张图像里面每一个segment
    for Segment in Segments:
        # 暂时对去重效果不好的数据不予处理
        if (len(Segment) != clusters):
            continue
        for segment_idx in range(len(Segment)):
            segment = Segment[segment_idx]
            start = (segment.ref_start, segment.read_start)
            end = (segment.ref_end, segment.read_end)
            coordinates[segment_idx].append(start)
            coordinates[segment_idx].append(end)
    centers = [[] for _ in range(clusters)]
    for i in range(clusters):
        points = np.array(coordinates[i])
        flitered_points = filter_points(points)
        if len(flitered_points) == 0:
            print(sv.id)
            raise ValueError(f"No points to cluster for segment with ID {sv.id}")

        kmeans = KMeans(2, n_init= 10).fit(flitered_points)

        # 绘制聚类中心,并转换为整数,按照起点坐标--终点坐标的方式进行排序
        centers_int = kmeans.cluster_centers_.astype(int)
        # 使用argsort根据第一个维度进行排序的索引
        sorted_indices = np.argsort(centers_int[:, 0])
        # 使用排序后的索引来排序聚类中心
        centers[i] = centers_int[sorted_indices]
    return centers

def filter_lines_by_length(lines, shrink_ratio=None):
    #高度方向，即行， 对应read；宽度方向，即列，对应ref
    read_shrink_ratio, ref_shrink_ratio = shrink_ratio
    lines_over_50 = []
    ##获取检测到的line的实际长度，并对长度在50以下的进行过滤
    for line in lines:
        # 根据高度和宽度的压缩比计算实际的线段长度
        actual_length = np.sqrt((line.length_read * read_shrink_ratio) ** 2 + (line.length_ref * ref_shrink_ratio) ** 2)
        line.set_actual_length(actual_length)

        if actual_length > 50:
            lines_over_50.append(line)

    # 使用列表推导式从每个对象中提取length属性
    lengths = [line.actual_length for line in lines_over_50]
    ## 如果 lengths的长度为2， 则说明该图像内仅仅有1条或者2条线，认为这种情况是去噪成功的
    if len(lengths) <= 2:
        return lines_over_50, shrink_ratio
    else:
        ##以下内容期待对lengths不为2的情况进行，即模型去噪失败的单张图像中的线段进行聚类, 通过聚类，找到最长的一组
        ## 以下代码是为了找到最佳的聚类数
        best_num_clusters = 2
        max_clusters = min(4, len(lengths) - 1)  # 确保最大聚类数量不超过样本数减1
        best_silhouette = float('-inf')
        for n_clusters in range(2, max_clusters + 1):
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            labels = kmeans.fit_predict(np.array(lengths).reshape(-1, 1))
            silhouette_avg = silhouette_score(np.array(lengths).reshape(-1, 1), labels)
            if silhouette_avg > best_silhouette:
                best_silhouette = silhouette_avg
                best_num_clusters = n_clusters

        ##以下内容根据最佳聚类数量，进行聚类，然后取长度最长的
        # K-means聚类
        kmeans = KMeans(n_clusters=best_num_clusters, random_state=42, n_init=10)
        labels = kmeans.fit_predict(np.array(lengths).reshape(-1, 1))

        # 选择要保留的聚类类别
        avg_lengths = []
        for i in range(best_num_clusters):
            avg_length = np.mean([lengths[j] for j in range(len(lengths)) if labels[j] == i])
            avg_lengths.append(avg_length)

        # 假设我们选择平均长度最大的那个类别
        target_cluster = np.argmax(avg_lengths)

        # 保留目标类别的线段
        filtered_lines = [lines_over_50[i] for i in range(len(lines_over_50)) if labels[i] == target_cluster]

        return filtered_lines, shrink_ratio

    # if len(lengths) == 2:
    #     return lines_over_50, shrink_ratio
    #
    # else:
    #     median_length = np.median(lengths)
    #     for line_over_50 in lines_over_50:
    #
    #         if line_over_50.actual_length >= 0.5 * median_length:
    #
    #             filtered_Line.append(line_over_50)
    #
    #     return filtered_Line, shrink_ratio

def filter_lines_by_segments(Lines, shrink_ratios):
    slope_threshold = 0.1
    intercept_threshold = 10.0
    filtered_lines = []
    filtered_shrink_ratios = []
    # 统计各个长度的出现次数
    counts = Counter(len(lines) for lines in Lines)

    # 如果只有长度为0或1的线段，将其从统计中移除
    if 0 in counts and counts[0] > 0:
        del counts[0]
    if 1 in counts and counts[1] > 0:
        del counts[1]

    # 查找出现最频繁的线段数量（排除了0和1后的最大值）
    most_frequent = max(counts, key=counts.get)
    if most_frequent != 2:
        freq_of_two = counts.get(2, 0)
        if abs(counts.get(most_frequent, 0) - freq_of_two) <= 5 and freq_of_two > 0:
            most_frequent = 2

    # 筛选出线段数量等于most_frequent的所有Lines
    for lines_idx in range(len(Lines)):
        if len(Lines[lines_idx]) == most_frequent:
            # 对于经过过滤后的lines进行再次过滤，如果所有直线之间没有明显的差异，则对该样本进行抛弃，避免正常图像内容影响 SV图像的识别
            processed = set()
            all_similar = True  # 初始化为True，假设所有线段都相似
            for line in Lines[lines_idx]:
                if len(processed) == 0:
                    processed.add(line)
                    continue

                slope1, intercept1 = line.slope, line.intercept  # 获取当前线段的斜率和截距
                # 遍历所有已处理的线段，比较斜率和截距
                is_line_similar = False  # 判断当前线段是否和所有已处理的线段相似
                for processed_line in processed:
                    slope2, intercept2 = processed_line.slope, processed_line.intercept
                    # 判断斜率和截距是否相似
                    if abs(slope1 - slope2) < slope_threshold and abs(intercept1 - intercept2) < intercept_threshold:
                        is_line_similar = True
                        break  # 如果当前线段相似，则跳出循环

                # 如果找到相似的线段，继续下一个循环；否则，将该线段加入processed集合
                if not is_line_similar:
                    processed.add(line)
                    all_similar = False  # 如果有不相似的线段，标记为False，表示该样本不可抛弃
                    break #有线段不相似，则结束整体的循环

            # 如果所有线段都很相似，则抛弃该样本
            if not all_similar:
                filtered_lines.append(Lines[lines_idx])
                filtered_shrink_ratios.append(shrink_ratios[lines_idx])

    if len(filtered_lines) == 0:
        print("all lines are similar, without lines are retained")
        return most_frequent, [Lines[i] for i in range(len(Lines)) if len(Lines[i]) == most_frequent], [shrink_ratios[i] for i in range(len(shrink_ratios)) if len(Lines[i]) == most_frequent]
    else:
        return most_frequent, filtered_lines, filtered_shrink_ratios

def delete_matrix_dotplots(files_lists):
    for file_list in files_lists:
        for file_path in file_list:
            try:
                file_path = file_path.strip()  # 清理路径上的多余空格或换行符
                if os.path.exists(file_path):  # 检查文件是否存在
                    os.unlink(file_path)  # 删除文件
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")


def kmeans_select_majority(sv_positions, n_clusters=2):
    """
    使用 KMeans 将 SV 的起点或终点聚为两类，筛选出数量较多的一类并计算均值。

    参数:
        sv_positions (list): SV 起点或终点的列表（一维数据）。
        n_clusters (int): KMeans 的聚类数，默认 2。

    返回:
        float: 数量较多一类的均值，作为最终的起点或终点。
        dict: 聚类结果，包含每一类的点和对应的均值。
    """
    sv_positions = np.array(sv_positions).reshape(-1, 1)  # 转换为二维数组

    # 使用 KMeans 进行聚类
    kmeans = KMeans(n_clusters=n_clusters, random_state=42).fit(sv_positions)
    labels = kmeans.labels_  # 获取聚类标签
    cluster_centers = kmeans.cluster_centers_  # 获取聚类中心

    # 分类点到两个簇中
    clusters = {}
    for cluster_id in range(n_clusters):
        clusters[cluster_id] = sv_positions[labels == cluster_id].flatten().tolist()

    # 找到数量较多的簇
    majority_cluster_id = max(clusters, key=lambda x: len(clusters[x]))
    majority_cluster_points = clusters[majority_cluster_id]

    # 计算数量较多一类的均值,并取整
    majority_mean = round(np.mean(majority_cluster_points))

    # 返回结果
    return majority_mean



def line2Segments(sv, result_path, options=None):

    ref2read_dotplots = sv.ref2read_dotplots
    Segments = []
    shrink_ratios = []
    segments = None
    Lines = []
    origin_matrix_path = []
    origin_dotplot_path = []
    clear_matrix_path = []
    clear_dotplot_path = []
    ref_start = []
    ref_end = []
    try:
        # 每一次循环处理一张图像，每一个segments里面存放的是单张图像的待处理segments
        for dotplot in ref2read_dotplots:
            ref_start.append(dotplot.ref_start)
            ref_end.append(dotplot.ref_end)
            origin_matrix_path.append(dotplot.origin_dotplot_file)

            # get clear dotplots path
            clear_path = dotplot.clear_matrix_file
            clear_matrix_path.append(clear_path)

            if options.visual:
                origin_dotplot_path.append(dotplot.origin_dotplot_file)
                clear_dotplot_path.append(dotplot.clear_dotplot_file)

            # line detection
            if sv.origin_type == "INV":
                origin_path = dotplot.origin_dotplot_file
                lines = sorted(lsd_line_detection_for_inv(origin_path), key=lambda x: x.read_start)
                print("This SV Type is INV:")
                print(sv.to_string())
            else:
                lines = lsd_line_detection(clear_path)

                if len(lines) == 0:
                    continue

                lines = sorted(lines, key=lambda x: x.read_start)

            # 根据每一张dotplot中线段的长度进行过滤

            filtered_lines, shrink_ratio = filter_lines_by_length(lines, dotplot.shrink_ratio)

            if len(filtered_lines) > 0:
                Lines.append(filtered_lines)
                shrink_ratios.append(shrink_ratio)


        #判断检测到线段最多的数目，并将该数目对应的dotplots进行聚类
        clusters, filtered_lines, filtered_shrink_ratios = filter_lines_by_segments(Lines, shrink_ratios)

        # convert line to segments
        for lines_idx in range(len(filtered_lines)):
            segments = sorted(convert_lines_to_segments(filtered_lines[lines_idx], filtered_shrink_ratios[lines_idx]),
                                      key=lambda x: x.read_start)

            Segments.append(segments)

        # 对所有的segments进行聚类，确定当前sv的字段特征，可以在提高处理的效率
        centers = segments_cluster(Segments, clusters, sv)

        ref_start_clustered = kmeans_select_majority(ref_start)
        ref_end_clustered = kmeans_select_majority(ref_end)
        # sv_type = parse_segments(centers)
        en_ref, en_read, sv_type, sv_region = search_segments(centers, ref_start_clustered, ref_end_clustered)

        if len(sv_type) == 0:
            raise ValueError("No sv operation is valid")

        print(sv.to_string() + "----->" + ",".join(sv_type) + ";" + f"region:{sv_region}")

        if any(sv.origin_type == sv_tp[0:3].upper() or (sv.origin_type == "INS" and sv_tp[0:3].upper() == "DUP") for sv_tp in sv_type):
            with open(result_path, 'a') as f:
                f.write(f'chorm: {sv.origin_chrom}, sv_id: {sv.id}, origin_sv_type: {sv.origin_type}, reptype: {sv.reptype}, origin_dotplot_num: {len(ref2read_dotplots)}, '
                        f'final_sv_type: {sv_type}, sv_region: {sv_region}\n'
                        f'encoded ref: {en_ref} \n'
                        f'encoded read: {en_read} \n')

            return False
        else:
            if options.visual:
                file_lists = [origin_matrix_path, origin_dotplot_path, clear_matrix_path, clear_dotplot_path]
            else:
                file_lists = [origin_matrix_path, clear_matrix_path]
            delete_matrix_dotplots(file_lists)

            return True

    except Exception as e:
        print(f"An error occurred: {e}")
        # 打印完整的异常堆栈信息
        traceback.print_exc()

        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # 准备要记录的信息
        error_info = f"Time: {current_time}, Error: {e}"
        error_traceback = traceback.format_exc()
        error_log_path = './error_log/' + 'error_log_' + result_path.split('/')[-1]
        # 使用'a'模式打开文件以追加内容
        with open(error_log_path, 'a') as file:
            file.write(error_info + "\n")
            file.write("Traceback:\n")
            sv_information = sv.to_string()
            file.write(sv_information + "\n")
            file.write(error_traceback + "\n")

        print(sv.to_string)

        if options.visual:
            file_lists = [origin_matrix_path, origin_dotplot_path, clear_matrix_path, clear_dotplot_path]
        else:
            file_lists = [origin_matrix_path, clear_matrix_path]
        delete_matrix_dotplots(file_lists)

        return True

    finally:

        ##运行完成后手动回收变量
        del ref2read_dotplots, Lines, segments, shrink_ratios
        gc.collect()
        # return False
