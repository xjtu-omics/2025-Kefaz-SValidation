from collections import deque

import numpy as np


# def move2zeros(segements, horzontal, vertical):

def parse_ref_read_by_overlap_gap(ref, read):
    ref_overlaps_idx = []
    ref_overlaps_length = []
    read_overlaps_idx = []
    read_overlaps_length = []
    ref_gaps_idx = []
    ref_gaps_length = []
    read_gaps_idx = []
    read_gaps_length = []
    ref_overlap_cnt = 0
    ref_gap_cnt = 0
    read_overlap_cnt = 0
    read_gap_cnt = 0
    Counter = []
    for i in range(1, len(ref) - 1, 2):
        ref_current_end = ref[i]
        ref_next_start = ref[i + 1]
        #ref方向存在overlap
        if(ref_current_end - ref_next_start) >= 30:
            ref_overlaps_idx.append(i)
            ref_overlap_cnt += 1
            ref_overlaps_length.append(ref_current_end - ref_next_start)
        elif(ref_next_start - ref_current_end) >= 30:
            ref_gaps_idx.append(i)
            ref_gap_cnt += 1
            ref_gaps_length.append(ref_next_start - ref_current_end)

    for i in range(1, len(read) - 1, 2):
        read_current_end = read[i]
        read_next_start = read[i + 1]
        #read方向存在overlap
        if(read_current_end - read_next_start) >=30:
            read_overlaps_idx.append(i)
            read_overlap_cnt += 1
            read_overlaps_length.append(read_current_end - read_next_start)
        elif(read_next_start - read_current_end) >=30:
            read_gaps_idx.append(i)
            read_gap_cnt += 1
            read_gaps_length.append(read_next_start - read_current_end)

    condition_a = False
    condition_b = False
    condition_ab = (ref_gap_cnt == 1  and ref_overlap_cnt == 0) or (read_gap_cnt == 1 and read_overlap_cnt == 0)
    if condition_ab:
        if (ref_gap_cnt == 1  and ref_overlap_cnt == 0 and read_gap_cnt == 0 and read_overlap_cnt == 0):
            condition_a = True
        elif (ref_gap_cnt == 0  and ref_overlap_cnt == 0 and read_gap_cnt == 1 and read_overlap_cnt == 0):
            condition_b = True
        else:
            #二者均不是，则比较gap的大小
            for i in range(max(len(ref_gaps_idx), len(read_gaps_idx))):
                ref_length = ref_gaps_length[i] if i < len(ref_gaps_length) else 0
                read_length = read_gaps_length[i] if i < len(read_gaps_length) else 0

                if (ref_length > read_length) and read_length > 0:
                    condition_a = True
                elif (ref_length < read_length) and ref_length > 0:
                    condition_b = True
    # condition_c = False
    condition_de = ref_gap_cnt == 0  and ref_overlap_cnt == 1 and read_gap_cnt <= 1 and read_overlap_cnt == 0
    condition_fg = ref_gap_cnt == 0  and ref_overlap_cnt == 1 and read_gap_cnt == 0 and read_overlap_cnt == 1
    condition_h = ref_gap_cnt == 1  and ref_overlap_cnt == 0 and read_gap_cnt == 0 and read_overlap_cnt == 1
    ins_of_dup = ref_overlap_cnt >= 2 and ref_gap_cnt ==0  and read_gap_cnt ==0
    del_of_dup = ref_gap_cnt == 0  and ref_overlap_cnt == 0 and read_gap_cnt == 0 and read_overlap_cnt == 1

    # print(condition_a)
    # print(condition_b)
    # print(condition_de)
    # print(condition_fg)
    # print(condition_h)
    # print(ins_of_dup)
    # print(condition_fg2)
    if condition_fg:
        #沿着ref的方向遍历，若后面的read坐标大于前面的read坐标，则该情况为del of dup
        for i in range(1, len(read) - 1, 2):
            read_current_end = read[i]
            read_next_start = read[i + 1]
            if read_current_end - read_next_start >= 30:
                del_of_dup = True
    # print(del_of_dup)
    if condition_a or condition_h or del_of_dup:
        return "DEL"
    elif condition_b:
        return "INS"
    elif condition_de or condition_fg:
        return "DUP"
    elif ins_of_dup:
        return "INS_of_DUP"
    else:
        error_message = ("No sv is recognized, gap_on_ref: {}, overlap_on_ref: {}, gap_on_read: {}, overlap_on_read: {}"
                         .format(ref_gaps_length, ref_overlaps_length, read_gaps_length, read_overlaps_length))
        raise ValueError(error_message)

def filter_merge_segments(segs):
    filtered_segs = []
    for seg in segs:
        # 对于前面的步骤过滤失败的长度小于20的噪声序列，再次进行过滤
        distances = np.linalg.norm(seg[1:] - seg[0:])
        if distances > 40:
            filtered_segs.append(seg)

    if len(filtered_segs) > 2:
        merged = [filtered_segs[0]]
        for current in filtered_segs[1:]:
            last = merged[-1]
            # 计算起始坐标和终止坐标的横纵坐标差距
            start_diff_x = abs(last[0][0] - current[0][0])
            start_diff_y = abs(last[0][1] - current[0][1])
            end_diff_x = abs(last[1][0] - current[1][0])
            end_diff_y = abs(last[1][1] - current[1][1])
            # 检查横坐标是否重叠
            if start_diff_x <= 10 and start_diff_y <= 10 and end_diff_x <= 10 and end_diff_y <= 10:
                # 合并坐标
                new_start_x = min(last[0][0], current[0][0])
                new_start_y = min(last[0][1], current[0][1])
                new_end_x = max(last[1][0], current[1][0])
                new_end_y = max(last[1][1], current[1][1])
                merged[-1] = ([new_start_x, new_start_y], [new_end_x, new_end_y])

            else:
                merged.append(current)
        return merged
    else:
        return filtered_segs



def parse_segments(segments):
    # 对坐标进行平移，使之从0，0开始
    if not np.array_equal(segments[0][0, :], np.array([0, 0])):
        translation = -segments[0][0, :]
        new_segments = [seg + translation for seg in segments]
    else:
        new_segments = segments
    # print(segments)
    # print(new_segments)
    ref = []
    read = []
    new_segments = filter_merge_segments(new_segments)
    for seg in new_segments:
        ref.append(seg[0][0])
        read.append(seg[0][1])
        ref.append(seg[1][0])
        read.append(seg[1][1])
    sv_type = parse_ref_read_by_overlap_gap(ref, read)
    return sv_type




'''def bfs_rearrangement(target, reference):
    """
    使用广度优先搜索找到从参考序列到目标序列的最短操作序列。

    :param target: 目标序列
    :param reference: 参考序列
    :return: 描述重排事件的操作序列
    """
    # 初始化队列，每个元素是一个元组：(当前序列, 操作历史)
    queue = deque([(reference, [])])
    visited = set()  # 记录已访问的序列，避免重复处理

    while queue:
        current, history = queue.popleft()

        # 检查是否达到目标序列
        if current == target:
            return history

        # 防止无限循环
        if current in visited:
            continue
        visited.add(current)

        # 尝试所有可能的操作
        if len(current) >= 2:  # 删除和倒置至少需要2个字符
            for i in range(len(current)):
                for j in range(i + 2, len(current) + 1):
                    # 删除操作
                    new_seq = current[:i] + current[j:]
                    queue.append((new_seq, history + [f'deletion:{current[i:j]}']))

                    # 倒置操作
                    inverted_seq = current[:i] + current[i:j][::-1] + current[j:]
                    queue.append((inverted_seq, history + [f'inversion:{current[i:j]}']))

        # 复制操作
        for i in range(len(current)):
            for j in range(i + 1, len(current) + 1):
                duplicated_seq = current[:j] + current[i:j] + current[j:]
                queue.append((duplicated_seq, history + [f'duplication:{current[i:j]}']))

    # 如果未找到路径，则返回空列表
    return []


# 示例
reference = "ABCD"
target = "AD"
operations = bfs_rearrangement(target, reference)
print("操作序列:", operations)'''
