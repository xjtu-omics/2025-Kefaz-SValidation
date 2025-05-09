import gc
from collections import deque, Counter
import numpy as np
from src.valid.encoder import encoder_seg
import psutil


def bfs_rearrangement_by_array(target, reference):
    # 计算 target 和 reference 的字符频率
    target_counter = Counter(target)
    reference_counter = Counter(reference)

    # 给 target 和 reference 去除均只出现一次的字符
    filtered_target = ''.join(
        [char for char in target if not (target_counter[char] == 1 and reference_counter[char] == 1)])
    filtered_reference = ''.join(
        [char for char in reference if not (reference_counter[char] == 1 and target_counter[char] == 1)])
    # 将输入字符串转换为字符列表
    target = list(filtered_target)
    reference = list(filtered_reference)
    queue = deque([(reference, [])])
    visited = {tuple(reference): []}

    # 获取当前进程
    process = psutil.Process()
    m = 0
    while queue:
        # print(m)
        m += 1
        # print(process.memory_info().rss / (1024 * 1024 * 1024)) # 将字节转换为GB)
        memory_usage = process.memory_info().rss / (1024 * 1024 * 1024)  # 将字节转换为GB
        if memory_usage > 100:
            raise MemoryError("Momery usage exceeded 100 GB")
        current, history = queue.popleft()

        if current == target:
            result = history
            del queue
            del visited
            #返回结果之前手动进行回收
            gc.collect()
            return result

        for i in range(len(current)):
            # 删除操作，注意使用[:]创建列表的副本以避免直接修改原列表
            new_seq = current[:i] + current[i + 1:]
            # if tuple(new_seq) == "ABDEGD":
            #     break
            new_history = history + [f'deletion:{current[i]}']
            if tuple(new_seq) not in visited or len(new_history) < len(visited[tuple(new_seq)]):
                queue.append((new_seq, new_history))
                visited[tuple(new_seq)] = new_history

        for i in range(len(target)):
            if target[i] not in current:
                for j in range(len(current) + 1):
                    # 插入操作
                    new_seq = current[:j] + [target[i]] + current[j:]
                    new_history = history + [f'insertion:{target[i]}']
                    # 检查新状态是否在已访问的序列中或找到更短的路径
                    if tuple(new_seq) not in visited or len(new_history) < len(visited[tuple(new_seq)]):
                        queue.append((new_seq, new_history))
                        visited[tuple(new_seq)] = new_history

        for i in range(len(current)):
            for j in range(i + 1, len(current) + 1):
                # 复制操作
                duplicated_seq = current[:j] + current[i:j] + current[j:]
                new_history = history + [f'duplication:{"".join(current[i:j])}']
                # 检查新状态是否在已访问的序列中或找到更短的路径
                if tuple(duplicated_seq) not in visited or len(new_history) < len(visited[tuple(duplicated_seq)]):
                    queue.append((duplicated_seq, new_history))
                    visited[tuple(duplicated_seq)] = new_history

    ##返回空之前进行一次垃圾回收
    gc.collect()
    return []

def move2zeros(segments):
    # 对坐标进行平移，使之从0，0开始
    if not np.array_equal(np.array(segments[0][0]), np.array([0, 0])):
        translation = -segments[0][0]
        new_segments = [seg + translation for seg in segments]
    else:
        new_segments = segments
        translation = np.array([0, 0])

    return new_segments, translation

def get_sv_region(ref_start, ref_end, ref, read, sv_type):
    # sv_region = None
    start = ref_start
    end = 0
    #只涉及到一种类型的变异
    if len(sv_type) == 1:
        sv_operation, sv_segment = sv_type[0].split(":")  # 解构操作和段落
        if sv_operation == "deletion":
            for key, value in ref.items():
                if value == sv_segment:
                    start += key[0]
                    end = start + abs(key[1] - key[0])
                    break
        elif sv_operation == "insertion":
            for key, value in read.items():
                if value == sv_segment:
                    start += key[0]
                    end = start + abs(key[1] - key[0])
                    break
        elif sv_operation == "duplication":
            for key, value in read.items():
                if value == sv_segment:
                    start += key[0]
                    end = start + abs(key[1] - key[0])
                    break
    #涉及到多种类型的变异的时候，第一个片段更新start，后面的片段只叠加更新end即可
    else:
        for idx, sv in enumerate(sv_type):
            sv_operation, sv_segment = sv.split(":")  # 解构操作和段落
            #pre_key, pre_value = None, None
            if sv_operation == "deletion":
                for key, value in ref.items():
                    if value == sv_segment:
                        if idx == 0:
                            start += key[0]
                            end = start + abs(key[1] - key[0])
                            break
                        else:
                            #如果不是第一段变异，则只更新end
                            end += abs(key[1] - key[0])
                            break

            elif sv_operation == "insertion":
                for key, value in read.items():
                    if value == sv_segment:
                        if idx == 0:
                            start += key[0]
                            end = start + abs(key[1] - key[0])
                            break
                        else:
                            #如果不是第一段变异，则只更新end
                            end += abs(key[1] - key[0])
                            break

            elif sv_operation == "duplication":
                for key, value in read.items():
                    if value == sv_segment:
                        if idx == 0:
                            start += key[0]
                            end = start + abs(key[1] - key[0])
                            break
                        else:
                            #如果不是第一段变异，则只更新end
                            end += abs(key[1] - key[0])
                            break

    sv_region = (start, end)
    return sv_region

def search_segments(segments, ref_start, ref_end):
    new_segments, translation = move2zeros(segments)
    ref, read = encoder_seg(new_segments)
    ref_letters = "".join(ref.values())
    read_letters = "".join(read.values())
    print(ref_letters)
    print(read_letters)
    sv_type = bfs_rearrangement_by_array(read_letters, ref_letters)

    if len(sv_type) > 0:
        ref_start = ref_start - translation[0]
        ref_end = ref_end - translation[0]
        sv_region = get_sv_region(ref_start, ref_end, ref, read, sv_type)
    else:
        sv_region = None
    return ref, read, sv_type, sv_region



if __name__ == "__main__":
    # # ref = [(0, 569),(700, 1114)]
    # ref = [(0, 873),(656, 1099)]
    # en_ref = encode_segements(ref)
    # # read = [(0, 568),(552, 966)]
    # read = [(0, 871),(886, 1331)]
    # en_read = encode_segements(read)
    # sv_type = bfs_rearrangement(en_read, en_ref)
    # segments = [[(0, 0), (670, 670)], [(579, 521), (1038, 980)]]
    # segments = [[(0, 0), (658, 655)], [(550, 638), (753, 842)],[(765, 869), (1050, 1152)]]
    segments = [[(0, 0), (250, 250)], [(100, 150), (350, 400)]]
    en_ref, en_read, sv_type = search_segments(segments)
    sv_type = bfs_rearrangement_by_array(en_read, en_ref )
    print(sv_type)