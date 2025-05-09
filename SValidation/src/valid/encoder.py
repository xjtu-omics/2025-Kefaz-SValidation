import string


def create_uppercase_letters_A_to_J():
    all_uppercase_letters = string.ascii_uppercase
    uppercase_letters_A_to_J = list(all_uppercase_letters[:10])
    return uppercase_letters_A_to_J


def calculate_read_coord(current, ref_points):
    eps = 1e-5
    #斜率
    m = (current[0][1] - current[1][1]) / (current[0][0] - current[1][0] + eps)
    #截距
    b = current[1][1] - current[1][0] * m

    return int(m * ref_points + b)


def calculate_ref_coord(current, read_points):
    eps = 1e-5
    #斜率
    m = (current[0][1] - current[1][1]) / (current[0][0] - current[1][0] + eps)
    #截距
    b = current[1][1] - current[1][0] * m

    return int((read_points - b) / m)
def is_fig_f_or_g(current, next):
    "return 0 --> fig.f, return 1 --> fig.g"
    eps = 1e-5
    m1 = (current[0][1] - current[1][1] + eps) / (current[0][0] - current[1][0] + eps)
    m2 = (next[0][1] - next[1][1] + eps) / (next[0][0] - next[1][0] + eps)
    b1 = current[1][1] - current[1][0] * m1
    b2 = next[1][1] - next[1][0] * m2
    if b1 < b2:
        return 1
    elif b1 > b2:
        return 0
    else:
        return -1

def contains(segments, new_seg, min_gap=50):
    "判断当前片段是否存在，若存在则返回False，若不存在，则返回True"
    if abs(new_seg[1] - new_seg[0]) < min_gap:
        return False, len(segments.keys())
    else:
        for i in range(len(segments.keys())):
            seg = list(segments.keys())[i]
            #如果发生重叠，则判定为存在
            if (abs(seg[0] - new_seg[0]) <= min_gap ) and (abs(seg[1] - new_seg[1]) <= min_gap):
                return False, i
        return True, len(segments.keys())
def parse_ref_read_from_encodered_segs(new_segments, min_gap=50):
    uppercase_letters_A_to_J = create_uppercase_letters_A_to_J()
    ref = {}
    read = {}
    for i in range(len(new_segments)):

        if len(ref) == 0 and len(read) == 0:
            ref_seg = (new_segments[0][0][0], new_segments[0][1][0])
            ref [ref_seg] = uppercase_letters_A_to_J[i]
            read_seg = (new_segments[0][0][1], new_segments[0][1][1])
            read [read_seg] = uppercase_letters_A_to_J[i]
        else:
            ref2add = (new_segments[i][0][0], new_segments[i][1][0])
            read2add = (new_segments[i][0][1], new_segments[i][1][1])
            ref_flag, ref_idx = contains(ref, ref2add, min_gap)
            read_flag, read_idx = contains(read, read2add, min_gap)
            #ref和 read均不存在，则直接进行添加
            if ref_flag and read_flag:
                ref [ref2add] = uppercase_letters_A_to_J[i]
                read [read2add] = uppercase_letters_A_to_J[i]
            # ref方向不存在，read方向长度不够或者已经存在
            elif ref_flag and (not read_flag):
                ref [ref2add] = uppercase_letters_A_to_J[read_idx]
            # read方向不存在，但是ref方向已经存在或者长度不够
            elif read_flag and (not ref_flag):
                read[read2add] = uppercase_letters_A_to_J[ref_idx]
            else:
                continue

    return ref, read


def encoder_seg(segments):
    min_gaps = 25
    #初始原点坐标进行添加
    #segments中，每一个list记录了一条线的起点坐标和终点坐标
    new_segments = [segments[0]]
    idx = 1
    while(idx < len(segments)):
        current = new_segments[-1]
        next = segments[idx]
        #沿着ref方向遍历，寻找gaps, next[0]表示线段的起点坐标，next[0][0]表示起点坐标的横坐标
        #ref存在，但是read方向不存在,即严格对应于fig（a）
        if (next[0][0] - current[1][0] >= min_gaps) and abs(next[0][1] - current[1][1]) < min_gaps:
            #将next的起点坐标和current的终点坐标进行添加，作为一条新的平行于ref的线
            new_segments.append([current[1],next[0]])
            idx -= 1
        #read find gaps, but ref don't，即严格对应于fig（b）
        elif (next[0][1] - current[1][1] >= min_gaps) and abs(next[0][0] - current[1][0]) < min_gaps:
            new_segments.append([current[1], next[0]])
            idx -= 1
        # both ref and read find gaps，即不严格对应于fig(a) or fig(b)
        elif (next[0][0] - current[1][0] >= min_gaps) and (next[0][1] - current[1][1] >= min_gaps):
            ref_gap = next[0][0] - current[1][0]
            read_gap = next[0][1] - current[1][1]
            if ref_gap > read_gap:
                new_segments.append([current[1], (calculate_ref_coord(next, current[1][1]), current[1][1])])
                new_segments.append([(calculate_ref_coord(next, current[1][1]), current[1][1]), next[0]])
            elif ref_gap < read_gap:
                new_segments.append([current[1], (current[1][0], calculate_read_coord(next, current[1][0]))])
                new_segments.append([(current[1][0], calculate_read_coord(next, current[1][0])), next[0]])
            else:
                middle = min(next[0][0], next[0][1])
                new_segments.append([current[1], (middle, middle)])
                new_segments.append([(middle, middle), next[0]])
            idx -= 1
        # only ref find overlap，read don't, 即严格对应于fig（d，e）
        elif (current[1][0] - next[0][0] > min_gaps) and (abs(current[1][1] - next[0][1]) <= min_gaps):
            new_segments[-1] = [current[0],(next[0][0], calculate_read_coord(current, next[0][0]))]
            #添加第二条线是进行判断，以区分串联重复和分散重复
            #通过当前seg的终点横坐标和next seg的终点横坐标进行判断
            #当前终点的横坐标大于next终点的横坐标，则对应于分散情形, fig d
            if current[1][0] > next[1][0]:
                new_segments.append([(next[0][0], calculate_read_coord(current, next[0][0])), (next[1][0], calculate_read_coord(current, next[1][0]))])
                new_segments.append([(next[1][0], calculate_read_coord(current, next[1][0])), current[1]])
                new_segments.append(next)
            #else情况，则当前终点的横坐标小于next终点的横坐标
            else:
                new_segments.append([(next[0][0], calculate_read_coord(current, next[0][0])),current[1]])
                new_segments.append([next[0], (current[1][0], calculate_read_coord(next, current[1][0]))])
                new_segments.append([(current[1][0], calculate_read_coord(next, current[1][0])), next[1]])
        ## ref find overlaps but read find gaps
        elif (current[1][0] - next[0][0] > min_gaps) and (next[0][1] - current[1][1] > min_gaps):
            ##首先添加第一条线
            new_segments[-1] = [current[0], (next[0][0], calculate_read_coord(current, next[0][0]))]
            ## 分散情形
            if current[1][0] > next[1][0]:
                new_segments.append([(next[0][0], calculate_read_coord(current, next[0][0])), (next[1][0], calculate_read_coord(current, next[1][0]))])
                new_segments.append([(next[1][0], calculate_read_coord(current, next[1][0])), current[1]])
                new_segments.append([current[1], (current[1][0], next[0][1])])  ## add gaps
                new_segments.append(next)

            ## 串联情形
            else:
                new_segments.append([(next[0][0], calculate_read_coord(current, next[0][0])), current[1]])
                new_segments.append([current[1], (current[1][0], next[0][1])]) ## add gaps
                new_segments.append([next[0], (current[1][0], calculate_read_coord(next, current[1][0]))])
                new_segments.append([(current[1][0], calculate_read_coord(next, current[1][0])), next[1]])

        # only read find overlaps, and ref find gaps, fig h
        elif (current[1][1] - next[0][1] > min_gaps) and (next[0][0] - current[1][0] > min_gaps):
            new_segments[-1] = [current[0],(calculate_ref_coord(current, next[0][1]), next[0][1])]
            new_segments.append([(calculate_ref_coord(current, next[0][1]), next[0][1]), current[1]])
            new_segments.append([current[1], (next[0][0], current[1][1])])
            new_segments.append([next[0], (calculate_ref_coord(next, current[1][1]), current[1][1])])
            new_segments.append([(calculate_ref_coord(next, current[1][1]), current[1][1]), next[1]])
        # ref and read find overlaps, fig.f
        elif  (current[1][1] - next[0][1] > min_gaps) and is_fig_f_or_g(current, next) == 0:
        # elif (current[1][0] - next[0][0] > 50) and is_fig_f_or_g(current, next) == 0:
            new_segments[-1] = [current[0], (calculate_ref_coord(current, next[0][1]), next[0][1])]
            new_segments.append([(calculate_ref_coord(current, next[0][1]), next[0][1]), (next[0][0], calculate_read_coord(current, next[0][0]))])
            new_segments.append([(next[0][0], calculate_read_coord(current, next[0][0])), current[1]])
            new_segments.append([current[1], (calculate_ref_coord(next, current[1][1]), current[1][1])])
            new_segments.append([(calculate_ref_coord(next, current[1][1]), current[1][1]), next[1]])
            # new_segments.append([(calculate_ref_coord(current, next[0][1]), next[0][1]),current[1]])
            # new_segments.append([next[0], (calculate_ref_coord(next, current[1][1]), current[1][1])])
            # new_segments.append([(calculate_ref_coord(next, current[1][1]), current[1][1]), next[1]])
        #
        elif (current[1][0] - next[0][0] > min_gaps) and (current[1][1] - next[0][1] > min_gaps) and is_fig_f_or_g(current, next) == 1:
            new_segments[-1] = [current[0], (next[0][0], calculate_read_coord(current, next[0][0]))]
            new_segments.append([(next[0][0], calculate_read_coord(current, next[0][0])), (calculate_ref_coord(current, next[0][1]), next[0][1])])
            new_segments.append([(calculate_ref_coord(current, next[0][1]), next[0][1]), current[1]])
            new_segments.append([current[1], (current[1][0],calculate_read_coord(next, current[1][0]))])
            new_segments.append([(current[1][0],calculate_read_coord(next, current[1][0])), next[1]])
            # new_segments.append([(next[0][0], calculate_read_coord(current, next[0][0])), current[1]])
            # new_segments.append([next[0], (current[1][0], calculate_read_coord(next, current[1][0]))])
            # new_segments.append([(current[1][0], calculate_read_coord(next, current[1][0])), next[1]])
        else:
            new_segments.append(next)

        idx += 1


    ref, read = parse_ref_read_from_encodered_segs(new_segments, min_gap=min_gaps)
    print(ref)
    print(read)
    return ref, read




if __name__ == "__main__":
    segments = [[(0, 0), (250, 250)], [(100, 150),(350, 400)]]
    encoder_seg(segments)
