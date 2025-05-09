import os.path
import h5py
import cv2
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import KMeans


class Line:
    def __init__(self, id, ref_start, ref_end, read_start, read_end):

        self.id = id
        self.read_start = read_start
        self.read_end = read_end
        self.ref_start = ref_start
        self.ref_end = ref_end

        self.length_read = self.read_end - self.read_start
        self.length_ref = self.ref_end - self.ref_start
        self.actual_length = None
        # # generate fitted line's info
        self.poly_line = np.polyfit((self.read_start, self.read_end), (self.ref_start, self.ref_end), deg=1)

        self.slope, self.intercept = self.poly_line[0], self.poly_line[1]

    def set_actual_length(self, actual_length):
        self.actual_length = actual_length

    def update_by_linear_line(self, target_line):
        """
        update info by given a taget-linear line, the cords extend

        """
        self.read_start = min(self.read_start, target_line.read_start)
        self.read_end = max(self.read_end, target_line.read_end)

        self.ref_start = min(self.ref_start, target_line.ref_start)
        self.ref_end = max(self.ref_end, target_line.ref_end)

        self.length_read = self.read_end - self.read_start
        self.length_ref = self.ref_end - self.ref_start

        # # re-generate fitted line's info
        self.poly_line = np.polyfit((self.read_start, self.read_end), (self.ref_start, self.ref_end), deg=1)

        self.slope, self.intercept = self.poly_line[0], self.poly_line[1]

    def is_linear_with(self, target_line):
        """
        use slope and intercept to calculate if the two lines are linear

        """

        # # STEP: generate a new line with both self and target line's cords
        new_poly_line = np.polyfit((self.read_start, self.read_end, target_line.read_start, target_line.read_end), (self.ref_start, self.ref_end, target_line.ref_start, target_line.ref_end), deg=1)

        new_slope, new_intercept = new_poly_line[0], new_poly_line[1]

        # # STEP: calculate the alteration ratio

        if new_slope - self.slope > 0.2:
            return False, -1, -1

        slop_alt_ratio = abs(new_slope - self.slope) / self.slope

        intercept_alt_ratio = abs((new_intercept - self.intercept) / self.intercept)

        if slop_alt_ratio > 0.1:
            return False, -1, -1

        if intercept_alt_ratio > 0.2:
            return False, -1, -1

        return True, slop_alt_ratio, intercept_alt_ratio

    def to_string(self):
        return "Line, ref: {}-{}, read: {}-{}, length: {}, slope: {}, intercept: {}".format(self.ref_start, self.ref_end, self.read_start, self.read_end, self.length, round(self.slope, 2), round(self.intercept, 2))


def lsd_line_detection_for_inv(img_file, min_line_length=10, min_line_gap=10):
    """
        Detect line via LineSegmentDetector
        """

    # # STEP: read img
    img = cv2.imread(img_file)
    img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # # Apply sharpening filter
    # kernel = np.array([[-1, -1, -1], [-1, 9, -1], [-1, -1, -1]])
    # img_sharpened = cv2.filter2D(img_gray, -1, kernel)

    # Apply edge enhancement
    img_enhanced = cv2.GaussianBlur(img_gray, (3, 3), 0)
    img_enhanced = cv2.addWeighted(img_enhanced, 1, img_enhanced, -0.1, 0)

    # Create default parametrization LSD
    lsd = cv2.createLineSegmentDetector(0)

    # Detect lines in the image
    lines = lsd.detect(img_enhanced)[0]  # Position 0 of the returned tuple are the detected lines

    # # STEP: convert results to objects
    line_objs = []

    line_id = 0
    if lines is not None:
        for line in lines:
            for ref_start, read_start, ref_end, read_end in line:
                if abs(read_end - read_start) > min_line_length and abs(ref_end - ref_start) > min_line_length:
                    line_objs.append(Line(line_id, int(ref_start), int(ref_end), int(read_start), int(read_end)))

                    line_id += 1

        return line_objs
    else:
        return line_objs


def merge_lines(line_objs):
    """
    判断并合并线段，如果线段数量大于2，对其进行判断并合并。
    返回合并后的线段列表。
    """
    if len(line_objs) < 2:
        # 如果线段数量小于2，即只有1条线，则直接返回，不做任何合并
        return line_objs

    merged_lines = []
    used_lines = set()  # 用于记录已经合并过的线段索引

    for i in range(len(line_objs)):
        if i in used_lines:
            continue  # 如果这个线段已经被合并过，就跳过

        base_line = line_objs[i]
        for j in range(i + 1, len(line_objs)):
            if j in used_lines:
                continue

            target_line = line_objs[j]
            can_merge, _, _ = base_line.is_linear_with(target_line)

            if can_merge:
                distance_read = abs(base_line.read_end - target_line.read_start)
                distance_ref = abs(base_line.ref_end - target_line.ref_start)
                if distance_read <= 15 and distance_ref <=15:  ##对两条共线的线段之间的距离进行判断，离得近就合并，
                    base_line.update_by_linear_line(target_line)
                    used_lines.add(j)

        merged_lines.append(base_line)

    return merged_lines

def lsd_line_detection(img_file, min_line_length=10, min_line_gap=10):
    """
    Detect line via LineSegmentDetector
    """

    # # STEP: read img
    line_objs = []

    if img_file is None:
        print("Error: img_file is None. Skipping line detection.")
        return line_objs # Return an empty list if img_file is None

    with h5py.File(img_file, 'r') as f:
        original_data = f['matrix'][:]
    original_data = (original_data.squeeze() * 255).astype(np.uint8)

    # Create default parametrization LSD
    lsd = cv2.createLineSegmentDetector(0)

    # Detect lines in the image
    lines = lsd.detect(original_data)[0]  # Position 0 of the returned tuple are the detected lines

    # # STEP: convert results to objects

    line_id = 0
    if lines is not None:
        for line in lines:
            for ref_start, read_start, ref_end, read_end in line:
                if read_end - read_start > min_line_length and ref_end - ref_start > min_line_length:
                    line_objs.append(Line(line_id, int(ref_start), int(ref_end), int(read_start), int(read_end)))
                    line_id += 1

        line_objs_merged = merge_lines(line_objs)

        return line_objs_merged
    else:
        return line_objs

def detect_lines_from_dotplots(sv_out_path):
    """
    detect lines from dotplots for each SV

    """

    cleared_img_path = os.path.join(sv_out_path, "predict/images")

    for img_file in os.listdir(cleared_img_path):
        if "dotplot-outputs" not in img_file:
            continue

        img_lines = lsd_line_detection(img_file)

        print(len(img_lines))

        for line in img_lines:
            print(line.to_string())


if __name__ == '__main__':
    img_file = "/mnt/d/data/validation/test_res_chr20/chr20-12622967-12625701-tDUP+DEL/predict/images/chr20-12620233-12628435-S1_2837.ref2readdotplot-outputs.png"

    lines = lsd_line_detection(img_file)

    for line in lines:
        print(line.to_string())