import os

import numpy as np
import tensorflow as tf
#import tensorflow as tf
tf.compat.v1.disable_eager_execution()
from tensorflow import image



def create_op(func, **placeholders):
    op = func(**placeholders)

    def f(**kwargs):
        feed_dict = {}
        for argname, argvalue in kwargs.items():
            placeholder = placeholders[argname]
            feed_dict[placeholder] = argvalue
        return tf.get_default_session().run(op, feed_dict=feed_dict)

    return f

pad = create_op(
    func=tf.image.pad_to_bounding_box,
    image=tf.compat.v1.placeholder(tf.float32),
    #image=tf.placeholder(tf.float32),
    offset_height=tf.compat.v1.placeholder(tf.int32, []),
    offset_width=tf.compat.v1.placeholder(tf.int32, []),
    target_height=tf.compat.v1.placeholder(tf.int32, []),
    target_width=tf.compat.v1.placeholder(tf.int32, []),
)
crop = create_op(
    func=tf.image.crop_to_bounding_box,
    image=tf.compat.v1.placeholder(tf.float32),
    offset_height=tf.compat.v1.placeholder(tf.int32, []),
    offset_width=tf.compat.v1.placeholder(tf.int32, []),
    target_height=tf.compat.v1.placeholder(tf.int32, []),
    target_width=tf.compat.v1.placeholder(tf.int32, []),
)

grayscale_to_rgb = create_op(
    func=tf.image.grayscale_to_rgb,
    images=tf.compat.v1.placeholder(tf.float32),
)

downscale = create_op(
    func=tf.image.resize,
    images=tf.compat.v1.placeholder(tf.float32, [None, None, None]),
    size=tf.compat.v1.placeholder(tf.int32, [2]),
    method=tf.image.ResizeMethod.AREA,
)

upscale = create_op(
    func=tf.image.resize,
    images=tf.compat.v1.placeholder(tf.float32, [None, None, None]),
    size=tf.compat.v1.placeholder(tf.int32, [2]),
    method=tf.image.ResizeMethod.BICUBIC,
)

def custom_resize(src, re_size):

    target_img = np.ones((re_size, re_size, 3)) * 255

    target_img_size = np.shape(target_img)
    base_img_size = np.shape(src)

    target_img[0: min(target_img_size[0], base_img_size[0]), 0: min(target_img_size[1], base_img_size[1]), :] = src[0: min(target_img_size[0], base_img_size[0]), 0: min(target_img_size[1], base_img_size[1]), :]

    return target_img


def resize(src, re_size=256, is_pad=False):

    height, width, _ = src.shape
    dst = src
    if height != width:
        if is_pad:
            size = max(height, width)
            # pad to correct ratio
            oh = (size - height) // 2
            ow = (size - width) // 2
            dst = pad(image=dst, offset_height=oh, offset_width=ow, target_height=size, target_width=size)
        else:
            # crop to correct ratio
            size = min(height, width)
            oh = (height - size) // 2
            ow = (width - size) // 2
            dst = crop(image=dst, offset_height=oh, offset_width=ow, target_height=size, target_width=size)

    assert(dst.shape[0] == dst.shape[1])

    size, _, _ = dst.shape
    if size > re_size:
        dst = downscale(images=dst, size=[re_size, re_size])
    elif size < re_size:
        dst = upscale(images=dst, size=[re_size, re_size])

    return dst


def combine(src, src_pair):

    # make sure that dimensions are correct
    height, width, _ = src.shape
    if height != src_pair.shape[0] or width != src_pair.shape[1]:
        raise Exception("differing sizes")

    # convert both images to RGB if necessary
    if src.shape[2] == 1:
        src = grayscale_to_rgb(images=src)

    if src_pair.shape[2] == 1:
        src_pair = grayscale_to_rgb(images=src_pair)

    # remove alpha channel
    if src.shape[2] == 4:
        src = src[:, :, :3]

    if src_pair.shape[2] == 4:
        src_pair = src_pair[:, :, :3]

    dst = np.concatenate([src, src_pair], axis=1)

    return dst
