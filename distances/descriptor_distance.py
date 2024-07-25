import math

from skimage.metrics import structural_similarity, normalized_root_mse
import numpy as np


def split_image(image, split_size=14):
    h, w = image.shape[0], image.shape[1]
    col_count = int(math.ceil(h / split_size))
    row_count = int(math.ceil(w / split_size))
    tiles_count = col_count * row_count
    tiles = np.zeros((tiles_count, split_size, split_size))
    for y in range(col_count):
        for x in range(row_count):
            ind = x + (y * row_count)
            tiles[ind:(ind + 1), :, :] = image[split_size * y: (y + 1) * split_size,
                                         split_size * x:(x + 1) * split_size]

    return tiles


def get_descriptor(patch, bin_bounds):
    descriptor = np.zeros(len(bin_bounds))

    for index, bin_point in enumerate(bin_bounds):
        if index < len(bin_bounds) - 1:
            low = bin_bounds[index]
            high = bin_bounds[index + 1]
        else:
            low = bin_bounds[index]
            high = np.inf

        descriptor[index] = ((low <= patch) & (patch < high)).sum()
    descriptor = descriptor / np.sum(descriptor)
    descriptor = list(descriptor)
    return descriptor


def descriptor_distance(img1, img2, m=2, bins_bound=None):
    p1 = split_image(img1, 2 ** m)
    p2 = split_image(img2, 2 ** m)

    sub_matrices = p1.shape[0]

    vec1 = []
    vec2 = []
    for i in range(sub_matrices):
        vec1 += get_descriptor(patch=p1[i], bin_bounds=bins_bound)
        vec2 += get_descriptor(patch=p2[i], bin_bounds=bins_bound)

    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)

    denom_1 = np.sqrt(np.mean((vec1 * vec1), dtype=np.float64))
    denom_2 = np.sqrt(np.mean((vec2 * vec2), dtype=np.float64))
    if denom_1 > denom_2:
        distance = normalized_root_mse(vec1, vec2, normalization='euclidean')
    else:
        distance = normalized_root_mse(vec2, vec1, normalization='euclidean')

    # distance = normalized_root_mse(vec1, vec2, normalization='euclidean')
    # distance = np.sqrt(np.sum((vec1 - vec2) ** 2))
    return distance


if __name__ == '__main__':
    img = np.arange(0, 64).reshape(8, 8)
    img2 = np.arange(4, 68).reshape(8, 8)
    dist = descriptor_distance(img, img2, 2, [0, 8, 16])
    print(dist)
