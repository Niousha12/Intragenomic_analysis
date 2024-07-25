import numpy as np


def calculate_global_distance(cgr_1, cgr_2):
    # Calculate normalization factor
    p = cgr_1.reshape(-1)
    q = cgr_2.reshape(-1)

    nw = np.sum(p * q)

    # Calculate means
    p_bar = np.sum(p ** 2 * q) / nw
    q_bar = np.sum(q ** 2 * p) / nw

    # Calculate variances
    sp = np.sum((p - p_bar) ** 2 * p * q) / nw
    sq = np.sum((q - q_bar) ** 2 * p * q) / nw

    # Calculate modified Pearson correlation coefficient
    rwpq = np.sum(((p - p_bar) / np.sqrt(sp)) * ((q - q_bar) / np.sqrt(sq)) * p * q) / nw

    # Calculate global distance
    d = 1 - rwpq

    return d
