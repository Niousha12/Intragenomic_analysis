import os

import numpy as np
import torch
import lpips
from skimage.metrics import structural_similarity, normalized_root_mse
from scipy import stats

from distances.descriptor_distance import descriptor_distance
from distances.global_distance import calculate_global_distance
from distances.ml_distance import get_triplet_distance, lpips_dist, lpips_fine_tuned

from chaos_game_representation import CGR

# epoch = 150
#
# root_folder = os.path.dirname(__file__)
# checkpoint_path_lpips = os.path.join(root_folder, "checkpoints", f"model_{epoch}.pth")
# lpip_model = lpips.LPIPS(net='vgg')
# lpip_model.load_state_dict(torch.load(checkpoint_path_lpips, map_location=torch.device('cpu')))

DISTANCE_PARAM_DICT = {"Euclidean": {'prob': False, 'require_norm': True},
                       "Cosine": {'prob': False, 'require_norm': True},
                       "Manhattan": {'prob': False, 'require_norm': True},
                       "Descriptor": {'prob': False, 'require_norm': True},
                       "DSSIM": {'prob': False, 'require_norm': True},
                       "K-S": {'prob': True, 'require_norm': False},
                       "Wasserstein": {'prob': True, 'require_norm': False},
                       "LPIPS": {'prob': False, 'require_norm': True}}


def get_dist_ensemble(cgr_1, cgr_2, dist_list=None, weight=None):
    assert dist_list is not None
    assert weight is not None
    assert len(weight) == len(dist_list)
    weight = np.asarray(weight)
    weight = weight / np.sum(weight)
    overall_dist = 0
    for i, dist in enumerate(dist_list):
        overall_dist += weight[i] * get_dist(cgr_1, cgr_2, dist)
    return overall_dist


def get_dist(cgr_img1, cgr_img2, dist_m="DSSIM"):
    if DISTANCE_PARAM_DICT[dist_m]['prob']:
        cgr_img1 = cgr_img1 / np.sum(cgr_img1)
        cgr_img2 = cgr_img2 / np.sum(cgr_img2)
    if DISTANCE_PARAM_DICT[dist_m]['require_norm']:
        cgr_img1 = CGR.normalize(cgr_img1)
        cgr_img2 = CGR.normalize(cgr_img2)

    distance = None
    if dist_m.lower() == "cosine":
        if cgr_img2.sum() == 0.0:
            cgr_img2 = cgr_img2 + 0.000000000000001
        distance = 1 - (np.dot(cgr_img2.reshape(-1), cgr_img1.reshape(-1)) / (
                np.linalg.norm(cgr_img2.reshape(-1)) * np.linalg.norm(cgr_img1.reshape(-1))))

    elif (dist_m.lower() == "euclidean") or (dist_m == "normalized_euclidean"):
        denom_1 = np.sqrt(np.mean((cgr_img1 * cgr_img1), dtype=np.float64))
        denom_2 = np.sqrt(np.mean((cgr_img2 * cgr_img2), dtype=np.float64))
        if denom_1 > denom_2:
            distance = normalized_root_mse(cgr_img1, cgr_img2, normalization='euclidean')
        else:
            distance = normalized_root_mse(cgr_img2, cgr_img1, normalization='euclidean')

    elif dist_m.lower() == "dssim":
        (similarity_score, differ) = structural_similarity(cgr_img1, cgr_img2, full=True, data_range=1)
        distance = 1 - similarity_score

    elif dist_m.lower() == "manhattan":
        distance = np.sum(np.abs(cgr_img2.reshape(-1) - cgr_img1.reshape(-1)))

    elif (dist_m.lower() == "k-s") or (dist_m.lower() == "ks"):
        temp = stats.kstest(cgr_img2.reshape(-1), cgr_img1.reshape(-1))
        distance = temp.statistic

    elif dist_m.lower() == "descriptor":
        bound = [0.0, 0.0001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5]
        distance = descriptor_distance(cgr_img1, cgr_img2, m=4, bins_bound=bound)

    elif dist_m.lower() == "wasserstein":
        points = np.array([i for i in range(cgr_img1.reshape(-1).shape[0])])
        points = points / np.max(points)
        distance = stats.wasserstein_distance(points, points, u_weights=cgr_img1.reshape(-1),
                                              v_weights=cgr_img2.reshape(-1))

    elif dist_m.lower() == "aid":
        aid_f1 = np.sum(cgr_img1[cgr_img1 > 0])
        aid_f2 = np.sum(cgr_img2[cgr_img2 > 0])
        aid_f_sum = np.sum((cgr_img1 + cgr_img2)[(cgr_img1 + cgr_img2) > 0])
        distance = 2 - ((aid_f1 + aid_f2) / aid_f_sum)

    elif dist_m.lower() == "global":
        distance = calculate_global_distance(cgr_img1, cgr_img2)

    elif dist_m.lower() == "lpips":
        distance = lpips_dist(cgr_img1, cgr_img2)

    elif dist_m.lower() == "lpips_trained":
        distance = lpips_fine_tuned(cgr_img1, cgr_img2, lpip_model)

    elif dist_m.lower() == "ensemble":
        distance = get_dist_ensemble(cgr_img1, cgr_img2, dist_list=["DSSIM", "Descriptor"], weight=[1, 1])

    elif dist_m.lower() == "triplet":
        distance = get_triplet_distance(cgr_img1, cgr_img2, triplet_model)

    return distance
