import numpy as np
import torch
from torchvision import transforms
from matplotlib import pyplot as plt
from PIL import Image
import lpips


def get_triplet_distance(cgr_1, cgr_2, model):
    transform = transforms.Compose([
        transforms.Resize((128, 128)),
        transforms.Grayscale(num_output_channels=1),
        transforms.ToTensor()
    ])

    saving_path_1 = f"../temp/Triplet/cgr_1.png"
    plt.imsave(saving_path_1, cgr_1, cmap='Blues')

    saving_path_2 = f"../temp/Triplet/cgr_2.png"
    plt.imsave(saving_path_2, cgr_2, cmap='Blues')

    image_1 = Image.open(saving_path_1).convert('RGB')
    image_2 = Image.open(saving_path_2).convert('RGB')

    x1 = transform(image_1).unsqueeze(0)
    x2 = transform(image_2).unsqueeze(0)

    with torch.no_grad():
        emb1 = model(x1)
        emb2 = model(x2)

    distance = (emb1 - emb2).pow(2).sum(1)
    return distance.item()


# LPIPS original
def prepare_image_for_lpips(image):
    # Convert image to 3 channels and normalize
    image = np.stack([image, image, image], axis=-1)
    max_value = np.max(image) / 2
    # image = (image / max_value - 1.0)  # convert image between -1 and 1
    # Convert to tensor and add batch dimension
    image = torch.tensor(image).permute(2, 0, 1).unsqueeze(0).float()
    return image


def lpips_dist(cgr_1, cgr_2):
    image1 = prepare_image_for_lpips(cgr_1)
    image2 = prepare_image_for_lpips(cgr_2)

    loss_fn_alex = lpips.LPIPS(net='alex')
    d = loss_fn_alex(image1, image2).squeeze()
    d = d.item()
    return d
    # print(f"lpips distance: {d}")


# LPIPS fine tuned
def lpips_fine_tuned(cgr_1, cgr_2, lpip_model):
    transform = transforms.Compose([
        transforms.Resize((256, 256)),
        transforms.ToTensor()
    ])

    saving_path_1 = f"../temp/Triplet/cgr_1.png"
    plt.imsave(saving_path_1, cgr_1, cmap='Blues')

    saving_path_2 = f"../temp/Triplet/cgr_2.png"
    plt.imsave(saving_path_2, cgr_2, cmap='Blues')

    image_1 = Image.open(saving_path_1).convert('RGB')
    image_2 = Image.open(saving_path_2).convert('RGB')

    x1 = transform(image_1).unsqueeze(0)
    x2 = transform(image_2).unsqueeze(0)

    lpips_distance = lpip_model(x1, x2).squeeze()
    return lpips_distance.item()
