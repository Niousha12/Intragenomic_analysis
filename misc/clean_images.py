from PIL import Image
import os


def process_images(directory, crop_margins):  # , output_size):
    """
    Process images in the specified directory: crop them based on given margins and resize them.

    Args:
    directory (str): Path to the directory containing the images.
    crop_margins (dict): A dictionary with keys 'left', 'right', 'top', 'bottom' indicating the percentage
                         to crop from each side.
    output_size (tuple): The desired (width, height) to resize the images to after cropping.
    """
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp')):
            file_path = os.path.join(directory, filename)
            with Image.open(file_path) as img:
                width, height = img.size

                # Calculate crop dimensions
                left_crop = width * crop_margins['left']
                right_crop = width * crop_margins['right']
                top_crop = height * crop_margins['top']
                bottom_crop = height * crop_margins['bottom']

                # Calculate coordinates to crop the image
                left = left_crop
                top = top_crop
                right = width - right_crop
                bottom = height - bottom_crop

                # Crop the image
                img_cropped = img.crop((left, top, right, bottom))

                # Resize the image
                # img_resized = img_cropped.resize(output_size, Image.ANTIALIAS)

                # Save the image back to the directory or another directory
                process_path = os.path.join(directory, f"processed")
                if not os.path.exists(process_path):
                    os.makedirs(process_path)
                img_cropped.save(os.path.join(process_path, f"processed_{filename}"))


# Specify the directory containing the images
image_directory = "../Figures/Representative/maize/Different_centroids"


# Specify the crop margins as percentages of each side
crop_margins = {
    'left': 0.12,  # 5% from the left
    'right': 0.095,  # 10% from the right
    'top': 0.11,  # 5% from the top
    'bottom': 0.095  # 10% from the bottom
}

# Specify the desired size (width, height) after cropping and resizing
# desired_size = (189, 100)  # example size: 200x200 pixels

process_images(image_directory, crop_margins)  # , desired_size)
