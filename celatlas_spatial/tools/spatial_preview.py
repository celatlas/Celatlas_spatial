import numpy as np


def _resize_mask_nearest(mask, target_shape):
    target_height, target_width = target_shape
    source_height, source_width = mask.shape[:2]
    row_indices = np.floor(np.arange(target_height) * source_height / target_height).astype(int)
    column_indices = np.floor(np.arange(target_width) * source_width / target_width).astype(int)
    return mask[row_indices[:, None], column_indices]


def apply_tissue_mask(image, mask):
    """Return an image with pixels outside the tissue mask set to black."""
    if mask.shape[:2] != image.shape[:2]:
        mask = _resize_mask_nearest(mask, image.shape[:2])
    masked_image = image.copy()
    masked_image[mask <= 0] = 0
    return masked_image


def blend_heatmap_with_mask(background, heatmap, mask, background_weight=0.6):
    """Blend a heatmap only inside the final tissue mask."""
    if background.shape[:2] != heatmap.shape[:2]:
        raise ValueError("background and heatmap must use the same canvas")
    if mask.shape[:2] != background.shape[:2]:
        mask = _resize_mask_nearest(mask, background.shape[:2])

    background_weight = float(background_weight)
    if not 0.0 <= background_weight <= 1.0:
        raise ValueError("background_weight must be between 0 and 1")

    output = background.copy()
    keep = mask > 0
    if not np.any(keep):
        return output
    blended = (
        background[keep].astype(np.float32) * background_weight
        + heatmap[keep].astype(np.float32) * (1.0 - background_weight)
    )
    output[keep] = np.clip(np.rint(blended), 0, 255).astype(np.uint8)
    return output
