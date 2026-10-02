"""Persist the actual HE-to-canvas affine map, independently of log formatting."""
import hashlib
from pathlib import Path
import numpy as np


def image_signature(path):
    path = Path(path).resolve(strict=True)
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            digest.update(block)
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def registration_manifest(warp, resize_xy, raw_size_xy, registered_size_xy, source_signature, registered_image):
    """Compose OpenCV centre-aware resize with the FINAL (possibly refined) warp."""
    warp = np.asarray(warp, dtype=np.float64)
    if warp.shape != (2, 3) or not np.isfinite(warp).all():
        raise ValueError('Invalid final HE registration matrix')
    fx, fy = map(float, resize_xy)
    if min(fx, fy) <= 0:
        raise ValueError('Invalid HE resize factor')
    resize = np.array([[fx, 0, (fx-1)/2], [0, fy, (fy-1)/2], [0, 0, 1.]])
    matrix = np.vstack([warp, [0., 0., 1.]]) @ resize
    inverse = np.linalg.inv(matrix)
    return {
        'schema': 'umap_he_transform_v2',
        'provenance': 'pipeline_final_matrix',
        'coordinate_convention': 'zero-based pixel centres; x right, y down',
        'raw_he_size_xy': list(raw_size_xy),
        'registered_canvas_size_xy': list(registered_size_xy),
        'resize_factor_xy': [fx, fy],
        'resized_to_registered': np.vstack([warp, [0., 0., 1.]]).tolist(),
        'raw_he_to_registered': matrix.tolist(),
        'registered_to_raw_he': inverse.tolist(),
        'source_he': source_signature,
        'registered_image': image_signature(registered_image),
        'registration_accuracy': 'not measured by matrix serialization; validate anatomical landmarks',
    }
