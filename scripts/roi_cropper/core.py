"""Standalone, native-pixel scientific-image rotation/cropping. No Celatlas imports."""
from __future__ import annotations

import csv
import io
import json
import math
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import cv2
import numpy as np
import tifffile
from PIL import Image


class TiffSource:
    """Read supported scientific images through one region-source interface.

    TIFF is memory-mapped for large-image region access. PNG/JPEG are decoded
    once into memory because their compression does not support native random
    region reads without an additional pyramid backend.
    """
    def __init__(self, path):
        self.path = Path(path).expanduser().resolve(strict=True)
        suffix = self.path.suffix.lower()
        if suffix in ('.png', '.jpg', '.jpeg'):
            self._open_raster()
        elif suffix in ('.tif', '.tiff'):
            self._open_tiff()
        else:
            raise ValueError('支持 TIFF/BigTIFF、PNG 和 JPEG；其他格式请先用 Bio-Formats 转换。')
        self.signature = self.current_signature()

    def _open_raster(self):
        with Image.open(self.path) as image:
            if getattr(image, 'n_frames', 1) != 1:
                raise ValueError('仅支持单帧 PNG/JPEG。')
            orientation = image.getexif().get(274, 1)
            if int(orientation) != 1:
                raise ValueError('图像 Orientation 必须为 1；请先应用方向标记，避免坐标含义不明确。')
            if image.mode not in ('L', 'I;16', 'I;16L', 'I;16B', 'RGB', 'RGBA'):
                raise ValueError(f'不支持 {image.mode} 模式；PNG/JPEG 请使用 8/16-bit 灰度或 RGB(A)。')
            channels = len(image.getbands())
            itemsize = 2 if image.mode.startswith('I;16') else 1
            decoded_bytes = image.width * image.height * channels * itemsize
            if decoded_bytes > 2 * 1024**3:
                raise ValueError('压缩图解码后超过 2 GiB；请转换为未压缩 TIFF/BigTIFF，以便区域读取。')
            array = np.asarray(image)
            if array.dtype.kind != 'u' or array.dtype.itemsize not in (1, 2):
                raise ValueError('PNG/JPEG 解码后不是 uint8/uint16，拒绝隐式转换位深。')
            self.array = np.ascontiguousarray(array)
            self.height, self.width = self.array.shape[:2]
            self.extrasamples = (2,) if self.array.ndim == 3 and self.array.shape[2] == 4 else ()
            self.resolution_tags = {}
            if image.info.get('dpi'):
                self.resolution_tags['DPI'] = str(tuple(image.info['dpi']))
            self.storage_format = image.format or self.path.suffix.lstrip('.').upper()

    def _open_tiff(self):
        with tifffile.TiffFile(self.path) as tif:
            if len(tif.series) != 1 or len(tif.pages) != 1:
                raise ValueError('仅支持单张 TIFF；请先选定通道及 Z 投影，不接受多页/金字塔作为单张图。')
            p = tif.pages[0]
            if p.axes not in ('YX', 'YXS') or p.dtype not in (np.dtype('uint8'), np.dtype('uint16')):
                raise ValueError('支持 uint8/uint16 灰度或 RGB(A) TIFF。')
            if p.axes == 'YXS' and p.shape[2] not in (3, 4):
                raise ValueError('RGB 通道数必须是 3 或 4。')
            if ((p.axes == 'YX' and int(p.photometric) != 1)
                    or (p.axes == 'YXS' and int(p.photometric) != 2)):
                raise ValueError('仅支持黑底灰度或 RGB；不支持反相灰度、调色板及其他颜色空间。')
            self.extrasamples = tuple(int(value) for value in p.extrasamples)
            orientation = p.tags.get('Orientation')
            if orientation and int(orientation.value) != 1:
                raise ValueError('TIFF Orientation 必须为 1；请先统一方向，避免坐标含义不明确。')
            if (p.is_tiled or int(p.compression) != 1 or int(p.predictor) != 1
                    or (p.planarconfig is not None and int(p.planarconfig) != 1)):
                raise ValueError('快速模式要求未压缩、连续条带 TIFF。请用 Bio-Formats 导出未压缩 TIFF/BigTIFF。')
            dtype = np.dtype(p.dtype).newbyteorder(tif.byteorder)
            bits = p.bitspersample if isinstance(p.bitspersample, tuple) else (p.bitspersample,)
            if any(int(b) != dtype.itemsize * 8 for b in bits):
                raise ValueError('不支持打包位深 TIFF。')
            offsets, counts = np.asarray(p.dataoffsets), np.asarray(p.databytecounts)
            if len(offsets) == 0 or np.any(offsets[1:] != offsets[:-1] + counts[:-1]):
                raise ValueError('条带不连续，无法安全地快速区域读取；请重新导出未压缩 TIFF。')
            self.height, self.width = p.shape[:2]
            channels = p.shape[2] if p.axes == 'YXS' else 1
            row_bytes = self.width * channels * dtype.itemsize
            if (np.any(counts[:-1] != int(p.rowsperstrip) * row_bytes)
                    or counts[-1] < (self.height - (len(counts)-1)*int(p.rowsperstrip))*row_bytes):
                raise ValueError('TIFF 条带像素布局不一致。')
            total = int(counts.sum())
            if total % row_bytes or total // row_bytes < self.height:
                raise ValueError('TIFF 条带长度与画布不一致。')
            if int(offsets[-1] + counts[-1]) > self.path.stat().st_size:
                raise ValueError('TIFF 文件不完整。')
            shape = (total // row_bytes, self.width) + ((channels,) if p.axes == 'YXS' else ())
            self.array = np.memmap(self.path, dtype=dtype, mode='r', offset=int(offsets[0]), shape=shape)[:self.height]
            self.resolution_tags = {k: str(p.tags[k].value) for k in
                                    ('XResolution', 'YResolution', 'ResolutionUnit') if k in p.tags}
            self.storage_format = 'BigTIFF' if tif.is_bigtiff else 'TIFF'

    def current_signature(self):
        stat = self.path.stat()
        # Nanosecond timestamps exceed JavaScript's exact integer range.
        return {'path': str(self.path), 'bytes': stat.st_size, 'mtime_ns': str(stat.st_mtime_ns),
                'width': self.width, 'height': self.height, 'dtype': str(self.array.dtype),
                'channels': self.array.shape[2] if self.array.ndim == 3 else 1}

    def check_unchanged(self):
        if self.current_signature() != self.signature:
            raise ValueError('源图在打开后被修改，请重新打开后再操作。')

    def thumbnail(self, max_dim=2200, box=None):
        self.check_unchanged()
        x0, y0, x1, y1 = box or (0, 0, self.width, self.height)
        if not (0 <= x0 < x1 <= self.width and 0 <= y0 < y1 <= self.height):
            raise ValueError('预览区域超出原图。')
        scale = min(1.0, max_dim / max(x1-x0, y1-y0))
        nx, ny = max(1, round((x1-x0)*scale)), max(1, round((y1-y0)*scale))
        # Do not use np.ix_ advanced indexing on a memmap here.  On some
        # network filesystems that turns a 2k preview into millions of tiny
        # random reads.  Read decimated *contiguous rows* and resize once.
        row_step = max(1, int(math.ceil((y1-y0) / max(ny, 1))))
        sampled = np.asarray(self.array[y0:y1:row_step, x0:x1])
        if sampled.shape[0] != ny or sampled.shape[1] != nx:
            interpolation = cv2.INTER_AREA if sampled.shape[0] >= ny or sampled.shape[1] >= nx else cv2.INTER_LINEAR
            sampled = cv2.resize(sampled, (nx, ny), interpolation=interpolation)
        return display_rgb(np.ascontiguousarray(sampled))


def display_rgb(array):
    if array.dtype != np.uint8:
        lo, hi = np.percentile(array, [.1, 99.9])
        if hi <= lo:
            lo, hi = float(array.min()), float(array.max())
        array = np.clip((array.astype(np.float32)-lo) * (255/max(hi-lo, 1)), 0, 255).astype(np.uint8)
    if array.ndim == 2:
        return np.repeat(array[..., None], 3, axis=2)
    # The preview is annotated with OpenCV before JPEG export; force a
    # writable buffer when the source is a read-only memmap.
    return np.array(array[..., :3], copy=True, order='C')


def png_bytes(rgb):
    out = io.BytesIO()
    Image.fromarray(rgb).save(out, format='PNG')
    return out.getvalue()


def calibration_from_points(points, length_um):
    pts = np.asarray(points, dtype=float)
    length_um = float(length_um)
    if pts.shape != (2, 2) or not np.isfinite(pts).all() or not math.isfinite(length_um) or length_um <= 0:
        raise ValueError('需要两个有效标尺端点和正数长度（µm）。')
    distance = float(np.linalg.norm(pts[1]-pts[0]))
    if distance < 10:
        raise ValueError('标尺至少跨越 10 个原图像素；请放大后点击真正的标尺端点。')
    return {'method': 'scale_bar_two_points', 'points_source_px': pts.tolist(),
            'length_um': length_um, 'distance_px': distance, 'um_per_pixel': length_um/distance}


def make_recipe(source, params):
    """Coordinates are zero-based pixel centers; angles positive CCW on screen.

    Flips act on source deltas around the selected source center, before rotation.
    This keeps the selected raw region fixed while changing output orientation.
    """
    source.check_unchanged()
    keys = ('center_x', 'center_y', 'angle_deg', 'width_mm', 'height_mm', 'um_per_pixel')
    values = {k: float(params[k]) for k in keys}
    if not all(math.isfinite(v) for v in values.values()):
        raise ValueError('角度、位置、长度必须为有限数值。')
    if not params.get('calibration_confirmed'):
        raise ValueError('请先确认 µm/像素或完成标尺两点标定。')
    u = values['um_per_pixel']
    if not (0 < u <= 1000) or not (0 < values['width_mm'] <= 100 and 0 < values['height_mm'] <= 100):
        raise ValueError('像素标定或物理尺寸超出合理范围。')
    nx, ny = [int(math.floor(values[k]*1000/u + .5)) for k in ('width_mm', 'height_mm')]
    if min(nx, ny) < 1 or nx*ny > 2_000_000_000:
        raise ValueError('输出尺寸无效或超过 20 亿像素；请核对标尺单位。')
    angle = (values['angle_deg'] + 180) % 360 - 180
    flip_horizontal = bool(params.get('flip_horizontal', False))
    flip_vertical = bool(params.get('flip_vertical', False))
    image_role = str(params.get('image_role', 'unspecified')).lower()
    if image_role not in ('unspecified', 'dapi', 'fluorescence', 'he', 'other'):
        raise ValueError('图像类型必须为 DAPI、fluorescence、HE、other 或 unspecified。')
    a = math.radians(angle)
    r = np.array([[math.cos(a), math.sin(a)], [-math.sin(a), math.cos(a)]])
    flip = np.diag([-1.0 if flip_horizontal else 1.0,
                    -1.0 if flip_vertical else 1.0])
    linear = r @ flip
    center = np.array([values['center_x'], values['center_y']])
    matrix = np.eye(3)
    matrix[:2, :2] = linear
    matrix[:2, 2] = np.array([(nx-1)/2, (ny-1)/2]) - linear @ center
    inverse = np.linalg.inv(matrix)
    edges = np.array([[-.5, -.5, 1], [nx-.5, -.5, 1], [nx-.5, ny-.5, 1], [-.5, ny-.5, 1]])
    corners = (edges @ inverse.T)[:, :2]
    outside = bool(np.any(corners[:, 0] < -.5-1e-7) or np.any(corners[:, 0] > source.width-.5+1e-7)
                   or np.any(corners[:, 1] < -.5-1e-7) or np.any(corners[:, 1] > source.height-.5+1e-7))
    full_edges = np.array([[-.5,-.5], [source.width-.5,-.5],
                           [source.width-.5,source.height-.5], [-.5,source.height-.5]]) @ linear.T
    global_offset = -.5 - full_edges.min(axis=0)
    rotated_origin = linear @ center + global_offset - np.array([(nx-1)/2, (ny-1)/2])
    calibration = params.get('calibration') or {'method': 'manual', 'um_per_pixel': u}
    if not math.isclose(float(calibration.get('um_per_pixel', u)), u, rel_tol=1e-7):
        raise ValueError('标定记录与 µm/像素不一致，请重新确认。')
    if calibration.get('method') == 'scale_bar_two_points':
        measured = calibration_from_points(calibration['points_source_px'], calibration['length_um'])
        if not math.isclose(measured['um_per_pixel'], u, rel_tol=1e-7):
            raise ValueError('标尺端点/长度与像素标定不一致。')
    recipe = {
        'schema': 'native_roi_crop_v2', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'source': source.signature, 'calibration': calibration, 'calibration_confirmed': True,
        'image_role': image_role, 'input_storage_format': source.storage_format,
        'um_per_pixel': u, 'requested_size_mm': [values['width_mm'], values['height_mm']],
        'actual_size_mm': [nx*u/1000, ny*u/1000], 'output_size_px': [nx, ny],
        'center_source_px': center.tolist(), 'angle_deg_ccw': angle,
        'angle_input_deg_ccw': values['angle_deg'], 'rotation_pivot_source_px': center.tolist(),
        'flip_horizontal': flip_horizontal, 'flip_vertical': flip_vertical,
        'operation_order': 'source -> horizontal/vertical flip about selected source center -> visual CCW rotation -> crop canvas',
        'source_to_crop': matrix.tolist(), 'crop_to_source': inverse.tolist(),
        'crop_corners_source_pixel_edges': corners.tolist(),
        'expanded_rotated_canvas_size_px': np.ceil(full_edges.max(axis=0)-full_edges.min(axis=0)-1e-8).astype(int).tolist(),
        'crop_origin_expanded_rotated_canvas_px': rotated_origin.tolist(),
        'source_to_expanded_rotated_translation_px': global_offset.tolist(),
        'outside_source': outside, 'allow_black_padding': bool(params.get('allow_black_padding', False)),
        'coordinate_convention': 'x right, y down; zero-based source pixel centers; image edges -0.5..size-0.5; angle positive visual counterclockwise; horizontal flip mirrors left-right and vertical flip mirrors top-bottom within the selected raw region',
        'interpolation': 'bilinear, one inverse warp from original pixels; no additional downsample',
        'preview_only_adjustments': 'browser brightness and zoom are NOT applied to exported TIFF',
    }
    if params.get('imported_recipe_source'):
        recipe['imported_recipe_source'] = params['imported_recipe_source']
    return recipe


def export_crop(source, recipe, output_root, progress=lambda value: None, tile_size=1024):
    if not isinstance(tile_size, int) or not 1 <= tile_size <= 8192:
        raise ValueError('导出块大小须为 1–8192 的整数。')
    source.check_unchanged()
    if recipe['source'] != source.signature:
        raise ValueError('裁切参数与当前源图不匹配，请重新打开/确认。')
    if recipe['outside_source'] and not recipe['allow_black_padding']:
        raise ValueError('6×6 mm 框部分超出原图；请移动/旋转，或明确允许黑色补边。')
    nx, ny = recipe['output_size_px']
    shape = (ny, nx) + source.array.shape[2:]
    size_bytes = int(np.prod(shape)) * source.array.dtype.itemsize
    output_root = Path(output_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(output_root).free < size_bytes + 64*1024**2:
        raise ValueError('导出目录可用空间不足。')
    outdir = Path(tempfile.mkdtemp(prefix=datetime.now().strftime('crop_%Y%m%d_%H%M%S_'), dir=output_root))
    (outdir/'recipe.pending.json').write_text(json.dumps(recipe, ensure_ascii=False, indent=2))
    out = None
    try:
        dtype = source.array.dtype.newbyteorder('=')
        partial = outdir/'crop.partial.tif'
        out = tifffile.memmap(partial, shape=shape, dtype=dtype,
            photometric='rgb' if len(shape) == 3 else 'minisblack',
            extrasamples=source.extrasamples or None,
            bigtiff=size_bytes >= 4*1024**3-32*1024**2,
            resolution=(10000/recipe['um_per_pixel'], 10000/recipe['um_per_pixel']), resolutionunit='CENTIMETER',
            metadata={'axes': 'YXS' if len(shape) == 3 else 'YX', 'um_per_pixel': recipe['um_per_pixel']})
        inv = np.array(recipe['crop_to_source'], dtype=float)
        total = math.ceil(nx/tile_size)*math.ceil(ny/tile_size)
        done = 0
        for y in range(0, ny, tile_size):
            for x in range(0, nx, tile_size):
                tw, th = min(tile_size,nx-x), min(tile_size,ny-y)
                p = np.array([[x,y,1],[x+tw-1,y,1],[x+tw-1,y+th-1,1],[x,y+th-1,1]]) @ inv.T
                sx0,sy0 = np.maximum(0, np.floor(p[:,:2].min(axis=0)).astype(int)-2)
                sx1,sy1 = np.minimum([source.width,source.height], np.ceil(p[:,:2].max(axis=0)).astype(int)+3)
                if sx0 >= sx1 or sy0 >= sy1:
                    out[y:y+th,x:x+tw] = 0
                else:
                    patch = np.ascontiguousarray(source.array[sy0:sy1,sx0:sx1], dtype=dtype)
                    local = inv[:2].copy()
                    local[:,2] += inv[:2,:2] @ [x,y] - [sx0,sy0]
                    out[y:y+th,x:x+tw] = cv2.warpAffine(patch, local, (tw,th),
                        flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP, borderMode=cv2.BORDER_CONSTANT, borderValue=0)
                done += 1
                progress({'completed_tiles': done, 'total_tiles': total, 'percent': round(100*done/total, 1), 'output_dir': str(outdir)})
        out.flush()
        del out
        out = None
        source.check_unchanged()
        os.replace(partial, outdir/'crop.tif')
        with_source = TiffSource(outdir/'crop.tif')
        Image.fromarray(with_source.thumbnail(2000)).save(outdir/'crop_preview.jpg', quality=94)
        del with_source
        preview = source.thumbnail(2000)
        polygon = (np.array(recipe['crop_corners_source_pixel_edges'])+.5) * [preview.shape[1]/source.width,preview.shape[0]/source.height] - .5
        cv2.polylines(preview, [np.round(polygon).astype(np.int32)], True, (0,255,80), 2)
        Image.fromarray(preview).save(outdir/'selection_on_source.jpg', quality=94)
        for name in ('source_to_crop','crop_to_source'):
            with (outdir/f'{name}.csv').open('w', newline='') as handle:
                csv.writer(handle).writerows(recipe[name])
        recipe = dict(recipe, output={'directory': str(outdir), 'tiff': str(outdir/'crop.tif'), 'bytes': (outdir/'crop.tif').stat().st_size})
        (outdir/'recipe.json').write_text(json.dumps(recipe, ensure_ascii=False, indent=2))
        (outdir/'recipe.pending.json').unlink()
        return recipe
    except Exception as exc:
        if out is not None:
            out.flush()
        (outdir/'FAILED.txt').write_text(str(exc))
        raise


def params_from_recipe(recipe):
    if recipe.get('schema') not in ('native_roi_crop_v1', 'native_roi_crop_v2'):
        raise ValueError('不支持的裁切参数格式。')
    return dict(center_x=recipe['center_source_px'][0], center_y=recipe['center_source_px'][1],
                angle_deg=recipe['angle_deg_ccw'], width_mm=recipe['requested_size_mm'][0],
                height_mm=recipe['requested_size_mm'][1], um_per_pixel=recipe['um_per_pixel'],
                calibration=recipe['calibration'], calibration_confirmed=True,
                allow_black_padding=recipe['allow_black_padding'],
                flip_horizontal=bool(recipe.get('flip_horizontal', False)),
                flip_vertical=bool(recipe.get('flip_vertical', False)),
                image_role=recipe.get('image_role', 'unspecified'))
