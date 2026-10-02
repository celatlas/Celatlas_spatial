
import os
import re
import ast
import json
import tarfile
import argparse
import shutil
import gzip
import cv2
import skimage
import tifffile
import numpy as np
import pandas as pd
import SimpleITK as sitk

# Set matplotlib to non-interactive backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from tqdm import tqdm
from datetime import datetime
from joblib import Parallel, delayed
from scipy.io import mmwrite
from scipy import ndimage as ndi
from skimage import measure
from skimage.filters import threshold_multiotsu
from skimage.transform import estimate_transform

from celatlas_spatial.celatlas import ArgFormatter
from celatlas_spatial.tools import utils
from celatlas_spatial.tools import reference
from celatlas_spatial.tools.process import Processor
from celatlas_spatial.tools.image_seg import SwinChipCut
from celatlas_spatial.tools.spatial_preview import apply_tissue_mask, blend_heatmap_with_mask
from celatlas_spatial.tools.spatial_canvas import resolve_tissue_canvas_shape
from celatlas_spatial.tools.spatial_registration_metadata import image_signature, registration_manifest
from celatlas_spatial.tools.step import Step, s_common
from celatlas_spatial.tools.matrix import CountMatrix
from celatlas_spatial.tools.downsample import deterministic_downsample_metrics
from celatlas_spatial.rna.mkref import Mkref_rna
from celatlas_spatial.tools.plotly_plot import Table_plot
from celatlas_spatial.tools.__init__ import BARCODE_FILE_NAME, MATRIX_FILE_NAME, FEATURE_FILE_NAME
from celatlas_spatial.__init__ import __VERSION__, ROOT_PATH, HELP_DICT, genome


class SquareBin:
    def __init__(self, micron_bin, bin, tissue_bbox, bin_folder=None, gene_count_detail=None):
        self.micron_bin = micron_bin
        self.tissue_bbox = tissue_bbox
        self.bin = bin
        self.up_row = 0
        self.up_col = 0
        self.num_bits = 0
        self.bin_folder = bin_folder
        self.gene_count_detail = gene_count_detail

    def upsampling_tissue_matrix(self):
        self.up_row = np.ceil(self.tissue_bbox[3] / self.bin).astype(np.int64)
        self.up_col = np.ceil(self.tissue_bbox[2] / self.bin).astype(np.int64)
        self.num_bits = np.ceil(np.log2(self.up_row * self.up_col)).astype(np.int64)
        if self.num_bits % 2 == 1:
            self.num_bits += 1
            self.num_bits = np.int64(self.num_bits)
        mtx = np.zeros((self.up_row, self.up_col), dtype=object)
        return mtx


class BinSegment(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self._bin = 1
        self._parallel = True
        self._raw_show = True
        self._tmp = False
        self._test = self.debug
        self._extend_list = f'{ROOT_PATH}/data/bclist/bclist_V2.0'

        self.bin_exec = None
        self.pixel_size = self.args.pixel_size
        self.bin = np.ceil(np.array(self._micron_bin) / self.pixel_size)

        self.up_row = 0
        self.up_col = 0
        self.num_bits = 0
        self.base_mapping = np.array(['A', 'C', 'G', 'T'], dtype='U1')

        self.tissue_name = self.args.sample
        self.ref_genome = self.args.genomeDir
        self.input = utils.find_directory_by_name(self.args.input, self.args.sample)

        count_detail_arg = getattr(self.args, 'count_detail', None)
        if count_detail_arg:
            self.count_dir = os.path.dirname(os.path.normpath(count_detail_arg))
        else:
            outdir_parent = os.path.dirname(os.path.normpath(self.outdir))
            if os.path.basename(outdir_parent) == '06.segment':
                outdir_parent = os.path.dirname(outdir_parent)
            self.count_dir = os.path.join(outdir_parent, '05.count')
        self.exp_bs_path = self.outdir  # binSegment output folder path
        self.exp_bin_im = os.path.join(self.exp_bs_path, 'images')
        self.exp_sq_bin = os.path.join(self.exp_bs_path, 'square_bin')
        self.exp_tmp_path = os.path.join(self.exp_bs_path, 'tmp')
        segment_parent = os.path.dirname(os.path.normpath(self.exp_bs_path))
        if os.path.basename(segment_parent) == '06.segment':
            self.runtime_mask_dir = os.path.join(segment_parent, 'mask')
        else:
            self.runtime_mask_dir = os.path.abspath(os.path.join(self.exp_bs_path, '..', 'mask'))

        # ========================================
        # Mode Definition: Three execution modes supported
        # ========================================
        # 1. ssDNA mode: ssDNA image + gene expression data
        # 2. gene_expr mode: Gene expression data only, no images
        # 3. HE mode: H&E stained image + gene expression data

        self.method = self.args.method  # Raw parameter: "gene_expr", "ssDNA"/"image", or "HE"
        self.segment_type = self.args.segment_type

        # Mode: Will be set to 'ssDNA', 'gene_expr', or 'HE' in detect_mode()
        self.mode = None

        self.segment = self.args.segment
        self.count = self.args.count
        self.rectify = self.args.rectify
        self.extend = self.args.extend

        # Gene expression enhancement parameters
        enhance_method_raw = getattr(self.args, 'enhance_method', 'percentile')
        self.enhance_method = None if enhance_method_raw == 'none' else enhance_method_raw

        enhance_params_str = getattr(self.args, 'enhance_params', None)
        if enhance_params_str:
            try:
                self.enhance_params = json.loads(enhance_params_str)
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse enhance_params JSON: {e}")
                self.enhance_params = None
        else:
            self.enhance_params = None

        # Gene expression aggregation resolution (in microns, for generating gem heatmap)
        self.gem_bin_size = getattr(self.args, 'gem_bin_size', 10)

        # Minimum UMI threshold (for bottom noise pre-filtering)
        self.umi_min_threshold = getattr(self.args, 'umi_min_threshold', 'auto')

        # ssDNA tissue-mask controls. A lower threshold scale retains dim
        # peripheral fluorescence; optional expansion compensates for small
        # edge losses after morphology without changing the default behavior.
        self.ssdna_threshold_scale = float(getattr(self.args, 'ssdna_threshold_scale', 1.0))
        self.ssdna_mask_expand_pixels = int(getattr(self.args, 'ssdna_mask_expand_pixels', 0))
        if self.ssdna_threshold_scale <= 0:
            raise ValueError('--ssdna-threshold-scale must be greater than zero')
        if self.ssdna_mask_expand_pixels < 0:
            raise ValueError('--ssdna-mask-expand-pixels must be zero or greater')

        # Registration type (for aligning HE image with gene expression mask)
        self.registration_type = getattr(self.args, 'registration_type', 'similarity')
        self.fluorescence_background = bool(getattr(self.args, 'fluorescence_background', False)) or (
            os.environ.get('CELATLAS_FLUORESCENCE_BACKGROUND', '').lower() in {'1', 'true', 'yes'}
        )
        # Enable SimpleITK fine-tuning optimization
        use_sitk_str = getattr(self.args, 'use_sitk_refinement', 'false')
        self.use_sitk_refinement = (use_sitk_str.lower() == 'true')
        # Enable feature-based fine registration
        use_feature_str = getattr(self.args, 'use_feature_refinement', 'false')
        self.use_feature_refinement = (use_feature_str.lower() == 'true')

        # New parameters for improved HE mode
        # Multi-region threshold for HE segmentation
        self.he_multi_region_threshold = getattr(self.args, 'he_multi_region_threshold', 0.1)
        # Tissue mask generation strategy
        self.tissue_mask_strategy = getattr(self.args, 'tissue_mask_strategy', 'gem_primary')
        self.gene_mask_filter = bool(getattr(self.args, 'gene_mask_filter', False))
        self.tissue_image_mask = None
        self.he_image_mask = None
        self.background_image_mask = None
        self.background_image_mask_before_reg = None
        self.he_image = None
        self.gem_mask = None
        self.gem_manual_mask_path = None
        self.gem = None
        if self.input is None or 'FilterBarcodes' not in self.input:
            raise FileNotFoundError(
                f"FilterBarcodes file for sample {self.tissue_name} was not found under {self.args.input}"
            )
        self.barcodes_detail = pd.read_csv(self.input['FilterBarcodes'], header=None, names=['x', 'y', 'barcode'])
        self.tissue_bbox = self._resolve_tissue_canvas()
        self.tissue_image = np.zeros((self.tissue_bbox[2], self.tissue_bbox[3]), dtype=np.uint8)
        self.chip_shape = np.array((self.tissue_bbox[2], self.tissue_bbox[3]), dtype=np.int32)

        if self.segment:
            self.prompt = ast.literal_eval(self.args.prompt) if self.args.prompt is not None else self.args.prompt
            self.model = self.args.model
            # pre-load tissue image for segmentation
            self.im_path = self.args.tif
            if self.im_path and os.path.exists(self.im_path):
                # Support multiple image formats (TIFF, JPEG, PNG, etc.)
                file_ext = os.path.splitext(self.im_path)[1].lower()
                if file_ext in ['.tif', '.tiff']:
                    # Use tifffile for TIFF format (supports 16-bit and multi-channel)
                    self.tissue_image = tifffile.imread(self.im_path)
                else:
                    # Use OpenCV for other formats (JPEG, PNG, etc.)
                    self.tissue_image = cv2.imread(self.im_path, cv2.IMREAD_UNCHANGED)
                    if self.tissue_image is None:
                        raise FileNotFoundError(f"Failed to read image file: {self.im_path}")
                    # OpenCV loads as BGR, convert to RGB if color image
                    if self.tissue_image.ndim == 3:
                        self.tissue_image = cv2.cvtColor(self.tissue_image, cv2.COLOR_BGR2RGB)

                self.image_d = self.tissue_image.copy()
                if self.tissue_image.dtype != np.uint8:
                    self.tissue_image = utils.convert_8bit(self.tissue_image)
            else:
                self.image_d = np.zeros_like(self.tissue_image)
            self.counts_path = os.path.join(self.count_dir, f'{self.sample}_counts.txt')
            self.roi_rect = np.array([0] * 4)

        if self.count:
            self.count_detail = self.args.count_detail
            self.tb_position = None
            self.in_bc = None
            self.gene_count_detail = None
            self.raw_count_detail = None
            self.bs_mtx = None
            self.features = None

        if self.rectify:
            self.slice_px = (self.bin[-1] * 10).astype(np.int64)  # 1000 micron here

        if self.extend:
            self._bc_list = pd.read_csv(self._extend_list, header=None)
            self._barcode_extend()

    def _is_he_image(self, filename):
        """
        Strictly determine if filename represents an H&E image (must contain _he suffix)

        Parameters:
        -----------
        filename : str
            Filename (without path)

        Returns:
        --------
        bool : True if it's an H&E image

        Matching Rules:
        ---------------
        - Must contain '_he' suffix, followed by extension (.tif, .png, .jpg, .jpeg)
        - Case insensitive

        Examples:
        ---------
        ✅ Matches:
          - sample_he.tif
          - ST110250_A1_he.png
          - tissue_HE.jpg
          - image_he.jpeg

        ❌ Does NOT match:
          - ST110250_A1.tif (no _he suffix)
          - the_best.tif (contains 'he', but not as '_he' suffix)
          - somewhere.jpg (contains 'he' but not as suffix)
          - ssDNA.tif (does not contain at all)
        """
        filename_lower = filename.lower()

        # Strict match: *_he.{tif,tiff,png,jpg,jpeg}
        # Use regex to ensure _he is before extension
        he_pattern = r'_he\.(tif|tiff|png|jpg|jpeg)$'

        return bool(re.search(he_pattern, filename_lower))

    def detect_mode(self):
        """
        Auto-detect binSegment execution mode

        Returns:
        --------
        str : 'ssDNA', 'gene_expr', or 'HE'

        Mode Definitions:
        -----------------
        1. ssDNA mode: ssDNA image + gene expression data
        2. gene_expr mode: Gene expression data only, no images
        3. HE mode: H&E stained image + gene expression data

        Detection Logic:
        ----------------
        - If method == 'gene_expr' → gene_expr mode
        - If method in {'ssDNA', 'image'} and has image:
            - Use strict filename matching to determine if HE image
            - Yes → HE mode
            - No → ssDNA mode (default)
        - Priority: HE > ssDNA > gene_expr

        Prevent False Positives:
        ------------------------
        - Use regex word boundary matching
        - Avoid misclassifying words containing 'he' like 'the', 'other', 'somewhere'
        """

        if self.method == 'gene_expr':
            # Pure expression mode
            mode = 'gene_expr'
            print("[Mode Detection] gene_expr mode: expression data only, no images")

        elif self.method == 'HE':
            # HE mode: Must have HE image
            mode = 'HE'
            print("[Mode Detection] HE mode: gene expression + HE image registration")

            # Verify HE image is provided
            if hasattr(self, 'im_path') and self.im_path and os.path.exists(self.im_path):
                filename = os.path.basename(self.im_path)
                if not self._is_he_image(filename):
                    print(f"[Mode Warning] --method=HE but image filename '{filename}' doesn't match HE pattern (*_he.tif/png/jpg)")
                    print(f"[Mode Warning] Continuing anyway, but please use proper HE image naming convention")
            else:
                print(f"[Mode Warning] --method=HE but no --tif image provided!")
                print(f"[Mode Warning] Will fallback to gene_expr mode")
                mode = 'gene_expr'

        elif self.method in {'image', 'ssDNA'}:
            # Has image, need to distinguish HE vs ssDNA
            if hasattr(self, 'im_path') and self.im_path and os.path.exists(self.im_path):
                filename = os.path.basename(self.im_path)

                # Use strict filename matching
                if self._is_he_image(filename):
                    mode = 'HE'
                    print(f"[Mode Detection] HE mode: detected H&E staining image from filename: {filename}")
                else:
                    mode = 'ssDNA'
                    print(f"[Mode Detection] ssDNA mode: detected ssDNA/other image from filename: {filename}")
            else:
                # No image but method is 'image', fallback to gene_expr mode
                mode = 'gene_expr'
                print("[Mode Detection] gene_expr mode: ssDNA/image method but no image file found, fallback to gene_expr")
        else:
            # Unknown method, default to gene_expr
            mode = 'gene_expr'
            print(f"[Mode Detection] gene_expr mode: unknown method '{self.method}', fallback to gene_expr")

        self.mode = mode
        print(f"[Mode Detection] Final mode: {mode}")
        print(f"[Mode Validation] Mode is locked and cannot be changed during this run")
        print("=" * 60)

        return mode

    def _set_background_image_mask(self, mask):
        """Store a mode-agnostic background mask."""
        self.background_image_mask = mask
        self.he_image_mask = mask

    def _set_background_image_mask_before_reg(self, mask):
        self.background_image_mask_before_reg = mask
        self.he_image_mask_before_reg = mask

    def _get_background_image_mask(self):
        return self.background_image_mask if self.background_image_mask is not None else self.he_image_mask

    def _get_background_image_mask_before_reg(self):
        if self.background_image_mask_before_reg is not None:
            return self.background_image_mask_before_reg
        return getattr(self, 'he_image_mask_before_reg', None)

    def _read_tissue_bbox_metadata(self):
        bbox_file = self.input.get('TissueBbox') if self.input else None
        if not bbox_file or not os.path.exists(bbox_file):
            return None
        try:
            bbox_raw = str(pd.read_csv(bbox_file)['bbox'][0])
            bbox = np.array(bbox_raw.split('\t')).astype(np.int32)
            if len(bbox) < 4 or bbox[2] <= 0 or bbox[3] <= 0:
                print(f"[Canvas] Ignoring invalid tissue_bbox: {bbox_raw}")
                return None
            return bbox[:4]
        except Exception as e:
            print(f"[Canvas] Warning: failed to read tissue_bbox {bbox_file}: {e}")
            return None

    def _read_barcode_to_pos_shape(self):
        h5_file = self.input.get('BarcodeToPos') if self.input else None
        if not h5_file or not os.path.exists(h5_file):
            return None
        try:
            import h5py
            with h5py.File(h5_file, 'r') as handle:
                shapes = []

                def collect_shape(_, obj):
                    if hasattr(obj, 'shape') and len(obj.shape) >= 2:
                        shapes.append(tuple(int(x) for x in obj.shape[:2]))

                handle.visititems(collect_shape)
            if not shapes:
                print(f"[Canvas] Warning: no 2D dataset found in barcodeToPos h5: {h5_file}")
                return None
            # Prefer the largest matrix; barcodeToPos h5 normally stores one bpMatrix dataset.
            height, width = max(shapes, key=lambda shape: shape[0] * shape[1])
            if height > 0 and width > 0:
                return int(height), int(width)
        except Exception as e:
            print(f"[Canvas] Warning: failed to read barcodeToPos h5 {h5_file}: {e}")
        return None

    def _barcode_coordinate_shape(self):
        try:
            coords = self.barcodes_detail[['x', 'y']].apply(pd.to_numeric, errors='coerce').dropna()
            if coords.empty:
                return None
            max_x = int(np.ceil(coords['x'].max()))
            max_y = int(np.ceil(coords['y'].max()))
            if max_x < 0 or max_y < 0:
                return None
            return max_y + 1, max_x + 1
        except Exception as e:
            print(f"[Canvas] Warning: failed to infer canvas from FilterBarcodes: {e}")
            return None

    def _resolve_tissue_canvas(self):
        """Resolve full chip canvas from barcodeToPos and FilterBarcodes coordinates."""
        tissue_bbox_metadata = self._read_tissue_bbox_metadata()
        h5_shape = self._read_barcode_to_pos_shape()
        coord_shape = self._barcode_coordinate_shape()

        height, width, source, messages = resolve_tissue_canvas_shape(
            tissue_bbox_metadata=tissue_bbox_metadata,
            h5_shape=h5_shape,
            coord_shape=coord_shape,
        )
        for message in messages:
            print(f"[Canvas] {message}")

        tissue_bbox = np.array([0, 0, int(height), int(width)], dtype=np.int32)
        print(
            f"[Canvas] Resolved tissue canvas from {source}: "
            f"height={tissue_bbox[2]}, width={tissue_bbox[3]}"
        )
        return tissue_bbox

    @utils.add_log
    def prepare(self):
        # Detect and set mode
        self.detect_mode()

        # create output folder
        os.makedirs(self.exp_bs_path, exist_ok=True)
        os.makedirs(self.exp_bin_im, exist_ok=True)
        os.makedirs(self.exp_sq_bin, exist_ok=True)
        os.makedirs(self.exp_tmp_path, exist_ok=True)

        # prepare some required files
        if self.segment:
            if self.im_path and os.path.exists(self.im_path):
                os.system(f'cp -r {self.im_path} {os.path.join(self.exp_bin_im, f"{self.tissue_name}.tif")}')

        if self.count:
            bins = self.bin
            bins = np.append(bins, 1.) if self._raw_show else bins
            for i in range(len(bins)):
                if bins[i] != 1:
                    bin_path = os.path.join(self.exp_sq_bin, f'{self.tissue_name}_bin{self._micron_bin[i]}')
                else:
                    bin_path = os.path.join(self.exp_sq_bin, f'{self.tissue_name}_Raw')
                os.makedirs(bin_path, exist_ok=True)
                os.makedirs(os.path.join(bin_path, 'filtered_feature_bc_matrix'), exist_ok=True)
                os.makedirs(os.path.join(bin_path, 'spatial'), exist_ok=True)

                self.scalefactors_json(bin_path, bins[i])
    @staticmethod
    def read_reads_stat(stat_path):
        reads_stat = {}
        pattern = re.compile(r'(\d{1,3},)+\d{3}')

        with open(stat_path, 'r') as f:
            for line in f:
                key, value = line.strip().split(':')
                match = pattern.search(value)
                if match:
                    value = int(match.group().replace(',', ''))
                else:
                    value = value.strip()
                reads_stat[key.strip()] = value
        return reads_stat

    def scalefactors_json(self, output, bin_size):
        # 10X/Seurat expects tissue position columns 5/6 to be in full-resolution
        # image pixels. Scale factors therefore map full-resolution coordinates
        # to the stored hires/lowres PNGs and must not depend on square-bin size.
        fullres_max_dim = float(max(self.chip_shape))
        if fullres_max_dim <= 0:
            fullres_max_dim = 1.0
        hires_scale = 2000.0 / fullres_max_dim
        lowres_scale = 600.0 / fullres_max_dim

        factors = {
            "tissue_hires_scalef": hires_scale,
            "tissue_lowres_scalef": lowres_scale,
            "fiducial_diameter_fullres": bin_size,
            "spot_diameter_fullres": bin_size
        }
        json_data = json.dumps(factors, indent=4)
        with open(os.path.join(output, 'scalefactors_json.json'), 'w') as json_file:
            json_file.write(json_data)

    @staticmethod
    def _resize_longest_side(image, max_side, interpolation=cv2.INTER_AREA):
        """Resize image with aspect ratio preserved and longest side fixed."""
        height, width = image.shape[:2]
        max_dim = max(height, width)
        if max_dim <= 0:
            return image
        scale = float(max_side) / float(max_dim)
        target_width = max(1, int(round(width * scale)))
        target_height = max(1, int(round(height * scale)))
        return cv2.resize(image, (target_width, target_height), interpolation=interpolation)

    def _write_spatial_preview_images(self, image):
        lowres = self._resize_longest_side(image, 600, interpolation=cv2.INTER_AREA)
        hires = self._resize_longest_side(image, 2000, interpolation=cv2.INTER_LANCZOS4)
        if image.ndim == 3:
            lowres = cv2.cvtColor(lowres, cv2.COLOR_RGB2BGR)
            hires = cv2.cvtColor(hires, cv2.COLOR_RGB2BGR)
        cv2.imwrite(os.path.join(self.exp_bin_im, 'tissue_lowres_image.png'), lowres)
        cv2.imwrite(os.path.join(self.exp_bin_im, 'tissue_hires_image.png'), hires)
        print(
            "[Spatial Image] Saved aspect-preserving previews: "
            f"lowres={lowres.shape[1]}x{lowres.shape[0]}, "
            f"hires={hires.shape[1]}x{hires.shape[0]}"
        )

    def _refresh_ssdna_spatial_previews(self):
        """Mask report previews after the final ssDNA tissue region is known."""
        if self.mode != 'ssDNA' or self.tissue_image_mask is None:
            return

        masked_image = apply_tissue_mask(self.tissue_image, self.tissue_image_mask)
        self._write_spatial_preview_images(masked_image)
        print(
            "[ssDNA mode] Refreshed spatial preview images using the final tissue mask; "
            "the full registered TIFF remains unchanged"
        )

    @utils.add_log
    def slice_registration(self, model, mode="ssDNA", enhance_method='percentile', enhance_params=None):
        """
        Main slice registration workflow

        Parameters:
        -----------
        model : str
            Model path (for image method)
        mode : str
            Execution mode: 'ssDNA', 'gene_expr', or 'HE'
        enhance_method : str
            Gene expression enhancement method (for gene_expr/HE methods only)
        enhance_params : dict
            Enhancement method parameter dictionary (for gene_expr/HE methods only)
        """
        # Use configurable gem_bin_size (in microns) instead of hardcoded 10
        resolution = np.ceil(self.gem_bin_size / self.pixel_size).astype(int)
        print(f"Using gem_bin_size={self.gem_bin_size} microns (resolution={resolution} pixels) for gene expression aggregation")

        if mode == "gene_expr":
            # Pure gene expression mode: No HE registration
            self.gem_registration(resolution=resolution, enhance_method=enhance_method, enhance_params=enhance_params)
            # self.save_gene_expr_enhancement_comparison()  # Debug visualization, not needed in production

            # Use gem as tissue_image
            self.tissue_image = self.gem
            self._set_background_image_mask(self.gem_mask)
            print("[gene_expr mode] No HE registration, using gem as tissue_image")

        elif mode == "HE":
            # HE mode: Gene expression + HE registration
            self.gem_registration(resolution=resolution, enhance_method=enhance_method, enhance_params=enhance_params)

            # Must have HE image
            if self.args.tif and os.path.exists(self.args.tif):
                filename = os.path.basename(self.args.tif)
                print(f"[HE mode] Performing HE registration with: {filename}")
                self.he_registration(self.args.tif)
                self.save_background_registration_comparison()  # Generates overlays folder (debug plot disabled internally)
            else:
                print(f"[HE mode] ERROR: --tif parameter required for HE mode!")
                print(f"[HE mode] Fallback to gene_expr mode")
                self.tissue_image = self.gem
                self._set_background_image_mask(self.gem_mask)

            self.save_gene_expr_enhancement_comparison()  # Debug plot disabled internally

        elif mode == "ssDNA":
            # ssDNA image segmentation mode
            if self.gene_mask_filter:
                print("[ssDNA mode] gene-mask filter enabled; generating/loading GEM mask before ssDNA segmentation")
                self.gem_registration(resolution=resolution, enhance_method=enhance_method, enhance_params=enhance_params)
            self.chipset_registration(model=model)

        self.save_register_image()

    @utils.add_log
    def chipset_registration(self, model, image_size=224):
        # image-segmentation base on transformer-unet framework
        self.chipset_segmentation(model, image_size=image_size)

        # load chip image
        chip_masks = np.zeros_like(self.tissue_image)
        positions = self.barcodes_detail.to_numpy()
        chip_contours = np.array([
            np.min(positions[:, 0]), np.max(positions[:, 0]), np.min(positions[:, 1]), np.max(positions[:, 1])
        ])
        chip_masks[chip_contours[2]: chip_contours[3], chip_contours[0]: chip_contours[1]] = 255

        row_coords, col_coords = np.where(self.pred == 1)
        min_col, min_row, max_col, max_row = np.min(col_coords), np.min(row_coords), np.max(col_coords), np.max(row_coords)
        bbox_tissue = [min_col, min_row, max_col, max_row]
        bbox_tissue = np.array(
            [[bbox_tissue[0], bbox_tissue[1]],
             [bbox_tissue[2], bbox_tissue[1]],
             [bbox_tissue[2], bbox_tissue[3]],
             [bbox_tissue[0], bbox_tissue[3]]]
        )
        bbox_chip = np.array(
            [[chip_contours[0], chip_contours[2]],
             [chip_contours[1], chip_contours[2]],
             [chip_contours[1], chip_contours[3]],
             [chip_contours[0], chip_contours[3]]]
        )

        # option: ‘similarity’ or ‘affine’ or ‘projective’
        tform = estimate_transform('affine', bbox_tissue, bbox_chip)
        H = tform.params
        tissue_image = self.tissue_image
        tissue_image = cv2.warpPerspective(tissue_image, H, (self.tissue_bbox[3], self.tissue_bbox[2]))
        self.tissue_image = tissue_image  # update tissue image after registration

    @utils.add_log
    def chipset_segmentation(self, model_path, image_size):
        """
        Tissue slide segmentation

        Uses Swin-UNet deep learning model for tissue vs background segmentation.
        """
        im_chip_cut = self.tissue_image.copy()

        # Use Swin-UNet strategy
        chip_cut = SwinChipCut(image_size=image_size, model_path=model_path)
        self.pred = chip_cut.f_predict(self.tissue_image)

        im_chip_cut[self.pred == 0] = 0
        self.tissue_image = im_chip_cut

    def umi_bin_counts(self, counts_path, resolution):
        counts = pd.read_csv(counts_path, sep='\t', dtype={'Barcode': str,
                                                           'readcount': int,
                                                           'UMI2': int,
                                                           'UMI': int,
                                                           'geneID': int,
                                                           'mark': str})
        barcodes = self.barcodes_detail.copy()
        barcodes.columns = ['pxl_col_in_fullres', 'pxl_row_in_fullres', 'barcode']
        counts = counts.merge(barcodes, left_on='Barcode', right_on='barcode', how='inner')
        counts = counts.drop(columns=['barcode'])

        row_bins = np.arange(0, self.tissue_bbox[2] + resolution, resolution)
        col_bins = np.arange(0, self.tissue_bbox[3] + resolution, resolution)
        counts['row_bin'] = pd.cut(counts['pxl_row_in_fullres'], row_bins, labels=False, include_lowest=True, right=False)
        counts['col_bin'] = pd.cut(counts['pxl_col_in_fullres'], col_bins, labels=False, include_lowest=True, right=False)

        # get gem image
        num_bits = int(2 * np.ceil(np.ceil(np.log2(len(row_bins) * len(col_bins))) / 2))
        uuid = counts['row_bin'].values * len(col_bins) + counts['col_bin'].values
        uuid_binary = ((uuid[:, np.newaxis] & (1 << np.arange(num_bits))) > 0).astype(int)
        unicode = self.base_mapping[uuid_binary[:, ::-1].reshape(-1, 2).dot(np.array([2, 1]))]
        unicode = unicode.reshape(-1, num_bits // 2).view('U' + str(num_bits // 2)).ravel()
        counts['unicode'] = unicode

        counts = counts[['unicode', 'row_bin', 'col_bin', 'UMI']]
        counts = counts.rename(columns={'unicode': 'Barcode'})
        return counts, len(row_bins) - 1, len(col_bins) - 1

    def enhance_gene_expression_contrast(self, gem_data, method='percentile', **kwargs):
        """
        Enhance dynamic range of gene expression data to improve binary segmentation

        Parameters:
        -----------
        gem_data : numpy.ndarray
            Gene expression heatmap data
        method : str
            Enhancement method: 'percentile', 'clahe', 'gamma', 'bilateral_percentile'
        **kwargs : dict
            Method-specific parameters:
            - percentile: p_low=2, p_high=98
            - clahe: clip_limit=2.0, tile_grid_size=(8,8)
            - gamma: gamma=0.5
            - bilateral_percentile: d=9, sigma_color=75, sigma_space=75, p_low=2, p_high=98

        Returns:
        --------
        enhanced : numpy.ndarray (uint8)
            Enhanced gene expression data
        """
        gem_8bit = utils.convert_8bit(gem_data.copy())

        if method == 'percentile':
            # Percentile stretching: expand histogram, make high values higher and low values lower
            p_low = kwargs.get('p_low', 2)
            p_high = kwargs.get('p_high', 98)

            valid_data = gem_8bit[gem_8bit > 0]
            if len(valid_data) > 0:
                v_min, v_max = np.percentile(valid_data, (p_low, p_high))
                enhanced = np.clip((gem_8bit - v_min) / (v_max - v_min + 1e-8) * 255, 0, 255).astype(np.uint8)
                print(f"Percentile stretching: [{v_min:.2f}, {v_max:.2f}] -> [0, 255]")
            else:
                enhanced = gem_8bit

        elif method == 'clahe':
            # CLAHE: Contrast Limited Adaptive Histogram Equalization
            clip_limit = kwargs.get('clip_limit', 2.0)
            tile_grid_size = kwargs.get('tile_grid_size', (8, 8))

            clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=tile_grid_size)
            enhanced = clahe.apply(gem_8bit)
            print(f"CLAHE applied: clipLimit={clip_limit}, tileGridSize={tile_grid_size}")

        elif method == 'gamma':
            # Gamma correction: power transformation
            # gamma < 1: brighten low values
            # gamma > 1: darken low values, highlight high values
            gamma = kwargs.get('gamma', 0.5)

            normalized = gem_8bit / 255.0
            enhanced = np.power(normalized, gamma) * 255
            enhanced = np.clip(enhanced, 0, 255).astype(np.uint8)
            print(f"Gamma correction applied: gamma={gamma}")

        elif method == 'bilateral_percentile':
            # Combined strategy: bilateral filter for edge-preserving smoothing, then percentile stretch
            d = kwargs.get('d', 9)
            sigma_color = kwargs.get('sigma_color', 75)
            sigma_space = kwargs.get('sigma_space', 75)
            p_low = kwargs.get('p_low', 2)
            p_high = kwargs.get('p_high', 98)

            # Bilateral filtering
            smoothed = cv2.bilateralFilter(gem_8bit, d, sigma_color, sigma_space)

            # Percentile stretching
            valid_data = smoothed[smoothed > 0]
            if len(valid_data) > 0:
                v_min, v_max = np.percentile(valid_data, (p_low, p_high))
                enhanced = np.clip((smoothed - v_min) / (v_max - v_min + 1e-8) * 255, 0, 255).astype(np.uint8)
                print(f"Bilateral + Percentile: smoothed then stretched [{v_min:.2f}, {v_max:.2f}] -> [0, 255]")
            else:
                enhanced = smoothed
        else:
            print(f"Unknown enhancement method '{method}', using original data")
            enhanced = gem_8bit

        # Bottom noise suppression: set values below threshold to 0, enhance black/white contrast
        suppress_noise = kwargs.get('suppress_noise', True)
        noise_threshold = kwargs.get('noise_threshold', 'auto')  # 'auto' or specific value

        if suppress_noise:
            if noise_threshold == 'auto':
                # Auto-calculate noise threshold: use Otsu or percentile method
                valid_data = enhanced[enhanced > 0]
                if len(valid_data) > 0:
                    # Use low percentile as noise threshold
                    noise_thresh = np.percentile(valid_data, 10)  # 10th percentile as noise threshold
                else:
                    noise_thresh = 0
            else:
                noise_thresh = noise_threshold

            # Suppress values below threshold to 0
            enhanced[enhanced < noise_thresh] = 0
            print(f"Bottom noise suppression: threshold={noise_thresh:.2f}, pixels suppressed={np.sum(enhanced == 0)}")

        return enhanced

    @utils.add_log
    def gem_registration(self, resolution, enhance_method='percentile', enhance_params=None):
        """
        Gene expression data registration and binary mask generation

        Parameters:
        -----------
        resolution : int
            Resolution in pixels
        enhance_method : str
            Dynamic range enhancement method: 'percentile', 'clahe', 'gamma', 'bilateral_percentile', None
        enhance_params : dict
            Parameter dictionary for enhancement method
        """
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
        umi_counts, rows, cols = self.umi_bin_counts(self.counts_path, resolution)
        umi_counts = umi_counts.groupby(['Barcode', 'row_bin', 'col_bin']).agg({'UMI': 'sum'}).reset_index()
        self.gem = np.zeros((rows, cols), dtype=np.int64)
        self.gem[umi_counts['row_bin'].values, umi_counts['col_bin'].values] = umi_counts['UMI'].values

        if np.any(self.gem < 0):
            raise ValueError("GEM data contains negative values, which are not supported by threshold_multiotsu")

        # Bottom noise pre-filtering: remove low-UMI bins before Otsu thresholding to improve edge detection
        umi_min_threshold = getattr(self, 'umi_min_threshold', 'auto')
        if umi_min_threshold != 'none':
            if umi_min_threshold == 'auto':
                # Auto-calculate: use 5th percentile as noise threshold
                valid_gem = self.gem[self.gem > 0]
                if len(valid_gem) > 0:
                    umi_thresh = np.percentile(valid_gem, 5)
                    print(f"Auto UMI threshold (5th percentile): {umi_thresh:.1f}")
                else:
                    umi_thresh = 0
            else:
                umi_thresh = float(umi_min_threshold)
                print(f"Manual UMI threshold: {umi_thresh:.1f}")

            # Filter low-UMI bins
            bins_before = np.sum(self.gem > 0)
            self.gem[self.gem < umi_thresh] = 0
            bins_after = np.sum(self.gem > 0)
            print(f"UMI filtering: removed {bins_before - bins_after} low-expression bins ({bins_after}/{bins_before} retained)")

        # Dynamic range enhancement: optimize binary segmentation quality
        if enhance_method is not None:
            enhance_params = enhance_params or {}
            gem_enhanced = self.enhance_gene_expression_contrast(self.gem, method=enhance_method, **enhance_params)
            print(f"Gene expression contrast enhanced using '{enhance_method}' method")
        else:
            gem_enhanced = utils.convert_8bit(self.gem)
            print("No contrast enhancement applied")

        # Save original and enhanced data for comparison visualization
        gem_original_temp = utils.convert_8bit(self.gem.copy())

        # Perform threshold segmentation using enhanced data
        thresholds = threshold_multiotsu(gem_enhanced, classes=3)
        self.gem_mask = gem_enhanced >= thresholds[0]
        self.gem_mask = skimage.morphology.remove_small_objects(self.gem_mask, min_size=2).astype(np.uint8) * 255
        self.gem_mask = cv2.morphologyEx(self.gem_mask, cv2.MORPH_CLOSE, kernel, iterations=3)
        self.gem_mask = skimage.morphology.remove_small_objects(self.gem_mask > 0, min_size=10).astype(np.uint8) * 255
        self.gem_mask = cv2.dilate(self.gem_mask, kernel, iterations=6)
        self.gem_mask = cv2.morphologyEx(self.gem_mask, cv2.MORPH_CLOSE, kernel, iterations=5)
        self.gem_mask = utils.hole_fill(self.gem_mask)
        self.gem_mask = cv2.erode(self.gem_mask, kernel, iterations=3)
        self.gem_mask = cv2.resize(self.gem_mask, (self.tissue_bbox[3], self.tissue_bbox[2]), interpolation=cv2.INTER_AREA)

        # Check for manually edited GEM mask at 06.segment/mask/manual_mask.png.
        manual_mask_path = self._find_manual_mask('gem')
        if manual_mask_path:
            self.gem_manual_mask_path = manual_mask_path
            print(f"[Manual Mask] Found manual mask at: {manual_mask_path}")
            manual_mask = cv2.imread(manual_mask_path, cv2.IMREAD_GRAYSCALE)
            if manual_mask is not None:
                manual_mask_resized = self._resize_manual_mask(manual_mask, manual_mask_path)
                # Apply binary threshold to ensure 0/255 values
                _, manual_mask_binary = cv2.threshold(manual_mask_resized, 127, 255, cv2.THRESH_BINARY)
                self.gem_mask = manual_mask_binary
                print(f"[Manual Mask] Loaded and applied manual mask (size: {manual_mask_binary.shape})")
            else:
                print(f"[Manual Mask] Warning: Failed to load {manual_mask_path}, using auto-generated mask")
        else:
            self.gem_manual_mask_path = None
            print(f"[Manual Mask] No manual mask found, using auto-generated mask")
            print(f"[Manual Mask] To use manual mask: save edited mask as {self._manual_mask_hint('gem')}")

        # Save original and enhanced GEM for later visualization (both resized to tissue_bbox to ensure consistent dimensions)
        # Use INTER_NEAREST to preserve sharp dot patterns, avoid blurring
        self.gem_original = cv2.resize(gem_original_temp, (self.tissue_bbox[3], self.tissue_bbox[2]), interpolation=cv2.INTER_NEAREST)
        self.gem_enhanced = cv2.resize(gem_enhanced.copy(), (self.tissue_bbox[3], self.tissue_bbox[2]), interpolation=cv2.INTER_NEAREST)

        self.gem = utils.convert_8bit(self.gem)
        self.gem = cv2.resize(self.gem, (self.tissue_bbox[3], self.tissue_bbox[2]), interpolation=cv2.INTER_NEAREST)

    @utils.add_log
    def he_registration(self, tif_path):
        if tif_path and os.path.exists(tif_path):
            self.get_tissue_image(tif_path, resize=True)
        else:
            self.tissue_image = self.gem
            self._set_background_image_mask(self.gem_mask)

        if self.rectify:
            processor = Processor(output=self.exp_bin_im, image=self.he_image)
            # roi_mask = processor.detect_border(border_value=179, tolerance=5, debug=False)
            roi_mask, _ = processor.detect_polygons(border_value=179, tolerance=5, min_area=100)
            degree, slope, center = processor.find_top_line(roi_mask)

            self.tissue_image, new_center = processor.get_rotate_image(self.tissue_image, degree, center)
            he_mask, new_center = processor.get_rotate_image(roi_mask, degree, center)
            self._set_background_image_mask(he_mask)

            points = self.barcodes_detail[["x", "y"]].values
            rotated_points = processor.rotate_points(points, -degree, new_center).astype(np.int32)
            min_x, min_y = rotated_points.min(axis=0)
            # max_x, max_y = rotated_points.max(axis=0)
            dx = np.abs(min_x) if min_x < 0 else 0
            dy = np.abs(min_y) if min_y < 0 else 0
            offset = np.array([dx, dy])
            rotated_points += offset
            self.barcodes_detail[["x", "y"]] = rotated_points

            motion_m = np.float32([[1, 0, dx], [0, 1, dy]])
            self.tissue_image = cv2.warpAffine(self.tissue_image, motion_m, (self.tissue_image.shape[1] + dx, self.tissue_image.shape[0] + dy))
            he_mask = cv2.warpAffine(self.he_image_mask, motion_m, (self.he_image_mask.shape[1] + dx, self.he_image_mask.shape[0] + dy))
            self._set_background_image_mask(he_mask)
            self.roi_rect = processor.get_nonzero_bbox(self.he_image_mask)

    @staticmethod
    def register_mask_images_sitk(image_fixed, image_moving, registration_type='affine', verbose=True, use_contour_init=True, initial_transform_matrix=None):
        """
        Optimized image registration using SimpleITK

        Parameters:
        -----------
        image_fixed : ndarray
            Reference mask image (target coordinate system, usually gem_mask)
        image_moving : ndarray
            Moving mask image (source coordinate system, usually he_mask)
        registration_type : str
            Registration type: 'rigid' (rigid), 'similarity' (similarity), 'affine' (affine)
        verbose : bool
            Whether to print detailed information
        use_contour_init : bool
            Whether to use contour registration result as initial transform (recommended)
        initial_transform_matrix : ndarray (2x3), optional
            Externally provided initial transform matrix (takes priority over use_contour_init)

        Returns:
        --------
        transM : ndarray (2x3)
            Affine transformation matrix in OpenCV format
        """
        # Ensure consistent image dimensions
        if image_fixed.shape != image_moving.shape:
            image_moving = cv2.resize(image_moving, (image_fixed.shape[1], image_fixed.shape[0]),
                                     interpolation=cv2.INTER_NEAREST)

        # Convert to SimpleITK images
        fixed_image = sitk.GetImageFromArray(image_fixed.astype(np.float32))
        moving_image = sitk.GetImageFromArray(image_moving.astype(np.float32))

        # Create registration method
        registration_method = sitk.ImageRegistrationMethod()

        # Use multi-resolution pyramid strategy to improve accuracy and speed
        registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
        registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
        registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

        # Choose similarity metric - use Mattes mutual information, better for mask images
        registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
        registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
        registration_method.SetMetricSamplingPercentage(0.2)

        # Set interpolation method
        registration_method.SetInterpolator(sitk.sitkLinear)

        # Set optimizer - use gradient descent
        registration_method.SetOptimizerAsGradientDescent(
            learningRate=1.0,
            numberOfIterations=500,  # Increase iterations to improve accuracy
            convergenceMinimumValue=1e-6,
            convergenceWindowSize=10
        )
        registration_method.SetOptimizerScalesFromPhysicalShift()

        # Choose transformation based on registration type
        # Prioritize externally provided initial transform matrix
        if initial_transform_matrix is not None:
            if verbose:
                print("[SimpleITK Registration] Using provided initial transform matrix...")
            # Convert OpenCV 2x3 matrix to SimpleITK AffineTransform
            initial_transform = sitk.AffineTransform(2)
            # OpenCV format: [a b tx]
            #                [c d ty]
            # SimpleITK needs: matrix=[a,b,c,d], translation=[tx,ty]
            matrix_params = [initial_transform_matrix[0,0], initial_transform_matrix[0,1],
                           initial_transform_matrix[1,0], initial_transform_matrix[1,1]]
            translation_params = [initial_transform_matrix[0,2], initial_transform_matrix[1,2]]
            initial_transform.SetMatrix(matrix_params)
            initial_transform.SetTranslation(translation_params)

            # Set image center as rotation center
            center = [image_fixed.shape[1] / 2.0, image_fixed.shape[0] / 2.0]
            initial_transform.SetCenter(center)

        elif use_contour_init:
            # Use contour registration result as initial transform
            if verbose:
                print("[SimpleITK Registration] Using contour-based initialization...")

            contours_fixed, _ = cv2.findContours(image_fixed.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            contours_moving, _ = cv2.findContours(image_moving.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

            # Multi-region safe: use image moments for global centroid, merged points for angle
            M_fixed = cv2.moments(image_fixed.astype(np.uint8))
            M_moving = cv2.moments(image_moving.astype(np.uint8))
            center_fixed = (int(M_fixed['m10'] / M_fixed['m00']), int(M_fixed['m01'] / M_fixed['m00']))
            center_moving = (int(M_moving['m10'] / M_moving['m00']), int(M_moving['m01'] / M_moving['m00']))

            all_pts_fixed = np.vstack(contours_fixed)
            all_pts_moving = np.vstack(contours_moving)
            angle_fixed = cv2.fitEllipseAMS(all_pts_fixed)[2] if len(all_pts_fixed) >= 5 else 0.0
            angle_moving = cv2.fitEllipseAMS(all_pts_moving)[2] if len(all_pts_moving) >= 5 else 0.0

            rotation_angle = -(angle_fixed - angle_moving)
            area_fixed = float(np.sum(image_fixed > 0))
            area_moving = float(np.sum(image_moving > 0))
            scale = np.sqrt(area_fixed / area_moving)

            if verbose:
                print(f"[Contour Init] Rotation: {rotation_angle:.2f}°, Scale: {scale:.3f}")
                print(f"[Contour Init] Center fixed: {center_fixed} ({len(contours_fixed)} regions), Center moving: {center_moving} ({len(contours_moving)} regions)")

            # Create initial transform - use SimpleITK geometric center initialization, then fine-tune with contour parameters
            if registration_type == 'rigid':
                initial_transform = sitk.CenteredTransformInitializer(
                    fixed_image,
                    moving_image,
                    sitk.Euler2DTransform(),
                    sitk.CenteredTransformInitializerFilter.GEOMETRY
                )
                # Set angle calculated from contours
                initial_transform.SetAngle(np.radians(rotation_angle))
            elif registration_type == 'similarity':
                initial_transform = sitk.CenteredTransformInitializer(
                    fixed_image,
                    moving_image,
                    sitk.Similarity2DTransform(),
                    sitk.CenteredTransformInitializerFilter.GEOMETRY
                )
                # Set scale and angle calculated from contours
                initial_transform.SetScale(scale)
                initial_transform.SetAngle(np.radians(rotation_angle))
            else:  # affine
                initial_transform = sitk.CenteredTransformInitializer(
                    fixed_image,
                    moving_image,
                    sitk.AffineTransform(2),
                    sitk.CenteredTransformInitializerFilter.GEOMETRY
                )
                # For affine transform, set initial matrix
                angle_rad = np.radians(rotation_angle)
                cos_a = np.cos(angle_rad)
                sin_a = np.sin(angle_rad)
                initial_transform.SetMatrix([scale * cos_a, -scale * sin_a,
                                            scale * sin_a, scale * cos_a])
        else:
            # Use geometric center initialization
            if registration_type == 'rigid':
                # Rigid transform: only rotation and translation allowed
                initial_transform = sitk.CenteredTransformInitializer(
                    fixed_image,
                    moving_image,
                    sitk.Euler2DTransform(),
                    sitk.CenteredTransformInitializerFilter.GEOMETRY
                )
            elif registration_type == 'similarity':
                # Similarity transform: rotation, translation, uniform scaling
                initial_transform = sitk.CenteredTransformInitializer(
                    fixed_image,
                    moving_image,
                    sitk.Similarity2DTransform(),
                    sitk.CenteredTransformInitializerFilter.GEOMETRY
                )
            else:  # affine (default)
                # Affine transform: rotation, translation, scaling, shear
                initial_transform = sitk.CenteredTransformInitializer(
                    fixed_image,
                    moving_image,
                    sitk.AffineTransform(fixed_image.GetDimension()),
                    sitk.CenteredTransformInitializerFilter.GEOMETRY
                )

        registration_method.SetInitialTransform(initial_transform, inPlace=False)

        # Execute registration
        if verbose:
            print(f"[SimpleITK Registration] Starting {registration_type} registration...")

        final_transform = registration_method.Execute(fixed_image, moving_image)

        # Output registration quality information
        if verbose:
            print(f"[SimpleITK Registration] Optimizer stop: {registration_method.GetOptimizerStopConditionDescription()}")
            print(f"[SimpleITK Registration] Final metric value: {registration_method.GetMetricValue():.6f}")
            print(f"[SimpleITK Registration] Iterations: {registration_method.GetOptimizerIteration()}")

        # Handle CompositeTransform - extract actual transform
        if isinstance(final_transform, sitk.CompositeTransform):
            # CompositeTransform contains initial transform, we need to extract the optimized transform
            # Get the last transform (optimized result)
            actual_transform = final_transform.GetBackTransform()
        else:
            actual_transform = final_transform

        # Convert to OpenCV format 2x3 transformation matrix
        if isinstance(actual_transform, sitk.Euler2DTransform):
            # Rigid transform
            angle = actual_transform.GetAngle()
            translation = actual_transform.GetTranslation()
            center = actual_transform.GetCenter()

            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            tx = translation[0] + center[0] - cos_a * center[0] + sin_a * center[1]
            ty = translation[1] + center[1] - sin_a * center[0] - cos_a * center[1]

            transM = np.array([
                [cos_a, -sin_a, tx],
                [sin_a, cos_a, ty]
            ], dtype=np.float32)

            if verbose:
                print(f"[Registration] Rigid - Rotation: {np.degrees(angle):.2f}°, Translation: ({tx:.1f}, {ty:.1f})")

        elif isinstance(actual_transform, sitk.Similarity2DTransform):
            # Similarity transform
            scale = actual_transform.GetScale()
            angle = actual_transform.GetAngle()
            translation = actual_transform.GetTranslation()
            center = actual_transform.GetCenter()

            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            tx = translation[0] + center[0] - scale * cos_a * center[0] + scale * sin_a * center[1]
            ty = translation[1] + center[1] - scale * sin_a * center[0] - scale * cos_a * center[1]

            transM = np.array([
                [scale * cos_a, -scale * sin_a, tx],
                [scale * sin_a, scale * cos_a, ty]
            ], dtype=np.float32)

            if verbose:
                print(f"[Registration] Similarity - Scale: {scale:.3f}, Rotation: {np.degrees(angle):.2f}°, Translation: ({tx:.1f}, {ty:.1f})")

        else:  # AffineTransform
            # Affine transform
            matrix = actual_transform.GetMatrix()
            translation = actual_transform.GetTranslation()

            transM = np.array([
                [matrix[0], matrix[1], translation[0]],
                [matrix[2], matrix[3], translation[1]]
            ], dtype=np.float32)

            # Extract parameters from matrix for display
            scale_x = np.sqrt(matrix[0]**2 + matrix[2]**2)
            scale_y = np.sqrt(matrix[1]**2 + matrix[3]**2)
            angle = np.arctan2(matrix[2], matrix[0])
            shear = np.arctan2(matrix[1], matrix[3]) - angle

            if verbose:
                print(f"[Registration] Affine - Scale: ({scale_x:.3f}, {scale_y:.3f}), Rotation: {np.degrees(angle):.2f}°, "
                      f"Shear: {np.degrees(shear):.2f}°, Translation: ({translation[0]:.1f}, {translation[1]:.1f})")

        return transM

    @staticmethod
    def feature_based_registration(image_fixed, image_moving, he_image_gray, gem_heatmap,
                                   coarse_transform=None, use_roi=True, verbose=True):
        """
        Feature-based fine registration using SIFT/ORB

        Parameters:
        -----------
        image_fixed : ndarray
            Reference mask (gem_mask)
        image_moving : ndarray
            Moving mask (he_mask)
        he_image_gray : ndarray
            HE grayscale image (for feature extraction)
        gem_heatmap : ndarray
            Gene expression heatmap (for feature extraction)
        coarse_transform : ndarray (2x3), optional
            Coarse registration transform matrix, will be applied first if provided
        use_roi : bool
            Whether to use ROI cropping strategy to handle partial sequencing issues
        verbose : bool
            Whether to print detailed information

        Returns:
        --------
        transform_matrix : ndarray (2x3)
            Fine registration affine transformation matrix
        matched_points : int
            Number of successfully matched feature points
        """
        if verbose:
            print("[Feature Registration] Starting feature-based fine registration...")

        # Step 1: Apply coarse registration to HE image if available
        if coarse_transform is not None:
            he_image_warped = cv2.warpAffine(he_image_gray, coarse_transform,
                                            (image_fixed.shape[1], image_fixed.shape[0]))
            he_mask_warped = cv2.warpAffine(image_moving, coarse_transform,
                                           (image_fixed.shape[1], image_fixed.shape[0]))
            if verbose:
                print("[Feature Registration] Applied coarse transform to HE image")
        else:
            he_image_warped = cv2.resize(he_image_gray, (image_fixed.shape[1], image_fixed.shape[0]))
            he_mask_warped = cv2.resize(image_moving, (image_fixed.shape[1], image_fixed.shape[0]))

        # Step 2: ROI extraction - solve partial sequencing problem
        if use_roi:
            # Get gem_mask bounding box (this is the sequencing region)
            contours_fixed, _ = cv2.findContours(image_fixed.astype(np.uint8),
                                                cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if len(contours_fixed) == 0:
                if verbose:
                    print("[Feature Registration] WARNING: No contours in gem_mask, using full image")
                use_roi = False
            else:
                x, y, w, h = cv2.boundingRect(max(contours_fixed, key=cv2.contourArea))
                # Expand ROI boundary by 10% to ensure edge features are included
                margin = int(max(w, h) * 0.1)
                x_roi = max(0, x - margin)
                y_roi = max(0, y - margin)
                w_roi = min(image_fixed.shape[1] - x_roi, w + 2 * margin)
                h_roi = min(image_fixed.shape[0] - y_roi, h + 2 * margin)

                # Crop ROI region
                gem_roi = gem_heatmap[y_roi:y_roi+h_roi, x_roi:x_roi+w_roi]
                he_roi = he_image_warped[y_roi:y_roi+h_roi, x_roi:x_roi+w_roi]
                gem_mask_roi = image_fixed[y_roi:y_roi+h_roi, x_roi:x_roi+w_roi]

                if verbose:
                    print(f"[Feature Registration] ROI extracted: ({x_roi},{y_roi}) size={w_roi}x{h_roi}")
                    print(f"  Full image size: {image_fixed.shape}, ROI covers {(w_roi*h_roi)/(image_fixed.shape[0]*image_fixed.shape[1])*100:.1f}% of area")
        else:
            gem_roi = gem_heatmap
            he_roi = he_image_warped
            gem_mask_roi = image_fixed
            x_roi, y_roi = 0, 0

        # Step 3: Enhance image contrast to improve feature point quality
        gem_roi_enhanced = cv2.equalizeHist(gem_roi) if gem_roi.dtype == np.uint8 else cv2.equalizeHist(
            (gem_roi / gem_roi.max() * 255).astype(np.uint8))
        he_roi_enhanced = cv2.equalizeHist(he_roi) if he_roi.dtype == np.uint8 else cv2.equalizeHist(
            (he_roi / he_roi.max() * 255).astype(np.uint8))

        # Step 4: Feature point detection - prefer SIFT, fallback to ORB if unavailable
        try:
            # Try using SIFT (requires opencv-contrib-python)
            detector = cv2.SIFT_create(nfeatures=2000)
            method_name = "SIFT"
        except AttributeError:
            # If SIFT unavailable, use ORB
            detector = cv2.ORB_create(nfeatures=2000)
            method_name = "ORB"

        if verbose:
            print(f"[Feature Registration] Using {method_name} feature detector")

        # Detect keypoints and compute descriptors
        kp_gem, desc_gem = detector.detectAndCompute(gem_roi_enhanced, gem_mask_roi)
        kp_he, desc_he = detector.detectAndCompute(he_roi_enhanced, None)

        if desc_gem is None or desc_he is None or len(kp_gem) < 4 or len(kp_he) < 4:
            if verbose:
                print(f"[Feature Registration] WARNING: Insufficient features (gem={len(kp_gem) if kp_gem else 0}, he={len(kp_he) if kp_he else 0})")
                print("[Feature Registration] Falling back to coarse transform")
            return coarse_transform if coarse_transform is not None else np.eye(2, 3, dtype=np.float32), 0

        if verbose:
            print(f"[Feature Registration] Detected features: gem={len(kp_gem)}, he={len(kp_he)}")

        # Step 5: Feature point matching
        if method_name == "SIFT":
            # SIFT uses FLANN matcher
            FLANN_INDEX_KDTREE = 1
            index_params = dict(algorithm=FLANN_INDEX_KDTREE, trees=5)
            search_params = dict(checks=50)
            matcher = cv2.FlannBasedMatcher(index_params, search_params)
        else:
            # ORB uses BFMatcher
            matcher = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=False)

        matches = matcher.knnMatch(desc_gem, desc_he, k=2)

        # Lowe's ratio test to filter matching points
        good_matches = []
        for m_n in matches:
            if len(m_n) == 2:
                m, n = m_n
                if m.distance < 0.7 * n.distance:  # ratio test threshold
                    good_matches.append(m)

        if verbose:
            print(f"[Feature Registration] Good matches after ratio test: {len(good_matches)}")

        if len(good_matches) < 4:
            if verbose:
                print("[Feature Registration] WARNING: Insufficient good matches, falling back to coarse transform")
            return coarse_transform if coarse_transform is not None else np.eye(2, 3, dtype=np.float32), 0

        # Step 6: Extract matching point coordinates
        pts_gem = np.float32([kp_gem[m.queryIdx].pt for m in good_matches])
        pts_he = np.float32([kp_he[m.trainIdx].pt for m in good_matches])

        # Convert ROI coordinates back to full image coordinates
        if use_roi:
            pts_gem[:, 0] += x_roi
            pts_gem[:, 1] += y_roi
            pts_he[:, 0] += x_roi
            pts_he[:, 1] += y_roi

        # Step 7: Estimate affine transform using RANSAC (from HE to GEM)
        # Note: We transform from he coordinate system to gem coordinate system, so pts_he is source, pts_gem is target
        transform_refine, inliers = cv2.estimateAffinePartial2D(
            pts_he, pts_gem,
            method=cv2.RANSAC,
            ransacReprojThreshold=5.0,  # RANSAC threshold (pixels)
            maxIters=2000,
            confidence=0.99
        )

        if transform_refine is None:
            if verbose:
                print("[Feature Registration] WARNING: RANSAC failed, falling back to coarse transform")
            return coarse_transform if coarse_transform is not None else np.eye(2, 3, dtype=np.float32), 0

        inlier_count = np.sum(inliers) if inliers is not None else 0

        if verbose:
            print(f"[Feature Registration] RANSAC inliers: {inlier_count}/{len(good_matches)}")
            print(f"[Feature Registration] Transform matrix:\n{transform_refine}")

        # Step 8: Combine transformation matrices if coarse registration exists
        if coarse_transform is not None:
            # Convert two 2x3 matrices to 3x3 for matrix multiplication
            coarse_3x3 = np.vstack([coarse_transform, [0, 0, 1]])
            refine_3x3 = np.vstack([transform_refine, [0, 0, 1]])
            combined_3x3 = refine_3x3 @ coarse_3x3
            final_transform = combined_3x3[:2, :]
            if verbose:
                print("[Feature Registration] Combined coarse + refine transforms")
        else:
            final_transform = transform_refine

        return final_transform, inlier_count

    @staticmethod
    def contours_registration(image_fixed, image_moving):
        """
        Contour-based registration (supports images of different sizes and multiple regions)

        Parameters:
        -----------
        image_fixed : ndarray
            Reference mask image (target coordinate system, usually gem_mask)
        image_moving : ndarray
            Moving mask image (source coordinate system, usually he_mask)

        Returns:
        --------
        transform_matrix : ndarray (2x3)
            Affine transformation matrix, transforms moving coordinate system to fixed coordinate system
        """
        def get_global_features(mask_image, contours):
            """
            Calculate global center and angle from ALL contours (multi-region safe).
            Uses image moments to compute the overall centroid across all regions,
            and fitEllipseAMS on the merged contour points for angle estimation.
            """
            # Use image moments for global centroid (handles multi-region correctly)
            M = cv2.moments(mask_image)
            if M['m00'] == 0:
                raise ValueError("Mask image has zero area, cannot compute centroid")
            cx = int(M['m10'] / M['m00'])
            cy = int(M['m01'] / M['m00'])

            # For angle: merge all contour points and fit ellipse
            all_points = np.vstack(contours)
            if len(all_points) >= 5:  # fitEllipseAMS requires at least 5 points
                (x, y), (MA, ma), angle = cv2.fitEllipseAMS(all_points)
            else:
                angle = 0.0

            # Total area = sum of all contour areas (or use mask pixel count)
            total_area = float(np.sum(mask_image > 0))

            return (cx, cy), angle, total_area

        contours_fixed, _ = cv2.findContours(image_fixed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        contours_moving, _ = cv2.findContours(image_moving, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        center_fixed, angle_fixed, area_fixed = get_global_features(image_fixed, contours_fixed)
        center_moving, angle_moving, area_moving = get_global_features(image_moving, contours_moving)

        print(f"[Contour Features] Fixed  - center: {center_fixed}, angle: {angle_fixed:.2f}°, area: {area_fixed:.0f} px, regions: {len(contours_fixed)}")
        print(f"[Contour Features] Moving - center: {center_moving}, angle: {angle_moving:.2f}°, area: {area_moving:.0f} px, regions: {len(contours_moving)}")

        # Calculate scale based on total mask area ratio (multi-region safe)
        scale = np.sqrt(area_fixed / area_moving)

        # Global angle search: test 0° to 360° range, 1° step for each candidate
        # This finds optimal registration at any angle
        print("[Angle Search] Searching optimal rotation angle (0°-360°, step=1°)...")

        best_angle = 0
        best_score = -1
        angle_scores = []

        # Coarse search: test every 1° (360 angles total)
        coarse_angles = list(range(360))
        total_coarse = len(coarse_angles)
        print(f"[Coarse Search] Testing {total_coarse} angles...")

        for idx, test_angle in enumerate(coarse_angles):
            # Create test transformation
            test_matrix = cv2.getRotationMatrix2D(center_moving, test_angle, scale)
            test_center_transformed = test_matrix @ np.array([center_moving[0], center_moving[1], 1])
            test_tx = center_fixed[0] - test_center_transformed[0]
            test_ty = center_fixed[1] - test_center_transformed[1]
            test_matrix[0, 2] += test_tx
            test_matrix[1, 2] += test_ty

            # Apply transformation to moving mask
            moving_transformed = cv2.warpAffine(image_moving, test_matrix,
                                               (image_fixed.shape[1], image_fixed.shape[0]))

            # Calculate overlap (Dice coefficient)
            intersection = np.sum((image_fixed > 0) & (moving_transformed > 0))
            score = 2.0 * intersection / (np.sum(image_fixed > 0) + np.sum(moving_transformed > 0)) if intersection > 0 else 0

            angle_scores.append((test_angle, score))

            if score > best_score:
                best_score = score
                best_angle = test_angle

            # Progress display: every 10 angles or last one
            if (idx + 1) % 10 == 0 or (idx + 1) == total_coarse:
                progress = (idx + 1) / total_coarse * 100
                print(f"  Progress: {idx+1}/{total_coarse} ({progress:.1f}%) - Current best: {best_angle}° (Dice: {best_score:.5f})")

        print(f"[Coarse Search] Completed! Best angle: {best_angle}° (Dice: {best_score:.5f})")

        # Fine search: within ±1° of best angle, test every 0.1° (21 angles)
        fine_search_range = [best_angle + i * 0.1 for i in range(-10, 11)]
        total_fine = len(fine_search_range)
        print(f"[Fine Search] Refining around {best_angle:.1f}° (testing {total_fine} angles in ±1° range, step=0.1°)...")

        for idx, test_angle in enumerate(fine_search_range):
            # Create test transformation
            test_matrix = cv2.getRotationMatrix2D(center_moving, test_angle, scale)
            test_center_transformed = test_matrix @ np.array([center_moving[0], center_moving[1], 1])
            test_tx = center_fixed[0] - test_center_transformed[0]
            test_ty = center_fixed[1] - test_center_transformed[1]
            test_matrix[0, 2] += test_tx
            test_matrix[1, 2] += test_ty

            # Apply transformation to moving mask
            moving_transformed = cv2.warpAffine(image_moving, test_matrix,
                                               (image_fixed.shape[1], image_fixed.shape[0]))

            # Calculate overlap
            intersection = np.sum((image_fixed > 0) & (moving_transformed > 0))
            score = 2.0 * intersection / (np.sum(image_fixed > 0) + np.sum(moving_transformed > 0)) if intersection > 0 else 0

            if score > best_score:
                best_score = score
                best_angle = test_angle

            # Progress display
            progress = (idx + 1) / total_fine * 100
            print(f"  Progress: {idx+1}/{total_fine} ({progress:.1f}%) - Angle: {test_angle:.1f}° (Dice: {score:.5f})")

        # Use best angle
        rotation_angle = best_angle
        print(f"[Fine Search] Completed! Final best angle: {rotation_angle:.2f}° (Dice score: {best_score:.5f})")

        # Step 1: Rotate and scale around center_moving in moving coordinate system
        transform_matrix = cv2.getRotationMatrix2D(center_moving, rotation_angle, scale)

        # Step 2: Calculate new position of moving center after rotation+scale
        center_moving_transformed = transform_matrix @ np.array([center_moving[0], center_moving[1], 1])

        # Step 3: Calculate translation from transformed position to fixed center
        tx = center_fixed[0] - center_moving_transformed[0]
        ty = center_fixed[1] - center_moving_transformed[1]

        # Step 4: Add translation component
        transform_matrix[0, 2] += tx
        transform_matrix[1, 2] += ty

        print(f"[Registration] Final - Rotation: {rotation_angle:.2f}°, Scale: {scale:.3f}, Translation: ({tx:.1f}, {ty:.1f})")

        return transform_matrix

    @utils.add_log
    def splitplot_experiment(self, seg_type, mode):
        if seg_type == 'tissue':
            self.tissue_segmentation(mode)
        elif seg_type == 'cell':
            raise NotImplementedError("seg_type='cell' is not supported in the v1.8 binSegment flow")

    @utils.add_log
    def tissue_segmentation(self, mode="ssDNA"):
        def getArea(elem):
            return elem.area

        if mode == "ssDNA":
            # ssDNA mode: preserve all connected tissue regions from segmentation output.
            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (20, 20))
            otsu_threshold, _ = cv2.threshold(
                self.tissue_image,
                0,
                255,
                cv2.THRESH_BINARY | cv2.THRESH_OTSU,
            )
            threshold = int(round(otsu_threshold * self.ssdna_threshold_scale))
            threshold = max(0, min(255, threshold))
            _, tissue_mask = cv2.threshold(self.tissue_image, threshold, 255, cv2.THRESH_BINARY)
            if self.ssdna_threshold_scale != 1.0:
                print(
                    '[ssDNA mode] Otsu threshold adjusted: '
                    f'raw={otsu_threshold:.1f}, scale={self.ssdna_threshold_scale:.3f}, used={threshold}'
                )
            tissue_mask = cv2.morphologyEx(tissue_mask, cv2.MORPH_CLOSE, kernel, iterations=8)

            # Do not fill every enclosed hole: ssDNA slides commonly contain
            # genuine blank regions inside the tissue.  Fill only tiny holes
            # caused by pixel noise and leave larger internal cavities as
            # background so they are excluded from bin counts.
            tissue_bool = tissue_mask > 0
            hole_area = int(getattr(self.args, 'ssdna_min_hole_area', 5000))
            if hole_area > 0:
                inv = (~tissue_bool).astype(np.uint8)
                labels, num = ndi.label(inv)
                sizes = ndi.sum(inv, labels, index=np.arange(1, num + 1))
                border = set(np.unique(np.concatenate([
                    labels[0, :], labels[-1, :], labels[:, 0], labels[:, -1]
                ])))
                tissue_result_bool = tissue_bool.copy()
                for label_id, size in enumerate(sizes, start=1):
                    if label_id not in border and size < hole_area:
                        tissue_result_bool[labels == label_id] = True
                tissue_result = tissue_result_bool.astype(np.uint8) * 255
                print(f"[ssDNA mode] Filled only internal holes smaller than {hole_area} pixels")
            else:
                tissue_result = tissue_bool.astype(np.uint8) * 255

            # Additional background noise removal: remove small objects and erode slightly
            # This helps eliminate edge artifacts and isolated noise pixels
            small_kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
            tissue_result = skimage.morphology.remove_small_objects(tissue_result > 0, min_size=1000).astype(np.uint8) * 255
            # Slight erosion to remove noisy edges (adjust iterations to control strictness)
            tissue_result = cv2.erode(tissue_result, small_kernel, iterations=2)
            # Re-dilate to recover tissue area (but noise won't come back)
            tissue_result = cv2.dilate(tissue_result, small_kernel, iterations=2)
            if self.ssdna_mask_expand_pixels:
                diameter = 2 * self.ssdna_mask_expand_pixels + 1
                expand_kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (diameter, diameter))
                tissue_result = cv2.dilate(tissue_result, expand_kernel, iterations=1)
                print(
                    '[ssDNA mode] Expanded final tissue mask by '
                    f'{self.ssdna_mask_expand_pixels} pixel(s)'
                )

            self.tissue_image_mask = tissue_result
            self._set_background_image_mask(tissue_result)
            print("[ssDNA mode] Tissue mask generated from segmentation output with all connected regions preserved")
            if self.gene_mask_filter:
                self._apply_gene_mask_filter_to_tissue_mask("ssDNA")
        elif mode == "gene_expr":
            # Pure gene expression mode: use GEM mask
            self.tissue_image_mask = self.gem_mask
            self._set_background_image_mask(self.gem_mask)
            print("[gene_expr mode] Using GEM mask as tissue mask (expression-based)")
            if self.gene_mask_filter:
                print("[gene_expr mode] gene-mask filter enabled; GEM/manual mask is the final tissue region")

        elif mode == "HE":
            # HE mode: supported tissue mask strategies
            strategy = getattr(self, 'tissue_mask_strategy', 'gem_primary')
            background_mask = self._get_background_image_mask()

            print(f"\n[HE Mode] Tissue mask strategy: {strategy}")

            if strategy == 'gem_primary':
                self.tissue_image_mask = self._generate_gem_primary_tissue_mask()
                print("[HE mode] Using GEM-primary tissue mask (prevents phantom expression)")

            elif strategy == 'intersection':
                self.tissue_image_mask = cv2.bitwise_and(self.gem_mask, background_mask)
                kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
                self.tissue_image_mask = cv2.dilate(self.tissue_image_mask, kernel, iterations=2)
                print("[HE mode] Using GEM ∩ HE intersection (conservative, high confidence only)")

            elif strategy == 'union':
                self.tissue_image_mask = cv2.bitwise_or(self.gem_mask, background_mask)
                print("[HE mode] Using GEM ∪ HE union (liberal, may include false positives)")

            else:
                # Fallback to gem_primary
                print(f"[HE mode] Unknown strategy '{strategy}', falling back to 'gem_primary'")
                self.tissue_image_mask = self._generate_gem_primary_tissue_mask()

        tifffile.imwrite(os.path.join(self.exp_bin_im, f'{self.tissue_name}_tissue_cut.tif'), self.tissue_image_mask)
        self._refresh_ssdna_spatial_previews()
        if mode == "HE":
            # The initial registration overlay is written before the final tissue mask exists.
            # Refresh it now so the customer-facing GEM overlay excludes off-tissue noise.
            self.save_gem_background_image_overlay()

    def _apply_gene_mask_filter_to_tissue_mask(self, mode_label):
        if self.tissue_image_mask is None:
            print(f"[Gene Mask Filter] WARNING: tissue mask is not available in {mode_label} mode")
            return
        if self.gem_mask is None:
            print(f"[Gene Mask Filter] WARNING: GEM/manual mask is not available in {mode_label} mode; keeping original tissue mask")
            return

        gene_mask = self.gem_mask
        if gene_mask.shape[:2] != self.tissue_image_mask.shape[:2]:
            gene_mask = cv2.resize(
                gene_mask,
                (self.tissue_image_mask.shape[1], self.tissue_image_mask.shape[0]),
                interpolation=cv2.INTER_NEAREST,
            )

        tissue_before = int(np.sum(self.tissue_image_mask > 0))
        gene_area = int(np.sum(gene_mask > 0))
        if mode_label == "ssDNA" and self.gem_manual_mask_path:
            filtered = (gene_mask > 0).astype(np.uint8) * 255
            tissue_after = int(np.sum(filtered > 0))
            self.tissue_image_mask = filtered
            self._set_background_image_mask(filtered)
            print(
                f"[Gene Mask Filter] Applied manual mask override in ssDNA mode: "
                f"manual_mask={self.gem_manual_mask_path}, "
                f"ssDNA_tissue_area={tissue_before}, manual_mask_area={gene_area}, kept={tissue_after}"
            )
            return

        filtered = cv2.bitwise_and(
            (self.tissue_image_mask > 0).astype(np.uint8) * 255,
            (gene_mask > 0).astype(np.uint8) * 255,
        )
        tissue_after = int(np.sum(filtered > 0))
        self.tissue_image_mask = filtered
        print(
            f"[Gene Mask Filter] Applied in {mode_label} mode: "
            f"tissue_area={tissue_before}, gene_mask_area={gene_area}, kept_intersection={tissue_after}"
        )

    def _generate_gem_primary_tissue_mask(self):
        """
        Generate tissue mask with GEM as primary source, HE for boundary assistance and QC

        Design Philosophy:
        ------------------
        - Gene expression data (GEM) is ground truth for tissue location
        - HE image is used for visualization and optional boundary refinement
        - This prevents "phantom expression" in areas where HE shows tissue but no RNA exists

        Returns:
        --------
        tissue_mask : ndarray
            Final tissue mask based primarily on GEM expression
        """
        print("\n" + "="*70)
        print("[Tissue Mask Generation] GEM-Primary Strategy")
        print("="*70)

        # 1. Start with GEM mask as base
        tissue_mask = self.gem_mask.copy()
        gem_area = np.sum(tissue_mask > 0)
        print(f"[Step 1] Base GEM mask area: {gem_area:,} pixels")

        # 2. Check if HE mask is available for boundary assistance
        background_mask = self._get_background_image_mask()
        if background_mask is None:
            print("[Step 2] No background mask available, using pure GEM mask")
            return tissue_mask

        he_area = np.sum(background_mask > 0)
        print(f"[Step 2] HE mask area: {he_area:,} pixels")

        # 3. Calculate overlap metrics for QC
        intersection = cv2.bitwise_and(self.gem_mask, background_mask)
        union = cv2.bitwise_or(self.gem_mask, background_mask)
        intersection_area = np.sum(intersection > 0)
        union_area = np.sum(union > 0)

        # IoU (Intersection over Union)
        iou = intersection_area / union_area if union_area > 0 else 0

        # Overlap ratio relative to GEM
        gem_overlap_ratio = intersection_area / gem_area if gem_area > 0 else 0

        # Overlap ratio relative to HE
        he_overlap_ratio = intersection_area / he_area if he_area > 0 else 0

        print(f"[Step 3] Overlap Analysis:")
        print(f"  - Intersection area: {intersection_area:,} pixels")
        print(f"  - Union area: {union_area:,} pixels")
        print(f"  - IoU (Intersection/Union): {iou:.1%}")
        print(f"  - GEM overlap ratio: {gem_overlap_ratio:.1%}")
        print(f"  - HE overlap ratio: {he_overlap_ratio:.1%}")

        # 4. Decide whether to use HE for boundary assistance
        use_he_boundary = False

        if iou >= 0.75:
            # High overlap: HE and GEM are well aligned, can use HE for boundary refinement
            use_he_boundary = True
            print(f"[Step 4] High overlap (IoU={iou:.1%}), HE boundary will assist edge refinement")
        elif iou >= 0.5:
            # Medium overlap: Use with caution
            use_he_boundary = True
            print(f"[Step 4] Medium overlap (IoU={iou:.1%}), HE boundary used conservatively")
        else:
            # Low overlap: Don't trust HE boundary
            print(f"[Step 4] ⚠️  Low overlap (IoU={iou:.1%}), HE boundary NOT used")
            print(f"[Step 4] Possible causes: registration error, tissue damage, or necrosis")

        # 5. Apply HE boundary assistance if appropriate
        if use_he_boundary:
            # Strategy: Expand GEM mask slightly, then intersect with HE mask
            # This allows HE to smooth GEM boundaries without introducing false regions
            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))

            # Dilate GEM mask to create expansion zone
            gem_dilated = cv2.dilate(self.gem_mask, kernel, iterations=3)

            # Only keep expanded regions that are also in HE mask
            he_boundary_assist = cv2.bitwise_and(gem_dilated, background_mask)

            # Combine with original GEM mask
            tissue_mask_assisted = cv2.bitwise_or(tissue_mask, he_boundary_assist)

            added_area = np.sum(tissue_mask_assisted > 0) - gem_area
            print(f"[Step 5] HE boundary assistance added {added_area:,} pixels ({added_area/gem_area:.1%} of GEM area)")

            tissue_mask = tissue_mask_assisted
        else:
            print(f"[Step 5] Skipped HE boundary assistance")

        # 6. Morphological refinement
        kernel_small = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
        tissue_mask = cv2.morphologyEx(tissue_mask, cv2.MORPH_CLOSE, kernel_small, iterations=2)
        tissue_mask = cv2.morphologyEx(tissue_mask, cv2.MORPH_OPEN, kernel_small, iterations=1)

        final_area = np.sum(tissue_mask > 0)
        print(f"[Step 6] Morphological refinement complete, final area: {final_area:,} pixels")

        # 7. Quality Control: Detect HE-only regions (potential issues)
        he_only = cv2.bitwise_and(background_mask, cv2.bitwise_not(tissue_mask))
        he_only_area = np.sum(he_only > 0)

        if he_only_area > 0:
            he_only_ratio = he_only_area / he_area
            print(f"\n[QC Check] HE-only regions (no gene expression):")
            print(f"  - Area: {he_only_area:,} pixels ({he_only_ratio:.1%} of HE tissue)")

            if he_only_ratio > 0.2:
                print(f"  - ⚠️  WARNING: Large HE-only region detected!")
                print(f"  - Possible causes:")
                print(f"    • Necrotic tissue (visible in HE but RNA degraded)")
                print(f"    • Fibrotic/stromal tissue (low cell density)")
                print(f"    • Registration error (misalignment)")
                print(f"    • Blood vessels or empty spaces")

                # Save QC image
                qc_dir = os.path.join(self.outdir, 'images', 'qc')
                os.makedirs(qc_dir, exist_ok=True)
                qc_path = os.path.join(qc_dir, f'{self.tissue_name}_he_no_expression.png')
                cv2.imwrite(qc_path, he_only)
                print(f"  - QC image saved: {qc_path}")
            else:
                print(f"  - ✓ Acceptable level of HE-only regions")

        # 8. Summary
        print(f"\n[Summary] Tissue Mask Generation Complete:")
        print(f"  - Strategy: GEM-primary (expression-based)")
        print(f"  - Final mask area: {final_area:,} pixels")
        print(f"  - GEM coverage: {gem_area/final_area:.1%} of final mask")
        print(f"  - HE boundary assistance: {'Yes' if use_he_boundary else 'No'}")
        print("="*70 + "\n")

        return tissue_mask

    @staticmethod
    def get_offset(image):
        def get_boundary_point(image):
            mask = image.copy()
            mask[image != 0] = 255
            mask = cv2.dilate(mask, np.ones((3, 3), np.uint8), iterations=5)

            boundary_points = []
            contours, hierarchy = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            contours = sorted(contours, key=cv2.contourArea, reverse=True)[:4]
            for contour in contours:
                # use Douglas-Peucker algorithm to approximate the contour
                epsilon = 0.009 * cv2.arcLength(contour, True)
                approx = cv2.approxPolyDP(contour, epsilon, True)
                # if the approximated contour has four(or more) points, then assume that screen is found
                # make sure it is a square at least
                if len(approx) >= 4:
                    boundary_points.extend(approx)

            return boundary_points

        boundary_points = get_boundary_point(image)
        boundary_points = np.array(boundary_points)
        left_top_index = np.linalg.norm(boundary_points - np.array([0, 0]), axis=2).argmin()
        transform_dis = (boundary_points[left_top_index].squeeze() - 100).clip(min=0)
        return transform_dis.astype(np.int32)

    @utils.add_log
    def generate_tissue_barcode_position(self):
        bcd = self.barcodes_detail.to_numpy()
        locs = bcd[:, :2].astype(int)

        # write tissue barcode position
        in_tissue = np.where(self.segment, np.where(self.tissue_image_mask[locs[:, 1], locs[:, 0]] == 255, 1, 0), 1)
        fi_bc = pd.DataFrame({
            'barcode': bcd[:, 2],
            'in_tissue': in_tissue,
            'array_row': None,
            'array_col': None,
            'pxl_row_in_fullres': bcd[:, 1],
            'pxl_col_in_fullres': bcd[:, 0],
        })
        fi_bc = fi_bc.sort_values(by=['pxl_row_in_fullres', 'pxl_col_in_fullres']).reset_index(drop=True)
        fi_bc = fi_bc[['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']]
        fi_bc.to_csv(os.path.join(self.exp_bs_path, f'{self.tissue_name}_Barcodes_tissue_positions.csv'), index=False)
        self.tb_position = fi_bc

    @utils.add_log
    def generate_count_detail(self):
        barcode_columns = ['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
        self.in_bc = self.tb_position.loc[self.tb_position['in_tissue'] == 1].reset_index(drop=True)

        if self.rectify:
            dx, dy = self.roi_rect[:2]
            self.in_bc["pxl_col_in_fullres"] -= dx
            self.in_bc["pxl_row_in_fullres"] -= dy
            self.in_bc["pxl_col_in_fullres"] = self.in_bc["pxl_col_in_fullres"].clip(0, self.slice_px - 1)
            self.in_bc["pxl_row_in_fullres"] = self.in_bc["pxl_row_in_fullres"].clip(0, self.slice_px - 1)

            self.tissue_image = self.tissue_image[self.roi_rect[1]: self.roi_rect[1] + self.roi_rect[3], self.roi_rect[0]: self.roi_rect[0] + self.roi_rect[2]]
            background_mask = self._get_background_image_mask()
            background_mask = background_mask[self.roi_rect[1]: self.roi_rect[1] + self.roi_rect[3], self.roi_rect[0]: self.roi_rect[0] + self.roi_rect[2]]
            self.tissue_image = cv2.resize(self.tissue_image, (self.slice_px, self.slice_px))
            background_mask = cv2.resize(background_mask, (self.slice_px, self.slice_px))
            background_mask[background_mask > 0] = 255
            self._set_background_image_mask(background_mask)

            self.save_register_image(is_correct=True)

        self.gene_count_detail = pd.read_csv(self.count_detail, sep='\t', dtype={'Barcode': str, 'geneID': str, 'UMI': str, 'count': str})
        self.gene_count_detail = self.gene_count_detail.merge(self.in_bc[barcode_columns], left_on='Barcode', right_on='barcode', how='inner')
        self.gene_count_detail = self.gene_count_detail.drop(columns=['barcode'])
        if self._raw_show:
            self.raw_count_detail = self.gene_count_detail.copy(deep=True)
            self.raw_count_detail = self.raw_count_detail.rename(columns={'Barcode': 'unicode'})
            self.raw_count_detail = self.raw_count_detail.rename(columns={'pxl_row_in_fullres': 'row_in_bin'})
            self.raw_count_detail = self.raw_count_detail.rename(columns={'pxl_col_in_fullres': 'col_in_bin'})

    @utils.add_log
    def spatial_bin_generate(self):
        import os  # Import at function start to avoid scoping issues
        import gc
        mask_shape = self.tissue_image_mask.shape
        bbox = [0, 0, mask_shape[1], mask_shape[0]]

        def make_bin_args(*, copy_count_detail: bool) -> list[list[object]]:
            """Build bin handlers, optionally copying the large count table.

            A count-detail DataFrame can occupy many gigabytes.  In the normal
            service configuration (``BINSEGMENT_MAX_PARALLEL=1``), keeping a
            deep copy for every bin defeats the intended memory limit even
            though only one bin is processed at a time.  Sequential callers
            therefore create one copy lazily per iteration below.
            """
            args: list[list[object]] = []
            for i in range(len(self._micron_bin)):
                detail = self.gene_count_detail.copy(deep=True) if copy_count_detail else None
                args.append([
                    i,
                    SquareBin(
                        self._micron_bin[i],
                        self.bin[i],
                        bbox,
                        bin_folder=os.path.join(self.exp_sq_bin, f'{self.tissue_name}_bin{self._micron_bin[i]}'),
                        gene_count_detail=detail,
                    ),
                ])
            if self._raw_show:
                args.append([
                    None,
                    SquareBin(
                        'Raw',
                        self._bin,
                        bbox,
                        bin_folder=os.path.join(self.exp_sq_bin, f'{self.tissue_name}_Raw'),
                        gene_count_detail=self.raw_count_detail,
                    ),
                ])
            return args

        # With one worker, avoid retaining a deep copy for every bin.  Keep the
        # lightweight handlers for get_bin_summary(), but release each large
        # DataFrame as soon as its bin is complete.
        max_parallel = int(os.environ.get('BINSEGMENT_MAX_PARALLEL', '2')) if self._parallel else 1
        if max_parallel <= 1:
            bin_args: list[list[object]] = []
            for i in range(len(self._micron_bin)):
                handler = SquareBin(
                    self._micron_bin[i],
                    self.bin[i],
                    bbox,
                    bin_folder=os.path.join(self.exp_sq_bin, f'{self.tissue_name}_bin{self._micron_bin[i]}'),
                    gene_count_detail=self.gene_count_detail.copy(deep=True),
                )
                args = [i, handler]
                self.bin_files_execute(args)
                handler.gene_count_detail = None
                bin_args.append(args)
                gc.collect()
            if self._raw_show:
                handler = SquareBin(
                    'Raw', self._bin, bbox,
                    bin_folder=os.path.join(self.exp_sq_bin, f'{self.tissue_name}_Raw'),
                    gene_count_detail=self.raw_count_detail,
                )
                args = [None, handler]
                self.bin_files_execute(args)
                handler.gene_count_detail = None
                bin_args.append(args)
                gc.collect()
        else:
            bin_args = make_bin_args(copy_count_detail=True)
            # Memory-optimized parallel processing:
            # - Use threading backend to share memory instead of forking processes
            # - Limit n_jobs to 2 to reduce memory footprint (balance speed vs memory)
            # - Original: n_jobs=5 uses ~300GB RAM, crashes with OOM
            # - Optimized: n_jobs=2 uses ~80GB RAM, only 1.5x slower
            n_jobs = min(max_parallel, len(self._micron_bin))
            print(f"[Memory Optimization] Running bin generation with n_jobs={n_jobs} (threading backend)")
            print(f"  - Total bins to process: {len(bin_args)}")
            print(f"  - Expected memory usage: ~{n_jobs * 40}GB (vs ~{len(bin_args) * 60}GB if fully parallel)")
            Parallel(n_jobs=n_jobs, backend='threading', verbose=10)(
                map(delayed(self.bin_files_execute), bin_args)
            )
        self.bin_exec = bin_args

        del self.gene_count_detail  # release memory after bin segment
        del self.raw_count_detail

    @utils.add_log
    def bin_files_execute(self, args):
        import gc
        i, handler = args
        print(f"[Bin {handler.micron_bin}] Processing started...")

        is_raw = handler.micron_bin == 'Raw'
        total_reads_before_filter = None

        if is_raw:
            gene = self.raw_count_detail
        elif isinstance(handler.micron_bin, int):
            self.bin_segment(handler)
            gene = handler.gene_count_detail

            # Calculate total reads BEFORE filtering (for accurate Fraction Reads in Square calculation)
            # Must do this BEFORE filtering out background bins
            unique_gene_unfiltered = self.unique_gene_count(gene)
            total_reads_before_filter = unique_gene_unfiltered['count'].astype(int).sum()
            print(f"[Bin {handler.micron_bin}] Total reads before filtering: {total_reads_before_filter:,}")
            del unique_gene_unfiltered
            gc.collect()

            # Filter out bins that are not in tissue (background removal)
            if 'in_tissue_bin' in gene.columns:
                n_before = len(gene)
                gene = gene[gene['in_tissue_bin'] == True].copy()
                n_after = len(gene)
                print(f"[Bin {handler.micron_bin}] Filtered out {n_before - n_after} UMIs in background bins (kept {n_after}/{n_before} UMIs in tissue)")
                # Free memory from original unfiltered data
                del handler.gene_count_detail
                handler.gene_count_detail = gene
                gc.collect()
        else:
            raise ValueError("Please check the input parameters <bin>")

        # Generate tissue barcode positions first (this filters to in_tissue=1 bins only)
        tissue_barcode = self.tissue_barcode_position(gene)

        # Extract only barcodes that are in tissue
        in_tissue_barcodes = set(tissue_barcode['Barcode'].unique())

        # Generate unique gene count (now from filtered data)
        unique_gene = self.unique_gene_count(gene)

        # Filter unique_gene to only include in-tissue barcodes (double-check background removal)
        n_before = len(unique_gene)
        unique_gene = unique_gene[unique_gene['Barcode'].isin(in_tissue_barcodes)]
        n_after = len(unique_gene)
        if n_before > n_after:
            print(f"[Bin {handler.micron_bin}] Additional filtering: removed {n_before - n_after} UMI records from background bins in gene matrix")

        print(f"[Bin {handler.micron_bin}] Total reads after filtering: {unique_gene['count'].astype(int).sum():,}")
        self.gene_indicator_stats(handler, unique_gene, tissue_barcode, total_reads_before_filter)
        self.generate_gene_matrix(handler, unique_gene)

        base_dir = handler.bin_folder
        matrix_dir = os.path.join(base_dir, "filtered_feature_bc_matrix")
        spatial_dir = os.path.join(base_dir, "spatial")

        os.makedirs(matrix_dir, exist_ok=True)

        gene_tsv_path = os.path.join(handler.bin_folder, "matrix", "genes.tsv")
        print(f"Checking genes.tsv existence: {gene_tsv_path} -> exists: {os.path.exists(gene_tsv_path)}")
        if os.path.exists(gene_tsv_path):
            features_gz_path = os.path.join(matrix_dir, "features.tsv.gz")
            try:
                with open(gene_tsv_path, 'r') as f_in:
                    with gzip.open(features_gz_path, 'wt') as f_out:
                        for line in f_in:
                            f_out.write(line.strip() + "\tGene Expression\n")

                print(f"Generated features.tsv.gz at {features_gz_path}")
                os.remove(gene_tsv_path)
                print(f"Removed redundant genes.tsv at {gene_tsv_path}")
            except Exception as e:
                print(f"Error processing genes.tsv: {e}")

        for ext in ["barcodes.tsv", "features.tsv", "matrix.mtx"]:
            src_file = os.path.join(handler.bin_folder, "matrix", ext)
            dst_file = os.path.join(matrix_dir, ext)
            if os.path.exists(src_file):
                shutil.move(src_file, dst_file)
                gz_file = f"{dst_file}.gz"
                with open(dst_file, "rb") as f_in, gzip.open(gz_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(dst_file)
                print(f"Generated {gz_file}")


        tissue_positions_list_path = os.path.join(spatial_dir, "tissue_positions_list.csv")
        tissue_barcode.to_csv(tissue_positions_list_path, header=True, index=False)
        print(f"Generated tissue_positions_list.csv at {tissue_positions_list_path}")

        hires_img = os.path.join(self.exp_bin_im, "tissue_hires_image.png")
        lowres_img = os.path.join(self.exp_bin_im, "tissue_lowres_image.png")
        if os.path.exists(hires_img):
            shutil.copy(hires_img, os.path.join(spatial_dir, "tissue_hires_image.png"))
        if os.path.exists(lowres_img):
            shutil.copy(lowres_img, os.path.join(spatial_dir, "tissue_lowres_image.png"))


        scalefactor_file = os.path.join(handler.bin_folder, "scalefactors_json.json")
        if os.path.exists(scalefactor_file):
            shutil.move(scalefactor_file, os.path.join(spatial_dir, "scalefactors_json.json"))

        try:
            self.scale_tissue_positions(spatial_dir)
        except Exception as e:
            print(f"[ERROR] tissue_positions: {e}")


        count_detail_file = os.path.join(handler.bin_folder,
                                         f"{self.tissue_name}_Bin{handler.micron_bin}_count_detail.txt")
        if os.path.exists(count_detail_file):
            shutil.move(count_detail_file, os.path.join(base_dir, "stat.txt"))
             # Clean up temporary matrix folder after processing
        matrix_temp_dir = os.path.join(handler.bin_folder, 'matrix')
        if os.path.exists(matrix_temp_dir):
            shutil.rmtree(matrix_temp_dir)
            print(f"Cleaned up temporary matrix folder: {matrix_temp_dir}")

    @utils.add_log
    def scale_tissue_positions(self, spatial_dir):
        coords_file = os.path.join(spatial_dir, "tissue_positions_list.csv")
        scalefactors_file = os.path.join(spatial_dir, "scalefactors_json.json")
        output_file = os.path.join(spatial_dir, "tissue_positions_scaled.csv")
        output_csv_file = os.path.join(spatial_dir, "tissue_positions.csv")

        if not os.path.exists(coords_file):
            raise FileNotFoundError(f"The coords file was not found: {coords_file}")
        coords_data = pd.read_csv(coords_file)

        if not os.path.exists(scalefactors_file):
            raise FileNotFoundError(f"The scalefactors file was not found: {scalefactors_file}")
        with open(scalefactors_file) as f:
            scalefactors = json.load(f)
        spot_diameter = float(
            scalefactors.get(
                "spot_diameter_fullres",
                scalefactors.get("fiducial_diameter_fullres", 1.0),
            )
        )

        fullres_height, fullres_width = self._fullres_spatial_shape()
        hires_height, hires_width = self._spatial_image_shape(os.path.join(spatial_dir, "tissue_hires_image.png"), (2000, 2000))
        lowres_height, lowres_width = self._spatial_image_shape(os.path.join(spatial_dir, "tissue_lowres_image.png"), (600, 600))

        row_values = pd.to_numeric(coords_data["row_in_bin"], errors="coerce").fillna(0)
        col_values = pd.to_numeric(coords_data["col_in_bin"], errors="coerce").fillna(0)

        if spot_diameter <= 1.5:
            coords_data["pxl_row_in_fullres"] = row_values
            coords_data["pxl_col_in_fullres"] = col_values
        else:
            coords_data["pxl_row_in_fullres"] = (row_values + 0.5) * spot_diameter
            coords_data["pxl_col_in_fullres"] = (col_values + 0.5) * spot_diameter

        coords_data["pxl_row_in_fullres"] = coords_data["pxl_row_in_fullres"].clip(0, fullres_height - 1)
        coords_data["pxl_col_in_fullres"] = coords_data["pxl_col_in_fullres"].clip(0, fullres_width - 1)

        scalefactors["tissue_hires_scalef"] = float(max(hires_height, hires_width)) / float(max(fullres_height, fullres_width))
        scalefactors["tissue_lowres_scalef"] = float(max(lowres_height, lowres_width)) / float(max(fullres_height, fullres_width))
        scalefactors["fiducial_diameter_fullres"] = spot_diameter
        scalefactors["spot_diameter_fullres"] = spot_diameter
        with open(scalefactors_file, "w") as handle:
            json.dump(scalefactors, handle, indent=4)

        cols = ["Barcode", "in_tissue", "row_in_bin", "col_in_bin", "pxl_row_in_fullres", "pxl_col_in_fullres"]
        coords_data = coords_data[cols]
        coords_data.rename(
            columns={
                "Barcode": "barcode",
                "row_in_bin": "array_row",
                "col_in_bin": "array_col",
            },
            inplace=True,
        )

        coords_data.to_csv(output_file, index=False, header=False)
        coords_data.to_csv(output_csv_file, index=False, header=True)

        if os.path.exists(coords_file):
            os.remove(coords_file)
        os.rename(output_file, coords_file)
        print(f"The original has been updated tissue_positions_list.csv: {coords_file}")

    def _fullres_spatial_shape(self):
        if getattr(self, "tissue_image_mask", None) is not None:
            return int(self.tissue_image_mask.shape[0]), int(self.tissue_image_mask.shape[1])
        if getattr(self, "tissue_image", None) is not None and getattr(self.tissue_image, "shape", None) is not None:
            return int(self.tissue_image.shape[0]), int(self.tissue_image.shape[1])
        return int(self.chip_shape[0]), int(self.chip_shape[1])

    @staticmethod
    def _spatial_image_shape(image_path, fallback):
        if not os.path.exists(image_path):
            return fallback
        try:
            image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
            if image is None:
                return fallback
            return int(image.shape[0]), int(image.shape[1])
        except Exception:
            return fallback


    def bin_segment(self, bin_handler):
        import gc
        self.bs_mtx = bin_handler.upsampling_tissue_matrix()  # generate gene matrix
        up_row = bin_handler.up_row
        up_col = bin_handler.up_col
        bin = bin_handler.bin
        num_bits = bin_handler.num_bits

        # broadcast variables
        row_bins = np.arange(0, self.tissue_image_mask.shape[0] + bin, bin)
        col_bins = np.arange(0, self.tissue_image_mask.shape[1] + bin, bin)
        assert len(row_bins) == up_row + 1 and len(col_bins) == up_col + 1  # check bins

        # Memory optimization: calculate bin indices
        bin_handler.gene_count_detail['row_in_bin'] = pd.cut(bin_handler.gene_count_detail['pxl_row_in_fullres'], bins=row_bins, labels=False, include_lowest=True, right=False)
        bin_handler.gene_count_detail['col_in_bin'] = pd.cut(bin_handler.gene_count_detail['pxl_col_in_fullres'], bins=col_bins, labels=False, include_lowest=True, right=False)

        # Calculate unicode values more memory-efficiently
        unicode_values = bin_handler.gene_count_detail['row_in_bin'].values * up_col + bin_handler.gene_count_detail['col_in_bin'].values

        # Memory optimization: process unicode_values_bi in-place where possible
        unicode_values_bi = np.zeros((len(unicode_values), num_bits), dtype=np.uint8)
        for i in range(num_bits):
            unicode_values_bi[:, i] = (unicode_values >> i) & 1

        # Convert to unicode barcode string
        unicode = self.base_mapping[unicode_values_bi[:, ::-1].reshape(-1, 2).dot(np.array([2, 1]))]
        unicode = unicode.reshape(-1, num_bits // 2).view('U' + str(num_bits // 2)).ravel()
        bin_handler.gene_count_detail['unicode'] = unicode

        # Free temporary arrays immediately
        del unicode_values, unicode_values_bi, unicode
        gc.collect()

        # Calculate bin center coordinates and check if they are in tissue mask
        bin_handler.gene_count_detail['bin_center_row'] = (bin_handler.gene_count_detail['row_in_bin'] + 0.5) * bin
        bin_handler.gene_count_detail['bin_center_col'] = (bin_handler.gene_count_detail['col_in_bin'] + 0.5) * bin
        bin_handler.gene_count_detail['bin_center_row'] = bin_handler.gene_count_detail['bin_center_row'].clip(0, self.tissue_image_mask.shape[0] - 1).astype(int)
        bin_handler.gene_count_detail['bin_center_col'] = bin_handler.gene_count_detail['bin_center_col'].clip(0, self.tissue_image_mask.shape[1] - 1).astype(int)

        # Mark bins that are in tissue based on mask
        bin_handler.gene_count_detail['in_tissue_bin'] = self.tissue_image_mask[
            bin_handler.gene_count_detail['bin_center_row'].values,
            bin_handler.gene_count_detail['bin_center_col'].values
        ] == 255

        # Drop temporary columns
        bin_handler.gene_count_detail = bin_handler.gene_count_detail.drop(columns=['bin_center_row', 'bin_center_col'])

    @staticmethod
    def unique_gene_count(gene_detail):
        columns_extract = ['unicode', 'geneID', 'UMI', 'count']
        if 'Barcode' in gene_detail.columns:
            columns_extract.insert(1, 'Barcode')
        unique_gene_count_detail = gene_detail[columns_extract]
        unique_gene_count_detail = unique_gene_count_detail.rename(
            columns={'unicode': 'Barcode', 'Barcode': 'OriginalBarcode'}
        )
        return unique_gene_count_detail

    @staticmethod
    def tissue_barcode_position(gene_detail):
        # Check if in_tissue_bin column exists (from bin_segment method)
        if 'in_tissue_bin' in gene_detail.columns:
            columns_extract = ['unicode', 'row_in_bin', 'col_in_bin', 'in_tissue_bin']
            count_detail = gene_detail[columns_extract]
            count_detail = count_detail.rename(columns={'unicode': 'Barcode'})
            # Aggregate in_tissue_bin: if any UMI in a bin is in tissue, mark bin as in_tissue
            count_detail = count_detail.groupby(['Barcode', 'row_in_bin', 'col_in_bin']).agg({'in_tissue_bin': 'any'}).reset_index()
            count_detail['in_tissue'] = count_detail['in_tissue_bin'].astype(int)
            count_detail = count_detail.drop(columns=['in_tissue_bin'])
            # IMPORTANT: Only keep bins that are in tissue (background removal)
            count_detail = count_detail[count_detail['in_tissue'] == 1]
            print(f"Kept {len(count_detail)} bins in tissue after background filtering")
        else:
            # Fallback to old behavior if in_tissue_bin is not available (e.g., for raw bin)
            columns_extract = ['unicode', 'row_in_bin', 'col_in_bin']
            count_detail = gene_detail[columns_extract]
            count_detail = count_detail.rename(columns={'unicode': 'Barcode'})
            count_detail.insert(1, 'in_tissue', 1)
            count_detail = count_detail.drop_duplicates(subset=['Barcode'])
        return count_detail

    @staticmethod
    def generate_bin_downsample(bin_handler, count_detail):
        """Generate downsampling data for a specific bin size.

        This creates a downsample.tsv file for the bin with proper median genes calculation
        based on square bins (not raw barcodes).

        Args:
            bin_handler: BinHandler object for the current bin size
            count_detail: DataFrame with columns (Barcode, geneID, UMI, count)
        """
        try:
            # Skip for Raw bin
            if bin_handler.micron_bin == 'Raw':
                return

            df_downsample = deterministic_downsample_metrics(count_detail)

            # Save to file
            downsample_file = os.path.join(bin_handler.bin_folder, 'downsample.tsv')
            df_downsample.to_csv(downsample_file, index=False, sep='\t')

            print(f"[Bin {bin_handler.micron_bin}] Generated downsample data: median genes {df_downsample['median_gene_number'].iloc[-1]:.0f}, saturation {df_downsample['umi_saturation'].iloc[-1]:.1f}%")

        except Exception as e:
            print(f"[Bin {bin_handler.micron_bin}] Warning: Failed to generate downsample data: {e}")

    @staticmethod
    def gene_indicator_stats(bin_handler, count_detail, bin_position, total_reads_before_filter=None):
        def num_gt2(x):
            return (x > 1).sum()

        count_detail['count'] = count_detail['count'].astype(int)
        bin_count_detail = count_detail.groupby('Barcode').agg({
            'count': ['sum', num_gt2],
            'UMI': 'count',
            'geneID': 'nunique'
        })
        bin_count_detail.columns = ['readcount', 'UMI2', 'UMI', 'geneID']
        bin_count_detail.sort_values(by='UMI', ascending=False)

        estimated_cells = len(bin_position)
        total_genes = count_detail['geneID'].nunique()

        # Calculate Fraction Reads in Square
        # Total reads in tissue squares (after filtering)
        total_reads_in_bins = bin_count_detail['readcount'].sum()
        # Total reads from ALL barcodes (before any tissue filtering)
        # Use the pre-filter total if provided, otherwise fall back to current count_detail
        if total_reads_before_filter is not None:
            total_reads_all = total_reads_before_filter
        else:
            total_reads_all = count_detail['count'].sum()

        # Calculate fraction as percentage
        if total_reads_all > 0:
            fraction_reads_in_square = round(float(total_reads_in_bins) / total_reads_all * 100, 2)
        else:
            fraction_reads_in_square = 0.0

        # Division by zero protection: if no cells, set statistics to 0
        if estimated_cells == 0 or pd.isna(estimated_cells):
            mean_reads_per_cell = 0
            median_reads_per_cell = 0
            mean_umi_per_cell = 0
            median_umi_per_cell = 0
            mean_genes_per_cell = 0
            median_genes_per_cell = 0
        else:
            mean_reads_per_cell = int(bin_count_detail['readcount'].sum() / estimated_cells)
            median_reads_per_cell = int(bin_count_detail['readcount'].median())
            mean_umi_per_cell = int(bin_count_detail['UMI'].mean())
            median_umi_per_cell = int(bin_count_detail['UMI'].median())
            mean_genes_per_cell = int(bin_count_detail['geneID'].mean())
            median_genes_per_cell = int(bin_count_detail['geneID'].median())

        with open(os.path.join(bin_handler.bin_folder, f'stat.txt'), 'w') as f:
            f.write(f"Total Genes: {total_genes: ,.0f}\n")
            f.write(f"Estimated Number of square bin: {estimated_cells: ,.0f}\n")
            f.write(f"Fraction Reads in Square: {fraction_reads_in_square}%\n")
            f.write(f"Mean Reads per square bin: {mean_reads_per_cell: ,.0f}\n")
            f.write(f"Median Reads per square bin: {median_reads_per_cell: ,.0f}\n")
            f.write(f"Mean UMI per square bin: {mean_umi_per_cell: ,.0f}\n")
            f.write(f"Median UMI per square bin: {median_umi_per_cell: ,.0f}\n")
            f.write(f"Mean Genes per square bin: {mean_genes_per_cell: ,.0f}\n")
            f.write(f"Median Genes per square bin: {median_genes_per_cell: ,.0f}\n")

        # Generate downsample data for this bin size (for accurate saturation curves)
        BinSegment.generate_bin_downsample(bin_handler, count_detail)

    def get_bin_summary(self, use_raw=False):
        # get bin statistics
        def parse_stat_file(content):
            stats = {}
            for line in content.split('\n'):
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    try:
                        value = float(value.replace(',', ''))
                        if value.is_integer():
                            value = int(value)
                    except ValueError:
                        pass
                    stats[key] = value
            return stats

        bin_statistics = {}
        if self.bin_exec is not None:
            for i, bh in self.bin_exec:
                if not use_raw and i is None:
                    continue
                if bh is not None:
                    with open(os.path.join(bh.bin_folder, f'stat.txt'), 'r') as f:
                        bin_statistics[bh.micron_bin] = parse_stat_file(f.read())

        if bin_statistics:
            table_plot = Table_plot(data=bin_statistics, index_name='Bin Size')
            table_plot.figure.write_html(os.path.join(self.outdir, 'bin_statistics.html'))
        else:
            table_plot = None
        self.add_data(table_plot=table_plot.to_dict()) if table_plot else self.add_data(table_plot=None)

    @staticmethod
    def get_gtf(genome_dir):
        gtf_file = Mkref_rna.parse_genomeDir(genome_dir)['gtf']
        gp = reference.GtfParser(gtf_file)
        gp.get_id_name()
        features = gp.get_features()
        return features

    def generate_gene_matrix(self, bin_handler, gene_count):
        matrix_dir = os.path.join(bin_handler.bin_folder, 'matrix')
        count_matrix = CountMatrix.from_dataframe(gene_count, self.features, value="UMI")
        utils.check_mkdir(dir_name=matrix_dir)
        count_matrix.get_features().to_tsv(os.path.join(matrix_dir, FEATURE_FILE_NAME))
        pd.Series(count_matrix.get_barcodes()).to_csv(os.path.join(matrix_dir, BARCODE_FILE_NAME),
                                                      index=False, sep='\t', header=False)
        matrix_path = os.path.join(matrix_dir, MATRIX_FILE_NAME)
        mmwrite(matrix_path, count_matrix.get_matrix())  # use scipy to write matrix file

    def gzip_bin_files(self, bin_handler):

        gzip_path = os.path.join(
            self.exp_tmp_path,
            f"{os.path.basename(bin_handler.bin_folder)}_{datetime.now().strftime('%Y%m%d%H%M')}.tar.gz"
        )
        with tarfile.open(gzip_path, 'w:gz') as tar:
            tar.add(bin_handler.bin_folder, arcname=os.path.basename(bin_handler.bin_folder))
    def _barcode_extend(self):
        ref = {}
        for bc in self._bc_list.to_numpy():
            bc = bc[0]  # trap here
            ref[bc[:4]] = bc
        for row, item in tqdm(self.barcodes_detail.iterrows(), desc='ex_barcode', total=self.barcodes_detail.shape[0]):
            barcode = item[2]
            cur_barcodes = ''.join([ref[barcode[4 * i: 4 * i + 4]] for i in range(len(barcode) // 4)
                                    if barcode[4 * i: 4 * i + 4] in reference])
            # only if barcode length is 18, then we can use it
            if len(cur_barcodes) == 18:
                self.barcodes_detail.at[row, 2] = cur_barcodes
        self.barcodes_detail = self.barcodes_detail[self.barcodes_detail[2].str.len() == 18]

    def test_get_bin_summary(self):
        self.bin_exec = [
            [i, SquareBin(self._micron_bin[i],
                            self.bin[i],
                            self.tissue_bbox,
                            bin_folder=os.path.join(self.exp_sq_bin,
                                                    f'{self.tissue_name}_bin{self._micron_bin[i]}')
                            )] for i in range(len(self._micron_bin))
        ]
        self.bin_exec.append(
            [None, SquareBin('Raw',
                               self._bin,
                               self.tissue_bbox,
                               bin_folder=os.path.join(self.exp_sq_bin, f'{self.tissue_name}_Raw')
                               )]
        )
        self.get_bin_summary()

    def get_tissue_image(self, tif_path, resize=True, use_adaptive_thresh=False):
        """
        Extract tissue contour from HE image

        Parameters:
        -----------
        tif_path : str
            HE image path (supports TIFF, JPEG, PNG formats)
        resize : bool
            Whether to resize (auto-calculate appropriate scale rather than fixed 1/3)
        use_adaptive_thresh : bool
            Whether to use adaptive threshold (more robust for low-contrast images)
        """
        # Support multiple image formats
        file_ext = os.path.splitext(tif_path)[1].lower()
        if file_ext in ['.tif', '.tiff']:
            # Use tifffile for TIFF format (supports 16-bit and multi-channel)
            self.he_image, _ = utils.read_tiff_with_metadata(tif_path)
        else:
            # Use OpenCV for other formats (JPEG, PNG, etc.)
            self.he_image = cv2.imread(tif_path, cv2.IMREAD_UNCHANGED)
            if self.he_image is None:
                raise FileNotFoundError(f"Failed to read HE image: {tif_path}")
            # OpenCV loads as BGR, convert to RGB to preserve original colors
            if self.he_image.ndim == 3:
                self.he_image = cv2.cvtColor(self.he_image, cv2.COLOR_BGR2RGB)

        print(f"[HE Image] Original size: {self.he_image.shape}")
        self.registration_source_signature = image_signature(tif_path)
        self.registration_raw_size_xy = [int(self.he_image.shape[1]), int(self.he_image.shape[0])]
        self.registration_resize_xy = [1.0, 1.0]
        self.image_d = self.he_image.copy()

        # Adaptive resize: determine scale based on image size and tissue_bbox size
        if self.he_image.ndim == 3:
            self.he_image = cv2.cvtColor(self.he_image, cv2.COLOR_RGB2GRAY)
            self.he_image = cv2.bitwise_not(self.he_image)

            if resize:
                # Calculate appropriate scale: make HE image approximately 0.5-1x the tissue_bbox size
                # This balances computation cost with precision
                he_h, he_w = self.he_image.shape[:2]
                target_h, target_w = self.tissue_bbox[2], self.tissue_bbox[3]
                scale_h = target_h / he_h
                scale_w = target_w / he_w
                target_scale = min(scale_h, scale_w)  # Use smaller scale

                # Only resize if HE image is much larger than tissue_bbox
                if target_scale < 0.8:
                    resize_scale = max(target_scale, 0.3)  # Minimum scale 30%, avoid too small
                    self.registration_resize_xy = [float(resize_scale), float(resize_scale)]
                    self.he_image = cv2.resize(self.he_image, None, fx=resize_scale, fy=resize_scale, interpolation=cv2.INTER_AREA)
                    self.image_d = cv2.resize(self.image_d, None, fx=resize_scale, fy=resize_scale, interpolation=cv2.INTER_AREA)
                    print(f"[HE Image] Resized to {self.he_image.shape} (scale={resize_scale:.3f})")
                else:
                    print(f"[HE Image] No resize needed (HE size ~ tissue_bbox size)")
            else:
                self.image_d = self.image_d

        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))

        # Check for a manually edited HE mask in 06.segment/mask.
        manual_he_mask_path = self._find_manual_mask('he')
        require_manual_he = getattr(getattr(self, 'args', None), 'require_manual_he_mask', False)
        if require_manual_he and not manual_he_mask_path:
            raise ValueError('DAPI ROI registration requires a manual image mask: ' + self._manual_mask_hint('he'))
        manual_mask_loaded = False

        if manual_he_mask_path:
            print(f"[Manual HE Mask] Found manual HE mask at: {manual_he_mask_path}")
            manual_he_mask = cv2.imread(manual_he_mask_path, cv2.IMREAD_GRAYSCALE)
            if manual_he_mask is not None:
                # Resize manual mask to match HE image dimensions
                manual_he_mask_resized = cv2.resize(manual_he_mask, (self.he_image.shape[1], self.he_image.shape[0]), interpolation=cv2.INTER_AREA)
                # Apply binary threshold to ensure 0/255 values
                _, manual_he_mask_binary = cv2.threshold(manual_he_mask_resized, 127, 255, cv2.THRESH_BINARY)
                filled_mask = manual_he_mask_binary
                manual_mask_loaded = True
                print(f"[Manual HE Mask] Loaded and applied manual HE mask (size: {manual_he_mask_binary.shape})")
            else:
                print(f"[Manual HE Mask] Warning: Failed to load {manual_he_mask_path}, using auto-segmentation")
                if require_manual_he:
                    raise ValueError(f'Cannot read manual DAPI image mask: {manual_he_mask_path}')
        else:
            print(f"[Manual HE Mask] No manual HE mask found, using auto-segmentation")
            print(f"[Manual HE Mask] To use manual mask: save edited mask as {self._manual_mask_hint('he')}")

        # Auto-segmentation (only if no manual mask was loaded)
        if not manual_mask_loaded:
            # Improved threshold segmentation
            if use_adaptive_thresh:
                # Adaptive threshold (more robust for low-contrast images)
                binary = cv2.adaptiveThreshold(self.he_image, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                              cv2.THRESH_BINARY, 11, 2)
                print("Using adaptive threshold for HE tissue segmentation")
            else:
                # Otsu threshold
                _, binary = cv2.threshold(self.he_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
                print("Using Otsu threshold for HE tissue segmentation")

            binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel, iterations=5)
            contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

            if contours:
                # Support multiple tissue regions (not just the largest one)
                # Calculate area threshold: keep regions >= 10% of largest region
                areas = [cv2.contourArea(c) for c in contours]
                max_area = max(areas)

                # Configurable threshold (can be adjusted via parameter)
                area_threshold_ratio = getattr(self, 'he_multi_region_threshold', 0.1)
                area_threshold = max_area * area_threshold_ratio

                # Filter contours by area threshold
                valid_contours = [c for c in contours if cv2.contourArea(c) >= area_threshold]

                # Draw all valid contours
                mask = np.zeros(self.he_image.shape[:2], dtype=np.uint8)
                cv2.drawContours(mask, valid_contours, -1, 255, -1)
                filled_mask = mask

                print(f"[HE Segmentation] Multi-region support:")
                print(f"  - Total contours found: {len(contours)}")
                print(f"  - Largest region area: {max_area:.0f} pixels")
                print(f"  - Area threshold: {area_threshold:.0f} pixels ({area_threshold_ratio:.0%} of largest)")
                print(f"  - Valid regions kept: {len(valid_contours)}")

                if len(valid_contours) > 1:
                    print(f"  - ✓ Multiple tissue regions detected and preserved")
                    for i, c in enumerate(valid_contours):
                        area = cv2.contourArea(c)
                        print(f"    Region {i+1}: {area:.0f} pixels ({area/max_area:.1%} of largest)")
            else:
                filled_mask = np.zeros(self.he_image.shape[:2], dtype=np.uint8)
                print("WARNING: No contours found in HE image! Registration may fail.")

        self._set_background_image_mask(filled_mask)

        # Save mask before registration for visualization comparison
        self._set_background_image_mask_before_reg(filled_mask.copy())

        # Multi-stage registration: align he_image_mask to gem_mask
        # Stage 1: Coarse alignment using contour features
        print("[Stage 1] Coarse registration using contour features...")
        trans_matrix_coarse = self.contours_registration(self.gem_mask, self.he_image_mask)

        # Calculate Dice score after coarse registration (for display)
        he_mask_after_coarse = cv2.warpAffine(self.he_image_mask, trans_matrix_coarse,
                                             (self.tissue_bbox[3], self.tissue_bbox[2]))
        coarse_intersection = np.sum((self.gem_mask > 0) & (he_mask_after_coarse > 0))
        coarse_dice = 2.0 * coarse_intersection / (np.sum(self.gem_mask > 0) + np.sum(he_mask_after_coarse > 0)) if (np.sum(self.gem_mask > 0) + np.sum(he_mask_after_coarse > 0)) > 0 else 0

        # Apply coarse registration transform
        he_image_coarse = cv2.warpAffine(self.he_image, trans_matrix_coarse,
                                         (self.tissue_bbox[3], self.tissue_bbox[2]))
        he_mask_coarse = cv2.warpAffine(self.he_image_mask, trans_matrix_coarse,
                                        (self.tissue_bbox[3], self.tissue_bbox[2]))
        image_d_coarse = cv2.warpAffine(self.image_d, trans_matrix_coarse,
                                        (self.tissue_bbox[3], self.tissue_bbox[2]))

        # Stage 2: Fine registration options
        # Priority: feature-based > SimpleITK > coarse only
        use_feature_refinement = getattr(self, 'use_feature_refinement', False)
        use_sitk_refinement = getattr(self, 'use_sitk_refinement', False)

        if use_feature_refinement:
            # Feature-based fine registration
            print("[Stage 2] Fine registration using feature matching (SIFT/ORB)...")

            # Prepare GEM heatmap for feature extraction
            gem_heatmap_for_features = self.gem_enhanced if hasattr(self, 'gem_enhanced') else self.gem

            trans_matrix_fine, inlier_count = self.feature_based_registration(
                image_fixed=self.gem_mask,
                image_moving=self.he_image_mask,
                he_image_gray=self.he_image,
                gem_heatmap=gem_heatmap_for_features,
                coarse_transform=None,  # Don't pass coarse, already applied
                use_roi=True,  # Enable ROI cropping to handle partial sequencing
                verbose=True
            )

            # Apply fine registration transform
            self.he_image = cv2.warpAffine(he_image_coarse, trans_matrix_fine,
                                          (self.tissue_bbox[3], self.tissue_bbox[2]))
            he_mask_refined = cv2.warpAffine(he_mask_coarse, trans_matrix_fine,
                                             (self.tissue_bbox[3], self.tissue_bbox[2]))
            self._set_background_image_mask(he_mask_refined)
            self.image_d = cv2.warpAffine(image_d_coarse, trans_matrix_fine,
                                         (self.tissue_bbox[3], self.tissue_bbox[2]))

            # Combine two transformation matrices
            coarse_3x3 = np.vstack([trans_matrix_coarse, [0, 0, 1]])
            fine_3x3 = np.vstack([trans_matrix_fine, [0, 0, 1]])
            total_3x3 = fine_3x3 @ coarse_3x3
            trans_matrix = total_3x3[:2, :]
            print(f"[Registration] Two-stage registration completed (contour + features, {inlier_count} inliers)")

        elif use_sitk_refinement:
            # SimpleITK fine registration - continue optimizing from coarse result
            print("[Stage 2] Fine registration using SimpleITK optimization...")
            print(f"  Starting from coarse registration result (Dice={coarse_dice:.5f})")
            registration_type = getattr(self, 'registration_type', 'similarity')

            # IMPORTANT: Pass coarse matrix as initial value, let SimpleITK optimize from this good result
            trans_matrix_refined = self.register_mask_images_sitk(
                self.gem_mask,
                self.he_image_mask,  # Use original HE mask
                registration_type=registration_type,
                verbose=True,
                use_contour_init=False,  # Don't use internal contour initialization
                initial_transform_matrix=trans_matrix_coarse  # Pass coarse result as starting point!
            )

            # Apply refined transform directly (already includes coarse registration)
            self.he_image = cv2.warpAffine(self.he_image, trans_matrix_refined,
                                          (self.tissue_bbox[3], self.tissue_bbox[2]))
            he_mask_refined = cv2.warpAffine(self.he_image_mask, trans_matrix_refined,
                                             (self.tissue_bbox[3], self.tissue_bbox[2]))
            self._set_background_image_mask(he_mask_refined)
            self.image_d = cv2.warpAffine(self.image_d, trans_matrix_refined,
                                         (self.tissue_bbox[3], self.tissue_bbox[2]))

            trans_matrix = trans_matrix_refined
            print(f"[Registration] Two-stage registration completed (contour + SimpleITK)")
        else:
            # Use coarse registration result only
            self.he_image = he_image_coarse
            self._set_background_image_mask(he_mask_coarse)
            self.image_d = image_d_coarse
            trans_matrix = trans_matrix_coarse
            print(f"[Registration] Coarse registration completed (refinement disabled)")

        # Explicit fluorescence rendering preserves black outside the tissue ROI.
        background_mask = self._get_background_image_mask()
        self.image_d[background_mask == 0] = self._registration_background_value()
        self.tissue_image = self.image_d

        # Save transformation matrix for visualization
        self.registration_transform = trans_matrix

    def _registration_background_value(self):
        return 0 if getattr(self, 'fluorescence_background', False) else 255

    @utils.add_log
    def save_register_image(self, is_correct=False):
        # Determine which image to use based on mode
        if self.mode == 'HE':
            # HE mode: always use registered HE image (self.tissue_image = self.image_d)
            scanpy_image = self.tissue_image.copy()
            print(f"[HE Mode] Using registered HE image for tissue visualization (RGB, shape={scanpy_image.shape})")
        elif self.mode == 'gene_expr' and hasattr(self, 'gem_enhanced') and self.gem_enhanced is not None:
            # gene_expr mode: use enhanced gem image
            scanpy_image = self.gem_enhanced.copy()
            print(f"[gene_expr Mode] Using enhanced gene expression image for tissue visualization (grayscale, shape={scanpy_image.shape})")
        else:
            # ssDNA mode or gene_expr without enhancement: use original tissue_image
            scanpy_image = self.tissue_image.copy()
            print(f"[{self.mode} Mode] Using tissue_image for visualization (shape={scanpy_image.shape})")

        # shift = self.get_offset(scanpy_image)
        # transM = np.float32([[1, 0, -shift[0]], [0, 1, -shift[1]]])
        # scanpy_image = cv2.warpAffine(scanpy_image, transM, (scanpy_image.shape[1], scanpy_image.shape[0]))
        scanpy_image_600 = self._resize_longest_side(scanpy_image, 600, interpolation=cv2.INTER_AREA)
        scanpy_image_2000 = self._resize_longest_side(scanpy_image, 2000, interpolation=cv2.INTER_LANCZOS4)

        if is_correct:
            fine_tune = np.float32([[1, 0, -self.bin[-1] // 2], [0, 1, -self.bin[-1] // 2]])
            scanpy_image_600 = cv2.warpAffine(
                scanpy_image_600, fine_tune,
                (scanpy_image_600.shape[1], scanpy_image_600.shape[0])
            )
            scanpy_image_2000 = cv2.warpAffine(
                scanpy_image_2000, fine_tune,
                (scanpy_image_2000.shape[1], scanpy_image_2000.shape[0])
            )

        if scanpy_image.ndim == 3:
            scanpy_image_600 = cv2.cvtColor(scanpy_image_600, cv2.COLOR_RGB2BGR)
            scanpy_image_2000 = cv2.cvtColor(scanpy_image_2000, cv2.COLOR_RGB2BGR)

        registered_path = os.path.join(self.exp_bin_im, f'{self.tissue_name}_regist.tif')
        tifffile.imwrite(registered_path, scanpy_image)
        registered_stat = os.stat(registered_path)
        metadata = {
            'method': self.mode,
            'sample': self.tissue_name,
            'fluorescence_background': bool(getattr(self, 'fluorescence_background', False)),
            'shape': list(scanpy_image.shape),
            'image_bytes': registered_stat.st_size,
            'image_mtime_ns': registered_stat.st_mtime_ns,
        }
        if self.mode == 'HE' and hasattr(self, 'registration_source_signature'):
            if getattr(self, 'rectify', False) or is_correct:
                metadata['coordinate_mapping_status'] = 'unavailable: rectify adds unsaved crop/rotation transforms'
            else:
                metadata.update(registration_manifest(
                    self.registration_transform, self.registration_resize_xy,
                    self.registration_raw_size_xy,
                    [int(scanpy_image.shape[1]), int(scanpy_image.shape[0])],
                    self.registration_source_signature, registered_path,
                ))
        with open(registered_path + '.registration.json', 'w', encoding='utf-8') as handle:
            json.dump(metadata, handle, indent=2, allow_nan=False)
        cv2.imwrite(os.path.join(self.exp_bin_im, f'tissue_lowres_image.png'), scanpy_image_600)
        cv2.imwrite(os.path.join(self.exp_bin_im, f'tissue_hires_image.png'), scanpy_image_2000)
        print(
            "[Spatial Image] Saved aspect-preserving previews: "
            f"lowres={scanpy_image_600.shape[1]}x{scanpy_image_600.shape[0]}, "
            f"hires={scanpy_image_2000.shape[1]}x{scanpy_image_2000.shape[0]}"
        )

    @utils.add_log
    def save_gene_expr_enhancement_comparison(self):
        """
        Save comparison plot of gene expression before and after enhancement to evaluate optimization effect
        """
        # DISABLED: Debug comparison plot not needed for customers
        return

        if not hasattr(self, 'gem_original') or not hasattr(self, 'gem_enhanced'):
            print("No enhancement data available for comparison")
            return

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        # Original gene expression heatmap (grayscale)
        im1 = axes[0, 0].imshow(self.gem_original, cmap='gray', interpolation='nearest')
        axes[0, 0].set_title('Original Gene Expression')
        axes[0, 0].axis('off')
        plt.colorbar(im1, ax=axes[0, 0], fraction=0.046)

        # Enhanced gene expression heatmap (grayscale)
        im2 = axes[0, 1].imshow(self.gem_enhanced, cmap='gray', interpolation='nearest')
        axes[0, 1].set_title('Enhanced Gene Expression')
        axes[0, 1].axis('off')
        plt.colorbar(im2, ax=axes[0, 1], fraction=0.046)

        # Binary mask
        im3 = axes[0, 2].imshow(self.gem_mask, cmap='gray', interpolation='nearest')
        axes[0, 2].set_title('Binary Mask (After Enhancement)')
        axes[0, 2].axis('off')
        plt.colorbar(im3, ax=axes[0, 2], fraction=0.046)

        # Original data histogram
        valid_original = self.gem_original[self.gem_original > 0].flatten()
        axes[1, 0].hist(valid_original, bins=100, color='blue', alpha=0.7, edgecolor='black')
        axes[1, 0].set_title('Histogram: Original')
        axes[1, 0].set_xlabel('Intensity')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_yscale('log')

        # Enhanced data histogram
        valid_enhanced = self.gem_enhanced[self.gem_enhanced > 0].flatten()
        axes[1, 1].hist(valid_enhanced, bins=100, color='green', alpha=0.7, edgecolor='black')
        axes[1, 1].set_title('Histogram: Enhanced')
        axes[1, 1].set_xlabel('Intensity')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].set_yscale('log')

        # Overlay comparison (grayscale)
        overlay = cv2.addWeighted(self.gem_original, 0.5, self.gem_enhanced, 0.5, 0)
        im6 = axes[1, 2].imshow(overlay, cmap='gray', interpolation='nearest')
        axes[1, 2].set_title('Overlay Comparison (50/50)')
        axes[1, 2].axis('off')
        plt.colorbar(im6, ax=axes[1, 2], fraction=0.046)

        plt.tight_layout()

        # Save comparison plot
        comparison_path = os.path.join(self.exp_bin_im, 'gene_expr_enhancement_comparison.png')
        plt.savefig(comparison_path, dpi=150, bbox_inches='tight')

        plt.close(fig)  # Explicitly close figure after all saves
        print(f"Enhancement comparison saved to: {comparison_path}")

    @utils.add_log
    def save_he_registration_comparison(self):
        """Backward-compatible wrapper for background registration visualization."""
        self.save_background_registration_comparison()

    @utils.add_log
    def save_background_registration_comparison(self):
        """
        Save background-image and gene_expr registration comparison plot to evaluate registration quality
        """
        # DISABLED: Debug comparison plot not needed for customers
        # Still call overlay generation to produce report assets.
        self.save_gem_background_image_overlay()
        return

        background_mask = self._get_background_image_mask()
        if background_mask is None:
            print("No background image registration data available")
            return

        # Check if gem_mask is available
        if not hasattr(self, 'gem_mask') or self.gem_mask is None:
            print("No gem_mask available for comparison")
            return

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        # First row: gem_mask, background mask before registration, background mask after registration
        axes[0, 0].imshow(self.gem_mask, cmap='gray')
        axes[0, 0].set_title('GEM Mask (Reference)')
        axes[0, 0].axis('off')

        # Background mask before registration (need to save original)
        background_mask_before_reg = self._get_background_image_mask_before_reg()
        if background_mask_before_reg is not None:
            axes[0, 1].imshow(background_mask_before_reg, cmap='gray')
            axes[0, 1].set_title('Background Mask (Before Registration)')
        else:
            axes[0, 1].text(0.5, 0.5, 'Not Available', ha='center', va='center')
            axes[0, 1].set_title('Background Mask (Before Registration)')
        axes[0, 1].axis('off')

        # Background mask after registration
        axes[0, 2].imshow(background_mask, cmap='gray')
        axes[0, 2].set_title('Background Mask (After Registration)')
        axes[0, 2].axis('off')

        # Second row: overlay comparison, Dice coefficient, contour overlay
        # 1. Green-red overlay (gem=green, background=red)
        overlay_rgb = np.zeros((*self.gem_mask.shape, 3), dtype=np.uint8)
        overlay_rgb[:, :, 1] = self.gem_mask  # gem_mask -> green channel
        overlay_rgb[:, :, 0] = background_mask  # background mask -> red channel
        axes[1, 0].imshow(overlay_rgb)
        axes[1, 0].set_title('Overlay (GEM=Green, Background=Red)')
        axes[1, 0].axis('off')

        # 2. Dice coefficient and IoU
        intersection = np.logical_and(self.gem_mask > 0, background_mask > 0).sum()
        union = np.logical_or(self.gem_mask > 0, background_mask > 0).sum()
        gem_area = np.sum(self.gem_mask > 0)
        he_area = np.sum(background_mask > 0)

        dice = 2 * intersection / (gem_area + he_area) if (gem_area + he_area) > 0 else 0
        iou = intersection / union if union > 0 else 0

        # Calculate contour alignment - use contour center distance rather than edge pixel overlap
        # This avoids the impact of edge style differences
        contour_alignment_score = 0

        # Extract contours from both masks
        gem_contours, _ = cv2.findContours(self.gem_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        he_contours, _ = cv2.findContours(background_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        if len(gem_contours) > 0 and len(he_contours) > 0:
            # Take largest contour
            gem_contour = max(gem_contours, key=cv2.contourArea)
            he_contour = max(he_contours, key=cv2.contourArea)

            # Method 1: Calculate average distance of contour points (simplified Hausdorff distance)
            # Sample points from gem contour, calculate nearest distance to he contour
            gem_pts = gem_contour.reshape(-1, 2)
            he_pts = he_contour.reshape(-1, 2)

            # Downsample to speed up calculation
            sample_rate = max(1, len(gem_pts) // 1000)
            gem_pts_sampled = gem_pts[::sample_rate]
            he_pts_sampled = he_pts[::sample_rate]

            # Calculate bidirectional average nearest distance
            from scipy.spatial.distance import cdist
            if len(gem_pts_sampled) > 0 and len(he_pts_sampled) > 0:
                # GEM to HE distance
                dist_gem_to_he = cdist(gem_pts_sampled, he_pts_sampled, metric='euclidean')
                avg_dist_gem_to_he = np.mean(np.min(dist_gem_to_he, axis=1))

                # HE to GEM distance
                dist_he_to_gem = cdist(he_pts_sampled, gem_pts_sampled, metric='euclidean')
                avg_dist_he_to_gem = np.mean(np.min(dist_he_to_gem, axis=1))

                # Bidirectional average distance
                avg_contour_distance = (avg_dist_gem_to_he + avg_dist_he_to_gem) / 2

                # Convert to alignment score (smaller distance = higher score)
                # Assuming image size as reference, distance <50 pixels is excellent
                image_size = max(self.gem_mask.shape)
                contour_alignment_score = max(0, 1 - (avg_contour_distance / (image_size * 0.05)))

                print(f"[Contour Alignment] Avg distance: {avg_contour_distance:.1f} pixels, Score: {contour_alignment_score:.3f}")
                edge_dice_label = f'Contour Alignment: {contour_alignment_score:.3f} (0=bad, 1=perfect)'
            else:
                edge_dice_label = 'Contour Alignment: N/A'
        else:
            edge_dice_label = 'Contour Alignment: N/A'

        edge_dice = contour_alignment_score

        axes[1, 1].text(0.1, 0.88, f'Dice Score: {dice:.5f}', fontsize=14, transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.1, 0.78, f'IoU: {iou:.5f}', fontsize=14, transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.1, 0.68, edge_dice_label, fontsize=11, transform=axes[1, 1].transAxes, color='blue')
        axes[1, 1].text(0.1, 0.58, f'Overlap: {intersection} pixels', fontsize=10, transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.05, 0.45, 'Contour Alignment measures', fontsize=8, transform=axes[1, 1].transAxes, style='italic')
        axes[1, 1].text(0.05, 0.38, 'avg distance between contours', fontsize=8, transform=axes[1, 1].transAxes, style='italic')
        axes[1, 1].text(0.05, 0.28, 'Good: >0.8, Excellent: >0.9', fontsize=8, transform=axes[1, 1].transAxes, style='italic', color='gray')
        axes[1, 1].set_title('Registration Quality')
        axes[1, 1].axis('off')

        # 3. Contour overlay (draw contours on original image)
        if hasattr(self, 'tissue_image') and self.tissue_image is not None:
            # Use tissue_image as background
            if self.tissue_image.ndim == 2:
                contour_img = cv2.cvtColor(self.tissue_image, cv2.COLOR_GRAY2BGR)
            else:
                contour_img = self.tissue_image.copy()

            # Extract contours
            gem_contours, _ = cv2.findContours(self.gem_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            he_contours, _ = cv2.findContours(background_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

            # Draw contours: gem=green, background=red
            cv2.drawContours(contour_img, gem_contours, -1, (0, 255, 0), 3)  # green
            cv2.drawContours(contour_img, he_contours, -1, (255, 0, 0), 3)  # red

            axes[1, 2].imshow(cv2.cvtColor(contour_img, cv2.COLOR_BGR2RGB))
            axes[1, 2].set_title('Contour Overlay (GEM=Green, Background=Red)')
        else:
            axes[1, 2].text(0.5, 0.5, 'Tissue image not available', ha='center', va='center')
            axes[1, 2].set_title('Contour Overlay')
        axes[1, 2].axis('off')

        plt.tight_layout()

        # Save comparison plot
        comparison_path = os.path.join(self.exp_bin_im, 'he_registration_comparison.png')
        plt.savefig(comparison_path, dpi=150, bbox_inches='tight')

        plt.close(fig)  # Explicitly close figure after all saves
        print(f"Background registration comparison saved to: {comparison_path}")
        print(f"Registration quality: Dice={dice:.5f}, IoU={iou:.5f}, Contour-Alignment={edge_dice:.3f} (visual quality)")

        # Additional save: gem enhanced image and background image overlay comparison - needed for Image Alignment tab
        self.save_gem_background_image_overlay()  # ENABLED: Generate overlays folder with background alignment images

    @utils.add_log
    def save_gem_he_image_overlay(self):
        """Backward-compatible wrapper for background overlay generation."""
        self.save_gem_background_image_overlay()

    @utils.add_log
    def save_gem_background_image_overlay(self):
        """Save report overlays without constructing full-resolution RGB heatmaps.

        Registered slides are commonly larger than 12k x 12k.  Matplotlib's
        colormaps expand those images to float64 RGBA arrays, so one temporary
        can exceed 5 GB.  The customer-facing files are capped at 4000 pixels
        anyway; resize the source planes first and do all colour work on that
        preview canvas.
        """
        if not hasattr(self, 'gem_enhanced') or not hasattr(self, 'image_d'):
            print("Warning: gem_enhanced or image_d not available, skipping image overlay")
            return

        gem_img = self.gem_enhanced if self.gem_enhanced is not None else self.gem
        he_img = self.image_d if hasattr(self, 'image_d') else self.he_image

        try:
            max_dimension = int(os.environ.get('CELATLAS_OVERLAY_MAX_DIMENSION', '4000'))
            if max_dimension <= 0:
                raise ValueError
        except ValueError:
            max_dimension = 4000
            print("[Image Compression] Invalid CELATLAS_OVERLAY_MAX_DIMENSION; using 4000")

        source_h, source_w = he_img.shape[:2]
        scale = min(1.0, max_dimension / max(source_h, source_w))
        preview_size = (
            max(1, int(round(source_w * scale))),
            max(1, int(round(source_h * scale))),
        )
        if scale < 1.0:
            print(
                f"[Image Compression] Resizing sources from {source_w}x{source_h} "
                f"to {preview_size[0]}x{preview_size[1]} before colour rendering "
                f"(scale={scale:.3f})"
            )

        def preview_image(image, interpolation=cv2.INTER_AREA):
            if image.shape[:2] != (preview_size[1], preview_size[0]):
                return cv2.resize(image, preview_size, interpolation=interpolation)
            return image

        def preview_mask(mask):
            if mask is None:
                return None
            mask_uint8 = np.asarray(mask, dtype=np.uint8)
            resized = preview_image(mask_uint8, interpolation=cv2.INTER_NEAREST)
            return (resized > 0).astype(np.uint8) * 255

        gem_preview = preview_image(gem_img)
        if gem_preview.ndim == 3:
            gem_preview = cv2.cvtColor(gem_preview, cv2.COLOR_RGB2GRAY)
        gem_normalized = cv2.normalize(
            gem_preview, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U
        )
        del gem_preview

        he_preview = preview_image(he_img)
        if he_preview.ndim == 2:
            he_preview = cv2.cvtColor(he_preview, cv2.COLOR_GRAY2RGB)
        elif he_preview.shape[2] == 4:
            he_preview = cv2.cvtColor(he_preview, cv2.COLOR_RGBA2RGB)
        if he_preview.dtype == np.uint8:
            he_normalized = he_preview
        else:
            he_normalized = cv2.normalize(
                he_preview, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U
            )

        background_mask = preview_mask(self._get_background_image_mask())
        gem_mask = preview_mask(getattr(self, 'gem_mask', None))
        final_mask = preview_mask(getattr(self, 'tissue_image_mask', None))
        self.save_individual_overlay_images(
            gem_normalized,
            he_normalized,
            gem_mask=gem_mask,
            background_mask=background_mask,
            final_mask=final_mask,
        )

    def save_individual_overlay_images(
        self,
        gem_normalized,
        he_normalized,
        *,
        gem_mask=None,
        background_mask=None,
        final_mask=None,
    ):
        """
        Save separate comparison image files (supports three modes)

        File naming depends on self.mode:
        - HE mode: 2_he_registered.jpg, 5_gem_heatmap_on_he.png
        - ssDNA mode: 2_ssdna_registered.jpg, 5_gem_heatmap_on_ssdna.png
        - gene_expr mode: Only save GEM expression image, no background
        """
        images_dir = os.path.join(self.outdir, 'images', 'overlays')
        os.makedirs(images_dir, exist_ok=True)

        # Define compression parameters - to reduce image file size
        # PNG compression level: 0-9, higher value = higher compression
        png_compression = [cv2.IMWRITE_PNG_COMPRESSION, 9]
        # JPEG quality: 0-100, higher value = better quality but larger file
        jpeg_quality = [cv2.IMWRITE_JPEG_QUALITY, 85]

        # ========================================
        # 1. GEM expression heatmap (needed for all modes)
        # ========================================
        gem_expression_path = os.path.join(images_dir, '1_gem_expression.png')
        cv2.imwrite(gem_expression_path,
                   cv2.applyColorMap(gem_normalized, cv2.COLORMAP_VIRIDIS), png_compression)
        self._sync_gem_expression_to_runtime_mask_dir(gem_expression_path)
        print(f"[{self.mode} Mode] 1_gem_expression.png saved")

        # OpenCV colormaps stay uint8.  Convert once to RGB because image_d is RGB.
        gem_heatmap = cv2.cvtColor(
            cv2.applyColorMap(gem_normalized, cv2.COLORMAP_JET),
            cv2.COLOR_BGR2RGB,
        )

        # ========================================
        # 2. Background image (filename varies by mode)
        # ========================================
        if self.mode == 'HE':
            # HE mode: save HE image
            bg_filename = '2_he_registered.jpg'
            bg_label = 'HE'
        elif self.mode == 'ssDNA':
            # ssDNA mode: save ssDNA image
            bg_filename = '2_ssdna_registered.jpg'
            bg_label = 'ssDNA'
        else:  # gene_expr mode
            # gene_expr mode: no background image
            bg_filename = None
            bg_label = None

        if bg_filename:
            background_value = self._registration_background_value()
            background_name = 'black' if background_value == 0 else 'white'
            he_with_bg = he_normalized.copy()
            if background_mask is not None:
                he_with_bg[background_mask == 0] = background_value

            # Apply the fluorescence option consistently to tissue overlays.
            if getattr(self, 'fluorescence_background', False):
                he_normalized = he_with_bg

            cv2.imwrite(os.path.join(images_dir, bg_filename),
                       cv2.cvtColor(he_with_bg, cv2.COLOR_RGB2BGR), jpeg_quality)
            print(f"[{self.mode} Mode] {bg_filename} saved as JPEG with quality={jpeg_quality[1]} ({background_name} background)")

        # ========================================
        # 3. Tissue segmentation mask (mode isolation)
        # ========================================
        if self.mode in ['HE', 'ssDNA']:
            if gem_mask is not None and background_mask is not None:
                union_mask = cv2.bitwise_or(gem_mask, background_mask)

                # Save pure mask (compressed)
                cv2.imwrite(os.path.join(images_dir, '3_tissue_segmentation_mask.png'), union_mask, png_compression)
                print(f"[{self.mode} Mode] 3_tissue_segmentation_mask.png saved (union of GEM and {bg_label})")

                # ========================================
                # 3a. Pure HE/ssDNA mask (separate from GEM mask)
                # ========================================
                cv2.imwrite(os.path.join(images_dir, f'3a_{bg_label.lower()}_mask.png'), background_mask, png_compression)
                print(f"[{self.mode} Mode] 3a_{bg_label.lower()}_mask.png saved (pure {bg_label} mask only)")

                # ========================================
                # 3b. GEM Heatmap overlay on background (mask region only)
                # NOTE: This image is used in the Interactive Image Alignment Viewer
                # ========================================
                he_with_gem_overlay = np.full_like(he_normalized, self._registration_background_value())
                mask_bool = union_mask > 0
                if mask_bool.any():
                    blended = cv2.addWeighted(he_normalized, 0.7, gem_heatmap, 0.3, 0)
                    he_with_gem_overlay[mask_bool] = blended[mask_bool]
                    del blended
                cv2.imwrite(os.path.join(images_dir, '3b_tissue_segmentation_gem_heatmap.png'),
                           cv2.cvtColor(he_with_gem_overlay, cv2.COLOR_RGB2BGR), png_compression)
                del he_with_gem_overlay, mask_bool
                print(f"[{self.mode} Mode] 3b_tissue_segmentation_gem_heatmap.png saved (30% GEM on {bg_label}, {background_name} background) - Used in Interactive Viewer")

                # ========================================
                # 3c. Pure GEM Heatmap (mask region only, no background)
                # ========================================
                gem_mask_bool = gem_mask > 0
                if gem_mask_bool.any():
                    mask_values = gem_normalized[gem_mask_bool]
                    gem_normalized_pure = np.zeros_like(gem_normalized)
                    value_min = mask_values.min()
                    value_max = mask_values.max()
                    if value_max > value_min:
                        gem_normalized_pure[gem_mask_bool] = (
                            (mask_values.astype(np.float32) - value_min)
                            * (255.0 / (value_max - value_min))
                        ).astype(np.uint8)
                    else:
                        gem_normalized_pure[gem_mask_bool] = mask_values
                    gem_heatmap_pure_rgb = cv2.cvtColor(
                        cv2.applyColorMap(gem_normalized_pure, cv2.COLORMAP_JET),
                        cv2.COLOR_BGR2RGB,
                    )
                    gem_heatmap_pure = np.full_like(gem_heatmap, 255)
                    gem_heatmap_pure[gem_mask_bool] = gem_heatmap_pure_rgb[gem_mask_bool]
                    del gem_normalized_pure, gem_heatmap_pure_rgb, mask_values
                else:
                    gem_heatmap_pure = np.full_like(gem_heatmap, 255)
                cv2.imwrite(os.path.join(images_dir, '3c_gem_heatmap_only.png'),
                           cv2.cvtColor(gem_heatmap_pure, cv2.COLOR_RGB2BGR), png_compression)
                del gem_heatmap_pure, gem_mask_bool
                print(f"[{self.mode} Mode] 3c_gem_heatmap_only.png saved (pure GEM expression heatmap, white background)")

                # ========================================
                # 4. Segmentation region visualization (filename varies by mode)
                # ========================================
                if self.mode == 'HE':
                    seg_vis_filename = '4_tissue_segmentation_on_he.jpg'
                else:  # ssDNA
                    seg_vis_filename = '4_tissue_segmentation_on_ssdna.jpg'

                he_with_seg = he_normalized.copy()
                green_overlay = np.zeros_like(he_with_seg)
                green_overlay[:, :, 1] = 180
                tinted = cv2.addWeighted(he_with_seg, 0.7, green_overlay, 0.3, 0)
                union_bool = union_mask > 0
                he_with_seg[union_bool] = tinted[union_bool]
                del green_overlay, tinted, union_bool
                contours, _ = cv2.findContours(union_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(he_with_seg, contours, -1, (0, 255, 0), 2)

                cv2.imwrite(os.path.join(images_dir, seg_vis_filename),
                           cv2.cvtColor(he_with_seg, cv2.COLOR_RGB2BGR), jpeg_quality)
                del he_with_seg, contours
                print(f"[{self.mode} Mode] {seg_vis_filename} saved")

                # ========================================
                # 5. GEM Heatmap overlay (full image, filename varies by mode)
                # ========================================
                if self.mode == 'HE':
                    overlay_filename = '5_gem_heatmap_on_he.png'
                else:  # ssDNA
                    overlay_filename = '5_gem_heatmap_on_ssdna.png'

                if final_mask is None:
                    final_mask = union_mask
                expression_mask = (
                    (final_mask > 0) & (gem_normalized > 0)
                )
                if background_mask is not None:
                    expression_mask &= background_mask > 0
                expression_mask = expression_mask.astype(np.uint8) * 255
                masked_overlay = blend_heatmap_with_mask(
                    he_normalized,
                    gem_heatmap,
                    expression_mask,
                    background_weight=0.6,
                )
                cv2.imwrite(os.path.join(images_dir, overlay_filename),
                           cv2.cvtColor(masked_overlay, cv2.COLOR_RGB2BGR), png_compression)
                del masked_overlay, expression_mask
                print(
                    f"[{self.mode} Mode] {overlay_filename} saved "
                    "(40% nonzero GEM inside final tissue mask)"
                )

                full_canvas_filename = overlay_filename.replace('.png', '.full_canvas_qc.png')
                overlay_alpha = cv2.addWeighted(he_normalized, 0.6, gem_heatmap, 0.4, 0)
                cv2.imwrite(os.path.join(images_dir, full_canvas_filename),
                           cv2.cvtColor(overlay_alpha, cv2.COLOR_RGB2BGR), png_compression)
                del overlay_alpha
                print(f"[{self.mode} Mode] {full_canvas_filename} saved (unfiltered full-canvas QC)")

                # Calculate statistics
                source_gem_mask = getattr(self, 'gem_mask', None)
                source_background_mask = self._get_background_image_mask()
                gem_area = np.count_nonzero(source_gem_mask)
                bg_area = np.count_nonzero(source_background_mask)
                overlap_area = 0
                # Count the full-resolution intersection in row chunks so the
                # statistics stay exact without another slide-sized boolean array.
                for row_start in range(0, source_gem_mask.shape[0], 512):
                    row_stop = min(row_start + 512, source_gem_mask.shape[0])
                    overlap_area += np.count_nonzero(
                        (source_gem_mask[row_start:row_stop] > 0) &
                        (source_background_mask[row_start:row_stop] > 0)
                    )
                union_area = gem_area + bg_area - overlap_area

                print(f"[Tissue Segmentation] GEM mask area: {gem_area} pixels")
                print(f"[Tissue Segmentation] {bg_label} mask area: {bg_area} pixels")
                print(f"[Tissue Segmentation] Union area: {union_area} pixels")
                print(f"[Tissue Segmentation] Coverage: {union_area} pixels total from GEM + {bg_label}")

        elif self.mode == 'gene_expr':
            if gem_mask is not None:
                cv2.imwrite(os.path.join(images_dir, '3_tissue_segmentation_mask.png'), gem_mask, png_compression)
                print(f"[{self.mode} Mode] 3_tissue_segmentation_mask.png saved (GEM only)")

                gem_area = np.count_nonzero(self.gem_mask)
                print(f"[Tissue Segmentation] GEM mask area: {gem_area} pixels")

        # ========================================
        # Print summary information (by mode)
        # ========================================
        print(f"\n[{self.mode} Mode] Individual overlay images saved to: {images_dir}/")
        print("  - 1_gem_expression.png: GEM heatmap")

        if self.mode == 'HE':
            print("  - 2_he_registered.jpg: HE staining image (JPEG compressed, original colors)")
            print("  - 3_tissue_segmentation_mask.png: Tissue segmentation mask (GEM∪HE union)")
            print("  - 3b_tissue_segmentation_gem_heatmap.png: GEM heatmap overlaid on HE (30% transparency, mask region only)")
            print("  - 3c_gem_heatmap_only.png: Pure GEM heatmap (mask region only, transparent background)")
            print("  - 4_tissue_segmentation_on_he.jpg: Segmentation region overlaid on HE")
            print("  - 5_gem_heatmap_on_he.png: Nonzero GEM inside final tissue mask (40% transparency)")
            print("  - 5_gem_heatmap_on_he.full_canvas_qc.png: Unfiltered full-canvas GEM overlay for QC")
        elif self.mode == 'ssDNA':
            print("  - 2_ssdna_registered.jpg: ssDNA image (JPEG compressed)")
            print("  - 3_tissue_segmentation_mask.png: Tissue segmentation mask (GEM∪ssDNA union)")
            print("  - 3b_tissue_segmentation_gem_heatmap.png: GEM heatmap overlaid on ssDNA (30% transparency, mask region only)")
            print("  - 3c_gem_heatmap_only.png: Pure GEM heatmap (mask region only, transparent background)")
            print("  - 4_tissue_segmentation_on_ssdna.jpg: Segmentation region overlaid on ssDNA")
            print("  - 5_gem_heatmap_on_ssdna.png: Nonzero GEM inside final tissue mask (40% transparency)")
            print("  - 5_gem_heatmap_on_ssdna.full_canvas_qc.png: Unfiltered full-canvas GEM overlay for QC")
        else:  # gene_expr
            print("  - 3_tissue_segmentation_mask.png: Tissue segmentation mask (GEM only)")
            print("  [Note] gene_expr mode: no background image, only gene expression data")

        print("\n[Downstream Analysis] These images can be used for:")
        print("  1. Scanpy/Squidpy spatial analysis: Load with tissue_positions.csv")
        print("  2. Cell type deconvolution: Overlay cell types on expression maps")
        print("  3. Spatial domain detection: Visualize clusters")
        print("  4. Gene expression patterns: Analyze spatial distributions")
        if self.mode in ['HE', 'ssDNA']:
            print(f"  5. Interactive viewers: Overlay GEM on {bg_label} background")

    @utils.add_log
    def clean_redundant_folders(self):
        """Clean redundant folders: matrix, images, tmp."""
        redundant_folders = [
            os.path.join(self.exp_bs_path, 'matrix'),
            os.path.join(self.exp_bs_path, 'tmp')
        ]

        for folder in redundant_folders:
            if os.path.exists(folder):
                shutil.rmtree(folder)
                self.clean_redundant_folders.logger.info(f"Deleted redundant folder: {folder}")
            else:
                self.clean_redundant_folders.logger.info(f"Folder does not exist: {folder}")


    @utils.add_log
    def clean_specific_folder(self):
        """Clean the specific folder named '{sample_name}_Bin[10, 20, 50, 100]'."""
        target_folder = os.path.join(self.exp_sq_bin, f"{self.tissue_name}_bin[10, 20, 50, 100]")
        if os.path.exists(target_folder):
            shutil.rmtree(target_folder)
            self.clean_specific_folder.logger.info(f"Deleted specific folder: {target_folder}")
        else:
            self.clean_specific_folder.logger.info(f"Specific folder does not exist: {target_folder}")

    @utils.add_log
    def run(self):
        self.prepare()
        if self.segment:
            self.slice_registration(self.model, self.mode,
                                   enhance_method=self.enhance_method,
                                   enhance_params=self.enhance_params)
            self.splitplot_experiment(self.segment_type, self.mode)
        if self.count:
            self.features = self.get_gtf(self.ref_genome)
            if self.segment_type == 'tissue':
                self.generate_tissue_barcode_position()
                self.generate_count_detail()
                self.spatial_bin_generate()
            elif self.segment_type == 'cell':
                pass
        self.get_bin_summary()

        barcode_csv = os.path.join(self.exp_bs_path, f"{self.tissue_name}_Barcodes_tissue_positions.csv")
        target_csv = os.path.join(os.path.dirname(self.exp_bs_path),
                                  f"{self.tissue_name}_Barcodes_tissue_positions.csv")
        if os.path.exists(barcode_csv):
            shutil.move(barcode_csv, target_csv)

        bin_size = self._micron_bin
        bin_folder = os.path.join(self.exp_sq_bin, f"{self.tissue_name}_bin{bin_size}")
        spatial_dir = os.path.join(bin_folder, 'spatial')
        tissue_positions_path = os.path.join(spatial_dir, 'tissue_positions.csv')

        os.makedirs(spatial_dir, exist_ok=True)

        if list(self.barcodes_detail.columns) != ['barcode', 'x', 'y']:
            self.barcodes_detail.columns = ['barcode', 'x', 'y']

        self.barcodes_detail.to_csv(tissue_positions_path, index=False)
        print(f"Generated tissue_positions.csv at {tissue_positions_path}")

        self.clean_redundant_folders()
        self.clean_specific_folder()

    def _manual_mask_candidates(self, kind):
        os.makedirs(self.runtime_mask_dir, exist_ok=True)
        if kind == 'gem':
            names = [
                'manual_mask.png',
                f'{self.tissue_name}_manual_mask.png',
                f'{self.tissue_name}_mask_manual.png',
                'mask_manual.png',
            ]
        elif kind == 'he':
            names = [
                'he_manual_mask.png',
                f'{self.tissue_name}_he_manual_mask.png',
                f'{self.tissue_name}_he_mask_manual.png',
                'he_mask_manual.png',
            ]
        else:
            raise ValueError(f"Unknown manual mask kind: {kind}")

        return [os.path.join(self.runtime_mask_dir, name) for name in names]

    def _find_manual_mask(self, kind):
        for path in self._manual_mask_candidates(kind):
            if os.path.exists(path):
                return path
        return None

    def _manual_mask_hint(self, kind):
        return self._manual_mask_candidates(kind)[0]

    def _resize_manual_mask(self, manual_mask, manual_mask_path):
        target_height, target_width = int(self.tissue_bbox[2]), int(self.tissue_bbox[3])
        source_height, source_width = manual_mask.shape[:2]
        if (source_height, source_width) == (target_height, target_width):
            return manual_mask

        source_ratio = source_width / source_height
        target_ratio = target_width / target_height
        ratio_error = abs(source_ratio - target_ratio) / target_ratio
        if ratio_error > 0.002:
            raise ValueError(
                "Manual GEM mask has a different coordinate canvas: "
                f"{manual_mask_path} is {source_width}x{source_height}, "
                f"but the expression canvas is {target_width}x{target_height}. "
                "Redraw the mask on the current 1_gem_expression.png canvas "
                "without cropping or changing its aspect ratio."
            )

        print(
            "[Manual Mask] Scaling mask from "
            f"{source_width}x{source_height} to {target_width}x{target_height} "
            "with nearest-neighbor interpolation"
        )
        return cv2.resize(manual_mask, (target_width, target_height), interpolation=cv2.INTER_NEAREST)

    def _sync_gem_expression_to_runtime_mask_dir(self, overlay_path):
        os.makedirs(self.runtime_mask_dir, exist_ok=True)
        target = os.path.join(self.runtime_mask_dir, '1_gem_expression.png')
        try:
            shutil.copy2(overlay_path, target)
            print(f"[Manual Mask] Copied GEM expression reference to: {target}")
        except OSError as exc:
            print(f"[Manual Mask] Warning: failed to copy GEM expression reference to {target}: {exc}")


@utils.add_log
def validate_parameters(args):
    # if not args.sample.startswith("ST") and not args.sample.startswith("OST"):
    #     raise ValueError('--sample must start with "ST" or "OST"')

    # Auto-detect HE image if not provided
    if not args.tif or not os.path.exists(args.tif):
        # Search for HE image files with pattern: {sample}_he.png or {sample}_he.jpg
        possible_locations = []

        # Check in input directory
        if args.input and os.path.exists(args.input):
            possible_locations.append(args.input)
            # Also check for sample subdirectory in input
            sample_subdir = os.path.join(args.input, args.sample)
            if os.path.exists(sample_subdir):
                possible_locations.append(sample_subdir)

        # Check in parent directories of input and common image locations
        if args.input:
            parent_dir = os.path.dirname(args.input)
            if os.path.exists(parent_dir):
                possible_locations.append(parent_dir)
                # Check for common image directories.
                for subdir in ['images', 'image', 'data']:
                    img_dir = os.path.join(parent_dir, subdir)
                    if os.path.exists(img_dir):
                        possible_locations.append(img_dir)
                        # Also check for sample subdirectory within these
                        sample_img_dir = os.path.join(img_dir, args.sample)
                        if os.path.exists(sample_img_dir):
                            possible_locations.append(sample_img_dir)

                # Check grandparent directory structure (common in workflow setups)
                grandparent_dir = os.path.dirname(parent_dir)
                if os.path.exists(grandparent_dir):
                    for subdir in ['images', 'image', 'data']:
                        img_dir = os.path.join(grandparent_dir, subdir)
                        if os.path.exists(img_dir):
                            possible_locations.append(img_dir)
                            sample_img_dir = os.path.join(img_dir, args.sample)
                            if os.path.exists(sample_img_dir):
                                possible_locations.append(sample_img_dir)

        # Search for HE image files
        he_image_found = None
        for location in possible_locations:
            for ext in ['.png', '.jpg', '.jpeg']:
                he_path = os.path.join(location, f"{args.sample}_he{ext}")
                if os.path.exists(he_path):
                    he_image_found = he_path
                    print(f"Auto-detected HE image: {he_path}")
                    break
            if he_image_found:
                break

        if he_image_found:
            args.tif = he_image_found
        elif args.method in ['image', 'ssDNA'] and args.segment:
            print(f"Warning: tissue image not found for {args.method} mode. Searched {len(possible_locations)} candidate locations")

    if args.rectify:
        assert os.path.exists(args.tif), "HE image is required for rectification! Please check the path."
    if args.segment:
        segment_required_args = ['input', 'outdir', 'model']
        # Image-guided and HE modes require a tissue image file.
        if args.method in ['image', 'ssDNA', 'HE']:
            segment_required_args.append('tif')

        if not all(getattr(args, arg, None) for arg in segment_required_args):
            missing_args = [arg for arg in segment_required_args if not getattr(args, arg, None)]
            raise ValueError(f'The following arguments are required when using --segment: {", ".join(missing_args)} '
                             f'(method: {args.method})')
        else:
            # Both image and HE modes need to verify tif file exists
            if args.tif and not os.path.exists(args.tif) and args.method in ['image', 'ssDNA', 'HE']:
                raise FileNotFoundError(f"tif file {args.tif} not found, please check it again")
            if args.model and not os.path.exists(args.model):
                raise FileNotFoundError(f"model file {args.model} not found, please check it again")
    if args.count:
        count_required_args = ['input', 'outdir', 'count_detail', 'genomeDir']
        if not all(getattr(args, arg, None) for arg in count_required_args):
            missing_args = [arg for arg in count_required_args if not getattr(args, arg, None)]
            raise ValueError(f'These arguments are required when using --count: {", ".join(missing_args)}')
        else:
            if not os.path.exists(args.count_detail):
                raise FileNotFoundError(f"count detail file {args.count_detail} not found, please check it again")
            if not os.path.exists(args.genomeDir):
                raise FileNotFoundError(f"genomeDir file {args.genomeDir} not found, please check it again")

    standing_parameters = {
        'sample': args.sample,
        'input': args.input,
        'outdir': args.outdir,
        'model': args.model,
        'prompt': args.prompt,
        'pixel_size': args.pixel_size,
        'genomeDir': args.genomeDir,
        'tif': args.tif,
        'count_detail': args.count_detail,
        'segment': args.segment,
        'count': args.count,
        'extend': args.extend,
    }

    sample = standing_parameters['sample']
    gene_type = standing_parameters['genomeDir']
    count_detail = standing_parameters['count_detail']
    pixel_size = standing_parameters['pixel_size']

    if args.count:
        if sample not in count_detail:
            raise ValueError(f"sample {sample} not fit in count detail path, check it again")
    if pixel_size <= 0:
        raise ValueError(f"pixel size {pixel_size} is invalid, please check it again")
    if isinstance(gene_type, list):
        if gene_type[0] not in genome:
            raise ValueError(f"reference genome {gene_type[0]} not found, please check it again")
    print("these parameters are valid.")

@utils.add_log
def binSegment(args):
    validate_parameters(args)
    with BinSegment(args, display_title="SquareBin") as runner:
        runner.run()
        # runner.test_get_bin_summary()

def get_opts_binSegment(parser, sub_program):
    parser.add_argument(
        '--input',
        help='data(Barcodes & Bbox) path of the current sample',
        type=str
    )
    parser.add_argument(
        '--model',
        help='image segmentation model path',
        type=str
    )
    parser.add_argument(
        '--segment-type',
        help='run segmentation mode, tissue or cell',
        type=str,
        default='tissue'
    )
    parser.add_argument(
        '--method',
        help='image segmentation method',
        type=str,
        default='gene_expr'
    )
    parser.add_argument(
        '--pixel-size',
        help='Each pixel corresponds to a certain number of micrometers.',
        type=float,

    )
    parser.add_argument(
        '--prompt',
        help="""The coordinates of the positions anchored on the chip are utilized for vision transformer model segmentation.
The coordinates are manually selected. (This argument is deprecated and will be removed in future versions)""",
        type=str,
    )
    parser.add_argument(
        '--tif',
        help='Fluorescence Tissue Imaging under Optical Microscopy',
        type=str,
    )
    parser.add_argument(
        '--genomeDir',
        help=HELP_DICT['genomeDir'],
    )
    parser.add_argument(
        '--count_detail',
        help="""The detailed information of barcode counts, including the sequence of each barcode, gene ID, UMI, and count, 
will be used for data integration this step.""",
        type=str,
    )
    parser.add_argument(
        '--enhance-method',
        help="""Dynamic range enhancement method for gene expression data (gene_expr method only).
Options: 'percentile' (default), 'clahe', 'gamma', 'bilateral_percentile', None.
This enhances the contrast of the gene expression heatmap to improve binary segmentation.""",
        type=str,
        default='percentile',
        choices=['percentile', 'clahe', 'gamma', 'bilateral_percentile', 'none']
    )
    parser.add_argument(
        '--enhance-params',
        help="""Parameters for enhancement method in JSON format.
Examples:
- percentile: '{"p_low":2, "p_high":98}'
- clahe: '{"clip_limit":2.0, "tile_grid_size":[8,8]}'
- gamma: '{"gamma":0.5}'
- bilateral_percentile: '{"d":9, "sigma_color":75, "sigma_space":75, "p_low":2, "p_high":98}'
- Add suppress_noise=false to disable bottom noise suppression
- Add noise_threshold=30 to set custom noise threshold
""",
        type=str,
        default=None
    )
    parser.add_argument(
        '--gem-bin-size',
        help="""Bin size in microns for gene expression aggregation (gene_expr method only).
This controls the resolution at which UMI counts are aggregated to generate the gene expression heatmap.
Smaller values (e.g., 10, 20) provide finer detail but may be noisier.
Larger values (e.g., 50, 100) provide smoother results but less spatial resolution.
Default: 10 microns.""",
        type=int,
        default=10
    )
    parser.add_argument(
        '--umi-min-threshold',
        help="""Minimum UMI threshold for bottom noise pre-filtering (gene_expr method only).
Bins with UMI counts below this threshold will be set to 0 BEFORE Otsu segmentation.
This helps improve edge detection by removing background noise.
Options:
  - 'auto' (default): Automatically use 5th percentile of non-zero UMI counts
  - 'none': Disable UMI filtering
  - Numeric value (e.g., 10, 20, 50): Use specific UMI threshold
Recommended: Start with 'auto', then adjust if edges are still noisy.""",
        type=str,
        default='auto'
    )
    parser.add_argument(
        '--ssdna-threshold-scale',
        type=float,
        default=1.0,
        help=(
            'Scale the Otsu threshold for ssDNA tissue masking. Values below 1.0 retain dim '
            'tissue edges; values above 1.0 are more conservative. Default: 1.0.'
        ),
    )
    parser.add_argument(
        '--ssdna-mask-expand-pixels',
        type=int,
        default=0,
        help=(
            'Expand the final ssDNA tissue mask outward by this many registered-image pixels. '
            'Use a small value only to compensate for edge loss. Default: 0.'
        ),
    )
    parser.add_argument(
        '--ssdna-min-hole-area', type=int, default=5000,
        help='Fill only enclosed ssDNA mask holes smaller than this area (pixels); larger internal blank regions remain background.',
    )
    parser.add_argument(
        '--registration-type',
        help="""Image registration method for aligning HE image to gene expression mask (gene_expr method only).
Two-stage registration: Stage 1 uses contour-based coarse alignment, Stage 2 (optional) uses SimpleITK refinement.
Options for Stage 2 refinement (if enabled):
  - 'similarity' (default): Similarity transform with uniform scaling, rotation, and translation
                  Good when tissue shape is preserved but size may differ
  - 'affine': Full affine transform with rotation, translation, scaling, and shear
                        Best for general cases with potential tissue distortion
  - 'rigid': Rigid transform with only rotation and translation (no scaling/shear)
             Use when tissue shapes are identical but just need alignment
Recommended: 'similarity' for most cases.""",
        type=str,
        default='similarity',
        choices=['affine', 'similarity', 'rigid']
    )
    parser.add_argument(
        '--use-sitk-refinement',
        help="""Enable SimpleITK refinement after coarse contour-based registration (gene_expr method only).
  - 'true': Enable two-stage registration (contour coarse + SimpleITK refinement)
  - 'false' (default): Use only contour-based registration (faster, usually sufficient)
Recommended: Start with 'false', enable if coarse registration is insufficient.""",
        type=str,
        default='false',
        choices=['true', 'false']
    )
    parser.add_argument(
        '--use-feature-refinement',
        help="""Enable feature-based refinement (SIFT/ORB) after coarse contour-based registration (gene_expr method only).
This uses a hybrid approach: coarse alignment via contour matching, then fine alignment via feature point matching.
  - 'true': Enable feature-based fine registration (RECOMMENDED for best accuracy)
           - Automatically handles partial sequencing (only HE subset has data)
           - Uses ROI extraction to avoid whole-image mismatches
           - SIFT/ORB features for sub-pixel accuracy
  - 'false' (default): Use only contour-based registration

Priority: if both --use-feature-refinement and --use-sitk-refinement are 'true', feature-based takes precedence.

Recommended: 'true' for high-quality HE images, 'false' if HE image quality is poor or lacks texture.""",
        type=str,
        default='false',
        choices=['true', 'false']
    )
    parser.add_argument('--fluorescence-background', action='store_true',
                        help='Use black outside the tissue ROI for fluorescence images registered through HE mode; does not change registration or the StarDist model.')
    parser.add_argument('--require-manual-he-mask', action='store_true', help='Require an image-side ROI mask for DAPI registration through HE mode.')
    parser.add_argument(
        '--gene-mask-filter',
        action='store_true',
        help="""Restrict final tissue/bin mask to the GEM/manual gene mask.
For gene_expr this makes the GEM/manual mask explicit. For image/ssDNA, a GEM/manual mask
is generated or loaded and intersected with the ssDNA tissue mask."""
    )
    if sub_program:
        parser.add_argument('--segment', action='store_true', help='segment tissue image by pre-trained model')
        parser.add_argument('--count', action='store_true', help='count of in-tissue image')
        parser.add_argument('--rectify', action='store_true', help='whether to rectify the tissue by the HE image roi')
        parser.add_argument('--extend', action='store_true', help='extend barcode from 12 to 18 (deprecated)')
        parser = s_common(parser)
    return parser

def main():
    parser = argparse.ArgumentParser(description='Celatlas Spatial', formatter_class=ArgFormatter)
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers(dest='subparser_assay')
    subparser_1st = subparsers.add_parser('rna')
    subparser_2nd = subparser_1st.add_subparsers()
    subparser_bin_segment = subparser_2nd.add_parser('binSegment_analysis', formatter_class=ArgFormatter)
    sample_default = 'OST110014'
    workspace = os.environ.get('CELATLAS_WORKSPACE', os.path.join(os.path.expanduser('~'), 'celatlas_spatial'))
    reference_root = os.environ.get('CELATLAS_REFERENCE_DIR', os.path.join(workspace, 'reference'))
    results_root = os.environ.get('CELATLAS_RESULTS_ROOT', os.path.join(workspace, 'results'))
    src_dir = os.environ.get('CELATLAS_SRC_DIR', os.path.join(workspace, 'src'))
    demo_result_dir = os.path.join(results_root, 'demo', sample_default)
    demo_mask_dir = os.path.join(demo_result_dir, '06.segment', 'mask')
    subparser_bin_segment.add_argument('--model', type=str, default=os.path.join(src_dir, 'swin_tiny.pth'), help='path to model')
    subparser_bin_segment.add_argument('--thread', type=int, default=8, help='number of threads')
    subparser_bin_segment.add_argument('--pixel-size', type=float, default=0.5, help='pixel size of microscope image')
    subparser_bin_segment.add_argument('--sample', type=str, default=sample_default, help='tissue name')
    subparser_bin_segment.add_argument('--prompt', type=str, default='[[5168,2928],[6384,4416]]', help='tissue position prompt')
    subparser_bin_segment.add_argument('--genomeDir', type=str, default=os.path.join(reference_root, 'Mus_musculus'), help='reference genome type')
    subparser_bin_segment.add_argument('--input', type=str, default=demo_mask_dir, help='path to runtime mask/input files')
    subparser_bin_segment.add_argument(
        '--outdir',
        type=str,
        default=os.path.join(demo_result_dir, '06.segment', '01.binsegment'),
        help='output path'
    )
    subparser_bin_segment.add_argument('--tif', type=str, default=os.path.join(demo_mask_dir, f'{sample_default}.tif'), help='tissue microscope image path')
    subparser_bin_segment.add_argument('--count_detail', type=str, default=os.path.join(demo_result_dir, '05.count', f'{sample_default}_count_detail.txt'), help='reads path')
    subparser_bin_segment.add_argument('--segment-type', type=str, default='tissue', help='run segmentation mode, tissue or cell')
    subparser_bin_segment.add_argument(
        '--method',
        type=str,
        default='gene_expr',
        choices=['gene_expr', 'ssDNA', 'image', 'HE'],
        help="""Segmentation method:
  - 'gene_expr': Gene expression only (no images)
  - 'ssDNA' / 'image': Image-based segmentation (ssDNA images)
  - 'HE': Gene expression + HE image registration (requires *_he.tif/png/jpg)"""
    )

    # New parameters for improved HE mode
    subparser_bin_segment.add_argument(
        '--he-multi-region-threshold',
        type=float,
        default=0.1,
        help='Area threshold ratio for keeping multiple HE tissue regions (default: 0.1 = 10%% of largest region)'
    )

    subparser_bin_segment.add_argument(
        '--tissue-mask-strategy',
        type=str,
        default='gem_primary',
        choices=['gem_primary', 'intersection', 'union'],
        help="""Tissue mask generation strategy (HE mode only):
  - 'gem_primary': Use GEM as primary source, HE for boundary assistance (recommended, prevents phantom expression)
  - 'intersection': Use GEM ∩ HE intersection (conservative, high confidence regions only)
  - 'union': Use GEM ∪ HE union (liberal, may include false positives)"""
    )

    subparser_bin_segment.add_argument('--rectify', action='store_true', help='rectify tissue image by HE image roi')
    subparser_bin_segment.add_argument('--segment', action='store_true', help='segment tissue image')
    subparser_bin_segment.add_argument('--count', action='store_true', help='gene count')
    subparser_bin_segment.add_argument('--gene-mask-filter', action='store_true', help='restrict final tissue mask to GEM/manual gene mask')
    subparser_bin_segment.add_argument('--fluorescence-background', action='store_true', help='use black outside the fluorescence tissue ROI')
    subparser_bin_segment.add_argument('--require-manual-he-mask', action='store_true')
    subparser_bin_segment.add_argument('--extend', action='store_true', help='extend barcode from 12 to 18')
    subparser_bin_segment.add_argument('--debug', action='store_true', help='debug mode')
    subparser_bin_segment.add_argument('--omics', type=str, default='rna', help='omics type')
    args = parser.parse_args()

    binSegment(args)  # process tissue count


if __name__ == '__main__':
    main()
