
import os
import re
import ast
import json
import tarfile
import argparse
import shutil
import gzip
import cv2
import torch
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
from skimage import measure
from skimage.filters import threshold_multiotsu
from skimage.transform import estimate_transform
from segment_anything import SamPredictor, sam_model_registry

from celatlas_spatial.celatlas import ArgFormatter
from celatlas_spatial.tools import utils
from celatlas_spatial.tools import reference
from celatlas_spatial.tools.process import Processor
from celatlas_spatial.tools.image_seg import SwinChipCut
from celatlas_spatial.tools.step import Step, s_common
from celatlas_spatial.tools.matrix import CountMatrix
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

        self.count_dir = os.path.join(os.path.dirname(os.path.normpath(self.outdir)), '05.count')
        self.exp_bs_path = self.outdir  # 06.binSegment folder path
        self.exp_bin_im = os.path.join(self.exp_bs_path, 'images')
        self.exp_sq_bin = os.path.join(self.exp_bs_path, 'square_bin')
        self.exp_tmp_path = os.path.join(self.exp_bs_path, 'tmp')

        # ========================================
        # Mode Definition: Three execution modes supported
        # ========================================
        # 1. ssDNA mode: ssDNA image + gene expression data (legacy 'image' mode variant)
        # 2. gene_expr mode: Gene expression data only, no images
        # 3. HE mode: H&E stained image + gene expression data (legacy 'image' mode variant)

        self.method = self.args.method  # Raw parameter: "gene_expr", "image", or "HE"
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

        # Registration type (for aligning HE image with gene expression mask)
        self.registration_type = getattr(self.args, 'registration_type', 'similarity')
        # Enable SimpleITK fine-tuning optimization
        use_sitk_str = getattr(self.args, 'use_sitk_refinement', 'false')
        self.use_sitk_refinement = (use_sitk_str.lower() == 'true')
        # Enable feature-based fine registration
        use_feature_str = getattr(self.args, 'use_feature_refinement', 'false')
        self.use_feature_refinement = (use_feature_str.lower() == 'true')

        self.tissue_image_mask = None
        self.he_image_mask = None
        self.he_image = None
        self.gem_mask = None
        self.gem = None
        self.tissue_bbox = np.array(pd.read_csv(self.input['TissueBbox'])['bbox'][0].split('\t')).astype(np.int32)
        self.tissue_image = np.zeros((self.tissue_bbox[2], self.tissue_bbox[3]), dtype=np.uint8)
        self.chip_shape = np.array((self.tissue_bbox[2], self.tissue_bbox[3]), dtype=np.int32)
        self.barcodes_detail = pd.read_csv(self.input['FilterBarcodes'], header=None, names=['x', 'y', 'barcode'])

        if self.segment:
            self.bs_out = os.path.join(self.args.bs_out, self.tissue_name) if self.args.bs_out is not None else None
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
            if self.bs_out:
                os.makedirs(self.bs_out, exist_ok=True)
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
        - If method == 'image' and has image:
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

        elif self.method == 'image':
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
                print("[Mode Detection] gene_expr mode: method=image but no image file found, fallback to gene_expr")
        else:
            # Unknown method, default to gene_expr
            mode = 'gene_expr'
            print(f"[Mode Detection] gene_expr mode: unknown method '{self.method}', fallback to gene_expr")

        self.mode = mode
        print(f"[Mode Detection] Final mode: {mode}")
        print(f"[Mode Validation] Mode is locked and cannot be changed during this run")
        print("=" * 60)

        return mode

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
                if self.bs_out:
                    os.system(f'cp -r {self.im_path} {os.path.join(self.bs_out, f"{self.tissue_name}.tif")}')

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
        hires_scale = 2000 / int(max(self.chip_shape) / bin_size)
        lowres_scale = 600 / int(max(self.chip_shape) / bin_size)

        factors = {
            "tissue_hires_scalef": hires_scale,
            "tissue_lowres_scalef": lowres_scale,
            "fiducial_diameter_fullres": bin_size,
            "spot_diameter_fullres": bin_size
        }
        json_data = json.dumps(factors, indent=4)
        with open(os.path.join(output, 'scalefactors_json.json'), 'w') as json_file:
            json_file.write(json_data)

    @utils.deprecated
    def segment_anything_registration(self, model_path):
        # load sam model
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        sam = sam_model_registry["vit_h"](checkpoint=model_path)
        sam = sam.to(device)
        predictor = SamPredictor(sam)

        # load image and convert to 8bit make sure it can be processed by sam model
        tissue_image = self.tissue_image
        tissue_image = cv2.cvtColor(tissue_image, cv2.COLOR_GRAY2BGR)
        predictor.set_image(tissue_image)

        manual_prompt = np.array(self.prompt)  # manually select prompt
        tissue_masks, scores, logits = predictor.predict(
            point_labels=np.full(manual_prompt.shape[0], 1),
            point_coords=manual_prompt,
            multimask_output=True
        )
        tissue_masks = np.transpose(tissue_masks, (1, 2, 0)).astype(np.uint8)
        tissue_masks = cv2.cvtColor(tissue_masks, cv2.COLOR_BGR2GRAY) * 255
        tissue_image = cv2.cvtColor(tissue_image, cv2.COLOR_BGR2GRAY)

        test_output = tissue_image.copy()
        test_output[tissue_masks == 0] = 0
        if self.bs_out:
            cv2.imwrite(os.path.join(self.bs_out, f'{self.tissue_name}_test.png'), test_output)

        # load chip image
        chip_image = np.zeros((self.tissue_bbox[2], self.tissue_bbox[3]), dtype=np.uint8)
        chip = self.barcodes_detail.to_numpy()
        chip_contours = np.array([np.min(chip[:, 0]), np.max(chip[:, 0]), np.min(chip[:, 1]), np.max(chip[:, 1])])
        chip_masks = np.zeros_like(chip_image)
        chip_masks[chip_contours[2]: chip_contours[3], chip_contours[0]: chip_contours[1]] = 255

        row_coords, col_coords = np.where(tissue_masks == 255)
        min_row, max_row = np.min(row_coords), np.max(row_coords)
        min_col, max_col = np.min(col_coords), np.max(col_coords)
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
        tissue_image = test_output
        tissue_image = cv2.warpPerspective(tissue_image, H, (chip_image.shape[1], chip_image.shape[0]))

        # scanpy image for visualization
        scanpy_image = tissue_image.copy()
        scanpy_image_600 = cv2.resize(scanpy_image, (600, 600), interpolation=cv2.INTER_LANCZOS4)
        scanpy_image_2000 = cv2.resize(scanpy_image, (2000, 2000), interpolation=cv2.INTER_LANCZOS4)

        self.tissue_image = tissue_image

        torch.cuda.empty_cache()  # release gpu memory after registration
        if self.bs_out:
            tifffile.imwrite(os.path.join(self.bs_out, f'{self.tissue_name}_regist.tif'), tissue_image)
        tifffile.imwrite(os.path.join(self.exp_bin_im, f'{self.tissue_name}_regist.tif'), tissue_image)
        cv2.imwrite(os.path.join(self.exp_bin_im, f'tissue_lowres_image.png'), scanpy_image_600)
        cv2.imwrite(os.path.join(self.exp_bin_im, f'tissue_hires_image.png'), scanpy_image_2000)

    @utils.add_log
    def slice_registration(self, model, method="image", enhance_method='percentile', enhance_params=None):
        """
        Main slice registration workflow

        Parameters:
        -----------
        model : str
            Model path (for image method)
        method : str
            Registration method: 'gene_expr', 'image', or 'HE'
        enhance_method : str
            Gene expression enhancement method (for gene_expr/HE methods only)
        enhance_params : dict
            Enhancement method parameter dictionary (for gene_expr/HE methods only)
        """
        # Use configurable gem_bin_size (in microns) instead of hardcoded 10
        resolution = np.ceil(self.gem_bin_size / self.pixel_size).astype(int)
        print(f"Using gem_bin_size={self.gem_bin_size} microns (resolution={resolution} pixels) for gene expression aggregation")

        if method == "gene_expr":
            # Pure gene expression mode: No HE registration
            self.gem_registration(resolution=resolution, enhance_method=enhance_method, enhance_params=enhance_params)
            # self.save_gene_expr_enhancement_comparison()  # Debug visualization, not needed in production

            # Use gem as tissue_image
            self.tissue_image = self.gem
            self.he_image_mask = self.gem_mask
            print("[gene_expr mode] No HE registration, using gem as tissue_image")

        elif method == "HE":
            # HE mode: Gene expression + HE registration
            self.gem_registration(resolution=resolution, enhance_method=enhance_method, enhance_params=enhance_params)

            # Must have HE image
            if self.args.tif and os.path.exists(self.args.tif):
                filename = os.path.basename(self.args.tif)
                print(f"[HE mode] Performing HE registration with: {filename}")
                self.he_registration(self.args.tif)
                self.save_he_registration_comparison()  # Generates overlays folder (debug plot disabled internally)
            else:
                print(f"[HE mode] ERROR: --tif parameter required for HE mode!")
                print(f"[HE mode] Fallback to gene_expr mode")
                self.tissue_image = self.gem
                self.he_image_mask = self.gem_mask

            self.save_gene_expr_enhancement_comparison()  # Debug plot disabled internally

        elif method == "image":
            # ssDNA image segmentation mode
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
        if self.bs_out:
            cv2.imwrite(os.path.join(self.bs_out, f'{self.tissue_name}_chip_cut.png'), im_chip_cut)

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
            self.he_image_mask = self.gem_mask

        if self.rectify:
            processor = Processor(output=self.exp_bin_im, image=self.he_image)
            # roi_mask = processor.detect_border(border_value=179, tolerance=5, debug=False)
            roi_mask, _ = processor.detect_polygons(border_value=179, tolerance=5, min_area=100)
            degree, slope, center = processor.find_top_line(roi_mask)

            self.tissue_image, new_center = processor.get_rotate_image(self.tissue_image, degree, center)
            self.he_image_mask, new_center = processor.get_rotate_image(roi_mask, degree, center)

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
            self.he_image_mask = cv2.warpAffine(self.he_image_mask, motion_m, (self.he_image_mask.shape[1] + dx, self.he_image_mask.shape[0] + dy))
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

            # Get initial transform matrix using contour method - internal function for calculation
            def get_contour_features(contour):
                M = cv2.moments(contour)
                cx = int(M['m10'] / M['m00'])
                cy = int(M['m01'] / M['m00'])
                (x, y), (MA, ma), angle = cv2.fitEllipseAMS(contour)
                return (cx, cy), angle

            contours_fixed, _ = cv2.findContours(image_fixed.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            contours_moving, _ = cv2.findContours(image_moving.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

            contour_fixed = max(contours_fixed, key=cv2.contourArea)
            contour_moving = max(contours_moving, key=cv2.contourArea)
            center_fixed, angle_fixed = get_contour_features(contour_fixed)
            center_moving, angle_moving = get_contour_features(contour_moving)

            rotation_angle = -(angle_fixed - angle_moving)
            area_fixed = cv2.contourArea(contour_fixed)
            area_moving = cv2.contourArea(contour_moving)
            scale = np.sqrt(area_fixed / area_moving)

            if verbose:
                print(f"[Contour Init] Rotation: {rotation_angle:.2f}°, Scale: {scale:.3f}")
                print(f"[Contour Init] Center fixed: {center_fixed}, Center moving: {center_moving}")

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
        Contour-based registration (supports images of different sizes)

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
        def get_contour_features(contour):
            # Calculate moments of the contour
            M = cv2.moments(contour)
            # Calculate the centroid
            cx = int(M['m10'] / M['m00'])
            cy = int(M['m01'] / M['m00'])
            # Calculate the angle of the contour
            (x, y), (MA, ma), angle = cv2.fitEllipseAMS(contour)
            return (cx, cy), angle

        contours_fixed, _ = cv2.findContours(image_fixed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        contours_moving, _ = cv2.findContours(image_moving, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        contour_fixed = contours_fixed[0]
        contour_moving = contours_moving[0]
        center_fixed, angle_fixed = get_contour_features(contour_fixed)
        center_moving, angle_moving = get_contour_features(contour_moving)

        print(f"[Contour Features] Fixed - center: {center_fixed}, angle: {angle_fixed:.2f}°")
        print(f"[Contour Features] Moving - center: {center_moving}, angle: {angle_moving:.2f}°")

        # Calculate scale based on contour area ratio (automatically adapt to different resolutions)
        area_fixed = cv2.contourArea(contour_fixed)
        area_moving = cv2.contourArea(contour_moving)
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
    def splitplot_experiment(self, seg_type, method):
        if seg_type == 'tissue':
            self.tissue_segmentation(method)
        elif seg_type == 'cell':
            self.cell_segmentation()  # future work

    @utils.add_log
    def tissue_segmentation(self, method="image"):
        def getArea(elem):
            return elem.area

        if method == "image":
            # tissue threshold segmentation
            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (20, 20))
            _, tissue_mask = cv2.threshold(self.tissue_image, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
            tissue_mask = cv2.morphologyEx(tissue_mask, cv2.MORPH_CLOSE, kernel, iterations=8)

            # choose tissue prop
            label_image = measure.label(tissue_mask, connectivity=2)
            props = measure.regionprops(label_image, intensity_image=tissue_mask)
            props.sort(key=getArea, reverse=True)

            areas = [p['area'] for p in props]
            label_num = len(areas) if np.std(areas) * 10 < np.mean(areas) else int(np.sum(areas >= np.mean(areas)))

            result = np.zeros_like(self.tissue_image).astype(np.uint8)
            for i in range(label_num):
                prop = props[i]
                result += np.where(label_image != prop.label, 0, 1).astype(np.uint8)

            # find seed point for flood fill
            tissue_result = utils.hole_fill(result)

            # Additional background noise removal: remove small objects and erode slightly
            # This helps eliminate edge artifacts and isolated noise pixels
            small_kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
            tissue_result = skimage.morphology.remove_small_objects(tissue_result > 0, min_size=1000).astype(np.uint8) * 255
            # Slight erosion to remove noisy edges (adjust iterations to control strictness)
            tissue_result = cv2.erode(tissue_result, small_kernel, iterations=2)
            # Re-dilate to recover tissue area (but noise won't come back)
            tissue_result = cv2.dilate(tissue_result, small_kernel, iterations=2)

            self.tissue_image_mask = tissue_result
            print(f"Tissue mask generated with enhanced background removal")
        elif method == "gene_expr" or method == "HE":
            # Both gene_expr and HE modes use he_image_mask
            self.tissue_image_mask = self.he_image_mask

        if self.bs_out:
            tifffile.imwrite(os.path.join(self.bs_out, f'{self.tissue_name}_tissue_cut.tif'), self.tissue_image_mask)
        tifffile.imwrite(os.path.join(self.exp_bin_im, f'{self.tissue_name}_tissue_cut.tif'), self.tissue_image_mask)

    @utils.add_log
    def cell_segmentation(self):
        """
        discussion: sooner we may use cellpose-3.0 to implement cell segmentation
        """
        pass

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
            self.he_image_mask = self.he_image_mask[self.roi_rect[1]: self.roi_rect[1] + self.roi_rect[3], self.roi_rect[0]: self.roi_rect[0] + self.roi_rect[2]]
            self.tissue_image = cv2.resize(self.tissue_image, (self.slice_px, self.slice_px))
            self.he_image_mask = cv2.resize(self.he_image_mask, (self.slice_px, self.slice_px))
            self.he_image_mask[self.he_image_mask > 0] = 255

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
        mask_shape = self.tissue_image_mask.shape
        bin_args = [[i, SquareBin(
            self._micron_bin[i], self.bin[i], [0, 0, mask_shape[1], mask_shape[0]],
            bin_folder=os.path.join(self.exp_sq_bin, f'{self.tissue_name}_bin{self._micron_bin[i]}'),
            gene_count_detail=self.gene_count_detail.copy(deep=True))] for i in range(len(self._micron_bin))
        ]
        if self._raw_show:
            bin_args.append(
                [None, SquareBin('Raw', self._bin, [0, 0, mask_shape[1], mask_shape[0]],
                                 bin_folder=os.path.join(self.exp_sq_bin, f'{self.tissue_name}_Raw'),
                                 gene_count_detail=self.raw_count_detail)]
            )
        if not self._parallel:
            for args in bin_args:
                self.bin_files_execute(args)
        else:
            # Memory-optimized parallel processing:
            # - Use threading backend to share memory instead of forking processes
            # - Limit n_jobs to 2 to reduce memory footprint (balance speed vs memory)
            # - Original: n_jobs=5 uses ~300GB RAM, crashes with OOM
            # - Optimized: n_jobs=2 uses ~80GB RAM, only 1.5x slower
            max_parallel = int(os.environ.get('BINSEGMENT_MAX_PARALLEL', '2'))
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
                os.system(f"gzip {dst_file}") 


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

        if not os.path.exists(coords_file):
            raise FileNotFoundError(f"The coords file was not found: {coords_file}")
        coords_data = pd.read_csv(coords_file)

        if not os.path.exists(scalefactors_file):
            raise FileNotFoundError(f"The scalefactors file was not found: {scalefactors_file}")
        with open(scalefactors_file) as f:
            scalefactors = json.load(f)
        hires_scalef = scalefactors["tissue_hires_scalef"]

        coords_data["pxl_row_in_fullres"] = coords_data["row_in_bin"] * hires_scalef
        coords_data["pxl_col_in_fullres"] = coords_data["col_in_bin"] * hires_scalef

        cols = ["Barcode", "in_tissue", "row_in_bin", "col_in_bin", "pxl_row_in_fullres", "pxl_col_in_fullres"]
        coords_data = coords_data[cols]
        coords_data.rename(columns={"Barcode": "barcode"}, inplace=True)

        coords_data.to_csv(output_file, index=False, header=False)

        if os.path.exists(coords_file):
            os.remove(coords_file)
        os.rename(output_file, coords_file)
        print(f"The original has been updated tissue_positions_list.csv: {coords_file}")


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
        unique_gene_count_detail = gene_detail[columns_extract]
        unique_gene_count_detail = unique_gene_count_detail.rename(columns={'unicode': 'Barcode'})
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
        import pandas as pd
        import numpy as np

        try:
            # Skip for Raw bin
            if bin_handler.micron_bin == 'Raw':
                return

            # Prepare data for downsampling
            # Create index array where each row is repeated by its 'count' value
            count_detail_array = count_detail.copy()
            count_detail_array['count'] = count_detail_array['count'].astype(int)
            bin_read_index = np.array(
                count_detail_array.index.repeat(count_detail_array['count']),
                dtype='int32'
            )
            np.random.shuffle(bin_read_index)

            # Initialize downsample results
            downsample_results = {
                'read_fraction': [0.0],
                'median_gene_number': [0.0],
                'umi_saturation': [0.0],
                'read_saturation': [0.0]
            }

            # Downsample at different fractions
            for fraction in np.arange(0.1, 1.1, 0.1):
                fraction = round(fraction, 1)

                # Calculate subsample size
                total_reads = len(bin_read_index)
                subsample_size = int(total_reads * fraction)
                subsample_index = bin_read_index[:subsample_size]

                # Get unique entries (deduplicated UMIs)
                unique_index, counts = np.unique(subsample_index, return_counts=True)
                n_deduped_reads = len(unique_index)

                # Calculate saturation
                saturation = round((1 - n_deduped_reads / subsample_size) * 100, 2)

                # Calculate median genes per bin
                df_subsample = count_detail_array.loc[unique_index, :]
                median_genes = float(
                    df_subsample.groupby('Barcode').agg({'geneID': 'nunique'}).median().iloc[0]
                )

                downsample_results['read_fraction'].append(fraction)
                downsample_results['median_gene_number'].append(median_genes)
                downsample_results['umi_saturation'].append(saturation)
                downsample_results['read_saturation'].append(saturation)

            # Save to file
            df_downsample = pd.DataFrame(downsample_results)
            downsample_file = os.path.join(bin_handler.bin_folder, 'downsample.tsv')
            df_downsample.to_csv(downsample_file, index=False, sep='\t')

            print(f"[Bin {bin_handler.micron_bin}] Generated downsample data: median genes {downsample_results['median_gene_number'][-1]:.0f}, saturation {downsample_results['umi_saturation'][-1]:.1f}%")

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
                    self.he_image = cv2.resize(self.he_image, None, fx=resize_scale, fy=resize_scale, interpolation=cv2.INTER_AREA)
                    self.image_d = cv2.resize(self.image_d, None, fx=resize_scale, fy=resize_scale, interpolation=cv2.INTER_AREA)
                    print(f"[HE Image] Resized to {self.he_image.shape} (scale={resize_scale:.3f})")
                else:
                    print(f"[HE Image] No resize needed (HE size ~ tissue_bbox size)")
            else:
                self.image_d = self.image_d

        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))

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
            # Use largest contour
            max_contour = max(contours, key=cv2.contourArea)
            mask = np.zeros(self.he_image.shape[:2], dtype=np.uint8)
            cv2.drawContours(mask, [max_contour], 0, 255, -1)
            filled_mask = mask
            print(f"HE tissue contour extracted: {len(contours)} contours found, largest area={cv2.contourArea(max_contour):.0f}")
        else:
            filled_mask = np.zeros(self.he_image.shape[:2], dtype=np.uint8)
            print("WARNING: No contours found in HE image! Registration may fail.")

        self.he_image_mask = filled_mask

        # Save mask before registration for visualization comparison
        self.he_image_mask_before_reg = filled_mask.copy()

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
            self.he_image_mask = cv2.warpAffine(he_mask_coarse, trans_matrix_fine,
                                               (self.tissue_bbox[3], self.tissue_bbox[2]))
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
            self.he_image_mask = cv2.warpAffine(self.he_image_mask, trans_matrix_refined,
                                               (self.tissue_bbox[3], self.tissue_bbox[2]))
            self.image_d = cv2.warpAffine(self.image_d, trans_matrix_refined,
                                         (self.tissue_bbox[3], self.tissue_bbox[2]))

            trans_matrix = trans_matrix_refined
            print(f"[Registration] Two-stage registration completed (contour + SimpleITK)")
        else:
            # Use coarse registration result only
            self.he_image = he_image_coarse
            self.he_image_mask = he_mask_coarse
            self.image_d = image_d_coarse
            trans_matrix = trans_matrix_coarse
            print(f"[Registration] Coarse registration completed (refinement disabled)")
        self.image_d[self.he_image_mask == 0] = 0
        self.tissue_image = self.image_d

        # Save transformation matrix for visualization
        self.registration_transform = trans_matrix

    @utils.add_log
    def save_register_image(self, is_correct=False):
        # For gene_expr method, use enhanced gem image as tissue image
        if hasattr(self, 'gem_enhanced') and self.gem_enhanced is not None:
            # Use enhanced gem image (already correct size and grayscale)
            scanpy_image = self.gem_enhanced.copy()
            print(f"Using enhanced gene expression image for tissue visualization (grayscale, shape={scanpy_image.shape})")
        else:
            # Use original tissue_image (image method or gene_expr without enhancement)
            scanpy_image = self.tissue_image.copy()

        # shift = self.get_offset(scanpy_image)
        # transM = np.float32([[1, 0, -shift[0]], [0, 1, -shift[1]]])
        # scanpy_image = cv2.warpAffine(scanpy_image, transM, (scanpy_image.shape[1], scanpy_image.shape[0]))
        scanpy_image_600 = cv2.resize(scanpy_image, (600, 600))
        scanpy_image_2000 = cv2.resize(scanpy_image, (2000, 2000))

        if is_correct:
            fine_tune = np.float32([[1, 0, -self.bin[-1] // 2], [0, 1, -self.bin[-1] // 2]])
            scanpy_image_600 = cv2.warpAffine(scanpy_image_600, fine_tune, (600, 600))
            scanpy_image_2000 = cv2.warpAffine(scanpy_image_2000, fine_tune, (2000, 2000))

        if scanpy_image.ndim == 3:
            scanpy_image_600 = cv2.cvtColor(scanpy_image_600, cv2.COLOR_RGB2BGR)
            scanpy_image_2000 = cv2.cvtColor(scanpy_image_2000, cv2.COLOR_RGB2BGR)

        if self.bs_out:
            tifffile.imwrite(os.path.join(self.bs_out, f'{self.tissue_name}_regist.tif'), scanpy_image)
        tifffile.imwrite(os.path.join(self.exp_bin_im, f'{self.tissue_name}_regist.tif'), scanpy_image)
        cv2.imwrite(os.path.join(self.exp_bin_im, f'tissue_lowres_image.png'), scanpy_image_600)
        cv2.imwrite(os.path.join(self.exp_bin_im, f'tissue_hires_image.png'), scanpy_image_2000)

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

        if self.bs_out:
            comparison_path_bs = os.path.join(self.bs_out, 'gene_expr_enhancement_comparison.png')
            plt.savefig(comparison_path_bs, dpi=150, bbox_inches='tight')

        plt.close(fig)  # Explicitly close figure after all saves
        print(f"Enhancement comparison saved to: {comparison_path}")

    @utils.add_log
    def save_he_registration_comparison(self):
        """
        Save HE image and gene_expr registration comparison plot to evaluate registration quality
        """
        # DISABLED: Debug comparison plot not needed for customers
        # Still call save_gem_he_image_overlay to generate necessary overlays
        self.save_gem_he_image_overlay()
        return

        if not hasattr(self, 'he_image_mask') or self.he_image_mask is None:
            print("No HE image registration data available")
            return

        # Check if gem_mask is available
        if not hasattr(self, 'gem_mask') or self.gem_mask is None:
            print("No gem_mask available for comparison")
            return

        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        # First row: gem_mask, he_image_mask (before registration), he_image_mask (after registration)
        axes[0, 0].imshow(self.gem_mask, cmap='gray')
        axes[0, 0].set_title('GEM Mask (Reference)')
        axes[0, 0].axis('off')

        # HE mask before registration (need to save original)
        if hasattr(self, 'he_image_mask_before_reg'):
            axes[0, 1].imshow(self.he_image_mask_before_reg, cmap='gray')
            axes[0, 1].set_title('HE Mask (Before Registration)')
        else:
            axes[0, 1].text(0.5, 0.5, 'Not Available', ha='center', va='center')
            axes[0, 1].set_title('HE Mask (Before Registration)')
        axes[0, 1].axis('off')

        # HE mask after registration
        axes[0, 2].imshow(self.he_image_mask, cmap='gray')
        axes[0, 2].set_title('HE Mask (After Registration)')
        axes[0, 2].axis('off')

        # Second row: overlay comparison, Dice coefficient, contour overlay
        # 1. Green-red overlay (gem=green, he=red)
        overlay_rgb = np.zeros((*self.gem_mask.shape, 3), dtype=np.uint8)
        overlay_rgb[:, :, 1] = self.gem_mask  # gem_mask -> green channel
        overlay_rgb[:, :, 0] = self.he_image_mask  # he_mask -> red channel
        axes[1, 0].imshow(overlay_rgb)
        axes[1, 0].set_title('Overlay (GEM=Green, HE=Red)')
        axes[1, 0].axis('off')

        # 2. Dice coefficient and IoU
        intersection = np.logical_and(self.gem_mask > 0, self.he_image_mask > 0).sum()
        union = np.logical_or(self.gem_mask > 0, self.he_image_mask > 0).sum()
        gem_area = np.sum(self.gem_mask > 0)
        he_area = np.sum(self.he_image_mask > 0)

        dice = 2 * intersection / (gem_area + he_area) if (gem_area + he_area) > 0 else 0
        iou = intersection / union if union > 0 else 0

        # Calculate contour alignment - use contour center distance rather than edge pixel overlap
        # This avoids the impact of edge style differences
        contour_alignment_score = 0

        # Extract contours from both masks
        gem_contours, _ = cv2.findContours(self.gem_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        he_contours, _ = cv2.findContours(self.he_image_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

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
            he_contours, _ = cv2.findContours(self.he_image_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

            # Draw contours: gem=green, he=red
            cv2.drawContours(contour_img, gem_contours, -1, (0, 255, 0), 3)  # green
            cv2.drawContours(contour_img, he_contours, -1, (255, 0, 0), 3)  # red

            axes[1, 2].imshow(cv2.cvtColor(contour_img, cv2.COLOR_BGR2RGB))
            axes[1, 2].set_title('Contour Overlay (GEM=Green, HE=Red)')
        else:
            axes[1, 2].text(0.5, 0.5, 'Tissue image not available', ha='center', va='center')
            axes[1, 2].set_title('Contour Overlay')
        axes[1, 2].axis('off')

        plt.tight_layout()

        # Save comparison plot
        comparison_path = os.path.join(self.exp_bin_im, 'he_registration_comparison.png')
        plt.savefig(comparison_path, dpi=150, bbox_inches='tight')

        if self.bs_out:
            comparison_path_bs = os.path.join(self.bs_out, 'he_registration_comparison.png')
            plt.savefig(comparison_path_bs, dpi=150, bbox_inches='tight')

        plt.close(fig)  # Explicitly close figure after all saves
        print(f"HE registration comparison saved to: {comparison_path}")
        print(f"Registration quality: Dice={dice:.5f}, IoU={iou:.5f}, Contour-Alignment={edge_dice:.3f} (visual quality)")

        # Additional save: gem enhanced image and HE original image overlay comparison - needed for Image Alignment tab
        self.save_gem_he_image_overlay()  # ENABLED: Generate overlays folder with HE alignment images

    @utils.add_log
    def save_gem_he_image_overlay(self):
        """Save overlay comparison of gem enhanced image and HE original image for intuitive registration quality assessment"""
        if not hasattr(self, 'gem_enhanced') or not hasattr(self, 'image_d'):
            print("Warning: gem_enhanced or image_d not available, skipping image overlay")
            return

        # DISABLED: Skip generating debug comparison plot (not needed for customers)
        # Only generate individual overlay images for report
        skip_debug_plot = True

        # Prepare images (needed for save_individual_overlay_images)
        gem_img = self.gem_enhanced if self.gem_enhanced is not None else self.gem

        # Use image_d (registered color HE image) instead of he_image
        he_img = self.image_d if hasattr(self, 'image_d') else self.he_image

        # Convert GEM to grayscale (if 3-channel)
        if gem_img.ndim == 3:
            gem_gray = cv2.cvtColor(gem_img, cv2.COLOR_RGB2GRAY)
        else:
            gem_gray = gem_img.copy()

        # Keep HE color or convert
        if he_img.ndim == 2:
            # If HE is grayscale, convert to pseudo-color RGB
            he_color = cv2.cvtColor(he_img, cv2.COLOR_GRAY2RGB)
        elif he_img.shape[2] == 4:
            # If RGBA, convert to RGB
            he_color = cv2.cvtColor(he_img, cv2.COLOR_RGBA2RGB)
        else:
            he_color = he_img.copy()

        # Normalize GEM to 0-255
        gem_normalized = ((gem_gray - gem_gray.min()) / (gem_gray.max() - gem_gray.min()) * 255).astype(np.uint8)

        # Normalize HE to 0-255 (keep color)
        he_normalized = he_color.copy()
        if he_normalized.dtype != np.uint8:
            he_normalized = ((he_normalized - he_normalized.min()) / (he_normalized.max() - he_normalized.min()) * 255).astype(np.uint8)

        # Generate heatmap and color images (needed for save_individual_overlay_images)
        gem_heatmap = plt.cm.jet(gem_normalized / 255.0)[:, :, :3]  # RGB
        gem_heatmap = (gem_heatmap * 255).astype(np.uint8)
        overlay_alpha = cv2.addWeighted(he_normalized, 0.6, gem_heatmap, 0.4, 0)
        gem_color = plt.cm.viridis(gem_normalized / 255.0)[:, :, :3]
        gem_color = (gem_color * 255).astype(np.uint8)
        he_with_contour = he_normalized.copy()
        gem_contours, _ = cv2.findContours((gem_normalized > 20).astype(np.uint8),
                                          cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if len(gem_contours) > 0:
            cv2.drawContours(he_with_contour, gem_contours, -1, (255, 0, 0), 3)

        if not skip_debug_plot:
            # DISABLED: Debug comparison plot generation (customers don't need this)
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))

            # 1. GEM enhanced image (pseudo-color)
            axes[0, 0].imshow(gem_normalized, cmap='viridis')
            axes[0, 0].set_title('GEM Expression (Enhanced)', fontsize=14)
            axes[0, 0].axis('off')

            # 2. HE image (color)
            axes[0, 1].imshow(he_normalized)
            axes[0, 1].set_title('HE Image (Registered, H&E stain)', fontsize=14)
            axes[0, 1].axis('off')

            # 3. GEM heatmap overlay on HE image (semi-transparent)
            axes[0, 2].imshow(overlay_alpha)
            axes[0, 2].set_title('GEM Heatmap on HE (40% opacity)', fontsize=14)
            axes[0, 2].axis('off')
            axes[0, 2].text(0.5, -0.05, 'Shows gene expression overlay on tissue',
                           ha='center', transform=axes[0, 2].transAxes, fontsize=10, style='italic')

            # 4. Checkerboard comparison (alternating GEM and HE in color)
            checkerboard = np.zeros((gem_gray.shape[0], gem_gray.shape[1], 3), dtype=np.uint8)
            block_size = 50

            for i in range(0, gem_gray.shape[0], block_size):
                for j in range(0, gem_gray.shape[1], block_size):
                    if ((i // block_size) + (j // block_size)) % 2 == 0:
                        checkerboard[i:i+block_size, j:j+block_size] = gem_color[i:i+block_size, j:j+block_size]
                    else:
                        checkerboard[i:i+block_size, j:j+block_size] = he_normalized[i:i+block_size, j:j+block_size]
            axes[1, 0].imshow(checkerboard)
            axes[1, 0].set_title(f'Checkerboard (GEM vs HE, block={block_size}px)', fontsize=14)
            axes[1, 0].axis('off')

            # 5. GEM contour boundary on HE image
            axes[1, 1].imshow(he_with_contour)
            axes[1, 1].set_title('GEM Coverage on HE (red outline)', fontsize=14)
            axes[1, 1].axis('off')

            # 6. Side-by-side comparison (Split view)
            split_view = np.zeros_like(he_normalized)
            mid = he_normalized.shape[1] // 2
            split_view[:, :mid] = he_normalized[:, :mid]  # Left=HE
            split_view[:, mid:] = gem_color[:, mid:]      # Right=GEM

            # Draw white dividing line
            cv2.line(split_view, (mid, 0), (mid, split_view.shape[0]), (255, 255, 255), 3)

            axes[1, 2].imshow(split_view)
            axes[1, 2].set_title('Split View (Left=HE, Right=GEM)', fontsize=14)
            axes[1, 2].axis('off')

            plt.tight_layout()

            # Save comprehensive comparison plot
            overlay_path = os.path.join(self.outdir, 'images', 'gem_he_image_overlay.png')
            plt.savefig(overlay_path, dpi=150, bbox_inches='tight')

            if self.bs_out:
                overlay_path_bs = os.path.join(self.bs_out, 'gem_he_image_overlay.png')
                plt.savefig(overlay_path_bs, dpi=150, bbox_inches='tight')

            plt.close(fig)  # Explicitly close figure after all saves
            print(f"GEM-HE image overlay saved to: {overlay_path}")

        # Additional save: save each subplot as separate file
        self.save_individual_overlay_images(gem_normalized, he_normalized, gem_heatmap,
                                           gem_color, overlay_alpha, he_with_contour)

    def save_individual_overlay_images(self, gem_normalized, he_normalized, gem_heatmap,
                                      gem_color, overlay_alpha, he_with_contour):
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

        # Calculate image scaling ratio - reduce if image is too large
        max_dimension = 4000  # Maximum dimension limit is 4000 pixels
        h, w = he_normalized.shape[:2]
        scale = 1.0
        if max(h, w) > max_dimension:
            scale = max_dimension / max(h, w)
            print(f"[Image Compression] Resizing images from {w}x{h} to {int(w*scale)}x{int(h*scale)} (scale={scale:.3f})")

        def resize_if_needed(img):
            """Resize image if needed"""
            if scale < 1.0:
                new_size = (int(img.shape[1] * scale), int(img.shape[0] * scale))
                return cv2.resize(img, new_size, interpolation=cv2.INTER_AREA)
            return img

        # ========================================
        # 1. GEM expression heatmap (needed for all modes)
        # ========================================
        gem_resized = resize_if_needed(gem_normalized)
        cv2.imwrite(os.path.join(images_dir, '1_gem_expression.png'),
                   cv2.applyColorMap(gem_resized, cv2.COLORMAP_VIRIDIS), png_compression)
        print(f"[{self.mode} Mode] 1_gem_expression.png saved")

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
            he_resized = resize_if_needed(he_normalized)

            # Set background (outside mask) to pure white
            he_with_white_bg = he_resized.copy()
            if hasattr(self, 'he_image_mask') and self.he_image_mask is not None:
                mask_resized = resize_if_needed(self.he_image_mask)
                he_with_white_bg[mask_resized == 0] = 255  # White background outside mask

            cv2.imwrite(os.path.join(images_dir, bg_filename),
                       cv2.cvtColor(he_with_white_bg, cv2.COLOR_RGB2BGR), jpeg_quality)
            print(f"[{self.mode} Mode] {bg_filename} saved as JPEG with quality={jpeg_quality[1]} (white background)")

        # ========================================
        # 3. Tissue segmentation mask (mode isolation)
        # ========================================
        if self.mode in ['HE', 'ssDNA']:
            # HE and ssDNA modes: have background, calculate union
            if hasattr(self, 'gem_mask') and hasattr(self, 'he_image_mask'):
                # Calculate union
                union_mask = np.logical_or(self.gem_mask > 0, self.he_image_mask > 0).astype(np.uint8) * 255

                # Save pure mask (compressed)
                mask_resized = resize_if_needed(union_mask)
                cv2.imwrite(os.path.join(images_dir, '3_tissue_segmentation_mask.png'), mask_resized, png_compression)
                print(f"[{self.mode} Mode] 3_tissue_segmentation_mask.png saved (union of GEM and {bg_label})")

                # ========================================
                # 3b. GEM Heatmap overlay on background (mask region only)
                # NOTE: This image is used in the Interactive Image Alignment Viewer
                # ========================================
                # Start with pure white background
                he_with_gem_overlay = np.full_like(he_normalized, 255)

                # Overlay GEM heatmap within mask region (30% transparency)
                mask_bool = union_mask > 0
                if mask_bool.any():
                    # GEM heatmap 30%, background 70%
                    he_with_gem_overlay[mask_bool] = cv2.addWeighted(
                        he_normalized[mask_bool], 0.7,
                        gem_heatmap[mask_bool], 0.3, 0
                    )

                gem_heatmap_overlay_resized = resize_if_needed(he_with_gem_overlay)
                cv2.imwrite(os.path.join(images_dir, '3b_tissue_segmentation_gem_heatmap.png'),
                           cv2.cvtColor(gem_heatmap_overlay_resized, cv2.COLOR_RGB2BGR), png_compression)
                print(f"[{self.mode} Mode] 3b_tissue_segmentation_gem_heatmap.png saved (30% GEM on {bg_label}, white background) - Used in Interactive Viewer")

                # ========================================
                # 3c. Pure GEM Heatmap (mask region only, no background)
                # ========================================
                # Regenerate pure GEM heatmap from original GEM data (not HE-processed)
                # Use gem_enhanced or original gem
                gem_source = self.gem_enhanced if hasattr(self, 'gem_enhanced') and self.gem_enhanced is not None else self.gem

                # Convert to grayscale if needed
                if gem_source.ndim == 3:
                    gem_gray_pure = cv2.cvtColor(gem_source, cv2.COLOR_RGB2GRAY)
                else:
                    gem_gray_pure = gem_source.copy()

                # Use pure GEM mask (not union_mask) to show only GEM expression region
                gem_mask_bool = self.gem_mask > 0 if hasattr(self, 'gem_mask') and self.gem_mask is not None else np.ones_like(gem_gray_pure, dtype=bool)

                # Normalize GEM data within GEM mask region only
                if gem_mask_bool.any():
                    gem_gray_masked = gem_gray_pure.copy()
                    mask_region_data = gem_gray_masked[gem_mask_bool]
                    if mask_region_data.max() > mask_region_data.min():
                        # Normalize only the mask region data
                        gem_normalized_pure = np.zeros_like(gem_gray_masked, dtype=np.uint8)
                        gem_normalized_pure[gem_mask_bool] = ((mask_region_data - mask_region_data.min()) /
                                                              (mask_region_data.max() - mask_region_data.min()) * 255).astype(np.uint8)
                    else:
                        gem_normalized_pure = gem_gray_masked.astype(np.uint8)

                    # Apply jet colormap to generate heatmap
                    import matplotlib.pyplot as plt
                    gem_heatmap_pure_rgb = plt.cm.jet(gem_normalized_pure / 255.0)[:, :, :3]  # RGB
                    gem_heatmap_pure_rgb = (gem_heatmap_pure_rgb * 255).astype(np.uint8)

                    # Create final image with white background
                    gem_heatmap_pure = np.full((*gem_heatmap_pure_rgb.shape[:2], 3), 255, dtype=np.uint8)
                    gem_heatmap_pure[gem_mask_bool] = gem_heatmap_pure_rgb[gem_mask_bool]
                else:
                    # No mask data, create white image
                    gem_heatmap_pure = np.full((*gem_gray_pure.shape[:2], 3), 255, dtype=np.uint8)

                gem_heatmap_pure_resized = resize_if_needed(gem_heatmap_pure)
                cv2.imwrite(os.path.join(images_dir, '3c_gem_heatmap_only.png'),
                           cv2.cvtColor(gem_heatmap_pure_resized, cv2.COLOR_RGB2BGR), png_compression)
                print(f"[{self.mode} Mode] 3c_gem_heatmap_only.png saved (pure GEM expression heatmap, white background)")

                # ========================================
                # 4. Segmentation region visualization (filename varies by mode)
                # ========================================
                if self.mode == 'HE':
                    seg_vis_filename = '4_tissue_segmentation_on_he.jpg'
                else:  # ssDNA
                    seg_vis_filename = '4_tissue_segmentation_on_ssdna.jpg'

                he_with_seg = he_normalized.copy()
                # Mark valid tissue regions with semi-transparent green
                green_overlay = np.zeros_like(he_with_seg)
                green_overlay[:, :, 1] = 180  # Green channel
                he_with_seg[union_mask > 0] = cv2.addWeighted(
                    he_with_seg[union_mask > 0], 0.7,
                    green_overlay[union_mask > 0], 0.3, 0
                )

                # Draw union contour
                contours, _ = cv2.findContours(union_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(he_with_seg, contours, -1, (0, 255, 0), 2)  # Green contour

                he_with_seg_resized = resize_if_needed(he_with_seg)
                cv2.imwrite(os.path.join(images_dir, seg_vis_filename),
                           cv2.cvtColor(he_with_seg_resized, cv2.COLOR_RGB2BGR), jpeg_quality)
                print(f"[{self.mode} Mode] {seg_vis_filename} saved")

                # ========================================
                # 5. GEM Heatmap overlay (full image, filename varies by mode)
                # ========================================
                if self.mode == 'HE':
                    overlay_filename = '5_gem_heatmap_on_he.png'
                else:  # ssDNA
                    overlay_filename = '5_gem_heatmap_on_ssdna.png'

                overlay_resized = resize_if_needed(overlay_alpha)
                cv2.imwrite(os.path.join(images_dir, overlay_filename),
                           cv2.cvtColor(overlay_resized, cv2.COLOR_RGB2BGR), png_compression)
                print(f"[{self.mode} Mode] {overlay_filename} saved (40% GEM on {bg_label})")

                # Calculate statistics
                gem_area = np.sum(self.gem_mask > 0)
                bg_area = np.sum(self.he_image_mask > 0)
                union_area = np.sum(union_mask > 0)

                print(f"[Tissue Segmentation] GEM mask area: {gem_area} pixels")
                print(f"[Tissue Segmentation] {bg_label} mask area: {bg_area} pixels")
                print(f"[Tissue Segmentation] Union area: {union_area} pixels")
                print(f"[Tissue Segmentation] Coverage: {union_area} pixels total from GEM + {bg_label}")

        elif self.mode == 'gene_expr':
            # gene_expr mode: only GEM, no background
            if hasattr(self, 'gem_mask'):
                # Only save GEM mask
                mask_resized = resize_if_needed(self.gem_mask.astype(np.uint8) * 255)
                cv2.imwrite(os.path.join(images_dir, '3_tissue_segmentation_mask.png'), mask_resized, png_compression)
                print(f"[{self.mode} Mode] 3_tissue_segmentation_mask.png saved (GEM only)")

                gem_area = np.sum(self.gem_mask > 0)
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
            print("  - 5_gem_heatmap_on_he.png: GEM heatmap overlaid on HE (full image, 40% transparency)")
        elif self.mode == 'ssDNA':
            print("  - 2_ssdna_registered.jpg: ssDNA image (JPEG compressed)")
            print("  - 3_tissue_segmentation_mask.png: Tissue segmentation mask (GEM∪ssDNA union)")
            print("  - 3b_tissue_segmentation_gem_heatmap.png: GEM heatmap overlaid on ssDNA (30% transparency, mask region only)")
            print("  - 3c_gem_heatmap_only.png: Pure GEM heatmap (mask region only, transparent background)")
            print("  - 4_tissue_segmentation_on_ssdna.jpg: Segmentation region overlaid on ssDNA")
            print("  - 5_gem_heatmap_on_ssdna.png: GEM heatmap overlaid on ssDNA (full image, 40% transparency)")
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
            self.slice_registration(self.model, self.method,
                                   enhance_method=self.enhance_method,
                                   enhance_params=self.enhance_params)
            self.splitplot_experiment(self.segment_type, self.method)
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
                # Check for common image directories (images, rawdata, etc.)
                for subdir in ['images', 'rawdata', 'image', 'data']:
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
                    for subdir in ['images', 'rawdata', 'image', 'data']:
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
        elif args.method == 'image' and args.segment:
            print(f"Warning: HE image not found. Searched for {args.sample}_he.(png|jpg|jpeg) in {len(possible_locations)} locations")

    if args.rectify:
        assert os.path.exists(args.tif), "HE image is required for rectification! Please check the path."
    if args.segment:
        segment_required_args = ['input', 'outdir', 'model']
        # Both image and HE modes require tif file
        if args.method in ['image', 'HE']:
            segment_required_args.append('tif')

        if not all(getattr(args, arg, None) for arg in segment_required_args):
            missing_args = [arg for arg in segment_required_args if not getattr(args, arg, None)]
            raise ValueError(f'The following arguments are required when using --segment: {", ".join(missing_args)} '
                             f'(method: {args.method})')
        else:
            # Both image and HE modes need to verify tif file exists
            if args.tif and not os.path.exists(args.tif) and args.method in ['image', 'HE']:
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
        '--bs_out',
        help="""We specify this folder to store images segmented by the our pre-trained segmentation model,
in order to check the segmentation effect. If not specified, the default folder is the same as the output folder.
If the parameter bs_out is correctly set, a folder named 'binSegment' will be automatically created
at the corresponding path.""",
        type=str,
        default=None
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
    subparser_bin_segment.add_argument('--model', type=str, default='/home/zhoumy/lizt/code/work_space/spatial_bin/swin_tiny.pth', help='path to model')
    subparser_bin_segment.add_argument('--thread', type=int, default=8, help='number of threads')
    subparser_bin_segment.add_argument('--pixel-size', type=float, default=0.5, help='pixel size of microscope image')
    subparser_bin_segment.add_argument('--sample', type=str, default='OST110014', help='tissue name')
    subparser_bin_segment.add_argument('--prompt', type=str, default='[[5168,2928],[6384,4416]]', help='tissue position prompt')
    subparser_bin_segment.add_argument('--genomeDir', type=str, default='/mnt/strna/work_project/rawdata/celatlas_spatial/reference/Mus_musculus', help='reference genome type')
    subparser_bin_segment.add_argument('--input', type=str, default='/mnt/chip_mapping/', help='path to chip mapping')
    subparser_bin_segment.add_argument('--outdir', type=str, default='/mnt/strna/work_project/pipeline/celatlas_spatial/OST110014/06.binSegment', help='output path')
    subparser_bin_segment.add_argument('--bs_out', type=str, default='/mnt/strna/work_project/rawdata/stomics/binSegment/', help='path to binSegment')
    # subparser_bin_segment.add_argument('--tif', type=str, default='/mnt/strna/work_project/rawdata/celatlas_spatial/images/ST110121_C1.tif', help='tissue microscope image path')
    subparser_bin_segment.add_argument('--tif', type=str, default='/mnt/strna/work_project/rawdata/stomics/binSegment/STProtein_test/protein_label_region.tif', help='tissue microscope image path')
    subparser_bin_segment.add_argument('--count_detail', type=str, default='/mnt/strna/work_project/pipeline/celatlas_spatial/OST110014/05.count/OST110014_count_detail.txt', help='reads path')
    subparser_bin_segment.add_argument('--segment-type', type=str, default='tissue', help='run segmentation mode, tissue or cell')
    subparser_bin_segment.add_argument(
        '--method',
        type=str,
        default='gene_expr',
        choices=['gene_expr', 'image', 'HE'],
        help="""Segmentation method:
  - 'gene_expr': Gene expression only (no images)
  - 'image': Image-based segmentation (ssDNA images)
  - 'HE': Gene expression + HE image registration (requires *_he.tif/png/jpg)"""
    )
    subparser_bin_segment.add_argument('--rectify', action='store_true', help='rectify tissue image by HE image roi')
    subparser_bin_segment.add_argument('--segment', action='store_true', help='segment tissue image')
    subparser_bin_segment.add_argument('--count', action='store_true', help='gene count')
    subparser_bin_segment.add_argument('--extend', action='store_true', help='extend barcode from 12 to 18')
    subparser_bin_segment.add_argument('--debug', action='store_true', help='debug mode')
    subparser_bin_segment.add_argument('--omics', type=str, default='rna', help='omics type')
    args = parser.parse_args()

    args.segment = True  # image segmentation test
    args.count = True  # count test
    args.rectify = False  # rectify test
    args.subparser_assay = 'rna'
    args.tif = None

    binSegment(args)  # process tissue count


if __name__ == '__main__':
    main()
