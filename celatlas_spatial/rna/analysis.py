import base64
import argparse
import os

import tifffile

from celatlas_spatial.tools import analysis_wrapper, utils
from celatlas_spatial.celatlas import ArgFormatter
from celatlas_spatial.tools.step import Step
from celatlas_spatial.__init__ import __VERSION__


class Analysis(Step):
    """
    ## Features
    - Cell clustering with Seurat.

    - Calculate the marker gene of each cluster.

    - Cell type annotation(optional). You can provide markers of known cell types and annotate cell types for each cluster.

    ## Output
    - `markers.tsv` Marker genes of each cluster.

    - `tsne_coord.tsv` t-SNE coordinates and clustering information.

    - `{sample}/07.analsis/{sample}_auto_assign/` This result will only be obtained when `--type_marker_tsv`
    parameter is provided. The result contains 3 files:
        - `{sample}_auto_cluster_type.tsv` The cell type of each cluster; if cell_type is "NA",
    it means that the given marker is not enough to identify the cluster.
        - `{sample}_png/{cluster}_pctdiff.png` Percentage of marker gene expression in this cluster - percentage in all other clusters.
        - `{sample}_png/{cluster}_logfc.png` log2 (average expression of marker gene in this cluster / average expression in all other clusters + 1)
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)

        # Detect assay type
        self.assay = getattr(args, 'assay', None)
        if self.assay is None:
            # Auto-detect based on parameters
            self.assay = 'scrna' if args.square_bin_dir is None else 'spatial'

        if self.assay == 'spatial':
            # Spatial mode: setup paths for spatial data
            self.args.micron_bin = self._micron_bin
            self._square_bin_path = self.args.square_bin_dir
            self._image_path = os.path.join(os.path.dirname(os.path.normpath(self._square_bin_path)), 'images')
            self._matrix = 'matrix'
            self._spatial = 'spatial'

            self.args.matrix_file, self.args.position_file = self.generate_file_paths(
                sample=self.sample,
                square_bin_path=self._square_bin_path,
                matrix=self._matrix,
                spatial=self._spatial,
                micron_bins=self._micron_bin,
                use_raw=False,
            )
            self.args.im_shape = tifffile.imread(os.path.join(self._image_path, f'{self.sample}_regist.tif')).shape
            self.args.hires_image = os.path.join(self._image_path, 'tissue_hires_image.png')

            self.bioinfo_data = os.path.join(self.outdir, 'Bioinfodata')
        else:
            # scRNA mode: no spatial data needed
            self._square_bin_path = None
            self._image_path = None
            self.bioinfo_data = None

        self.display_title = display_title
        self._scatter = 'scatter.png'
        self._violin = 'violin.png'
        self._heatmap = 'spatial_gene_expression_distribution.png'
        self._spatial_cluster = 'spatial_cluster.png'
        self._umap = 'umap_cluster.png'
        self._table_id = 'marker_genes'

    def _render_html(self):
        """Override parent's _render_html to skip HTML generation.
        Analysis step doesn't generate reports - that's done by the separate report step.
        """
        pass

    @staticmethod
    def generate_file_paths(sample: str, square_bin_path: str, matrix: str, spatial: str, micron_bins: list[int], use_raw: bool) -> \
            tuple[dict[str, tuple[str, str]], dict[str, tuple[str, str]]]:
        """
        Generate file paths for matrix and spatial data for each bin size.  Also include raw data if it exists.
        :param sample: sample name
        :param square_bin_path: path to the square bin directory
        :param matrix: string for matrix data
        :param spatial: string for spatial data
        :param micron_bins: list of bin sizes
        :param use_raw: boolean to include raw data
        :return: dictionary of matrix file paths and dictionary of spatial file paths
        """
        def path_builder(bin_name: str) -> tuple[str, str]:
            base_path = os.path.join(square_bin_path, f'{sample}_{bin_name}')
            return (
                os.path.join(base_path, matrix),
                os.path.join(base_path, spatial, f'{sample}_{bin_name}_tissue_positions.tsv')
            )

        matrix_files = {}
        position_files = {}

        for bin_info in [('bin' + str(i), f'Bin{i}') for i in micron_bins]:
            matrix_files[bin_info[0]], position_files[bin_info[0]] = path_builder(bin_info[1])

        if use_raw:
            raw_path = os.path.join(square_bin_path, f'{sample}_Raw')
            if os.path.exists(raw_path):
                matrix_files['raw'], position_files['raw'] = path_builder('Raw')

        return matrix_files, position_files

    def run(self):
        with analysis_wrapper.Scanpy_wrapper(self.args, display_title=self.display_title) as scanpy_wrapper:
            scanpy_wrapper.run()
            self.set_metric_list(metric_list=scanpy_wrapper.get_metric_list())

        # scRNA mode: basic analysis only, no spatial visualizations
        if self.assay == 'scrna':
            return

        # Spatial mode: continue with spatial-specific processing
        with analysis_wrapper.Report_runner(self.args, display_title=self.display_title) as report_runner:
            report_runner.add_marker_help()

        df_tsne, df_marker = scanpy_wrapper.get_df()

        b_value = self.args.bin
        heatmap_path = os.path.join(self.bioinfo_data, '02.SpatialQC', f'{self.sample}_{self._heatmap}')
        spatial_cluster_path = os.path.join(self.bioinfo_data, '03.SpatialCluster', f'{self.sample}_{self._spatial_cluster}')
        umap_path = os.path.join(self.bioinfo_data, '03.SpatialCluster', f'{self.sample}_{self._umap}')

        heatmap = base64.b64encode(open(heatmap_path, 'rb').read()).decode('ascii')
        spatial_cluster = base64.b64encode(open(spatial_cluster_path, 'rb').read()).decode('ascii')
        umap = base64.b64encode(open(umap_path, 'rb').read()).decode('ascii')

        self.add_metric(name='Bin', value=b_value, help_info=f'Bin size {b_value} μm')

        self.add_data(heatmap=heatmap)
        self.add_data(spatial_cluster=spatial_cluster)
        self.add_data(umap=umap)

        scatter, violin = {}, {}
        for b in self._micron_bin:
            key = f'bin{b}'
            bin_scatter = os.path.join(self.bioinfo_data, '02.SpatialQC', f'{self.sample}_{key}_{self._scatter}')
            bin_violin = os.path.join(self.bioinfo_data, '02.SpatialQC', f'{self.sample}_{key}_{self._violin}')
            scatter[key] = base64.b64encode(open(bin_scatter, 'rb').read()).decode('ascii')
            violin[key] = base64.b64encode(open(bin_violin, 'rb').read()).decode('ascii')

        self.add_data(scatter=scatter)
        self.add_data(violin=violin)

        table_dict = self.get_table_dict(
            title='Marker Genes by Cluster',
            table_id=self._table_id,
            df_table=df_marker,
        )
        self.add_data(table_dict=table_dict)

@utils.add_log
def analysis(args):
    with Analysis(args, display_title='Analysis') as runner:
        runner.run()


def get_opts_analysis(parser, sub_program):
    analysis_wrapper.get_opts_analysis(parser, sub_program)


def main():
    parser = argparse.ArgumentParser(description='Celatlas Spatial', formatter_class=ArgFormatter)
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers(dest='subparser_assay')
    subparser_1st = subparsers.add_parser('rna')
    subparser_2nd = subparser_1st.add_subparsers()
    subparser_ba = subparser_2nd.add_parser('Bioinformatics_analysis', formatter_class=ArgFormatter)
    subparser_ba.add_argument('--thread', type=int, default=128, help='number of threads')
    subparser_ba.add_argument('--sample', type=str, default='ST110121_C1', help='tissue name')
    subparser_ba.add_argument('--outdir', type=str, default='/mnt/strna/work_project/pipeline/celatlas_spatial/ST110121_C1/07_1.analysis', help='output directory')
    subparser_ba.add_argument('--genomeDir', type=str, default='/mnt/strna/work_project/rawdata/celatlas_spatial/reference/Mus_musculus', help='reference genome type')
    subparser_ba.add_argument('--square_bin_dir', type=str, default='/mnt/strna/work_project/pipeline/celatlas_spatial/ST110121_C1/06_1.binSegment/square_bin', help='square bin directory')
    subparser_ba.add_argument('--pixel-size', type=float, default=0.72, help='x.x um per pixel')
    subparser_ba.add_argument('--bin', type=int, default=100, help='bin size(micron)')
    subparser_ba.add_argument('--debug', action='store_true', help='debug mode')
    args = subparser_ba.parse_args()
    args.subparser_assay = 'rna'

    analysis(args)  # Run the analysis


if __name__ == '__main__':
    main()
