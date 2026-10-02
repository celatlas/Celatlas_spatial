import os
import pathlib
import re

from collections import defaultdict

import pandas as pd

from celatlas_spatial.tools import utils
from celatlas_spatial.rna.mkref import Mkref_rna
from celatlas_spatial.tools.step import Step, s_common
from celatlas_spatial.__init__ import HELP_DICT


GENE_ID_PATTERN = re.compile(r'(?:^|[;\s])gene_id(?:\s+|=)"?[^;\s"]+')
STRAND_PARAM_PATTERN = re.compile(r'(^|\s)(-s|--strandSpecific)(\s|$)')


def validate_featurecounts_gtf(gtf):
    """Reject annotations that cannot support Celatlas exon/gene metrics."""

    path = pathlib.Path(gtf)
    if not path.is_file():
        raise ValueError(f'featureCounts annotation file does not exist: {path}')

    feature_counts = {'exon': 0, 'gene': 0}
    gene_id_counts = {'exon': 0, 'gene': 0}
    with path.open(errors='replace') as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] not in feature_counts:
                continue
            feature_type = fields[2]
            feature_counts[feature_type] += 1
            if GENE_ID_PATTERN.search(fields[8]):
                gene_id_counts[feature_type] += 1

    problems = []
    for feature_type in ('exon', 'gene'):
        if feature_counts[feature_type] == 0:
            problems.append(f'no {feature_type} records')
        elif gene_id_counts[feature_type] == 0:
            problems.append(f'{feature_type} records have no gene_id attribute')
    if problems:
        detail = '; '.join(problems)
        raise ValueError(
            f'Invalid GTF for Celatlas featureCounts: {path}: {detail}. '
            'Normalize the source GFF3/GTF so exon and gene records contain gene_id.'
        )


def _effective_featurecounts_param(featurecounts_param: str | None) -> str:
    """
    Return a single featureCounts parameter string with exactly one strand setting.

    Celatlas defaults to stranded exon counting for legacy direct invocations, but
    workflow-specific calls may already inject `-s` or `--strandSpecific`. In that
    case we must not prepend the default again.
    """

    param = (featurecounts_param or "").strip()
    if not param:
        return "-s 1"
    if STRAND_PARAM_PATTERN.search(param):
        return param
    return f"-s 1 {param}"


class FeatureCounts(Step):
    """
    ## Features
    - Assigning uniquely mapped reads to genomic features with FeatureCounts.
    ## Output
    - `{sample}` Numbers of reads assigned to features (or meta-features).
    - `{sample}_summary` Stat info for the overall summrization results, including number of
    successfully assigned reads and number of reads that failed to be assigned due to
    various reasons (these reasons are included in the stat info).
    - `{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam` featureCounts output BAM,
    sorted by coordinates;BAM file contains tags as following(Software Version>=1.1.8):
        - CB cell barcode
        - UB UMI
        - GN gene name
        - GX gene id
    - `{sample}_name_sorted.bam` featureCounts output BAM, sorted by read name.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.thread = 64 if self.thread > 64 else self.thread

        # set
        self.gtf = Mkref_rna.parse_genomeDir(self.args.genomeDir)['gtf']
        validate_featurecounts_gtf(self.gtf)
        self.featureCounts_param = args.featureCounts_param

        # gtf_type
        self.gtf_types = ['exon', 'gene']

        # stats
        self.feature_log_dict = defaultdict(dict)

        # out
        input_basename = os.path.basename(self.args.input)
        self.featureCounts_bam = f'{self.outdir}/{input_basename}.featureCounts.bam'
        self.nameSorted_bam = f'{self.outdir}/{self.sample}_nameSorted.bam'

    @staticmethod
    def read_log(log_file):
        """
        Args:
            log_file: featureCounts log summary file
        Returns:
            log_dict: {'Assigned': 123, ...}
        """
        # skip first line
        df = pd.read_csv(log_file, sep='\t', header=None, names=['name', 'value'], skiprows=1)
        log_dict = df.set_index('name')['value'].to_dict()
        return log_dict

    @utils.add_log
    def run_featureCounts(self, outdir, gtf_type):
        featurecounts_param = _effective_featurecounts_param(self.featureCounts_param)
        cmd = (
            'featureCounts '
            f'{featurecounts_param} '
            f'--largestOverlap '
            f'-M '
            f'-a {self.gtf} '
            f'-o {outdir}/{self.sample} '
            '-R BAM '
            f'-T {self.thread} '
            f'-t {gtf_type} '
            f'{self.args.input} '
            '2>&1 '
        )
        self.debug_subprocess_call(cmd)

    def run_get_log(self):
        tmp_dir = f'{self.outdir}/tmp/'
        for gtf_type in self.gtf_types:
            outdir = f'{tmp_dir}/{gtf_type}'
            pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
            log_file = f'{outdir}/{self.sample}.summary'

            self.run_featureCounts(outdir, gtf_type)
            self.feature_log_dict[gtf_type] = FeatureCounts.read_log(log_file)

            if gtf_type == self.args.gtf_type:
                cmd = f'mv {outdir}/* {self.outdir}'
                self.debug_subprocess_call(cmd)

        cmd = f'rm -r {tmp_dir}'
        self.debug_subprocess_call(cmd)

    def run(self):
        self.run_get_log()
        self.add_metrics()
        utils.sort_bam(self.featureCounts_bam, self.nameSorted_bam, threads=self.thread, by='name')
        self.remove_temp_file()

    def remove_temp_file(self):
        os.remove(self.featureCounts_bam)

    @utils.add_log
    def add_metrics(self):
        total = sum(self.feature_log_dict['exon'].values())

        Assigned_exon = self.feature_log_dict['exon']['Assigned']
        Assigned_intergenic = self.feature_log_dict['gene']['Unassigned_NoFeatures']
        """
        https://academic.oup.com/nargab/article/2/3/lqaa073/5910008
        Approximately 15% of genes had exon counts that were greater than genebody counts (by a median value of eight counts). 
        This was due to our conservative approach of excluding reads that overlapped features in multiple genes during 
        the read summarization step by featureCounts using the argument allowMultiOverlap=FALSE. Under this strategy, 
        some reads were counted towards the exon count set but not the genebody count set. This happens when a read 
        overlaps the exon in one gene and the intron of another gene—it is counted towards exon counts but not genebody counts 
        due to its overlap of multiple genebodies but not multiple exons.
        """
        Unassigned_ambiguity = self.feature_log_dict['exon']['Unassigned_Ambiguity']
        Assigned_intron = total - Assigned_exon - Assigned_intergenic - Unassigned_ambiguity

        self.add_metric(
            name='Feature Type',
            value=self.args.gtf_type.capitalize(),
            help_info='Specified by `--gtf_type`. For snRNA-seq, you need to add `--gtf_type gene` to include reads mapped to intronic regions. Staring from Celatlas-spatial v0.1.3, the default value of gtf_type is changed from `exon` to `gene`.'
        )
        self.add_metric(
            name='Reads Assigned To Exonic Regions',
            value=Assigned_exon,
            total=total,
            help_info='Reads that can be successfully assigned to exonic regions'
        )
        self.add_metric(
            name='Reads Assigned To Intronic Regions',
            value=Assigned_intron,
            total=total,
            help_info='Reads that can be successfully assigned to intronic regions'
        )
        self.add_metric(
            name='Reads Assigned To Intergenic Regions',
            value=Assigned_intergenic,
            total=total,
            help_info='Reads that can be successfully assigned to intergenic regions'
        )
        self.add_metric(
            name='Reads Unassigned Ambiguity',
            value=Unassigned_ambiguity,
            total=total,
            help_info='Alignments that overlap two or more features'
        )


@utils.add_log
def featureCounts(args):
    with FeatureCounts(args) as runner:
        runner.run()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument(
        '--gtf_type',
        help='Specify feature type in GTF annotation',
        default='gene',
        choices=['exon', 'gene'],
    )
    parser.add_argument('--genomeDir', help=HELP_DICT['genomeDir'])
    parser.add_argument('--featureCounts_param', help=HELP_DICT['additional_param'], default="")

    if sub_program:
        parser.add_argument('--input', help='Required. BAM file path.', required=True)
        parser = s_common(parser)
    return parser
