from celatlas_spatial.tools import utils
from celatlas_spatial.__init__ import __VERSION__
from celatlas_spatial.tools.__init__ import PATTERN_DICT
from celatlas_spatial.tools.spatial_barcode import Chemistry
from celatlas_spatial.tools.step import Step, s_common


def add_kit_version(chemistry):
    kit_dict = {
        '0': 'testing version',  # BBV0: scRNA only testing version
        '1': 'no longer in use',
        '2': 'kit V1',
        '3': 'kit V2',
    }
    if chemistry.startswith('strnaV'):
        s = chemistry.replace('strnaV', '')
        chem_version = s[0]
        if chem_version in kit_dict:
            kit = kit_dict[chem_version]
            chemistry = f'{chemistry} ({kit})'
    elif chemistry.startswith('BBV'):
        s = chemistry.replace('BBV', '')
        chem_version = s[0]
        if chem_version in kit_dict:
            kit = kit_dict[chem_version]
            chemistry = f'{chemistry} ({kit})'

    return chemistry


class Sample(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.version = __VERSION__
        self.chemistry = args.chemistry
        
        # Get sample name and chip number from environment variables
        import os
        self.sample_name = os.environ.get('CELATLAS_SAMPLE_NAME', '')
        self.chip_number = os.environ.get('CELATLAS_CHIP_NUMBER', '')

    @utils.add_log
    def run(self):
        if self.chemistry == 'auto':
            fq1 = self.args.fq1
            ch = Chemistry(fq1)
            chemistry = ch.check_chemistry()
            chemistry = ",".join(set(chemistry))
        else:
            chemistry = self.chemistry

        # Display logic for Sample ID
        if self.sample_name:
            # If sample_name is provided, show it as Sample ID
            sample_id_display = self.sample_name
        else:
            # If only chip_number is provided, show as N/A
            sample_id_display = "N/A"

        self.add_metric(
            name='Sample ID',
            value=sample_id_display,
        )
        self.add_metric(
            name='Chemistry',
            value=chemistry,
            display=add_kit_version(chemistry),
            help_info='For more information, see...',
        )
        
        self.add_metric(
            name='Chip Number',
            value=self.chip_number if self.chip_number else "N/A",
        )
        
        # Determine Transcriptome based on species from args
        species = getattr(self.args, 'Species', getattr(self.args, 'species', ''))
        transcriptome_mapping = {
            'Mus_musculus': 'GRCm39',
            'mouse': 'GRCm39',
            'Homo_sapiens': 'GRCh38', 
            'human': 'GRCh38'
        }
        transcriptome = transcriptome_mapping.get(species, 'N/A')
        
        self.add_metric(
            name='Transcriptome',
            value=transcriptome,
        )
        
        self.add_metric(
            name='Pipeline Version',
            value=self.version,
        )
        self.add_metric(
            name='Image Alignment',
            value="ssDNA",  
        )
        
        self.add_metric(
            name='Probe Set Name',
            value="N/A",
        )
        self.add_metric(
            name='Filter Probes',
            value="N/A",
        )


@utils.add_log
def sample(args):
    with Sample(args) as runner:
        runner.run()

def get_opts_sample(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='read1 fq file')
    parser.add_argument('--chemistry', choices=list(PATTERN_DICT.keys()), help='chemistry version', default='auto')
    return parser
