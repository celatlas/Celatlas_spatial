import argparse

from celatlas_spatial.tools import utils

from celatlas_spatial.__init__ import __VERSION__, ASSAY_LIST

class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main():

    description = '''
Celatlas Spatial Transcriptomics Analysis Pipeline

Version: {}

Usage Examples:
  # RNA analysis steps:
  celatlas_spatial rna mkref --help
  celatlas_spatial rna sample --help
  celatlas_spatial rna barcode --help
  celatlas_spatial rna analysis --help

For complete workflow, use the Shell scripts:
  bash Celatlas.sh --help
  bash Celatlas_reanalysis.sh --help
'''.format(__VERSION__)

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=ArgFormatter,
        epilog='For detailed usage, see README.md or run: bash Celatlas.sh --help'
    )

    parser.add_argument('-v', '--version', action='version', version=__VERSION__)

    subparsers = parser.add_subparsers(
        dest='subparser_assay',
        title='Available assays',
        description='Use "celatlas_spatial {assay} {step} --help" for step-specific help'
    )

    for assay in ASSAY_LIST:

        subparser_1st = subparsers.add_parser(assay)

        subparser_2nd = subparser_1st.add_subparsers()

        init_module = utils.find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__

        for step in __STEPS__:
            step_module = utils.find_step_module(assay, step)
            func = getattr(step_module, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            parser_step = subparser_2nd.add_parser(step, formatter_class=ArgFormatter)
            func_opts(parser_step, sub_program=True)
            parser_step.set_defaults(func=func)

    args = parser.parse_args()
    if len(args.__dict__) <= 1:
        parser.print_help()
        parser.exit()
    else:
        args.func(args)


if __name__ == '__main__':
    main()
