import argparse
import importlib
import os
import sys

from celatlas_spatial.__init__ import __VERSION__, ASSAY_LIST


class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def find_assay_init(assay):
    return importlib.import_module(f"celatlas_spatial.{assay}.__init__")


def find_step_module(assay, step):
    root_path = os.path.dirname(__file__)
    init_module = find_assay_init(assay)
    assay_step_path = os.path.join(root_path, assay, f"{step}.py")
    tools_step_path = os.path.join(root_path, "tools", f"spatial_{step}.py")
    if os.path.exists(assay_step_path):
        return importlib.import_module(f"celatlas_spatial.{assay}.{step}")
    if hasattr(init_module, "IMPORT_DICT") and step in init_module.IMPORT_DICT:
        return importlib.import_module(f"{init_module.IMPORT_DICT[step]}.{step}")
    if os.path.exists(tools_step_path):
        return importlib.import_module(f"celatlas_spatial.tools.spatial_{step}")
    raise ModuleNotFoundError(f"No module found for {assay}.{step}")


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

    from celatlas_spatial.runner.cli import add_runner_subparsers

    add_runner_subparsers(subparsers)

    requested_assay = None
    requested_step = None
    argv = sys.argv[1:]
    if len(argv) >= 1 and not argv[0].startswith("-"):
        requested_assay = argv[0]
    if len(argv) >= 2 and not argv[1].startswith("-"):
        requested_step = argv[1]

    for assay in ASSAY_LIST:
        subparser_1st = subparsers.add_parser(assay)
        subparser_2nd = subparser_1st.add_subparsers()

        init_module = find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__
        if requested_assay == assay and requested_step in __STEPS__:
            steps_to_register = [requested_step]
            placeholder_steps = []
        elif requested_assay == assay and requested_step is None:
            steps_to_register = []
            placeholder_steps = __STEPS__
        elif requested_assay and requested_assay != assay:
            steps_to_register = []
            placeholder_steps = []
        elif requested_assay == assay and requested_step not in __STEPS__:
            steps_to_register = []
            placeholder_steps = __STEPS__
        else:
            steps_to_register = []
            placeholder_steps = __STEPS__

        for step in placeholder_steps:
            subparser_2nd.add_parser(step, formatter_class=ArgFormatter)

        for step in steps_to_register:
            step_module = find_step_module(assay, step)
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
