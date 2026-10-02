import os
import abc
import sys
import json
import numbers
import subprocess
from datetime import datetime

from celatlas_spatial.tools import utils
from celatlas_spatial.__init__ import HELP_DICT


def cap_str_except_preposition(my_string):
    prepositions = {"and", "or", "the", "a", "of", "in", "per", "after", 'with'}
    lowercase_words = my_string.split(" ")

    final_words = [word if word in prepositions else word[0].upper() + word[1:] for word in lowercase_words]
    final_words = " ".join(final_words)
    return final_words

def s_common(parser):
    """
    subparser common arguments
    """
    parser.add_argument('--outdir', help='Output diretory.', required=True)
    parser.add_argument('--sample', help='Sample name.', required=True)
    parser.add_argument('--thread', help=HELP_DICT['thread'], default=4)
    parser.add_argument('--debug', help=HELP_DICT['debug'], action='store_true')
    return parser


class Step:
    def __init__(self, args, display_title=None):

        print(f'Args: {args}')
        self.args = args  
        self.outdir = args.outdir  
        self.sample = args.sample  
        self.assay = args.subparser_assay  

        self.thread = int(args.thread) if args.thread != '-1' else os.cpu_count()
        self.debug = args.debug  
        self.out_prefix = f'{self.outdir}/{self.sample}'  
        self.display_title = display_title 
        self._micron_bin = [10, 20, 50, 100]  
        self._started_at = datetime.now()
        self._exit_info = (None, None, None)

        # important! make outdir before path_dict because path_dict use relative path.
        utils.check_mkdir(self.outdir)

        # set
        class_name = self.__class__.__name__
        if not display_title:
            self._display_title = class_name
        else:
            self._display_title = display_title
        self._step_name = class_name[0].lower() + class_name[1:]
        self.__slots = ['data', 'metrics']
        self._step_summary_name = f'{self._step_name}_summary'

        self.__metric_list = []
        self.__help_content = []
        self.__commands = []
        self._path_dict = {}
        for slot in self.__slots:
            self._path_dict[slot] = f'{self.outdir}/../.{slot}.json'

        self.__content_dict = {}
        for slot, path in self._path_dict.items():
            if not os.path.exists(path):
                self.__content_dict[slot] = {}
            else:
                with open(path) as f:
                    try:
                        self.__content_dict[slot] = json.load(f)
                    except ValueError:
                        print(f'WARNING: Decoding "{path}" as json has failed. Will create empty json file.')
                        self.__content_dict[slot] = {}
            # clear step_summary
            self.__content_dict[slot][self._step_summary_name] = {}

        # out file
        self.__step_log_file = f'{self.outdir}/step.log'

    def add_metric(self, name, value, total=None, help_info=None, display=None, show=True, print_log=True):
        """
        add metric to metric_list

        Args
            total: int or float, used to calculate fraction
            help_info: str, help info for metric in html report
            display: str, controls how to display the metric in HTML report.
            show: bool, whether to add to `.data.json` and the visible Metrics section in `step.log`.
            print_log: bool, whether to print metric to stdout
        """

        name = cap_str_except_preposition(name)
        if help_info:
            help_info = help_info[0].upper() + help_info[1:]
            if help_info[-1] != '.':
                help_info += '.'
        if not display:
            if isinstance(value, numbers.Number):
                display = str(format(value, ','))
            else:
                display = value
        fraction = None
        if total:
            fraction = round(value / total * 100, 2)
            display += f'({fraction}%)'
        self.__metric_list.append(
            {
                "name": name,
                "value": value,
                "total": total,
                "fraction": fraction,
                "display": display,
                "help_info": help_info,
                "show": show,
            }
        )

        if print_log:
            print(f'{name}: {display}')

    def _dump_content(self):
        """
        dump content to json file
        """
        for slot, path in self._path_dict.items():
            if self.__content_dict[slot]:
                with open(path, 'w') as f:
                    json.dump(self.__content_dict[slot], f, indent=4)

    @staticmethod
    def _format_log_value(value, max_length=500):
        if isinstance(value, str):
            text = value
        else:
            text = repr(value)

        text = text.replace('\n', '\\n')
        if len(text) > max_length:
            text = f"{text[:max_length]}... <truncated, {len(text)} chars>"
        return text

    @staticmethod
    def _format_file_size(size):
        units = ["B", "KB", "MB", "GB", "TB"]
        value = float(size)
        for unit in units:
            if value < 1024 or unit == units[-1]:
                if unit == "B":
                    return f"{int(value)} {unit}"
                return f"{value:.1f} {unit}"
            value /= 1024

    def _write_step_log(self):
        completed_at = datetime.now()
        elapsed = completed_at - self._started_at
        exc_type, exc_value, _traceback = self._exit_info
        status = "failed" if exc_type else "completed"
        legacy_report_html = f"{self.outdir}/../{self.sample}_report.html"

        with open(self.__step_log_file, 'w') as writer:
            writer.write("Celatlas Step Log\n")
            writer.write("=================\n")
            writer.write(f"Status: {status}\n")
            writer.write(f"Started At: {self._started_at.isoformat(timespec='seconds')}\n")
            writer.write(f"Completed At: {completed_at.isoformat(timespec='seconds')}\n")
            writer.write(f"Elapsed: {elapsed}\n")
            writer.write(f"Sample: {self.sample}\n")
            writer.write(f"Assay: {self.assay}\n")
            writer.write(f"Step: {self._step_name}\n")
            writer.write(f"Display Title: {self._display_title}\n")
            writer.write(f"Outdir: {self.outdir}\n")
            writer.write(f"Thread: {self.thread}\n")
            writer.write(f"Debug: {self.debug}\n")
            if exc_type:
                writer.write(f"Exception: {exc_type.__name__}: {exc_value}\n")
            writer.write(f"Legacy HTML Report: disabled ({legacy_report_html} not written)\n")

            writer.write("\nArguments\n")
            writer.write("---------\n")
            for key, value in sorted(vars(self.args).items()):
                writer.write(f"{key}: {self._format_log_value(value)}\n")

            writer.write("\nCommands\n")
            writer.write("--------\n")
            if not self.__commands:
                writer.write("No subprocess commands recorded by this step.\n")
            for index, cmd in enumerate(self.__commands, start=1):
                writer.write(f"{index}. {self._format_log_value(cmd, max_length=2000)}\n")

            writer.write("\nMetrics\n")
            writer.write("-------\n")
            visible_metrics = [metric for metric in self.__metric_list if metric.get('show')]
            if not visible_metrics:
                writer.write("No visible metrics recorded.\n")
            for metric in visible_metrics:
                writer.write(f"- {metric['name']}: {metric['display']}\n")
                writer.write(f"  value: {self._format_log_value(metric['value'])}\n")
                if metric.get('total') is not None:
                    writer.write(f"  total: {self._format_log_value(metric['total'])}\n")
                if metric.get('fraction') is not None:
                    writer.write(f"  fraction: {metric['fraction']}%\n")
                if metric.get('help_info'):
                    writer.write(f"  help: {metric['help_info']}\n")

            hidden_metrics = [metric for metric in self.__metric_list if not metric.get('show')]
            if hidden_metrics:
                writer.write("\nHidden Metrics\n")
                writer.write("--------------\n")
                for metric in hidden_metrics:
                    writer.write(f"- {metric['name']}: {metric['display']}\n")

            writer.write("\nHelp Content\n")
            writer.write("------------\n")
            if not self.__help_content:
                writer.write("No help content recorded.\n")
            for item in self.__help_content:
                writer.write(f"- {item.get('name', 'N/A')}: {item.get('content', '')}\n")

            writer.write("\nData Summary\n")
            writer.write("------------\n")
            step_data = self.__content_dict['data'].get(self._step_summary_name, {})
            data_items = [
                (key, value)
                for key, value in step_data.items()
                if key not in {'metric_list', 'help_content'}
            ]
            if not data_items:
                writer.write("No extra data entries recorded.\n")
            for key, value in data_items:
                writer.write(f"{key}: {self._format_log_value(value)}\n")

            writer.write("\nSummary Files\n")
            writer.write("-------------\n")
            writer.write(f"step_log: {self.__step_log_file}\n")
            for slot, path in sorted(self._path_dict.items()):
                writer.write(f"{slot}: {path}\n")

            writer.write("\nOutput Directory Entries\n")
            writer.write("------------------------\n")
            try:
                with os.scandir(self.outdir) as scan:
                    current_log = os.path.abspath(self.__step_log_file)
                    entries = sorted(
                        (entry for entry in scan if os.path.abspath(entry.path) != current_log),
                        key=lambda entry: entry.name
                    )
            except OSError as error:
                writer.write(f"Unable to list output directory: {error}\n")
                return

            if not entries:
                writer.write("No output entries found.\n")
            for entry in entries[:200]:
                try:
                    stat = entry.stat()
                    kind = "dir" if entry.is_dir() else "file"
                    size = "-" if entry.is_dir() else self._format_file_size(stat.st_size)
                    modified = datetime.fromtimestamp(stat.st_mtime).isoformat(timespec='seconds')
                    writer.write(f"{kind}\t{size}\t{modified}\t{entry.name}\n")
                except OSError as error:
                    writer.write(f"unknown\t-\t-\t{entry.name} ({error})\n")
            if len(entries) > 200:
                writer.write(f"... {len(entries) - 200} more entries omitted\n")

    def _render_html(self):
        """Legacy per-step HTML reports are intentionally disabled."""
        return None

    def _add_content_data(self):
        step_summary = {'display_title': self._display_title}
        metric_list = []
        for metric in self.__metric_list:
            if metric['show']:
                metric_list.append(metric)
        step_summary['metric_list'] = metric_list
        step_summary['help_content'] = self.__help_content
        self.__content_dict['data'][self._step_summary_name].update(step_summary)

    def _add_content_metric(self):
        metric_dict = dict()
        for metric in self.__metric_list:
            name = metric['name']
            value = metric['value']
            fraction = metric['fraction']
            metric_dict[name] = value
            if fraction:
                metric_dict[f'{name} Fraction'] = fraction

        self.__content_dict['metrics'][self._step_summary_name].update(metric_dict)

    def add_data(self, **kwargs):
        """
        add data(other than metrics) to self.content_dict['data']
        for example: add plots and tables
        """
        for key, value in kwargs.items():
            self.__content_dict['data'][self._step_summary_name][key] = value

    def add_help_content(self, name, content):
        """
        add help info before metrics' help_info
        """
        self.__help_content.append(
            {
                'name': name,
                'content': content
            }
        )

    @utils.add_log
    def get_slot_key(self, slot, step_name, key):
        """
        read slot from json file
        """
        try:
            return self.__content_dict[slot][step_name + '_summary'][key]
        except KeyError:
            self.get_slot_key.logger.warning(f'{key} not found in {step_name}_summary.{slot}')
            raise

    def get_table_dict(self, title, table_id, df_table):
        """
        table_dict {title: '', table_id: '', df_table: pd.DataFrame}
        """
        table_dict = {'title': title, 'table': df_table.to_html(
            escape=False,
            index=False,
            table_id=table_id,
            justify="center"), 'id': table_id}
        return table_dict

    @utils.add_log
    def _clean_up(self):
        self._add_content_data()
        self._add_content_metric()
        self._dump_content()
        self._write_step_log()

    @utils.add_log
    def debug_subprocess_call(self, cmd):
        """
        debug subprocess call
        """
        self.__commands.append(cmd)
        self.debug_subprocess_call.logger.info(cmd)
        if cmd.find('2>&1') == -1:
            cmd += ' 2>&1 '
        subprocess.check_call(cmd, shell=True)

    def get_metric_list(self):
        return self.__metric_list

    def set_metric_list(self, metric_list):
        self.__metric_list = metric_list

    @abc.abstractmethod
    def run(self):
        sys.exit('Please implement run() method.')

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self._exit_info = args
        self._clean_up()
