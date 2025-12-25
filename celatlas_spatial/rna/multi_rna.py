from celatlas_spatial.rna.__init__ import __ASSAY__
from celatlas_spatial.tools.multi import Multi


class Multi_rna(Multi):
    """
    ## Usage
    ```
        multi_rna\\
        --mapfile ./rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --thread 8\\
        --mod shell
    ```
    Work for both single cell RNA-Seq and single nuclei RNA-Seq.
    """

    def starsolo(self, sample):
        step = 'starsolo'
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)


def main():
    multi = Multi_rna(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    def common_args(parser):
        parser.add_argument('--outdir', type=str, default='/mnt/strna/work_project/pipeline/celatlas_spatial')
        parser.add_argument('--genomeDir', type=str, default='/mnt/strna/work_project/rawdata/celatlas_spatial/reference/Mus_musculus')
        parser.add_argument('--mapfile', type=str, default='/home/zhoumy/rna.mapfile')
        parser.add_argument('--mod', type=str, default='shell')
        parser.add_argument('--queue', type=bool, default=False)
        parser.add_argument('--thread', type=int, default=128)
        return parser


    multi = Multi_rna(__ASSAY__)
    multi.parser = common_args(multi.parser)
    multi.run()
