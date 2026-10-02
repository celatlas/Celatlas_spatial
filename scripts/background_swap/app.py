#!/usr/bin/env python3
"""Standalone CLI: stage display images, review, apply with backup, rollback."""
import argparse
from pathlib import Path
import cv2

try:
    from .core import prepare, apply, rollback
except ImportError:
    from core import prepare, apply, rollback


def main():
    parser=argparse.ArgumentParser(description='Celatlas v1.8 配准后显示通道替换（不重新分割/计数）')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare',help='验证原配准，生成新底图及所有替换预览，不改原结果')
    p.add_argument('--sample-dir',required=True)
    p.add_argument('--sample',required=True)
    p.add_argument('--channel',required=True,help='与原始 DAPI 共享裁切画布的新显示通道')
    p.add_argument('--dapi',help='默认 06.segment/mask/<sample>_he.tif')
    p.add_argument('--he-mask',help='默认 06.segment/mask/<sample>_he_mask_manual.png')
    p.add_argument('--registration-log',required=True,help='生成现有 regist.tif 的那次 runner stdout.log')
    p.add_argument('--output-root',required=True,help='独立预览工作目录，需预留至少 2 GB')
    p.add_argument('--confirm-shared-canvas',action='store_true')
    p=sub.add_parser('apply',help='应用已检查的预览，自动备份；失败自动恢复')
    p.add_argument('--run',required=True)
    p.add_argument('--yes',action='store_true')
    p=sub.add_parser('rollback',help='恢复备份中的原显示文件，科学数据从不进入替换列表')
    p.add_argument('--backup',required=True)
    p.add_argument('--yes',action='store_true')
    args=parser.parse_args();cv2.setNumThreads(2)
    if args.command=='prepare':
        root=Path(args.sample_dir);mask_dir=root/'06.segment/mask'
        result=prepare(root,args.sample,args.channel,args.dapi or mask_dir/f'{args.sample}_he.tif',
                       args.he_mask or mask_dir/f'{args.sample}_he_mask_manual.png',args.registration_log,
                       args.output_root,args.confirm_shared_canvas,progress=lambda msg:print(msg,flush=True))
    elif args.command=='apply': result=apply(args.run,args.yes)
    else: result=rollback(args.backup,args.yes)
    print(f'{args.command} 完成: {result}',flush=True)


if __name__=='__main__': main()
