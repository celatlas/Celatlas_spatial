"""Display-only channel replacement for completed Celatlas v1.8 results.

No pipeline/model imports. Legacy transforms are recovered only for the
coarse-mask registration branch, and must reproduce the existing DAPI exactly.
"""
from __future__ import annotations

import csv
from contextlib import contextmanager
import fcntl
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
from datetime import datetime, timezone

import cv2
import numpy as np
import tifffile
from PIL import Image, ImageDraw

VERSION = '0.1.0'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def signature(path):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtime_ns=str(stat.st_mtime_ns), sha256=digest(path))


def write_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name + '.writing')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    os.replace(temp, path)


def rgb_tiff(path):
    with tifffile.TiffFile(path) as tif:
        p = tif.pages[0]
        if len(tif.pages) != 1 or p.axes != 'YXS' or p.shape[2] != 3 or p.dtype != np.uint8:
            raise ValueError('首版换底图要求单页 RGB uint8 TIFF，不自动改变位深/颜色空间。')
        if int(p.photometric) != 2 or int(p.tags.get('Orientation').value if p.tags.get('Orientation') else 1) != 1:
            raise ValueError('要求 RGB、Orientation=1。')
    try:
        return tifffile.memmap(path, mode='r')
    except ValueError:
        return tifffile.imread(path)


def read_cv(path, gray=False):
    out = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE if gray else cv2.IMREAD_COLOR)
    if out is None:
        raise ValueError(f'无法完整解码图像（检查文件是否损坏）: {path}')
    return out if gray else cv2.cvtColor(out, cv2.COLOR_BGR2RGB)


def save_image(path, arr, quality=92):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(path), cv2.cvtColor(arr, cv2.COLOR_RGB2BGR) if arr.ndim == 3 else arr,
                       [cv2.IMWRITE_JPEG_QUALITY, quality] if path.suffix == '.jpg' else [cv2.IMWRITE_PNG_COMPRESSION, 6]):
        raise IOError(f'无法保存 {path}')


def recover_transform(log_text):
    """Recover full precision scale from integer mask areas, not rounded logs."""
    def one(pattern):
        matches = re.findall(pattern, log_text)
        if len(matches) != 1:
            raise ValueError(f'需要单次成功运行日志，字段缺失或重复: {pattern}')
        return matches[0]
    if 'rectify=False' not in log_text or '[Registration] Coarse registration completed (refinement disabled)' not in log_text:
        raise ValueError('仅支持 rectify=False、无特征/SimpleITK 精配准的旧日志；其他情况需精确变换文件。')
    if '[Stage 2]' in log_text:
        raise ValueError('存在精配准，禁止仅凭粗配准日志推断最终矩阵。')
    h, w = map(int, one(r'\[HE Image\] Original size: \((\d+), (\d+), 3\)'))
    th, tw = map(int, one(r'\[Canvas\] Resolved tissue canvas[^\n]*height=(\d+), width=(\d+)'))
    resized = re.findall(r'\[HE Image\] Resized to \((\d+), (\d+)\)', log_text)
    factor = max(min(th/h, tw/w), .3) if min(th/h, tw/w) < .8 else 1.0
    rh, rw = (map(int, resized[0]) if len(resized) == 1 else (h, w))
    if len(resized) > 1 or (factor != 1.0 and not resized) or (round(h*factor),round(w*factor)) != (rh,rw):
        raise ValueError('日志缩放尺寸与 v1.8 规则不一致。')
    fixed = one(r'\[Contour Features\] Fixed\s+- center: \((\d+), (\d+)\), angle: [^,]+, area: (\d+) px')
    moving = one(r'\[Contour Features\] Moving\s+- center: \((\d+), (\d+)\), angle: [^,]+, area: (\d+) px')
    fx,fy,fa = map(int,fixed); mx,my,ma = map(int,moving)
    if min(fa,ma) <= 0:
        raise ValueError('mask 面积无效。')
    angle = float(one(r'\[Fine Search\] Completed! Final best angle: ([\-\d.]+)°'))
    if abs(angle*10-round(angle*10)) > 1e-6:
        raise ValueError('角度不在 v1.8 的 0.1 度搜索网格上。')
    scale = math.sqrt(fa/ma)
    matrix = cv2.getRotationMatrix2D((mx,my), angle, scale)
    matrix[:,2] += np.array([fx,fy]) - matrix @ [mx,my,1]
    pre = np.array([[factor,0,(factor-1)/2],[0,factor,(factor-1)/2],[0,0,1]])
    combined = np.vstack([matrix,[0,0,1]]) @ pre
    return dict(schema='celatlas_display_transform_v1', source_shape=[h,w,3],
                output_shape=[th,tw,3], resized_shape=[rh,rw], resize_factor=factor,
                resize_interpolation='INTER_AREA', warp_interpolation='INTER_LINEAR',
                mask_resize_interpolation='INTER_AREA followed by >127',
                angle_deg_ccw=angle, mask_area_scale=scale, fixed_area=fa, moving_area=ma,
                fixed_center=[fx,fy], moving_center=[mx,my],
                resized_to_registered=matrix.tolist(), source_to_registered=combined.tolist(),
                registered_to_source=np.linalg.inv(combined).tolist(),
                rendering='resize then warp (two stages), mask==0 black; combined matrix is coordinate mapping only')


def prepare_mask(path, transform):
    mask = read_cv(path, gray=True)
    rh,rw = transform['resized_shape']
    mask = (cv2.resize(mask,(rw,rh),interpolation=cv2.INTER_AREA)>127).astype(np.uint8)*255
    moments = cv2.moments(mask)
    center = [int(moments['m10']/moments['m00']),int(moments['m01']/moments['m00'])] if moments['m00'] else None
    if center != transform['moving_center'] or np.count_nonzero(mask) != transform['moving_area']:
        raise ValueError('当前 HE mask 与当次配准日志不一致，禁止套用。')
    h,w,_ = transform['output_shape']
    return cv2.warpAffine(mask,np.asarray(transform['resized_to_registered']),(w,h))


def warp(image, transform, mask):
    if list(image.shape) != transform['source_shape']:
        raise ValueError('原图尺寸与配准源不一致，需先恢复相同裁切画布。')
    f = transform['resize_factor']
    small = cv2.resize(image,None,fx=f,fy=f,interpolation=cv2.INTER_AREA) if f != 1 else image
    h,w,_ = transform['output_shape']
    result = cv2.warpAffine(small,np.asarray(transform['resized_to_registered']),(w,h))
    result[mask==0] = 0
    return result


def compare_pixels(a,b):
    if a.shape != b.shape:
        raise ValueError('配准画布尺寸不匹配。')
    different = total = maxerr = 0
    for y in range(0,a.shape[0],256):
        delta = np.abs(a[y:y+256].astype(np.int16)-b[y:y+256].astype(np.int16))
        different += int(np.count_nonzero(delta)); total += int(delta.sum()); maxerr=max(maxerr,int(delta.max()))
    return dict(different_channel_values=different, total_channel_values=int(a.size),
                max_absolute_error=maxerr, mean_absolute_error=total/a.size, exact=different==0)


def sampled(arr, size):
    w,h=size; scale=min(1.,max(size)/max(arr.shape[:2]))
    ys=np.minimum((np.arange(h)/scale).astype(int),arr.shape[0]-1)
    xs=np.minimum((np.arange(w)/scale).astype(int),arr.shape[1]-1)
    return np.asarray(arr[np.ix_(ys,xs)])


def decode_gem(gem_rgb, old_bg, old_full_overlay):
    """Invert saved lossless VIRIDIS; disambiguate duplicate colours using QC PNG."""
    vir=cv2.cvtColor(cv2.applyColorMap(np.arange(256,dtype=np.uint8)[:,None],cv2.COLORMAP_VIRIDIS),cv2.COLOR_BGR2RGB).reshape(256,3)
    jet=cv2.cvtColor(cv2.applyColorMap(np.arange(256,dtype=np.uint8)[:,None],cv2.COLORMAP_JET),cv2.COLOR_BGR2RGB).reshape(256,3)
    def pack(a):
        return (a[...,0].astype(np.uint32)<<16)|(a[...,1].astype(np.uint32)<<8)|a[...,2].astype(np.uint32)
    codes=pack(vir); unique=np.unique(codes)
    mapping=np.array([np.nonzero(codes==v)[0][0] for v in unique],dtype=np.uint8)
    result=np.empty(gem_rgb.shape[:2],dtype=np.uint8)
    for y in range(0,len(result),128):
        code=pack(gem_rgb[y:y+128]); idx=np.searchsorted(unique,code)
        if np.any(idx==len(unique)) or np.any(unique[np.minimum(idx,len(unique)-1)]!=code):
            raise ValueError('GEM 图不是原始无损 VIRIDIS 图，无法安全重建显示层。')
        val=mapping[idx]; bg=old_bg[y:y+128]; ref=old_full_overlay[y:y+128]
        for c in unique:
            choices=np.nonzero(codes==c)[0]
            if len(choices)>1:
                m=code==c
                for alternative in choices[1:]:
                    candidate=cv2.addWeighted(bg,.6,np.broadcast_to(jet[alternative],bg.shape).copy(),.4,0)
                    val[m & np.all(candidate==ref,axis=2)] = alternative
        rebuilt=cv2.addWeighted(bg,.6,jet[val],.4,0)
        if not np.array_equal(rebuilt,ref):
            raise ValueError('现有 full-canvas GEM 叠图不能按原配方精确复原，拒绝改变表达显示层。')
        result[y:y+128]=val
    return result,jet[result]


def bin_previews(root, sample, staged, old_registered, new_registered, mask):
    images=root/'06.segment/01.binsegment/images'; out=staged/'06.segment/01.binsegment/images'
    out.mkdir(parents=True,exist_ok=True)
    for name,interp in [('tissue_hires_image.png',cv2.INTER_LANCZOS4),('tissue_lowres_image.png',cv2.INTER_AREA)]:
        with Image.open(images/name) as im: size=im.size
        preview=cv2.resize(new_registered,size,interpolation=interp)
        save_image(out/name,preview)
        for target in (root/'06.segment/01.binsegment/square_bin').glob('*/spatial/'+name):
            dst=staged/target.relative_to(root); dst.parent.mkdir(parents=True,exist_ok=True)
            with Image.open(target) as im:
                if im.size!=size: raise ValueError(f'空间预览尺寸不同: {target}')
            shutil.copy2(out/name,dst)
    overlays=images/'overlays'; dest=out/'overlays'
    gem=read_cv(overlays/'1_gem_expression.png'); size=(gem.shape[1],gem.shape[0])
    old_bg=cv2.resize(old_registered,size,interpolation=cv2.INTER_AREA)
    existing_mask=read_cv(overlays/'3a_he_mask.png',True)>0
    old_bg[~existing_mask]=0
    old_full=read_cv(overlays/'5_gem_heatmap_on_he.full_canvas_qc.png')
    norm,heatmap=decode_gem(gem,old_bg,old_full)
    bg=cv2.resize(new_registered,size,interpolation=cv2.INTER_AREA)
    m=cv2.resize(mask,size,interpolation=cv2.INTER_NEAREST)>0
    if not np.array_equal(m,existing_mask):
        raise ValueError('已保存 HE mask 预览与重建 mask 不一致。')
    bg[~m]=0
    save_image(dest/'2_he_registered.jpg',bg,85)
    union=read_cv(overlays/'3_tissue_segmentation_mask.png',True)>0
    blend=cv2.addWeighted(bg,.7,heatmap,.3,0); blend[~union]=0
    save_image(dest/'3b_tissue_segmentation_gem_heatmap.png',blend)
    green=np.zeros_like(bg); green[:,:,1]=180
    tinted=bg.copy(); tinted[union]=cv2.addWeighted(bg,.7,green,.3,0)[union]
    contours,_=cv2.findContours(union.astype(np.uint8)*255,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
    cv2.drawContours(tinted,contours,-1,(0,255,0),2)
    save_image(dest/'4_tissue_segmentation_on_he.jpg',tinted,85)
    tissue=tifffile.memmap(images/f'{sample}_tissue_cut.tif',mode='r')
    keep=(cv2.resize(tissue,size,interpolation=cv2.INTER_NEAREST)>0)&m&(norm>0)
    blend=bg.copy()
    blend[keep]=np.clip(np.rint(bg[keep].astype(np.float32)*.6+heatmap[keep].astype(np.float32)*.4),0,255).astype(np.uint8)
    save_image(dest/'5_gem_heatmap_on_he.png',blend)
    save_image(dest/'5_gem_heatmap_on_he.full_canvas_qc.png',cv2.addWeighted(bg,.6,heatmap,.4,0))


def expanded_sample(labels, distance, size, progress=print):
    """Same full-width slab expansion as v1.8; only keep preview sample pixels."""
    from skimage.segmentation import expand_labels
    w,h=size; scale=min(1.,max(size)/max(labels.shape)); radius=max(1,math.ceil(distance))
    ys=np.minimum((np.arange(h)/scale).astype(int),labels.shape[0]-1)
    xs=np.minimum((np.arange(w)/scale).astype(int),labels.shape[1]-1)
    out=np.zeros((h,w),np.uint32)
    if distance<=0: return np.asarray(labels[np.ix_(ys,xs)])
    for y0 in range(0,labels.shape[0],256):
        y1=min(labels.shape[0],y0+256); r0=max(0,y0-radius); r1=min(labels.shape[0],y1+radius)
        chosen=np.nonzero((ys>=y0)&(ys<y1))[0]
        if chosen.size:
            slab=expand_labels(np.asarray(labels[r0:r1]),distance=distance)
            out[chosen]=slab[np.ix_(ys[chosen]-r0,xs)]
        if y0%2048==0: progress(f'细胞边界预览 {y0}/{labels.shape[0]} 行')
    return out


def cell_previews(root,sample,staged,new_registered,progress=print):
    cell=root/'06.segment/02.cellsegment'; vis=cell/'visualization'
    out=staged/'06.segment/02.cellsegment/visualization';out.mkdir(parents=True,exist_ok=True)
    info=json.loads((cell/f'{sample}.stardist_assignment_summary.json').read_text())
    labels=tifffile.memmap(cell/'stardist_inference'/f'{sample}.stardist_labels.tif',mode='r')
    if list(labels.shape)!=info['label_shape'] or labels.shape!=new_registered.shape[:2]:
        raise ValueError('核 label、assignment 与配准画布不同，不可更换。')
    for suffix in ('barcode_assignment.preview.jpg','barcode_assignment.full_chip_qc.jpg'):
        with Image.open(vis/f'{sample}.{suffix}') as im: size=im.size
        canvas=Image.fromarray(sampled(new_registered,size));draw=ImageDraw.Draw(canvas,'RGBA')
        x0,y0,bw,bh=info['bbox']
        if bw<=0 or bh<=0: raise ValueError('barcode bbox 无效。')
        assigned=[]
        with (cell/f'{sample}.stardist_barcode_assignments.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                x=(float(row['raw_x'])-x0)*size[0]/bw; y=(float(row['raw_y'])-y0)*size[1]/bh
                if 0<=x<size[0] and 0<=y<size[1]:
                    if int(row['cell_id'])>0: assigned.append((int(x),int(y)))
                    else: draw.point((int(x),int(y)),fill=(255,145,0,105))
        for x,y in assigned: draw.ellipse((x-1,y-1,x+1,y+1),fill=(0,220,255,155))
        canvas.save(out/f'{sample}.{suffix}',quality=92)
    name=f'{sample}.qc_boundary_overlay.preview.jpg'
    with Image.open(vis/name) as im: size=im.size
    lab=expanded_sample(labels,float(info['expand_pixels']),size,progress)
    lookup=np.zeros(int(labels.max())+1,np.uint32)
    with (cell/f'{sample}.cell_id_map.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'): lookup[int(row['raw_cell_id'])]=int(row['new_cell_id'])
    from skimage.segmentation import find_boundaries
    canvas=sampled(new_registered,size).copy()
    canvas[find_boundaries(lookup[lab],mode='outer')]=[0,255,255]
    save_image(out/name,canvas)
    full=f'{sample}.qc_boundary_overlay.jpg'
    with Image.open(vis/full) as im:
        if im.size!=size: raise ValueError('当前只支持 slide-sized preview-only QC，无全分辨率局部裁图。')
    shutil.copy2(out/name,out/full)


def allowed(relative, sample):
    """Write whitelist: never touch inference inputs, masks, labels or counts."""
    p=Path(relative)
    if p.is_absolute() or '..' in p.parts: return False
    names={'tissue_hires_image.png','tissue_lowres_image.png'}
    base=Path('06.segment/01.binsegment')
    if p.parent==base/'images': return p.name in names|{f'{sample}_display_regist.tif',f'{sample}_display_regist.tif.transform.json'}
    if p.parent==base/'images/overlays':
        return p.name in {'2_he_registered.jpg','3b_tissue_segmentation_gem_heatmap.png','4_tissue_segmentation_on_he.jpg',
                          '5_gem_heatmap_on_he.png','5_gem_heatmap_on_he.full_canvas_qc.png'}
    if len(p.parts)==6 and p.parts[:3]==('06.segment','01.binsegment','square_bin') and p.parts[4]=='spatial':
        return p.name in names
    if p.parent==Path('06.segment/02.cellsegment/visualization'):
        return p.name in {f'{sample}.{s}' for s in ('barcode_assignment.preview.jpg','barcode_assignment.full_chip_qc.jpg',
                                                  'qc_boundary_overlay.jpg','qc_boundary_overlay.preview.jpg')}
    return False


def safe_target(root,relative,sample):
    if not allowed(relative,sample): raise ValueError(f'目标不在显示文件白名单: {relative}')
    target=root/relative
    if target.is_symlink() or not target.resolve().is_relative_to(root):
        raise ValueError(f'禁止通过软链接修改项目外数据: {target}')
    return target


def prepare(root,sample,channel,dapi,mask,log_path,output_root,confirm_shared_canvas=False,progress=print):
    if not re.fullmatch(r'[A-Za-z0-9_-]+',sample): raise ValueError('sample ID 格式无效。')
    if not confirm_shared_canvas: raise ValueError('请核实新通道与 DAPI 共享裁切画布后添加 --confirm-shared-canvas。')
    root=Path(root).resolve(strict=True); channel=Path(channel).resolve(strict=True)
    dapi=Path(dapi).resolve(strict=True); mask=Path(mask).resolve(strict=True); log_path=Path(log_path).resolve(strict=True)
    registered=root/f'06.segment/01.binsegment/images/{sample}_regist.tif'
    if registered.with_name(f'{sample}_display_regist.tif').exists():
        raise ValueError('该样本已存在换底图结果；请先 rollback 上次替换，再生成新的预览。')
    sidecar=json.loads(Path(str(registered)+'.registration.json').read_text())
    if not sidecar.get('fluorescence_background') or sidecar.get('method')!='HE':
        raise ValueError('首版仅处理 HE 模式配准的黑底荧光结果。')
    if sidecar.get('image_bytes')!=registered.stat().st_size or str(sidecar.get('image_mtime_ns'))!=str(registered.stat().st_mtime_ns):
        raise ValueError('配准图与原流程签名不一致。')
    log_text=log_path.read_text(errors='replace')
    if f'--sample {sample} ' not in log_text: raise ValueError('日志 sample 不匹配。')
    transform=recover_transform(log_text)
    progress('记录输入指纹；恢复原始配准变换…')
    inputs=[signature(p) for p in (dapi,channel,mask,log_path,registered,Path(str(registered)+'.registration.json'))]
    old=rgb_tiff(registered); mask_warp=prepare_mask(mask,transform)
    reconstruction=warp(rgb_tiff(dapi),transform,mask_warp)
    validation=compare_pixels(reconstruction,old); del reconstruction
    progress(f'DAPI 完整像素校验: {validation}')
    if not validation['exact']: raise ValueError('DAPI 重建不是逐像素相同，禁止替换。可能用了旧日志/新裁图/不同插值版本。')
    new=warp(rgb_tiff(channel),transform,mask_warp)
    output_root=Path(output_root).resolve();output_root.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(output_root).free < new.nbytes*3+512*1024**2: raise ValueError('预览工作目录磁盘空间不足。')
    run=Path(tempfile.mkdtemp(prefix='swap_'+datetime.now().strftime('%Y%m%d_%H%M%S_'),dir=output_root))
    staged=run/'staged'; staged.mkdir()
    protected=[]
    for folder in ('06.segment/01.binsegment','06.segment/02.cellsegment'):
        for path in sorted((root/folder).rglob('*')):
            if path.is_file() and not allowed(path.relative_to(root),sample):
                # Stat guard for every non-display result; full hashes for principal scientific products.
                stat=path.stat(); entry=dict(path=str(path),bytes=stat.st_size,mtime_ns=str(stat.st_mtime_ns))
                if path.suffix in ('.gz','.tsv','.csv') or path.name==f'{sample}.stardist_labels.tif' or path==registered:
                    entry['sha256']=digest(path)
                protected.append(entry)
    try:
        transform.update(software_version=VERSION,source_dapi=inputs[0],source_display=inputs[1],manual_mask=inputs[2],
                         registration_log=inputs[3],validation=validation,shared_canvas_user_confirmed=True)
        rel=Path(f'06.segment/01.binsegment/images/{sample}_display_regist.tif')
        dest=staged/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        tifffile.imwrite(dest,new,photometric='rgb',metadata={'axes':'YXS','purpose':'display_only'})
        write_json(str(dest)+'.transform.json',transform)
        progress('更新 bin 显示层（复用原 GEM 图、mask 和尺寸）…')
        bin_previews(root,sample,staged,old,new,mask_warp)
        progress('更新 cell 显示层（复用原 assignment 和 QC ID，不重新分割/计数）…')
        cell_previews(root,sample,staged,new,progress)
        preview_scale=min(1.0,1230/max(old.shape[:2]))
        preview_size=tuple(max(1,round(v*preview_scale)) for v in old.shape[1::-1])
        save_image(run/'dapi_registered.jpg',cv2.resize(old,preview_size,interpolation=cv2.INTER_AREA))
        save_image(run/'display_registered.jpg',cv2.resize(new,preview_size,interpolation=cv2.INTER_AREA))
        files=[]
        for path in sorted(staged.rglob('*')):
            if not path.is_file(): continue
            relative=str(path.relative_to(staged)); target=safe_target(root,relative,sample)
            files.append(dict(relative=relative,new_sha256=digest(path),old_sha256=digest(target) if target.exists() else None))
        manifest=dict(schema='celatlas_background_swap_v1',version=VERSION,state='prepared',sample=sample,root=str(root),
                      created_utc=datetime.now(timezone.utc).isoformat(),inputs=inputs,protected=protected,files=files,
                      transform=transform,validation=validation,
                      unchanged=['registered DAPI/inference caches','labels/QC IDs/UMI/matrices/barcodes/coordinates/masks','cell_umi_heatmap (no channel background)',
                                 '07.outs, h5ad embedded images and report HTML (not updated by this scoped tool)'])
        write_json(run/'manifest.json',manifest)
        rows=''.join(f'<li>{html.escape(f["relative"])}</li>' for f in files)
        (run/'review.html').write_text('<!doctype html><meta charset="utf-8"><title>换底图检查</title>'
            '<h1>显示通道替换预览</h1><p>DAPI 重建校验：全部像素相同。原 DAPI、核标签和表达矩阵保持不变。'
            '同片通道的原始偏移仍需人工核对。这是生成时的预览；当前应用状态请查看 manifest.json 的 state。</p>'
            '<img width="45%" src="dapi_registered.jpg"><img width="45%" src="display_registered.jpg">'
            '<p>QC 边界用青色显示，避免与红色荧光混淆。Cell UMI Heatmap 是纯表达热图，保持不变。'
            '仅修改 06.segment；07.outs/h5ad/HTML 中已嵌入的图片不会自动更新。</p><ul>'+rows+'</ul>',encoding='utf-8')
        check_inputs(manifest)
        progress(f'已准备 {len(files)} 个显示文件；尚未写入原结果。检查 {run}/review.html')
        return run
    except Exception as exc:
        (run/'FAILED.txt').write_text(str(exc),encoding='utf-8');raise


def check_inputs(manifest):
    for entry in manifest['inputs']+manifest['protected']:
        path=Path(entry['path']); stat=path.stat()
        if stat.st_size!=entry['bytes'] or str(stat.st_mtime_ns)!=entry['mtime_ns']:
            raise ValueError(f'输入/分析结果在预览后改变: {path}')
        if 'sha256' in entry and digest(path)!=entry['sha256']:
            raise ValueError(f'输入/分析结果内容改变: {path}')


def atomic_copy(source,target):
    target.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='.background_swap_',dir=target.parent);os.close(fd)
    try:
        shutil.copy2(source,name);os.replace(name,target)
    finally:
        if Path(name).exists(): Path(name).unlink()


@contextmanager
def project_lock(root):
    """Serialize tool writers. The analysis pipeline itself must remain stopped."""
    lock=Path(root)/'06.segment/.background_swap.lock'
    with lock.open('a') as handle:
        try:
            fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('该样本已有换底图/回滚任务运行，请等待完成。') from exc
        try:
            yield
        finally:
            fcntl.flock(handle,fcntl.LOCK_UN)


def apply(run,confirm=False):
    if not confirm: raise ValueError('应用将更换列出的显示文件；请检查 review.html 后使用 --yes。')
    run=Path(run).resolve(strict=True)
    manifest=json.loads((run/'manifest.json').read_text())
    with project_lock(Path(manifest['root']).resolve(strict=True)):
        return _apply_locked(run,confirm)


def _apply_locked(run,confirm=False):
    if not confirm: raise ValueError('应用将更换列出的显示文件；请检查 review.html 后使用 --yes。')
    run=Path(run).resolve(strict=True);manifest=json.loads((run/'manifest.json').read_text())
    if manifest['schema']!='celatlas_background_swap_v1' or manifest['state']!='prepared': raise ValueError('任务不在 prepared 状态。')
    root=Path(manifest['root']).resolve(strict=True);check_inputs(manifest)
    targets=[]
    for f in manifest['files']:
        target=safe_target(root,f['relative'],manifest['sample']); staged=run/'staged'/f['relative']
        if not staged.resolve().is_relative_to(run/'staged') or digest(staged)!=f['new_sha256']: raise ValueError('预览文件改变。')
        if (digest(target) if target.exists() else None)!=f['old_sha256']: raise ValueError(f'目标在预览后改变: {target}')
        targets.append((f,target,staged))
    # Backup resides in the writable project with its own complete manifest.
    backup=root/'06.segment/background_swap_backups'/run.name
    backup.mkdir(parents=True,exist_ok=False)
    write_json(backup/'manifest.json',manifest)
    for f,target,staged in targets:
        if target.exists():
            saved=backup/'originals'/f['relative'];saved.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(target,saved)
            if digest(saved)!=f['old_sha256']: raise IOError('备份校验失败，未执行替换。')
    journal=dict(state='applying',attempted=[],applied=[]);write_json(backup/'journal.json',journal)
    try:
        for f,target,staged in targets:
            journal['attempted'].append(f['relative']);write_json(backup/'journal.json',journal)
            atomic_copy(staged,target)
            journal['applied'].append(f['relative']);write_json(backup/'journal.json',journal)
        check_inputs(manifest)
        journal['state']='applied';write_json(backup/'journal.json',journal)
        manifest['state']='applied';manifest['backup']=str(backup);write_json(run/'manifest.json',manifest)
    except Exception:
        # Only restore files actually attempted; original scientific products never enter this list.
        for f,target,_ in reversed(targets):
            if f['relative'] not in journal['attempted']: continue
            if f['old_sha256'] is None: target.unlink(missing_ok=True)
            else: atomic_copy(backup/'originals'/f['relative'],target)
        journal['state']='rolled_back_after_error';write_json(backup/'journal.json',journal)
        raise
    return backup


def rollback(backup,confirm=False):
    if not confirm: raise ValueError('回滚需 --yes；将恢复备份并移除本次新建显示文件。')
    backup=Path(backup).resolve(strict=True)
    manifest=json.loads((backup/'manifest.json').read_text())
    with project_lock(Path(manifest['root']).resolve(strict=True)):
        return _rollback_locked(backup,confirm)


def _rollback_locked(backup,confirm=False):
    if not confirm: raise ValueError('回滚需 --yes；将恢复备份并移除本次新建显示文件。')
    backup=Path(backup).resolve(strict=True);m=json.loads((backup/'manifest.json').read_text())
    journal=json.loads((backup/'journal.json').read_text())
    if journal['state'] not in ('applied','applying','rolling_back'): raise ValueError('该备份不是可回滚状态。')
    root=Path(m['root']).resolve(strict=True);items=[]
    for f in m['files']:
        if f['relative'] not in journal['attempted']: continue
        target=safe_target(root,f['relative'],m['sample'])
        current=digest(target) if target.exists() else None
        if current not in (f['new_sha256'],f['old_sha256']): raise ValueError(f'应用后文件被修改，拒绝覆盖: {target}')
        saved=backup/'originals'/f['relative']
        if f['old_sha256'] is not None and digest(saved)!=f['old_sha256']: raise ValueError('备份损坏。')
        items.append((f,target,saved))
    journal['state']='rolling_back';write_json(backup/'journal.json',journal)
    for f,target,saved in items:
        if f['old_sha256'] is None: target.unlink(missing_ok=True)
        else: atomic_copy(saved,target)
    journal['state']='rolled_back';write_json(backup/'journal.json',journal)
    return backup
