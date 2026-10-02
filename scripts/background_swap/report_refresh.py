"""Transactional display-only report/07.outs/h5ad refresh. No analysis rerun."""
from __future__ import annotations
import argparse
import base64
import hashlib
import json
import re
import shutil
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
import tifffile
from PIL import Image
from .core import atomic_copy, check_inputs, digest, project_lock, signature, write_json

IMAGE_URI = re.compile(r'data:image/[\w.+-]+;base64,([A-Za-z0-9+/=]+)')


def replace_report_images(html, mapping):
    hits = Counter()
    def replace(m):
        key = hashlib.sha256(base64.b64decode(m[1], validate=True)).hexdigest()
        if key not in mapping: return m[0]
        hits[key] += 1
        return mapping[key]
    result = IMAGE_URI.sub(replace, html)
    if IMAGE_URI.sub('IMAGE', result) != IMAGE_URI.sub('IMAGE', html):
        raise ValueError('Non-image HTML changed')
    return result, dict(hits)


def image_uri(path):
    mime = 'image/jpeg' if Path(path).suffix.lower() in ('.jpg','.jpeg') else 'image/png'
    return 'data:'+mime+';base64,'+base64.b64encode(Path(path).read_bytes()).decode('ascii')


def image_datasets(handle):
    names = []
    if 'uns/spatial' in handle:
        for library in handle['uns/spatial'].values():
            if 'images' in library:
                for key, ds in library['images'].items():
                    if key not in ('hires','lowres') or not isinstance(ds,h5py.Dataset):
                        raise ValueError(f'Unsupported spatial image: {ds.name}')
                    names.append(ds.name.lstrip('/'))
    return sorted(names)


def hash_value(h, value):
    a = np.asarray(value)
    if a.dtype.names:
        for name in a.dtype.names: h.update(name.encode()); hash_value(h,a[name])
    elif a.dtype.kind in 'OUS':
        for v in a.reshape(-1):
            raw = v if isinstance(v,bytes) else str(v).encode()
            h.update(len(raw).to_bytes(8,'little')); h.update(raw)
    else: h.update(a.tobytes(order='C'))


def h5_fingerprint(path, excluded=()):
    """All datasets and attributes; omit only permitted image pixel data."""
    records = {}
    with h5py.File(path,'r') as f:
        def visit(name,obj):
            h = hashlib.sha256()
            for k in sorted(obj.attrs): h.update(k.encode()); hash_value(h,obj.attrs[k])
            record = dict(attributes=h.hexdigest(),kind='group')
            if isinstance(obj,h5py.Dataset):
                record.update(kind='dataset',shape=list(obj.shape),dtype=str(obj.dtype))
                if name not in excluded:
                    h = hashlib.sha256()
                    if obj.ndim == 0: hash_value(h,obj[()])
                    else:
                        row_bytes = max(1,int(np.prod(obj.shape[1:]))*obj.dtype.itemsize)
                        step = max(1,8*1024**2//row_bytes)
                        for i in range(0,obj.shape[0],step): hash_value(h,obj[i:i+step])
                    record['data'] = h.hexdigest()
            records[name] = record
        visit('',f); f.visititems(visit)
    return records


def converted_image(ds,path):
    image = np.asarray(Image.open(path).convert('RGB'))
    if image.shape != ds.shape: raise ValueError(f'Image shape changed: {ds.name}')
    if ds.dtype == np.uint8: return image
    if ds.dtype.kind == 'f':
        old = ds[()]
        if not np.isfinite(old).all() or old.min()<0 or old.max()>1:
            raise ValueError('Float spatial images must be in [0,1]')
        return image.astype(ds.dtype)/ds.dtype.type(255)
    raise ValueError(f'Unsupported image dtype: {ds.dtype}')


def stage_h5(source,target,image_paths):
    with h5py.File(source,'r') as f: names = image_datasets(f)
    if not names: return None
    before = h5_fingerprint(source,names)
    target.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(source,target)
    expected = {}
    with h5py.File(target,'r+') as f:
        for name in names:
            image = converted_image(f[name],image_paths[name.rsplit('/',1)[-1]])
            f[name][...] = image
            expected[name] = hashlib.sha256(image.tobytes()).hexdigest()
    if before != h5_fingerprint(target,names): raise ValueError(f'Non-image HDF5 changed: {source}')
    with h5py.File(target,'r') as f:
        for name,value in expected.items():
            if hashlib.sha256(f[name][()].tobytes()).hexdigest()!=value: raise ValueError('Image verification failed')
    return dict(image_datasets=names,image_hashes=expected,non_image_data_identical=True,
                non_image_fingerprint=hashlib.sha256(json.dumps(before,sort_keys=True).encode()).hexdigest())


def plot_data(path):
    from anndata import AnnData
    from anndata._io.specs import read_elem
    from scipy.sparse import csr_matrix
    with h5py.File(path,'r') as f:
        obs = read_elem(f['obs']); a = AnnData(X=csr_matrix((len(obs),0)),obs=obs)
        a.obsm['spatial'] = f['obsm/spatial'][()]
        for key in ('spatial','cluster_colors','rank_genes_groups'):
            if key in f['uns']: a.uns[key] = read_elem(f['uns'][key])
    if not np.isfinite(a.obsm['spatial']).all(): raise ValueError('Non-finite spatial coordinates')
    return a


def draw_plots(root,sample,bin_size,staged,display):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scanpy as sc
    from anndata import AnnData
    from anndata._io.specs import read_elem
    from anndata.experimental import sparse_dataset
    from celatlas_spatial.tools.analysis_wrapper import Scanpy_wrapper
    b = Path('07.outs/binned_outputs'); c = Path('07.outs/cellsegmented_outputs')
    source = root/b/'h5ad'/f'{sample}_bin{bin_size}.h5ad'; a = plot_data(source)
    if sample not in a.uns.get('spatial',{}) or len(a.uns['spatial'])!=1:
        raise ValueError('Require one matching spatial library')
    meta = a.uns['spatial'][sample]; hires = root/'06.segment/01.binsegment/images/tissue_hires_image.png'
    image = plt.imread(hires)
    if image.shape!=meta['images']['hires'].shape: raise ValueError('Plot image shape differs')
    meta['images']['hires'] = image
    if not np.array_equal(a.obs[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy(),a.obsm['spatial']):
        raise ValueError('Coordinates are not full-resolution pixels')
    # Bypass constructor: use drawing helpers only, never run pipeline methods.
    helper = object.__new__(Scanpy_wrapper)
    helper.args = SimpleNamespace(sample=sample,hires_image=str(hires),fullres_image=str(display))
    helper.bin = f'bin{bin_size}'; helper.adata = {helper.bin:a}; outputs = []
    def save(rel):
        p = staged/rel; p.parent.mkdir(parents=True,exist_ok=True)
        plt.savefig(p,bbox_inches='tight',dpi=300); plt.close('all'); outputs.append(rel)
    for fields,name in [(['total_counts'],'spatial_total_counts'),(['n_genes_by_counts'],'spatial_gene_counts'),
                        (['total_counts','n_genes_by_counts'],'spatial_gene_expression_distribution')]:
        plt.figure(figsize=helper._spatial_figsize(a,base=3.6,ncols=len(fields)),dpi=300)
        sc.pl.spatial(a,color=fields,cmap='Spectral_r',spot_size=helper._spatial_spot_size(a),
                      alpha_img=.75,use_raw=False,show=False)
        save(b/'qc'/f'{sample}_{name}.png')
    plt.figure(figsize=helper._spatial_figsize(a,base=3.6),dpi=300)
    sc.pl.spatial(a,color=['cluster'],palette=sc.pl.palettes.default_20,size=1,
                  spot_size=helper._spatial_spot_size(a),alpha_img=.75,show=False)
    save(b/'clustering'/f'{sample}_spatial_cluster.png')
    marker_rel = b/'markers'/f'{sample}_spatial_markers.png'
    if (root/marker_rel).exists():
        de = Scanpy_wrapper.get_de_results(a,10,.25)
        genes = [de[de['cluster']==cl]['names'].iloc[0] for cl in sorted(de['cluster'].unique())[:4]]
        if not genes: raise ValueError('No recoverable marker genes')
        # Original marker visualization used raw counts before normalization.
        with h5py.File(source,'r') as f:
            var = read_elem(f['var']); indices = var.index.get_indexer(genes)
            if np.any(indices<0): raise ValueError('Marker missing from raw counts')
            raw = sparse_dataset(f['layers/raw'])
            values = np.column_stack([raw[:,int(i)].toarray().ravel() for i in indices])
        marker = AnnData(values,obs=a.obs.copy()); marker.var_names = genes
        marker.obsm['spatial'] = a.obsm['spatial'].copy(); marker.uns['spatial'] = a.uns['spatial']
        plt.figure(figsize=helper._spatial_figsize(marker,base=3.4,ncols=min(len(genes),4)),dpi=300)
        plt.subplots_adjust(wspace=.05,hspace=.05)
        sc.pl.spatial(marker,color=list(dict.fromkeys(genes)),cmap='Spectral_r',ncols=4,
                      spot_size=helper._spatial_spot_size(marker),alpha_img=.75,show=False)
        save(marker_rel)
    helper.bin = 'cell'; helper.adata = {'cell':plot_data(root/c/'h5ad'/f'{sample}_cell.h5ad')}
    helper.cell_spatial_dir = str(staged/c/'spatial'); Path(helper.cell_spatial_dir).mkdir(parents=True,exist_ok=True)
    for field,suffix,cat in [('cluster','cluster',True),('total_counts','total_counts',False),
                             ('n_genes_by_counts','n_genes',False),('pct_counts_mt','pct_counts_mt',False)]:
        rel = c/'spatial'/f'{sample}_cell_spatial_{suffix}.png'
        if (root/rel).exists():
            helper._write_cell_spatial_plot(field,rel.name,categorical=cat)
            if not (staged/rel).exists(): raise ValueError(f'Missing plot: {rel}')
            outputs.append(rel)
    rel = b/'clustering/tissue_hires_image.png'; shutil.copy2(hires,staged/rel); outputs.append(rel)
    return outputs


def checked_target(root,relative,sample):
    rel = Path(relative)
    if rel.is_absolute() or '..' in rel.parts: raise ValueError('Unsafe target')
    b = '07.outs/binned_outputs'; c = '07.outs/cellsegmented_outputs'
    allowed = {f'{sample}_spatial_analysis_report.html',f'{b}/clustering/tissue_hires_image.png',
               f'{b}/clustering/{sample}_spatial_cluster.png',f'{b}/markers/{sample}_spatial_markers.png'}
    allowed |= {f'{b}/qc/{sample}_spatial_{s}.png' for s in ('total_counts','gene_counts','gene_expression_distribution')}
    allowed |= {f'{c}/spatial/{sample}_cell_spatial_{s}.png' for s in ('cluster','total_counts','n_genes','pct_counts_mt')}
    allowed |= {f'{b}/h5ad/{sample}_bin{n}.h5ad' for n in (10,20,50,100)}
    allowed.add(f'{c}/h5ad/{sample}_cell.h5ad')
    target = root/rel
    if str(rel) not in allowed or target.is_symlink() or not target.resolve().is_relative_to(root):
        raise ValueError(f'Target outside whitelist: {relative}')
    return target


def prepare(swap_run,output_root,bin_size=50):
    swap_run = Path(swap_run).resolve(strict=True); swap = json.loads((swap_run/'manifest.json').read_text())
    if swap.get('state')!='applied': raise ValueError('Apply background swap first')
    root = Path(swap['root']).resolve(strict=True); sample = swap['sample']; check_inputs(swap)
    for f in swap['files']:
        if digest(root/f['relative'])!=f['new_sha256']: raise ValueError('Applied swap changed')
    output_root = Path(output_root).resolve(); output_root.mkdir(parents=True,exist_ok=True)
    h5paths = sorted(root.glob('07.outs/*/h5ad/*.h5ad'))
    if shutil.disk_usage(output_root).free < sum(p.stat().st_size for p in h5paths)*2+1024**3:
        raise ValueError('Insufficient staging space')
    run = Path(tempfile.mkdtemp(prefix='refresh_'+datetime.now().strftime('%Y%m%d_%H%M%S_'),dir=output_root))
    staged = run/'staged'; staged.mkdir(); report_rel = Path(f'{sample}_spatial_analysis_report.html')
    print('记录结果及 h5ad 指纹…',flush=True)
    inputs = [signature(root/report_rel),signature(swap_run/'manifest.json')]
    inputs += [signature(p) for p in sorted((root/'07.outs').rglob('*')) if p.is_file()]
    inputs += [signature(root/f['relative']) for f in swap['files']]
    display = root/f'06.segment/01.binsegment/images/{sample}_display_regist.tif'
    with tifffile.TiffFile(display) as t: shape = list(t.pages[0].shape)
    if shape!=swap['transform']['output_shape']: raise ValueError('Display canvas changed')
    print(f'画布 {shape[1]}×{shape[0]}；只重绘现有空间结果…',flush=True)
    outputs = draw_plots(root,sample,bin_size,staged,display)
    image_paths = {k:root/f'06.segment/01.binsegment/images/tissue_{k}_image.png' for k in ('hires','lowres')}
    h5_checks = {}; without_images = []
    for path in h5paths:
        rel = path.relative_to(root); checked_target(root,rel,sample)
        print(f'校验并更新 h5ad 显示层: {path.name}',flush=True)
        result = stage_h5(path,staged/rel,image_paths)
        if result is None: without_images.append(str(rel))
        else: outputs.append(rel); h5_checks[str(rel)] = result
    mapping = {}; required = set()
    def add(old_hash,new_file):
        uri = image_uri(new_file)
        if old_hash in mapping and mapping[old_hash]!=uri: raise ValueError('Ambiguous image mapping')
        mapping[old_hash] = uri
    for rel in outputs:
        if rel.suffix=='.png': add(digest(root/rel),staged/rel)
    for f in swap['files']:
        if f['old_sha256'] and Path(f['relative']).suffix in ('.png','.jpg'):
            add(f['old_sha256'],root/f['relative'])
    for rel in [Path('07.outs/binned_outputs/clustering/tissue_hires_image.png'),
                Path(f'07.outs/binned_outputs/clustering/{sample}_spatial_cluster.png'),
                Path(f'07.outs/binned_outputs/qc/{sample}_spatial_gene_expression_distribution.png')]:
        required.add(digest(root/rel))
    for f in swap['files']:
        if f['relative'].endswith(('qc_boundary_overlay.preview.jpg','barcode_assignment.preview.jpg','overlays/2_he_registered.jpg')):
            required.add(f['old_sha256'])
    updated,hits = replace_report_images((root/report_rel).read_text(),mapping)
    if not required.issubset(hits): raise ValueError('Report images do not match current results/swap backup')
    (staged/report_rel).write_text(updated); outputs.append(report_rel)
    files = []
    for rel in outputs:
        target = checked_target(root,rel,sample)
        files.append(dict(relative=str(rel),old_sha256=digest(target),new_sha256=digest(staged/rel)))
    m = dict(schema='celatlas_report_refresh_v2',state='prepared',root=str(root),sample=sample,
             swap_run=str(swap_run),inputs=inputs,protected=[],files=files,h5_checks=h5_checks,
             h5_without_images=without_images,canvas_shape=shape,report_image_replacements=sum(hits.values()),
             html_non_image_content_identical=True,bin_size=bin_size)
    check_inputs(m); write_json(run/'manifest.json',m)
    print(f'预览完成：{len(files)} 个文件，HTML 替换 {sum(hits.values())} 处图片；{run}',flush=True)
    return run


def apply(run,yes=False):
    if not yes: raise ValueError('Apply requires --yes')
    run = Path(run).resolve(strict=True); m = json.loads((run/'manifest.json').read_text())
    if m.get('schema')!='celatlas_report_refresh_v2' or m['state']!='prepared': raise ValueError('Not prepared')
    root = Path(m['root']).resolve(strict=True)
    with project_lock(root):
        check_inputs(m)
        for f in m['files']:
            target = checked_target(root,f['relative'],m['sample']); src = run/'staged'/f['relative']
            if not src.resolve().is_relative_to(run/'staged') or digest(src)!=f['new_sha256']:
                raise ValueError('Staged file changed')
            if digest(target)!=f['old_sha256']: raise ValueError('Original changed')
        backup = root/'06.segment/background_swap_backups'/run.name
        need = sum((root/f['relative']).stat().st_size for f in m['files'])
        if shutil.disk_usage(root).free < need*2+1024**3: raise ValueError('Insufficient backup space')
        backup.mkdir(parents=True,exist_ok=False); write_json(backup/'manifest.json',m)
        for f in m['files']:
            saved = backup/'originals'/f['relative']; saved.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(root/f['relative'],saved)
            if digest(saved)!=f['old_sha256']: raise ValueError('Backup verification failed')
        check_inputs(m)
        journal = dict(state='applying',attempted=[]); write_json(backup/'journal.json',journal)
        try:
            for f in m['files']:
                journal['attempted'].append(f['relative']); write_json(backup/'journal.json',journal)
                atomic_copy(run/'staged'/f['relative'],root/f['relative'])
                if digest(root/f['relative'])!=f['new_sha256']: raise ValueError('Published file mismatch')
            targets = {str(root/f['relative']) for f in m['files']}
            check_inputs(dict(inputs=[v for v in m['inputs'] if v['path'] not in targets],protected=[]))
            journal['state']='applied'; write_json(backup/'journal.json',journal)
            m.update(state='applied',backup=str(backup)); write_json(run/'manifest.json',m)
        except Exception:
            for rel in reversed(journal['attempted']): atomic_copy(backup/'originals'/rel,root/rel)
            journal['state']='rolled_back_after_error'; write_json(backup/'journal.json',journal)
            raise
    return backup


def rollback(backup,yes=False):
    if not yes: raise ValueError('Rollback requires --yes')
    backup = Path(backup).resolve(strict=True); m = json.loads((backup/'manifest.json').read_text())
    if m.get('schema')!='celatlas_report_refresh_v2': raise ValueError('Not a v2 report backup')
    root = Path(m['root']).resolve(strict=True)
    with project_lock(root):
        journal = json.loads((backup/'journal.json').read_text())
        if journal['state'] not in ('applied','applying','rolling_back'): raise ValueError('Not rollbackable')
        files = [f for f in m['files'] if f['relative'] in journal['attempted']]
        for f in files:
            target = checked_target(root,f['relative'],m['sample'])
            if digest(target) not in (f['old_sha256'],f['new_sha256']): raise ValueError('Later edits detected')
            if digest(backup/'originals'/f['relative'])!=f['old_sha256']: raise ValueError('Backup damaged')
        journal['state']='rolling_back'; write_json(backup/'journal.json',journal)
        for f in reversed(files): atomic_copy(backup/'originals'/f['relative'],root/f['relative'])
        journal['state']='rolled_back'; write_json(backup/'journal.json',journal)
    return backup


def main():
    p = argparse.ArgumentParser(description='报告、07.outs、h5ad 显示层刷新，不重跑分析')
    sub = p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare'); q.add_argument('--swap-run',required=True); q.add_argument('--output-root',required=True)
    q.add_argument('--bin',type=int,default=50,choices=[10,20,50,100])
    q=sub.add_parser('apply'); q.add_argument('--run',required=True); q.add_argument('--yes',action='store_true')
    q=sub.add_parser('rollback'); q.add_argument('--backup',required=True); q.add_argument('--yes',action='store_true')
    a=p.parse_args()
    if a.command=='prepare': result=prepare(a.swap_run,a.output_root,a.bin)
    elif a.command=='apply': result=apply(a.run,a.yes)
    else: result=rollback(a.backup,a.yes)
    print(f'{a.command} 完成: {result}',flush=True)


if __name__=='__main__': main()
