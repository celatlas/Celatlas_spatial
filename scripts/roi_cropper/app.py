#!/usr/bin/env python3
"""Local browser UI and CLI for reproducible native-resolution ROI crops."""
from __future__ import annotations

import argparse
import base64
import json
import secrets
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import cv2

try:
    from .core import TiffSource, calibration_from_points, export_crop, make_recipe, params_from_recipe, png_bytes
except ImportError:
    from core import TiffSource, calibration_from_points, export_crop, make_recipe, params_from_recipe, png_bytes


class CropApp:
    def __init__(self, image, output_root, calibration=None):
        self.image = str(Path(image).expanduser().resolve()) if image else None
        self.output_root = str(Path(output_root).expanduser().resolve())
        self.calibration = json.loads(Path(calibration).read_text()) if calibration else None
        self.token = secrets.token_urlsafe(32)
        self.sessions, self.jobs = {}, {}
        self.lock = threading.Lock()
        self.pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix='roi-export')

    def open_image(self, path):
        if not isinstance(path, (str, Path)) or not str(path).strip():
            raise ValueError('请先输入服务器上的图像路径。')
        source = TiffSource(path)
        preview = source.thumbnail()
        session = uuid.uuid4().hex
        calibration = None
        if self.calibration:
            expected = self.calibration['for_image']
            if (source.path == Path(expected['path']).resolve()
                    and source.width == expected['width'] and source.height == expected['height']
                    and source.signature['bytes'] == expected['bytes']):
                calibration = dict(self.calibration['calibration'])
        with self.lock:
            if len(self.sessions) >= 12:
                raise ValueError('已打开 12 张图，请重启工具释放缓存。')
            self.sessions[session] = source
        return dict(session=session, source=source.signature, preview='data:image/png;base64,'+base64.b64encode(png_bytes(preview)).decode(),
                    preview_size=[preview.shape[1], preview.shape[0]], calibration=calibration,
                    resolution_tags=source.resolution_tags, storage_format=source.storage_format)

    def source(self, session):
        if session not in self.sessions:
            raise ValueError('图像会话已失效，请重新打开。')
        return self.sessions[session]

    def submit(self, session, params):
        source = self.source(session)
        recipe = make_recipe(source, params)
        if recipe['outside_source'] and not recipe['allow_black_padding']:
            raise ValueError('选框超出原图，请调整位置，或明确允许黑色补边。')
        with self.lock:
            if any(j['state'] in ('queued','running') for j in self.jobs.values()):
                raise ValueError('已有导出任务正在运行，请等待完成。')
            job = uuid.uuid4().hex
            self.jobs[job] = {'state': 'queued', 'percent': 0}
        def update(data):
            with self.lock:
                self.jobs[job].update(data)
        def run():
            update({'state': 'running'})
            try:
                result = export_crop(source, recipe, self.output_root, progress=update)
                update({'state': 'done', 'percent': 100, 'recipe': result,
                        'output_dir': result['output']['directory']})
            except Exception as exc:
                update({'state': 'error', 'error': str(exc)})
        self.pool.submit(run)
        return {'job': job, 'recipe': recipe}


def handler_factory(app):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt, *args):
            # Do not log the authentication token embedded in the launch URL.
            print(f'[{self.log_date_time_string()}] {self.command} {urlparse(self.path).path}', flush=True)

        def authorized(self):
            query = parse_qs(urlparse(self.path).query)
            token = self.headers.get('X-Crop-Token') or query.get('token', [''])[0]
            origin = self.headers.get('Origin')
            same_origin = not origin or urlparse(origin).netloc == self.headers.get('Host')
            return same_origin and secrets.compare_digest(token, app.token)

        def send(self, data, status=200, content_type='application/json; charset=utf-8'):
            if not isinstance(data, bytes):
                data = json.dumps(data, ensure_ascii=False, allow_nan=False).encode()
            self.send_response(status)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('Referrer-Policy', 'no-referrer')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            if not self.authorized():
                return self.send({'error':'请使用终端打印的带 token 链接访问。'}, 403)
            parsed = urlparse(self.path)
            if parsed.path == '/':
                return self.send(Path(__file__).with_name('ui.html').read_bytes(), content_type='text/html; charset=utf-8')
            if parsed.path == '/api/config':
                return self.send({'image':app.image, 'output_root':app.output_root})
            if parsed.path == '/api/status':
                job = parse_qs(parsed.query).get('job',[''])[0]
                with app.lock:
                    state = dict(app.jobs.get(job, {'error':'未知任务'}))
                return self.send(state)
            self.send({'error':'Not found'}, 404)

        def do_POST(self):
            if not self.authorized():
                return self.send({'error':'Unauthorized'}, 403)
            try:
                length = int(self.headers.get('Content-Length', 0))
                if not 0 < length <= 1_000_000:
                    raise ValueError('请求长度无效。')
                data = json.loads(self.rfile.read(length))
                path = urlparse(self.path).path
                if path == '/api/open':
                    return self.send(app.open_image(data['path']))
                if path == '/api/recipe':
                    return self.send(make_recipe(app.source(data['session']), data['params']))
                if path == '/api/calibrate':
                    source = app.source(data['session'])
                    for x,y in data['points']:
                        if not (-.5 <= x <= source.width-.5 and -.5 <= y <= source.height-.5):
                            raise ValueError('标尺端点超出原图。')
                    return self.send(calibration_from_points(data['points'], data['length_um']))
                if path == '/api/detail':
                    source = app.source(data['session'])
                    box = [int(v) for v in data['box']]
                    preview = source.thumbnail(1800, box=box)
                    return self.send({'box':box,'image':'data:image/png;base64,'+base64.b64encode(png_bytes(preview)).decode()})
                if path == '/api/export':
                    return self.send(app.submit(data['session'], data['params']))
                self.send({'error':'Not found'}, 404)
            except Exception as exc:
                self.send({'error':str(exc)}, 400)
    return Handler


def main():
    parser = argparse.ArgumentParser(description='原生像素旋转、6×6 mm 裁切与参数记录（独立于 Celatlas 流程）')
    sub = parser.add_subparsers(dest='command', required=True)
    serve = sub.add_parser('serve', help='启动本机网页选区界面')
    serve.add_argument('--image', help='可选：启动后自动打开的图片；不填写则打开空白页面')
    serve.add_argument('--output-root', default='crop_exports')
    serve.add_argument('--calibration', help='可选：已测量的标尺 JSON，仍需要在页面中确认')
    serve.add_argument('--port', type=int, default=8877)
    replay = sub.add_parser('replay', help='对同坐标原始通道套用保存的裁切参数')
    replay.add_argument('--image', required=True)
    replay.add_argument('--recipe', required=True)
    replay.add_argument('--output-root', default='crop_exports')
    replay.add_argument('--confirm-shared-canvas', action='store_true',
                        help='确认不同通道具有完全相同的原始画布/像素间距/方向/起点；尺寸相同本身并不能证明对齐')
    args = parser.parse_args()
    cv2.setNumThreads(2)
    if args.command == 'serve':
        # No network exposure by default. Access a server through SSH forwarding.
        if args.image:
            TiffSource(args.image)
        app = CropApp(args.image, args.output_root, args.calibration)
        server = ThreadingHTTPServer(('127.0.0.1', args.port), handler_factory(app))
        print(f'浏览器打开 http://127.0.0.1:{server.server_port}/?token={app.token}', flush=True)
        print(f'导出目录：{app.output_root}；不修改原图，不调用 Celatlas 分割/分析。', flush=True)
        if not args.image:
            print('空白启动：在网页中输入服务器图像路径，再点击“打开新图并重置”。', flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print('正在关闭；若有导出任务，等待文件写完。', flush=True)
        finally:
            server.server_close()
            app.pool.shutdown(wait=True)
    else:
        source = TiffSource(args.image)
        recipe = json.loads(Path(args.recipe).read_text())
        original = recipe['source']
        if [source.width,source.height] != [original['width'],original['height']]:
            raise ValueError('两个通道原图尺寸不同，不能直接套用同一裁切坐标；请先恢复共同画布或配准。')
        if source.signature != original and not args.confirm_shared_canvas:
            raise ValueError('源图不同；核实同坐标画布后才可添加 --confirm-shared-canvas。')
        updated = make_recipe(source, params_from_recipe(recipe))
        updated['replayed_from_recipe'] = str(Path(args.recipe).resolve())
        updated['shared_canvas_user_confirmed'] = bool(args.confirm_shared_canvas)
        result = export_crop(source, updated, args.output_root,
                             progress=lambda v: print(f"\r导出 {v['percent']}%", end='', flush=True))
        print('\n'+result['output']['directory'])


if __name__ == '__main__':
    main()
