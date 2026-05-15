# -*- coding: utf-8 -*-
"""
下载 GSE268807 在 NCBI 上的全部分发文件到 02_raw/GSE268807/。

优先使用 **HTTPS**（与浏览器一致），并调用系统 **curl** 支持断点续传（-C -）。
FTP 数据通道在部分 Windows 网络环境下易失败，故不再用 RETR。

基础 URL：
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268807/<sub>/<filename>

子目录：suppl / soft / matrix / miniml

用法：
  python scripts/08_download_geo_gse268807.py
  python scripts/08_download_geo_gse268807.py --dry-run
  python scripts/08_download_geo_gse268807.py --curl curl.exe
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from ftplib import FTP, error_perm
from pathlib import Path


FTP_HOST = "ftp.ncbi.nlm.nih.gov"
HTTPS_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268807"
REMOTE_BASE = "/geo/series/GSE268nnn/GSE268807"
SUBDIRS = ("suppl", "soft", "matrix", "miniml")


def _list_name_to_size(ftp: FTP) -> dict[str, int]:
    lines: list[str] = []
    ftp.retrlines("LIST", lines.append)
    out: dict[str, int] = {}
    for ln in lines:
        parts = ln.split()
        if len(parts) < 9:
            continue
        try:
            sz = int(parts[4])
        except ValueError:
            continue
        name = parts[8]
        if name not in (".", ".."):
            out[name] = sz
    return out


def _curl_download(
    curl_exe: str,
    url: str,
    dest: Path,
    resume: bool,
    curl_extra: list[str],
) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = [curl_exe, "-L", "-f", "--retry", "5", "--retry-delay", "3"]
    cmd.extend(curl_extra)
    if resume:
        cmd.append("-C")
        cmd.append("-")
    cmd.extend(["-o", str(dest), url])
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
        if r.returncode != 0:
            sys.stderr.write(r.stderr or r.stdout or f"curl exit {r.returncode}\n")
            return False
        return dest.is_file() and dest.stat().st_size > 0
    except (subprocess.TimeoutExpired, OSError, FileNotFoundError):
        return False


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--curl", default="curl", help="curl 可执行文件，Windows 可写 curl.exe")
    p.add_argument("--no-resume", action="store_true", help="不使用断点续传（覆盖重下）")
    p.add_argument(
        "--curl-extra",
        action="append",
        default=None,
        help="传给 curl 的额外参数（可多次）。Windows 遇 schannel 错误可试：--curl-extra --ssl-no-revoke",
    )
    args = p.parse_args()
    curl_extra = args.curl_extra if args.curl_extra else []

    root = Path(__file__).resolve().parent.parent
    out_root = root / "02_raw" / "GSE268807"
    out_root.mkdir(parents=True, exist_ok=True)

    (out_root / "README_FTP下载说明.txt").write_text(
        "来源：NCBI GEO（HTTPS 镜像，与 ftp.ncbi.nlm.nih.gov 同步）\n"
        f"基础路径：{HTTPS_BASE}/\n"
        "子目录：suppl（RNA MTX 三联）、soft、matrix、miniml\n"
        "下载脚本：scripts/08_download_geo_gse268807.py（curl 断点续传）\n"
        "ATAC / fragments 常在 SRA，请查 GEO GSE268807 页面。\n"
        "小样本流程：任选一套最小文库（如 G150_D）用 scripts/09_mtx_gz_to_h5ad.py 生成 h5ad，再跑 02→03。\n",
        encoding="utf-8",
    )

    def _p(*a, **k):
        k.setdefault("flush", True)
        print(*a, **k)

    # 每个子目录单独连 FTP 列目录后立即断开，避免长时间 curl 下载导致控制连接 421 超时
    _p("【列目录】FTP", FTP_HOST, "（每子目录独立连接，仅枚举文件名与大小）")

    stats = {"skipped": 0, "downloaded": 0, "failed": 0}
    resume = not args.no_resume

    for sub in SUBDIRS:
        remote_dir = f"{REMOTE_BASE}/{sub}"
        ftp = FTP(FTP_HOST, timeout=120)
        ftp.login()
        ftp.set_pasv(True)
        try:
            ftp.cwd(remote_dir)
            sizes = _list_name_to_size(ftp)
        except error_perm as e:
            _p(f"【警告】无法列出 {remote_dir}: {e}")
            try:
                ftp.quit()
            except Exception:
                pass
            continue
        try:
            ftp.quit()
        except Exception:
            pass

        names = sorted(sizes.keys())
        _p(f"\n【目录】{sub} 文件数={len(names)}")

        for name in names:
            url = f"{HTTPS_BASE}/{sub}/{name}"
            dest = out_root / sub / name
            sz = sizes[name]

            if args.dry_run:
                _p(f"  {name} ({sz / 1e6:.2f} MB) <- {url}")
                continue

            if args.no_resume and dest.is_file():
                dest.unlink()

            use_resume = False
            if dest.is_file():
                loc = dest.stat().st_size
                if loc == sz:
                    _p(f"  [跳过] {name}")
                    stats["skipped"] += 1
                    continue
                if loc > sz:
                    _p(f"  [清理] {name} 本地 {loc} > 远程 {sz}，删除后重下")
                    dest.unlink()
                elif resume and loc < sz:
                    use_resume = True

            _p(f"  [下载] {name} ({sz / 1e6:.2f} MB)")
            ok = _curl_download(
                args.curl, url, dest, resume=use_resume, curl_extra=curl_extra
            )
            if ok and dest.stat().st_size == sz:
                stats["downloaded"] += 1
            elif ok:
                _p(f"    【警告】大小不一致 本地={dest.stat().st_size} 远程={sz}，可删文件后重跑")
                stats["failed"] += 1
            else:
                _p(f"    【失败】{name}")
                stats["failed"] += 1

    _p("\n【汇总】", stats if not args.dry_run else "(dry-run)")
    return 1 if stats.get("failed") else 0


if __name__ == "__main__":
    sys.exit(main())
