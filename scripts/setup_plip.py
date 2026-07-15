#!/usr/bin/env python3
"""
Check if the PLIP conda environment exists. If not, create and install everything.

Usage:
    python scripts/setup_plip.py            # check/create
    python scripts/setup_plip.py --check    # only check, don't create
"""
import argparse
import os
import subprocess
import sys

ENV_NAME = "PLIP"
CONDA = os.environ.get("CONDA_EXE") or "conda"


def conda_envs():
    """Return set of existing conda environment names."""
    out = subprocess.check_output([CONDA, "env", "list"], text=True)
    names = set()
    for line in out.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if parts:
            names.add(parts[0])
    return names


def env_exists():
    return ENV_NAME in conda_envs()


def create_env():
    print(f"[setup] 创建 conda 环境 {ENV_NAME} (python 3.10) ...")
    cmds = [
        f"{CONDA} create -n {ENV_NAME} python=3.10 -y",
        f"{CONDA} run -n {ENV_NAME} conda install -c conda-forge imagemagick -y",
        f"{CONDA} run -n {ENV_NAME} conda install -c conda-forge pymol-open-source -y",
        f"{CONDA} run -n {ENV_NAME} pip install plip",
    ]
    # If pip install plip fails, try source install
    for cmd in cmds:
        print(f"  $ {cmd}")
        ret = os.system(cmd)
        if ret != 0 and "pip install plip" in cmd:
            print("  pip install plip 失败，尝试源码安装...")
            source_cmds = [
                f"{CONDA} run -n {ENV_NAME} git clone https://github.com/pharmai/plip.git /tmp/plip_src",
                f"{CONDA} run -n {ENV_NAME} pip install /tmp/plip_src",
            ]
            for sc in source_cmds:
                print(f"  $ {sc}")
                os.system(sc)
            break

    # verify
    ret = os.system(f"{CONDA} run -n {ENV_NAME} plip -h > /dev/null 2>&1")
    if ret == 0:
        print(f"\n[ok] {ENV_NAME} 环境创建成功，plip 可用。")
    else:
        print(f"\n[warn] 环境已创建但 plip 验证失败，请手动检查。")


def main():
    ap = argparse.ArgumentParser(description="Setup PLIP conda environment")
    ap.add_argument("--check", action="store_true", help="Only check, don't create")
    args = ap.parse_args()

    if env_exists():
        print(f"[ok] conda 环境 '{ENV_NAME}' 已存在，无需创建。")
        # quick verify plip
        ret = os.system(f"{CONDA} run -n {ENV_NAME} plip -h > /dev/null 2>&1")
        if ret == 0:
            print(f"[ok] plip 可用。")
        else:
            print(f"[warn] 环境存在但 plip 不可用，尝试安装...")
            if not args.check:
                os.system(f"{CONDA} run -n {ENV_NAME} pip install plip")
        return 0

    if args.check:
        print(f"[info] conda 环境 '{ENV_NAME}' 不存在。用 --check 模式不自动创建。")
        return 1

    create_env()
    return 0


if __name__ == "__main__":
    sys.exit(main())
