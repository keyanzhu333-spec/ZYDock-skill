#!/usr/bin/env python3
"""
Check that the NVIDIA DiffDock API is reachable and the key works.
Runs a tiny real docking (aspirin + a small test protein) to confirm end-to-end.

Usage:
    export NVIDIA_API_KEY="nvapi-xxxxxxxx"
    python scripts/check_api.py
"""
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diffdock_api import get_api_key, dock, NIM_URL  # noqa: E402

TEST_LIGAND = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin


def main():
    print("=" * 60)
    print("DiffDock_ZY  API 连通性检查")
    print("=" * 60)

    # 1. key present?
    print("\n[1/3] 检查 NVIDIA_API_KEY 环境变量...")
    try:
        key = get_api_key()
        print(f"  ✓ 已设置 (以 {key[:8]}... 开头)")
    except SystemExit as e:
        print(e)
        return 1

    # 2. fetch a small test protein
    print("\n[2/3] 下载测试蛋白 (RCSB 1IEP)...")
    try:
        with urllib.request.urlopen(
            "https://files.rcsb.org/download/1IEP.pdb", timeout=60
        ) as r:
            protein = r.read().decode()
        print(f"  ✓ 下载成功 ({len(protein)} 字符)")
    except Exception as e:
        print(f"  ✗ 下载失败: {e}")
        print("    (网络问题，但不影响你用本地已有的 PDB 文件对接)")
        return 1

    # 3. real API call
    print(f"\n[3/3] 调用 API 做一次测试对接 ...\n  端点: {NIM_URL}")
    try:
        data = dock(protein, TEST_LIGAND, num_poses=5, steps=18)
    except RuntimeError as e:
        print(f"  ✗ API 调用失败: {e}")
        if "401" in str(e) or "403" in str(e):
            print("    → key 无效或过期，请去 build.nvidia.com 重新生成")
        return 1

    status = data.get("status")
    poses = data.get("ligand_positions") or []
    conf = data.get("position_confidence") or []
    print(f"  ✓ 成功! status={status}, 返回 {len(poses)} 个位姿")
    if conf:
        print(f"    最高置信度: {max(conf):.3f}")

    print("\n" + "=" * 60)
    print("✓ 全部通过 —— API 可用，可以开始对接了。")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
