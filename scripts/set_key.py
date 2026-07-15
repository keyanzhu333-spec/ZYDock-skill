#!/usr/bin/env python3
"""
Manage the saved NVIDIA API key for the ZYDock skill.

The key is stored OUTSIDE the skill directory, in:
    ~/.config/ZYDock/credentials   (chmod 600, owner-only)
so the skill itself never contains any key and stays safe to share.

Usage:
    python scripts/set_key.py nvapi-xxxxxxxx      # save a key
    python scripts/set_key.py                     # prompt for a key, then save
    python scripts/set_key.py --show              # show status (masked)
    python scripts/set_key.py --clear             # delete the saved key
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diffdock_api import save_api_key, CRED_FILE, _read_cred_file  # noqa: E402


def mask(key: str) -> str:
    if len(key) <= 12:
        return key[:4] + "..."
    return f"{key[:8]}...{key[-4:]}"


def main():
    ap = argparse.ArgumentParser(description="Manage the NVIDIA API key for ZYDock")
    ap.add_argument("key", nargs="?", help="API key to save (nvapi-...)")
    ap.add_argument("--show", action="store_true", help="Show saved key status (masked)")
    ap.add_argument("--clear", action="store_true", help="Delete the saved key")
    args = ap.parse_args()

    if args.clear:
        if os.path.isfile(CRED_FILE):
            os.remove(CRED_FILE)
            print(f"✓ 已删除保存的 key: {CRED_FILE}")
        else:
            print("没有已保存的 key。")
        return 0

    if args.show:
        key = _read_cred_file()
        if key:
            print(f"已保存 key: {mask(key)}")
            print(f"位置: {CRED_FILE}")
        else:
            print("尚未保存 key。运行 `python scripts/set_key.py nvapi-...` 保存。")
        return 0

    key = args.key
    if not key:
        try:
            key = input("请粘贴你的 NVIDIA API key (nvapi-...): ").strip()
        except (EOFError, KeyboardInterrupt):
            key = ""
    if not key:
        sys.exit("未输入 key，已取消。")

    if not key.startswith("nvapi-"):
        print(f"⚠ 警告: key 通常以 'nvapi-' 开头，你输入的以 '{key[:6]}' 开头。仍继续保存。")

    path = save_api_key(key)
    print(f"✓ 已保存到 {path} (权限 600，仅本人可读)")
    print("  之后运行对接脚本会自动使用，无需再次输入。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
