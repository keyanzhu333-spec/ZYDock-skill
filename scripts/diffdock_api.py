#!/usr/bin/env python3
"""
Shared helper for calling the hosted NVIDIA NIM DiffDock API.

API key handling — the key is NEVER stored inside this skill directory.
Resolution order:
  1. NVIDIA_API_KEY environment variable (highest priority, for temporary override)
  2. User credentials file OUTSIDE the skill: ~/.config/ZYDock/credentials
  3. First-run interactive prompt (if a terminal is attached): asks the user,
     then saves the key to the credentials file (chmod 600) for future runs.

All docking scripts in this skill import from here.
"""
import json
import os
import stat
import sys
import urllib.request
import urllib.error

NIM_URL = "https://health.api.nvidia.com/v1/biology/mit/diffdock"

# Credentials live OUTSIDE the skill dir so the skill stays key-free and shareable.
CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".config", "ZYDock")
CRED_FILE = os.path.join(CONFIG_DIR, "credentials")


def _read_cred_file():
    """Return the saved key from the credentials file, or None."""
    try:
        with open(CRED_FILE) as f:
            key = f.read().strip()
        return key or None
    except OSError:
        return None


def save_api_key(key: str):
    """Save the key to the user credentials file (outside the skill), chmod 600."""
    key = key.strip()
    os.makedirs(CONFIG_DIR, exist_ok=True)
    with open(CRED_FILE, "w") as f:
        f.write(key + "\n")
    os.chmod(CRED_FILE, stat.S_IRUSR | stat.S_IWUSR)  # 0600, owner-only
    return CRED_FILE


def get_api_key(allow_prompt: bool = True) -> str:
    """
    Resolve the NVIDIA API key.
      1) env var  2) saved credentials file  3) first-run interactive prompt.
    The key is never stored inside the skill directory.
    """
    # 1) environment variable
    key = os.environ.get("NVIDIA_API_KEY")
    if key:
        return key.strip()

    # 2) saved credentials file (set on a previous first run)
    key = _read_cred_file()
    if key:
        return key

    # 3) first run: ask the user, then remember it
    if allow_prompt and sys.stdin.isatty():
        sys.stderr.write(
            "\n首次使用 ZYDock —— 需要 NVIDIA API key。\n"
            "  本技能全部计算在 NVIDIA 云端完成，需要一个 key。\n"
            "  获取地址: https://build.nvidia.com  (形如 nvapi-xxxx)\n"
            "  输入的 key 将保存到 " + CRED_FILE + " (权限 600)，之后自动使用。\n\n"
        )
        try:
            entered = input("请粘贴你的 NVIDIA API key: ").strip()
        except (EOFError, KeyboardInterrupt):
            entered = ""
        if entered:
            path = save_api_key(entered)
            sys.stderr.write(f"✓ 已保存到 {path}，后续无需再次输入。\n\n")
            return entered

    # nothing worked
    sys.exit(
        "ERROR: 未找到 NVIDIA API key。\n"
        "  三种提供方式(任选其一):\n"
        "    1) 交互式:直接运行脚本，按提示粘贴 key(会记住)\n"
        "    2) 命令保存: python scripts/set_key.py nvapi-xxxxxxxx\n"
        '    3) 环境变量: export NVIDIA_API_KEY="nvapi-xxxxxxxx"\n'
        "  获取 key: https://build.nvidia.com\n"
    )


def clean_protein_text(pdb_text: str) -> tuple:
    """
    Remove all HETATM records (waters, ions, co-crystallized ligands, crystallization
    additives) from protein PDB text. Returns (cleaned_text, n_removed).
    Should be called BEFORE docking so the receptor is pure protein only.
    """
    lines = []
    removed = 0
    for line in pdb_text.splitlines(keepends=True):
        if line.startswith("HETATM"):
            removed += 1
            continue
        if line.startswith("CONECT"):
            continue
        lines.append(line)
    return "".join(lines), removed


def read_protein(protein_path: str = None, protein_sequence: str = None, clean: bool = True):
    """
    Return (protein_text, is_sequence). Exactly one input must be given.
    If clean=True (default), removes HETATM records from PDB files before returning.
    """
    if protein_path and protein_sequence:
        sys.exit("ERROR: provide either --protein OR --protein_sequence, not both.")
    if protein_path:
        with open(protein_path) as f:
            text = f.read()
        if clean:
            text, n = clean_protein_text(text)
            if n > 0:
                print(f"[clean] 已去除受体中 {n} 个杂原子记录（水/离子/共结晶配体/添加剂），保留纯蛋白")
        return text, False
    if protein_sequence:
        return protein_sequence, True
    sys.exit("ERROR: no protein given. Use --protein <pdb> or --protein_sequence <seq>.")


def read_ligand(ligand: str):
    """
    Return (ligand_text, ligand_file_type).
    If `ligand` points to an existing .sdf/.mol2 file, read it; else treat as SMILES.
    """
    if os.path.isfile(ligand):
        ext = os.path.splitext(ligand)[1].lower()
        with open(ligand) as f:
            content = f.read()
        ftype = "sdf" if ext in (".sdf", ".mol", ".mol2") else "txt"
        return content, ftype
    # treat as a raw SMILES string
    return ligand, "txt"


def dock(protein_text, ligand_text, ligand_file_type="txt",
         num_poses=20, steps=18, time_divisions=20,
         is_sequence=False, timeout=300):
    """
    Call the DiffDock NIM API once and return the parsed JSON dict.
    Raises RuntimeError on HTTP error.
    """
    key = get_api_key()

    payload = {
        "ligand": ligand_text,
        "ligand_file_type": ligand_file_type,
        "num_poses": num_poses,
        "time_divisions": time_divisions,
        "steps": steps,
        "save_trajectory": False,
        "is_staged": False,
    }
    # protein by structure OR by sequence
    if is_sequence:
        payload["protein"] = ""
        payload["protein_sequence"] = protein_text
    else:
        payload["protein"] = protein_text

    req = urllib.request.Request(
        NIM_URL,
        data=json.dumps(payload).encode(),
        headers={
            "Authorization": f"Bearer {key}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode()
    except urllib.error.HTTPError as e:
        detail = e.read().decode()[:2000]
        raise RuntimeError(f"HTTP {e.code}: {detail}")
    except urllib.error.URLError as e:
        raise RuntimeError(f"Network error reaching {NIM_URL}: {e.reason}")

    return json.loads(body)


def save_poses(data: dict, out_dir: str):
    """
    Extract ranked poses + confidences from an API response and write:
      rank{N}_conf{score}.sdf  and  confidence_scores.txt  and  raw_response.json
    Returns list of (rank, pose_index, confidence).
    """
    os.makedirs(out_dir, exist_ok=True)
    poses = data.get("ligand_positions") or []
    conf = data.get("position_confidence") or []

    with open(os.path.join(out_dir, "raw_response.json"), "w") as f:
        json.dump(data, f, indent=2)

    if not poses:
        return []

    order = sorted(range(len(conf)), key=lambda i: conf[i], reverse=True) if conf \
        else list(range(len(poses)))

    ranking = []
    lines = []
    for rank, i in enumerate(order, 1):
        c = conf[i] if conf else float("nan")
        fname = f"rank{rank}_conf{c:.3f}.sdf" if conf else f"rank{rank}.sdf"
        with open(os.path.join(out_dir, fname), "w") as f:
            f.write(poses[i])
        level = confidence_level(c) if conf else "N/A"
        lines.append(f"rank{rank}\tpose{i+1}\t{c:.4f}\t{level}")
        ranking.append((rank, i + 1, c))

    with open(os.path.join(out_dir, "confidence_scores.txt"), "w") as f:
        f.write("rank\tpose\tconfidence\tlevel\n")
        f.write("\n".join(lines) + "\n")

    return ranking


def confidence_level(c: float) -> str:
    if c > 0:
        return "High"
    if c > -1.5:
        return "Moderate"
    return "Low"
