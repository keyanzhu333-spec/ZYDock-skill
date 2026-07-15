# NVIDIA NIM DiffDock API 参考

## 端点

```
POST https://health.api.nvidia.com/v1/biology/mit/diffdock
```

## 认证

请求头携带 API key（从环境变量 `NVIDIA_API_KEY` 读取）：

```
Authorization: Bearer $NVIDIA_API_KEY
Content-Type: application/json
Accept: application/json
```

获取 key：登录 https://build.nvidia.com → DiffDock 模型 → 生成 key。

## 请求体字段

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `protein` | string | — | **PDB 文件的文本内容本身**（不是路径）。用序列时留空 |
| `protein_sequence` | string | — | 氨基酸序列（不给 PDB 时使用，API 内部用 ESMFold 折叠） |
| `ligand` | string | — | 配体：SMILES 字符串 或 SDF 文本内容 |
| `ligand_file_type` | string | `"txt"` | `"txt"`=SMILES；`"sdf"`=SDF 内容 |
| `num_poses` | int | 20 | 生成位姿数量。难例可增至 30-40 |
| `steps` | int | 18-20 | 扩散逆过程采样步数，增大可能提精度但更慢 |
| `time_divisions` | int | 20 | 时间划分 |
| `save_trajectory` | bool | false | 是否返回扩散中间轨迹 |
| `is_staged` | bool | false | 分阶段模式开关 |

## 响应字段

| 字段 | 说明 |
|---|---|
| `status` | `"success"` 表示成功 |
| `ligand_positions` | list，每项是一个位姿的 **SDF 文本**（含 3D 坐标） |
| `position_confidence` | list，与 `ligand_positions` 一一对应的置信度分数 |
| `trajectory` | 扩散中间轨迹（save_trajectory=true 时） |
| `protein` | 输入蛋白原样回传 |
| `ligand` | 输入配体原样回传 |
| `details` | 状态说明文字 |

## curl 示例

```bash
# 1. 准备 PDB
curl -s https://files.rcsb.org/download/1IEP.pdb -o /tmp/1iep.pdb

# 2. 组请求体（protein 传文本内容）
python3 -c "
import json
json.dump({
    'protein': open('/tmp/1iep.pdb').read(),
    'ligand': 'CC(=O)Oc1ccccc1C(=O)O',
    'ligand_file_type': 'txt',
    'num_poses': 10, 'time_divisions': 20, 'steps': 18,
    'save_trajectory': False, 'is_staged': False
}, open('/tmp/req.json','w'))
"

# 3. 调用
curl -s -X POST https://health.api.nvidia.com/v1/biology/mit/diffdock \
  -H "Authorization: Bearer $NVIDIA_API_KEY" \
  -H "Content-Type: application/json" \
  -d @/tmp/req.json -o /tmp/out.json
```

## 置信度解读

| 分数 | 等级 | 含义 |
|---|---|---|
| > 0 | High | 强预测 |
| -1.5 ~ 0 | Moderate | 合理，需验证 |
| < -1.5 | Low | 不确定 |

置信度反映模型对**结构**的确定性，**不代表结合亲和力**。

## 常见 HTTP 错误

- **401 / 403**：key 无效或过期 → 重新生成
- **422**：请求字段格式错 → 检查 `ligand_file_type` 与配体是否匹配、字段名拼写
- **超时**：大蛋白上传慢，脚本默认 300s 超时
