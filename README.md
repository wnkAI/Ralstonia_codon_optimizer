# Codon Optimizer for *Ralstonia eutropha* H16

针对 *Ralstonia eutropha* H16（*Cupriavidus necator* H16）的反向翻译 + 密码子优化工具。

## 安装依赖

```bash
pip install biopython
```

## 使用方法

### 1. 交互模式（推荐）

直接跑，按提示选参数：

```bash
python E:/codon_optimizer/codon_optimizer.py
```

会依次询问：
1. 蛋白序列（FASTA 路径 或 直接粘贴）
2. 禁用酶切位点：`a` 常见克隆位点 / `b` 完整列表 / `c` 自定义 / `n` 不禁用
3. 是否补 ATG 起始 / 终止密码子
4. GC 上下界（默认 45%-75%）
5. 最大连续同碱基（默认 6）+ 是否避免 SD 序列（AGGAGG）
6. 输出目录

### 2. 命令行模式

```bash
python codon_optimizer.py -i protein.fasta --enzymes common --gc-low 45 --gc-high 75
```

主要参数：

| 参数 | 默认 | 说明 |
|---|---|---|
| `-i` | (必填) | FASTA 文件 或 直接传蛋白序列字符串 |
| `-o` | `E:/codon_optimizer/output` | 输出基目录 |
| `--enzymes` | `common` | `common` / `full` / `none` |
| `--gc-low` | 45 | GC 下界（%） |
| `--gc-high` | 75 | GC 上界（%） |
| `--max-homopolymer` | 6 | 允许的最大连续同碱基 |
| `--no-sd-check` | off | 关闭 AGGAGG 检查 |
| `--no-atg` | off | 不自动加 ATG |
| `--no-stop` | off | 不自动加终止密码子 |
| `--seed` | None | 固定随机种子（结果可复现） |

## 输出

每次跑完会在 `output/<name>/` 下生成：
- `<name>_optimized.fasta` —— 优化后的 DNA 序列
- `<name>_report.txt` —— GC%、翻译验证、移除的酶切位点列表、密码子使用统计

## 内置酶切列表

- `common`：13 个常见克隆酶（EcoRI / BamHI / HindIII / XhoI / NdeI / NcoI / NotI / XbaI / SpeI / PstI / KpnI / SacI / SalI）
- `full`：完整列表（用户提供的 200+ 酶名串），通过 Biopython `Restriction_Dictionary` 解析

### 实用建议

| 场景 | 推荐设置 |
|---|---|
| 日常基因（任意长度） | `--enzymes common`（13 个位点） |
| 短基因（< 100 aa）严格设计 | `--enzymes full`（默认 `--min-site-length 6`，~115 位点） |
| 长基因（> 300 aa）严格设计 | `--enzymes full --min-site-length 7`（~30 位点） |
| 优化失败时 | 提高 `--min-site-length` / 改用 `common` / 拉宽 `--gc-low/--gc-high` |

**关于 `full` 模式**：默认跳过 ≤5 nt 位点（含退化 5-cutter 如 CCNGG/CCWGG）。原因：高 GC 长序列下 4-cutter 平均每 256bp 出现一次，几乎不可能全部避开。需严格模式可加 `--min-site-length 1`，但仅推荐用于 < 100 aa 短肽。

## 算法

加权随机选择 + 局部回溯：
1. 按密码子使用频率加权随机抽
2. 每加一个密码子检查最近 50 nt 窗口：禁用酶切位点（双链）/ 同碱基长串 / SD 序列 / GC 滑窗
3. 失败回溯 1-3 个密码子重选；多次失败重启（最多 200 次）
4. 完成后做全局校验 + 反向翻译比对蛋白

## 密码子使用表

来自 Kazusa Codon Usage Database (taxid 381666) 的近似值。GC-rich 偏好显著（GCG/CGC/AAC/GAC/CTG/CCG 等占主导）。如需更新，编辑 `CODON_USAGE` 字典即可。
