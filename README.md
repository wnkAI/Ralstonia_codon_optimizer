# Ralstonia Codon Optimizer

A reverse-translation and codon-optimization tool tailored for ***Ralstonia eutropha* H16** (*Cupriavidus necator* H16).

Given a protein sequence, it produces a DNA sequence using the host's codon usage bias, while avoiding restriction sites, GC-content extremes, long homopolymers, and Shine-Dalgarno-like motifs.

## Requirements

```bash
pip install biopython
```

Python 3.8+.

## Usage

### Interactive mode (recommended)

```bash
python codon_optimizer.py
```

Walks you through every parameter:
1. Protein input (FASTA path or pasted sequence)
2. Forbidden enzymes: `a` common cloning sites / `b` full list / `c` custom names / `n` none
3. Auto-add ATG start / stop codon
4. GC bounds (default 40%–75%)
5. Max homopolymer length (default 6) + AGGAGG (Shine-Dalgarno) avoidance
6. Output folder

### CLI mode

```bash
python codon_optimizer.py -i protein.fasta --enzymes common --gc-low 45 --gc-high 75
```

Custom enzyme set:

```bash
python codon_optimizer.py -i protein.fasta --enzymes custom --custom-enzymes "BsaI,BbsI"
```

| Flag | Default | Meaning |
|---|---|---|
| `-i` | (required) | FASTA path or raw protein-sequence string |
| `-o` | `E:/codon_optimizer/output` | Output base folder |
| `--enzymes` | `common` | `common` / `full` / `custom` / `none` |
| `--custom-enzymes` | `""` | Comma-separated names, used with `--enzymes custom` |
| `--gc-low` | 40 | GC% lower bound |
| `--gc-high` | 75 | GC% upper bound |
| `--max-homopolymer` | 6 | Max consecutive identical nt allowed |
| `--no-sd-check` | off | Skip AGGAGG (Shine-Dalgarno) avoidance |
| `--no-atg` | off | Do not auto-add ATG start codon |
| `--no-stop` | off | Do not auto-add stop codon |
| `--seed` | None | RNG seed for reproducibility |
| `--min-site-length` | 6 | Skip enzyme sites shorter than this many nt |

## Output

Every run writes two files into `output/<name>/`:
- `<name>_optimized.fasta` — the optimized DNA sequence
- `<name>_report.txt` — GC%, translation check, removed-site list, codon-usage breakdown

## Built-in enzyme sets

- `common` — 13 frequent cloning enzymes: EcoRI, BamHI, HindIII, XhoI, NdeI, NcoI, NotI, XbaI, SpeI, PstI, KpnI, SacI, SalI
- `full` — 200+ enzymes parsed from a built-in concatenated name string, resolved via Biopython's `Restriction_Dictionary`
- `custom` — comma-separated names you provide
- `none` — disable site avoidance

### Practical guidance

| Scenario | Recommended setting |
|---|---|
| Everyday gene of any length | `--enzymes common` (13 sites) |
| Short gene (< 100 aa), strict design | `--enzymes full` (default `--min-site-length 6`, ~115 sites) |
| Long gene (> 300 aa), strict design | `--enzymes full --min-site-length 7` (~30 sites) |
| Optimization fails | Raise `--min-site-length`, switch to `common`, or widen the GC range |

**About `full` mode:** sites of length ≤ 5 (including degenerate 5-cutters such as `CCNGG`, `CCWGG`) are skipped by default. In a high-GC long sequence, 4-cutters appear once every ~256 bp on average, so avoiding all of them is mathematically infeasible. Pass `--min-site-length 1` for strict mode (only realistic for peptides < 100 aa).

## Algorithm

Weighted random codon picking with local backtracking:
1. For each amino acid, draw a candidate codon weighted by host usage frequency.
2. After each codon, check the trailing 50-nt window:
   - No forbidden restriction site (both strands, IUPAC-aware regex)
   - No homopolymer longer than the threshold
   - No AGGAGG
   - GC content within bounds
3. On failure, retry with the remaining codons; if all fail at this position, undo 2–8 codons and resume.
4. Restart from scratch up to 200 times if backtracking exhausts.
5. Final check: full-sequence scan for forbidden sites + reverse-translation match against the input protein.

## Codon usage table

Approximate fractions per amino acid derived from the Kazusa Codon Usage Database for *Cupriavidus necator* H16 (NCBI taxid 381666). High-GC bias is reflected in the dominance of GCG (Ala), CGC (Arg), AAC (Asn), GAC (Asp), GGC (Gly), CTG (Leu), CCG (Pro), and so on. To update with newer data, edit the `CODON_USAGE` dict at the top of `codon_optimizer.py`.

## Repository layout

```
codon_optimizer/
├── codon_optimizer.py   # main script (single file)
├── README.md            # this file
└── .gitignore
```

## License

MIT (or update as needed).
