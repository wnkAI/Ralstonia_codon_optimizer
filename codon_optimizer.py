#!/usr/bin/env python3
"""
Codon optimizer for Ralstonia eutropha H16 (Cupriavidus necator H16).

Usage:
    python codon_optimizer.py                 # interactive mode
    python codon_optimizer.py -i prot.fasta   # CLI mode
"""

import argparse
import random
import re
import sys
from pathlib import Path

try:
    from Bio.Restriction import Restriction_Dictionary
except ImportError:
    print("Error: biopython not installed. Run: pip install biopython")
    sys.exit(1)


# Codon usage for Ralstonia eutropha H16 (Cupriavidus necator H16, taxid 381666).
# Approximate fractions per amino acid derived from Kazusa Codon Usage Database.
# High-GC genome (~66%) so GC-rich codons dominate.
CODON_USAGE = {
    'A': [('GCG', 0.49), ('GCC', 0.31), ('GCT', 0.10), ('GCA', 0.10)],
    'R': [('CGC', 0.55), ('CGG', 0.18), ('CGT', 0.16), ('AGG', 0.04), ('CGA', 0.04), ('AGA', 0.03)],
    'N': [('AAC', 0.82), ('AAT', 0.18)],
    'D': [('GAC', 0.70), ('GAT', 0.30)],
    'C': [('TGC', 0.70), ('TGT', 0.30)],
    'Q': [('CAG', 0.80), ('CAA', 0.20)],
    'E': [('GAG', 0.68), ('GAA', 0.32)],
    'G': [('GGC', 0.50), ('GGT', 0.20), ('GGG', 0.20), ('GGA', 0.10)],
    'H': [('CAC', 0.70), ('CAT', 0.30)],
    'I': [('ATC', 0.70), ('ATT', 0.27), ('ATA', 0.03)],
    'L': [('CTG', 0.55), ('CTC', 0.20), ('CTT', 0.10), ('TTG', 0.10), ('CTA', 0.03), ('TTA', 0.02)],
    'K': [('AAG', 0.70), ('AAA', 0.30)],
    'M': [('ATG', 1.00)],
    'F': [('TTC', 0.70), ('TTT', 0.30)],
    'P': [('CCG', 0.53), ('CCC', 0.27), ('CCT', 0.10), ('CCA', 0.10)],
    'S': [('TCG', 0.30), ('AGC', 0.25), ('TCC', 0.20), ('TCT', 0.10), ('AGT', 0.10), ('TCA', 0.05)],
    'T': [('ACC', 0.45), ('ACG', 0.35), ('ACT', 0.10), ('ACA', 0.10)],
    'W': [('TGG', 1.00)],
    'Y': [('TAC', 0.70), ('TAT', 0.30)],
    'V': [('GTG', 0.50), ('GTC', 0.30), ('GTT', 0.15), ('GTA', 0.05)],
    '*': [('TGA', 0.60), ('TAA', 0.30), ('TAG', 0.10)],
}

TRANSLATION_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

COMMON_CLONING_SITES = [
    'EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NdeI', 'NcoI', 'NotI',
    'XbaI', 'SpeI', 'PstI', 'KpnI', 'SacI', 'SalI',
]

# Full enzyme list as concatenated string (user-supplied).
FULL_ENZYME_STRING = (
    "AatIIAcc65IAccB7IAccIAccIIAccIIIAciIAcyIAflIIAflIIIAgeIAhaIIIAluIAlw21I"
    "Alw26IAlw44IAlwIAlwNIAor51HIApaIApaLIApoIAscIAseIAsp700IAsuIAsuIIAvaI"
    "AvaIIAvrIIBalIBamHIBanIBanIIBbsIBbvIBclIBcnIBfaIBglIBglIIBlnIBpu1102I"
    "BsaAIBsaBIBsaHIBsaIBsaJIBsaWIBsh1236IBsiEIBsiHKAIBsiWIBsiYIBslIBsmAI"
    "BsmIBsp120IBsp1286IBspDIBspEIBspHIBspMIBspMIIBsrBIBsrDIBsrFIBsrGIBsrI"
    "BssHIIBstBIBstEIIBstNIBstOIBstUIBstXIBstYIBstZIBsu36ICac8ICfoICfr10I"
    "Cfr13ICfrIClaICpoICsp6ICspIDdeIDpnIDpnIIDraIDraIIDraIIIDrdIDsaIEaeI"
    "EagIEam1105IEarIEco47IEco47IIIEco52IEco72IEco81IEcoICRIEcoNIEcoO109I"
    "EcoRIEcoRIIEcoRVEsp3IEspIFnu4HIFnuDIIFseIFspIHaeIIHaeIIIHapIIHgaIHgiAI"
    "HhaIHincIIHindIIHindIIIHinfIHpaIHpaIIHphIKasIKpn2IKpnIMaeIMaeIIMaeIII"
    "MboIMboIIMcrIMflIMluIMnlIMroIMscIMseIMslIMspIMstIMunIMvaIMwoINaeINarI"
    "NciINcoINdeINdeIINheINlaIIINlaIVNotINruINspBIINspINspVPacIPaeR7IPflMI"
    "PinAIPleIPmaCIPmeIPmlIPpuMIPsp1406IPssIPstIPvuIPvuIIRsaIRsrIISacISacII"
    "SalISapISau3AISau96ISauIScaIScrFISduISecISfaNISfcISfiISfuISgfISmaI"
    "SnaBISpeISphISplISrfISse8387ISspISstIStuIStyITaqITfiITru9IVspIXbaI"
    "XcmIXhoIXhoIIXmaIXmaIIIXmnI"
)


# IUPAC degenerate base -> regex character class
IUPAC_TO_REGEX = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[CG]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]',
    'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
    'N': '[ACGT]',
}


def reverse_complement(seq: str) -> str:
    table = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN')
    return seq.translate(table)[::-1]


def site_to_regex(site: str) -> re.Pattern:
    site = site.upper()
    site = ''.join(c for c in site if c in IUPAC_TO_REGEX)
    if not site:
        return re.compile('(?!)')  # match nothing
    fwd = ''.join(IUPAC_TO_REGEX[b] for b in site)
    rev = ''.join(IUPAC_TO_REGEX[b] for b in reverse_complement(site))
    if fwd == rev:
        return re.compile(fwd)
    return re.compile(f'(?:{fwd})|(?:{rev})')


def split_concatenated_enzyme_names(s: str):
    s = s.strip().replace(' ', '').replace('\n', '')
    known = sorted(Restriction_Dictionary.rest_dict.keys(), key=len, reverse=True)
    out, i = [], 0
    skipped = 0
    while i < len(s):
        matched = None
        for name in known:
            if s.startswith(name, i):
                matched = name
                break
        if matched is None:
            skipped += 1
            i += 1
            continue
        out.append(matched)
        i += len(matched)
    if skipped:
        print(f"  Note: {skipped} characters could not be parsed (skipped).")
    return out


def get_forbidden_sites(enzyme_names, min_site_length: int = 6):
    """Map enzymes to recognition sites; skip sites shorter than min_site_length
    (default 6, since 4-/5-cutters and degenerate 5-cutters are usually
    unavoidable in long sequences)."""
    rest = Restriction_Dictionary.rest_dict
    sites = {}
    missing = []
    too_short = []
    for name in enzyme_names:
        if name not in rest:
            missing.append(name)
            continue
        site = rest[name].get('site')
        if not site or '?' in site:
            continue
        site = site.upper().replace('U', 'T')
        site = ''.join(c for c in site if c in IUPAC_TO_REGEX)
        if not site:
            continue
        if len(site) < min_site_length:
            too_short.append((name, site))
            continue
        sites.setdefault(site, []).append(name)
    return sites, missing, too_short


def has_homopolymer(seq: str, max_run: int) -> bool:
    return bool(re.search(rf'(.)\1{{{max_run},}}', seq))


def has_sd_sequence(seq: str) -> bool:
    return 'AGGAGG' in seq


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count('G') + seq.count('C')) / len(seq)


def gc_window_ok(seq: str, window: int, low: float, high: float) -> bool:
    """Check the latest `window`-nt block's GC; relaxed for short sequences."""
    if len(seq) < window:
        return True
    gc = gc_content(seq[-window:])
    return low <= gc <= high


class CodonOptimizer:
    def __init__(self, forbidden_sites, gc_low=0.40, gc_high=0.75,
                 max_homopolymer=6, avoid_sd=True, seed=None):
        self.forbidden_patterns = [site_to_regex(s) for s in forbidden_sites]
        self.gc_low = gc_low
        self.gc_high = gc_high
        self.max_homopolymer = max_homopolymer
        self.avoid_sd = avoid_sd
        self.rng = random.Random(seed)

    def _check_local(self, seq: str, lookback: int = 30) -> bool:
        check = seq[-(lookback + 20):] if len(seq) > lookback + 20 else seq
        for pat in self.forbidden_patterns:
            if pat.search(check):
                return False
        if has_homopolymer(check, self.max_homopolymer):
            return False
        if self.avoid_sd and has_sd_sequence(check):
            return False
        if not gc_window_ok(seq[-60:], window=50, low=self.gc_low, high=self.gc_high):
            return False
        return True

    def _weighted_choice(self, codons):
        c = [x for x, _ in codons]
        w = [x for _, x in codons]
        return self.rng.choices(c, weights=w, k=1)[0]

    def optimize(self, protein: str, max_attempts: int = 200) -> str:
        protein = protein.upper().strip()
        for ch in protein:
            if ch not in CODON_USAGE:
                raise ValueError(f"Unknown amino acid: '{ch}'")
        for _ in range(max_attempts):
            seq = self._try_optimize(protein)
            if seq is not None:
                # final full-sequence sanity check
                if all(not p.search(seq) for p in self.forbidden_patterns):
                    return seq
        raise RuntimeError(
            f"Failed to optimize after {max_attempts} attempts. "
            "Try one of: (1) raise --min-site-length to 7 or 8 to drop more "
            "short/degenerate cutters; (2) switch --enzymes to 'common'; "
            "(3) widen --gc-low / --gc-high; (4) raise --max-homopolymer to 7 "
            "or pass --no-sd-check."
        )

    def _try_optimize(self, protein: str):
        seq = ''
        i = 0
        backtracks = 0
        max_backtracks = max(len(protein) * 5, 50)

        while i < len(protein):
            aa = protein[i]
            codons = list(CODON_USAGE[aa])
            tried = set()
            placed = False
            for _ in range(len(codons)):
                remaining = [(c, w) for c, w in codons if c not in tried]
                if not remaining:
                    break
                codon = self._weighted_choice(remaining)
                tried.add(codon)
                if self._check_local(seq + codon):
                    seq += codon
                    placed = True
                    break
            if placed:
                i += 1
            else:
                undo = min(self.rng.randint(2, 8), i)
                if undo == 0:
                    return None
                seq = seq[:-3 * undo]
                i -= undo
                backtracks += 1
                if backtracks > max_backtracks:
                    return None
        return seq


def translate(dna: str) -> str:
    return ''.join(
        TRANSLATION_TABLE.get(dna[i:i + 3], 'X')
        for i in range(0, len(dna) - 2, 3)
    )


def find_all_forbidden_sites(seq: str, forbidden_sites: dict):
    """Locate all (overlapping) forbidden-site matches; report 1-based positions."""
    hits = []
    for site, enzymes in forbidden_sites.items():
        inner = site_to_regex(site).pattern
        overlap_pat = re.compile(f'(?=({inner}))')
        for m in overlap_pat.finditer(seq):
            matched = m.group(1)
            hits.append((m.start() + 1, matched, site, ','.join(enzymes)))
    return sorted(hits)


def codon_usage_summary(seq: str) -> dict:
    counts = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        counts[codon] = counts.get(codon, 0) + 1
    return counts


def read_protein(path_or_text: str):
    """Parse FASTA file, FASTA-formatted string, or raw protein string.
    Raises ValueError on multi-record FASTA, empty input, or invalid AA."""
    p = Path(path_or_text)
    try:
        is_file = p.exists() and p.is_file()
    except OSError:
        is_file = False
    if is_file:
        text = p.read_text()
        source_name = p.stem
    else:
        text = path_or_text
        source_name = 'protein'

    text = text.strip()
    if not text:
        raise ValueError("Empty protein input.")

    if text.startswith('>'):
        lines = text.splitlines()
        header_count = sum(1 for l in lines if l.lstrip().startswith('>'))
        if header_count > 1:
            raise ValueError(
                f"Multi-record FASTA detected ({header_count} records). "
                "Provide a single sequence per run."
            )
        name = lines[0].lstrip('>').strip().split()[0] if lines[0].lstrip('>').strip() else source_name
        seq_lines = [l for l in lines[1:] if not l.lstrip().startswith('>')]
        seq = ''.join(seq_lines)
    else:
        name = source_name
        seq = text

    seq = re.sub(r'\s+', '', seq).upper()
    if not seq:
        raise ValueError("No sequence content found after parsing.")
    invalid = sorted(set(seq) - set(CODON_USAGE.keys()))
    if invalid:
        raise ValueError(
            f"Invalid amino acid character(s) in input: {invalid}. "
            "Standard 20 AAs (and '*' for stop) are accepted."
        )
    return name, seq


def write_output(out_dir: Path, name: str, dna: str, protein: str,
                 forbidden_sites: dict, params: dict):
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = out_dir / f'{name}_optimized.fasta'
    with fasta_path.open('w', encoding='utf-8') as f:
        f.write(f'>{name}_optimized_REH16\n')
        for i in range(0, len(dna), 60):
            f.write(dna[i:i + 60] + '\n')

    report_path = out_dir / f'{name}_report.txt'
    translated = translate(dna)
    p_clean = protein.rstrip('*')
    t_clean = translated.rstrip('*')
    match = (t_clean == p_clean)
    hits = find_all_forbidden_sites(dna, forbidden_sites)
    gc = gc_content(dna)
    counts = codon_usage_summary(dna)

    with report_path.open('w', encoding='utf-8') as f:
        f.write("Codon Optimization Report\n")
        f.write("=========================\n\n")
        f.write(f"Sequence name : {name}\n")
        f.write("Host          : Ralstonia eutropha H16 (Cupriavidus necator H16)\n")
        f.write(f"Protein length: {len(protein)} aa\n")
        f.write(f"DNA length    : {len(dna)} nt\n")
        f.write(f"GC content    : {gc * 100:.2f}%\n")
        f.write(f"Translation   : {'MATCH' if match else 'MISMATCH (BUG!)'}\n\n")
        f.write("Parameters:\n")
        for k, v in params.items():
            f.write(f"  {k}: {v}\n")
        f.write("\n")
        f.write(f"Forbidden enzyme sites checked: {len(forbidden_sites)}\n")
        f.write(f"Forbidden sites in output     : {len(hits)}\n")
        if hits:
            f.write("  Found:\n")
            for pos, found, site, enzyme in hits:
                f.write(f"    pos={pos:6d}  matched={found}  site={site}  enzymes={enzyme}\n")
        f.write("\n")
        f.write("Codon usage in output (count):\n")
        by_aa = {}
        for codon, n in counts.items():
            aa = TRANSLATION_TABLE.get(codon, '?')
            by_aa.setdefault(aa, []).append((codon, n))
        for aa in sorted(by_aa):
            total = sum(n for _, n in by_aa[aa])
            row = '  '.join(f"{c}={n}({n / total * 100:.0f}%)" for c, n in sorted(by_aa[aa]))
            f.write(f"  {aa} (total={total}):  {row}\n")

    return fasta_path, report_path


def _ask(prompt, default=None):
    s = input(f"   {prompt}{(' [' + str(default) + ']') if default is not None else ''}: ").strip()
    if not s and default is not None:
        return str(default)
    return s


def interactive_mode():
    print("\n=== Codon Optimizer for Ralstonia eutropha H16 ===\n")

    print("1. Protein input")
    print("   Enter FASTA path, or paste sequence. To paste multi-line, press Enter")
    print("   on an empty prompt and then paste; finish with empty line.")
    inp = input("   > ").strip()
    if not inp:
        print("   (paste protein sequence; empty line to finish)")
        lines = []
        while True:
            try:
                l = input()
            except EOFError:
                break
            if not l:
                break
            lines.append(l)
        raw = '\n'.join(lines)
    else:
        raw = inp
    try:
        name, protein = read_protein(raw)
    except ValueError as e:
        print(f"   Error: {e}")
        return
    print(f"   -> name='{name}', length={len(protein)} aa\n")

    print("2. Forbidden restriction sites")
    print("   (a) Common cloning sites only (13 enzymes)")
    print("   (b) Full list from your enzyme string")
    print("   (c) Custom comma-separated names")
    print("   (n) None")
    choice = _ask("Choose", "a").lower()
    if choice == 'a':
        enzyme_names = list(COMMON_CLONING_SITES)
    elif choice == 'b':
        enzyme_names = split_concatenated_enzyme_names(FULL_ENZYME_STRING)
    elif choice == 'c':
        custom = _ask("Enzyme names (comma-separated)", "")
        enzyme_names = [x.strip() for x in custom.split(',') if x.strip()]
    else:
        enzyme_names = []
    min_len = int(_ask("Min site length to forbid (4-/5-cutters often unavoidable)", "6"))
    forbidden_sites, missing, too_short = get_forbidden_sites(enzyme_names, min_len)
    print(f"   -> {len(enzyme_names)} enzymes -> {len(forbidden_sites)} unique sites (>= {min_len} nt)")
    if too_short:
        print(f"   Skipped {len(too_short)} short-cut sites (< {min_len} nt): "
              f"{', '.join(n for n, _ in too_short[:8])}{'...' if len(too_short) > 8 else ''}")
    if missing:
        print(f"   Not in Biopython DB: {missing[:5]}{'...' if len(missing) > 5 else ''}")
    print()

    print("3. Codons")
    add_atg = _ask("Add ATG start codon if missing? [Y/n]", "Y").lower() != 'n'
    add_stop = _ask("Add stop codon (TAA) if missing? [Y/n]", "Y").lower() != 'n'
    print()

    print("4. GC content (genome ~66%)")
    gc_low = float(_ask("Lower bound %", "40")) / 100
    gc_high = float(_ask("Upper bound %", "75")) / 100
    print()

    print("5. Other constraints")
    max_run = int(_ask("Max homopolymer length", "6"))
    avoid_sd = _ask("Avoid SD sequence (AGGAGG)? [Y/n]", "Y").lower() != 'n'
    print()

    print("6. Output")
    default_out = f"./output/{name}"
    out_dir = Path(_ask("Output folder", default_out))
    print()

    if add_atg and not protein.startswith('M'):
        protein = 'M' + protein
    if add_stop and not protein.endswith('*'):
        protein = protein + '*'

    print("Optimizing...")
    opt = CodonOptimizer(
        forbidden_sites=list(forbidden_sites.keys()),
        gc_low=gc_low,
        gc_high=gc_high,
        max_homopolymer=max_run,
        avoid_sd=avoid_sd,
    )
    try:
        dna = opt.optimize(protein)
    except RuntimeError as e:
        print(f"Error: {e}")
        return

    params = {
        'gc_low': f'{gc_low * 100:.0f}%',
        'gc_high': f'{gc_high * 100:.0f}%',
        'max_homopolymer': max_run,
        'avoid_sd': avoid_sd,
        'enzymes_used': len(enzyme_names),
        'unique_sites': len(forbidden_sites),
    }
    fasta, report = write_output(out_dir, name, dna, protein, forbidden_sites, params)
    print("\nDone.")
    print(f"  FASTA : {fasta}")
    print(f"  Report: {report}")


def cli_main(args):
    try:
        name, protein = read_protein(args.input)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    if args.enzymes == 'common':
        enzyme_names = list(COMMON_CLONING_SITES)
    elif args.enzymes == 'full':
        enzyme_names = split_concatenated_enzyme_names(FULL_ENZYME_STRING)
    elif args.enzymes == 'custom':
        enzyme_names = [x.strip() for x in args.custom_enzymes.split(',') if x.strip()]
        if not enzyme_names:
            print("Error: --enzymes custom requires --custom-enzymes 'name1,name2,...'")
            sys.exit(1)
    else:
        enzyme_names = []
    forbidden_sites, missing, too_short = get_forbidden_sites(
        enzyme_names, min_site_length=args.min_site_length
    )
    if too_short:
        print(f"Skipped {len(too_short)} sites < {args.min_site_length} nt "
              f"(use --min-site-length 1 to include 4-cutters etc.)")
    if missing:
        print(f"Not in Biopython DB: {missing[:5]}{'...' if len(missing) > 5 else ''}")

    if not args.no_atg and not protein.startswith('M'):
        protein = 'M' + protein
    if not args.no_stop and not protein.endswith('*'):
        protein = protein + '*'

    opt = CodonOptimizer(
        forbidden_sites=list(forbidden_sites.keys()),
        gc_low=args.gc_low / 100,
        gc_high=args.gc_high / 100,
        max_homopolymer=args.max_homopolymer,
        avoid_sd=not args.no_sd_check,
        seed=args.seed,
    )
    dna = opt.optimize(protein)

    out_dir = Path(args.output) / name
    params = {
        'gc_low': f'{args.gc_low}%',
        'gc_high': f'{args.gc_high}%',
        'max_homopolymer': args.max_homopolymer,
        'avoid_sd': not args.no_sd_check,
        'enzymes_used': len(enzyme_names),
        'unique_sites': len(forbidden_sites),
        'seed': args.seed,
    }
    fasta, report = write_output(out_dir, name, dna, protein, forbidden_sites, params)
    print(f"FASTA : {fasta}")
    print(f"Report: {report}")


def main():
    parser = argparse.ArgumentParser(
        description='Codon optimizer for Ralstonia eutropha H16'
    )
    parser.add_argument('-i', '--input', help='FASTA file path or protein sequence string')
    parser.add_argument('-o', '--output', default='./output',
                        help='Output base folder (default: ./output, relative to cwd)')
    parser.add_argument('--enzymes', choices=['common', 'full', 'none', 'custom'], default='common',
                        help='Enzyme set to forbid (default: common). Use "custom" with --custom-enzymes.')
    parser.add_argument('--custom-enzymes', default='',
                        help='Comma-separated enzyme names (used with --enzymes custom). e.g. BsaI,BbsI')
    parser.add_argument('--gc-low', type=float, default=40.0, help='GC%% lower bound (default 40)')
    parser.add_argument('--gc-high', type=float, default=75.0, help='GC%% upper bound (default 75)')
    parser.add_argument('--max-homopolymer', type=int, default=6,
                        help='Max consecutive same nt allowed (default 6)')
    parser.add_argument('--no-sd-check', action='store_true', help='Do not avoid AGGAGG')
    parser.add_argument('--no-atg', action='store_true', help='Do not auto-add ATG start')
    parser.add_argument('--no-stop', action='store_true', help='Do not auto-add stop codon')
    parser.add_argument('--seed', type=int, default=None, help='RNG seed for reproducibility')
    parser.add_argument('--min-site-length', type=int, default=6,
                        help='Skip enzyme sites shorter than this many nt (default 6; '
                             '4/5-cutters and degenerate 5-cutters are usually '
                             'unavoidable). Set to 1 for strict mode.')
    args = parser.parse_args()

    if args.input is None:
        interactive_mode()
    else:
        cli_main(args)


if __name__ == '__main__':
    main()
