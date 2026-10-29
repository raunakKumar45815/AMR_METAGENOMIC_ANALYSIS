#!/usr/bin/env bash
#
# build_colocation_tables.sh
#
# Builds three ARG co-location tables by linking, on Read_ID, the DeepARG
# mapping report with DIAMOND searches of the same ARG-mapped reads against
# mobileOG-db (MGE) and VFDB (VF).
#
#   ARG-MGE.tsv      Sample Read_ID ARG_gene ARG_class MGE_gene MGE_category
#   ARG-VF.tsv       Sample Read_ID ARG_gene ARG_class VF_gene VF_category
#   ARG-MGE-VF.tsv   Sample Read_ID ARG_gene ARG_class MGE_gene MGE_category VF_gene VF_category
#   summary.tsv      per-sample counts at every stage
#
# ARG-MGE-VF is a strict three-way intersection: the same read must carry an
# ARG, an MGE and a VF hit. It is always the smallest of the three.
#
# ---------------------------------------------------------------- input formats
# ARG  (DeepARG .mapping.ARG, header line starting with #ARG)
#   1 #ARG  2 query-start  3 query-end  4 read_id  5 predicted_ARG-class
#   6 best-hit  7 probability  8 identity  9 alignment-length
#   10 alignment-bitscore  11 alignment-evalue  12 counts
#
# MGE  (DIAMOND vs mobileOG-db, no header, 10 columns)
#   1 qseqid  2 sseqid  3 pident  4 length  5 qstart  6 qend
#   7 sstart  8 send  9 evalue  10 bitscore
#   sseqid is pipe delimited:
#     mobileOG_id|Gene_Name|Best_Hit_Accession|Major_Category|Minor_Category|Source|Evidence
#   Gene_Name is frequently NA, so the gene falls back to the accession then the
#   mobileOG id (order set by --mge-gene-fields).
#
# VF  (DIAMOND vs VFDB, header line, 12 columns, already annotated)
#   1 qseqid  2 sseqid  3 Gene  4 Category  5 pident  6 length
#   7 qstart  8 qend  9 sstart  10 send  11 evalue  12 bitscore
#
# Every column index above is overridable; see --help.
#
# Usage:
#   ./build_colocation_tables.sh --argdir deeparg/ --mgedir diamond_mge/ \
#       --vfdir diamond_vf/ --outdir colocation/
#
set -euo pipefail

# ----------------------------------------------------------------- defaults
ARG_DIR=""; MGE_DIR=""; VF_DIR=""; OUT_DIR="colocation"; SAMPLE_LIST=""

ARG_SUFFIX=".mapping.ARG"
MGE_SUFFIX=".mge.diamond.tsv"
VF_SUFFIX=".vf.diamond.tsv"

# ARG columns
A_GENE_COL=1; A_READ_COL=4; A_CLASS_COL=5; A_PROB_COL=7
ARG_MIN_PROB=0          # DeepARG has already filtered; 0 = keep everything it reported

# MGE columns (10-column layout above)
M_Q_COL=1; M_S_COL=2; M_PID_COL=3; M_LEN_COL=4; M_EVAL_COL=9; M_BIT_COL=10
MGE_GENE_FIELDS="2,3,1"   # pipe fields tried in order for the gene name
MGE_CAT_FIELD=4

# VF columns (12-column annotated layout above)
V_Q_COL=1; V_S_COL=2; V_GENE_COL=3; V_CAT_COL=4
V_PID_COL=5; V_LEN_COL=6; V_EVAL_COL=11; V_BIT_COL=12
VF_GENE_SOURCE="gene"     # gene | accession | sseqid
MGE_COLS_SET=0; VF_COLS_SET=0; AUTO_COLS=1

# Hit filters. Your data sits at 72-80% amino acid identity, so the identity
# default is deliberately permissive - raise it consciously, and read the
# per-filter drop counts printed at the end before you trust any table.
MIN_ID=60
MIN_LEN=50
MAX_EVAL=1e-10

MGE_MAP=""; VF_MAP=""     # optional gene/accession <TAB> category fallback tables
UNKNOWN="unclassified"
NA_PATTERN="^(NA|na|N/A|-|\\.|)$"
INCLUDE_VF_GENE=1
KEEP_TEMP=0

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log() { printf '%s\n' "$*" >&2; }
usage() { sed -n '2,/^set -euo/p' "$0" | sed 's/^#//; s/^ //' | head -n -1; exit 0; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --argdir) ARG_DIR="$2"; shift 2 ;;
    --mgedir) MGE_DIR="$2"; shift 2 ;;
    --vfdir)  VF_DIR="$2"; shift 2 ;;
    --outdir) OUT_DIR="$2"; shift 2 ;;
    --samples) SAMPLE_LIST="$2"; shift 2 ;;
    --arg-suffix) ARG_SUFFIX="$2"; shift 2 ;;
    --mge-suffix) MGE_SUFFIX="$2"; shift 2 ;;
    --vf-suffix)  VF_SUFFIX="$2"; shift 2 ;;
    --min-id)   MIN_ID="$2"; shift 2 ;;
    --min-len)  MIN_LEN="$2"; shift 2 ;;
    --max-eval) MAX_EVAL="$2"; shift 2 ;;
    --arg-min-prob) ARG_MIN_PROB="$2"; shift 2 ;;
    --mge-gene-fields) MGE_GENE_FIELDS="$2"; shift 2 ;;
    --mge-cat-field)   MGE_CAT_FIELD="$2"; shift 2 ;;
    --mge-cols) IFS=, read -r M_Q_COL M_S_COL M_PID_COL M_LEN_COL M_EVAL_COL M_BIT_COL <<<"$2"; MGE_COLS_SET=1; shift 2 ;;
    --vf-cols)  IFS=, read -r V_Q_COL V_S_COL V_GENE_COL V_CAT_COL V_PID_COL V_LEN_COL V_EVAL_COL V_BIT_COL <<<"$2"; VF_COLS_SET=1; shift 2 ;;
    --no-auto-cols) AUTO_COLS=0; shift ;;
    --arg-cols) IFS=, read -r A_GENE_COL A_READ_COL A_CLASS_COL <<<"$2"; shift 2 ;;
    --vf-gene-source) VF_GENE_SOURCE="$2"; shift 2 ;;
    --mge-map) MGE_MAP="$2"; shift 2 ;;
    --vf-map)  VF_MAP="$2"; shift 2 ;;
    --no-vf-gene) INCLUDE_VF_GENE=0; shift ;;
    --keep-temp) KEEP_TEMP=1; shift ;;
    -h|--help) usage ;;
    *) die "unknown option: $1 (try --help)" ;;
  esac
done

[[ -n "$ARG_DIR" ]] || die "--argdir is required"
[[ -d "$ARG_DIR" ]] || die "not a directory: $ARG_DIR"
[[ -n "$MGE_DIR" || -n "$VF_DIR" ]] || die "at least one of --mgedir / --vfdir is required"
[[ -z "$MGE_DIR" || -d "$MGE_DIR" ]] || die "not a directory: $MGE_DIR"
[[ -z "$VF_DIR"  || -d "$VF_DIR"  ]] || die "not a directory: $VF_DIR"
[[ -z "$MGE_MAP" || -f "$MGE_MAP" ]] || die "no such file: $MGE_MAP"
[[ -z "$VF_MAP"  || -f "$VF_MAP"  ]] || die "no such file: $VF_MAP"
case "$VF_GENE_SOURCE" in gene|accession|sseqid) ;; *) die "--vf-gene-source must be gene, accession or sseqid" ;; esac

mkdir -p "$OUT_DIR"
TMP="$(mktemp -d)"
cleanup() { if [[ $KEEP_TEMP -eq 1 ]]; then log "temp kept: $TMP"; else rm -rf "$TMP"; fi; }
trap cleanup EXIT

ARG_MGE="$OUT_DIR/ARG-MGE.tsv"
ARG_VF="$OUT_DIR/ARG-VF.tsv"
TRIPLE="$OUT_DIR/ARG-MGE-VF.tsv"
SUMMARY="$OUT_DIR/summary.tsv"
STATS="$TMP/filter_stats.tsv"; : > "$STATS"

# ----------------------------------------------------------------- samples
if [[ -n "$SAMPLE_LIST" ]]; then
  [[ -f "$SAMPLE_LIST" ]] || die "no such file: $SAMPLE_LIST"
  mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLE_LIST" | tr -d '\r')
else
  mapfile -t SAMPLES < <(
    find "$ARG_DIR" -maxdepth 1 -type f -name "*${ARG_SUFFIX}" -printf '%f\n' \
      | sed "s/${ARG_SUFFIX//./\\.}\$//" | sort
  )
fi
[[ ${#SAMPLES[@]} -gt 0 ]] || die "no files matching *${ARG_SUFFIX} in $ARG_DIR"

# ----------------------------------------------------------------- preflight
# Count how many samples actually have an MGE / VF file. Getting the suffix
# wrong silently produces three empty tables, so refuse to run instead.
n_mge=0; n_vf=0; first_mge=""; first_vf=""
for sample in "${SAMPLES[@]}"; do
  if [[ -n "$MGE_DIR" && -f "$MGE_DIR/${sample}${MGE_SUFFIX}" ]]; then
    n_mge=$((n_mge + 1)); [[ -z "$first_mge" ]] && first_mge="$MGE_DIR/${sample}${MGE_SUFFIX}"
  fi
  if [[ -n "$VF_DIR" && -f "$VF_DIR/${sample}${VF_SUFFIX}" ]]; then
    n_vf=$((n_vf + 1)); [[ -z "$first_vf" ]] && first_vf="$VF_DIR/${sample}${VF_SUFFIX}"
  fi
done

suffix_help() {  # $1 = dir, $2 = suffix, $3 = flag name
  log ""
  log "None of the ${#SAMPLES[@]} samples has a file matching <sample>${2} in ${1}"
  log "Files actually in that directory:"
  find "$1" -maxdepth 1 -type f -printf '  %f\n' 2>/dev/null | head -3 >&2
  log ""
  log "Set the suffix to match. If your files are just <sample>.tsv, use:  ${3} .tsv"
}
if [[ -n "$MGE_DIR" && $n_mge -eq 0 ]]; then suffix_help "$MGE_DIR" "$MGE_SUFFIX" "--mge-suffix"; die "no MGE files matched"; fi
if [[ -n "$VF_DIR"  && $n_vf  -eq 0 ]]; then suffix_help "$VF_DIR"  "$VF_SUFFIX"  "--vf-suffix";  die "no VF files matched"; fi

# ----------------------------------------------------------------- layout detection
# Works out column positions from the file itself, so a headerless outfmt 6 file
# and an annotated file with extra Gene/Category columns both parse correctly.
# Rule: e-value and bitscore are the last two columns; pident is the first
# numeric field at or after 3 that is followed by another number (the alignment
# length); anything between sseqid and pident is annotation.
detect_layout() {
  awk -F'\t' '
    /^#/ { next }
    tolower($1) ~ /^qseqid$/ || tolower($1) == "query" { next }
    NF >= 4 {
      pid = 0
      for (i = 3; i < NF - 1; i++)
        if ($i ~ /^[0-9]+(\.[0-9]+)?$/ && $(i+1) ~ /^[0-9]+$/) { pid = i; break }
      if (pid == 0) pid = 3
      printf "%d %d %d %d %d %d\n", pid, pid+1, NF-1, NF, (pid > 3 ? 3 : 0), (pid > 4 ? 4 : 0)
      exit
    }
  ' "$1"
}

if [[ $AUTO_COLS -eq 1 && $MGE_COLS_SET -eq 0 && -n "$first_mge" ]]; then
  read -r p l e b _g _c < <(detect_layout "$first_mge") || true
  if [[ -n "${b:-}" && "$b" -gt 3 ]]; then
    M_PID_COL=$p; M_LEN_COL=$l; M_EVAL_COL=$e; M_BIT_COL=$b
    log "MGE layout         : ${b} cols -> pident $p, length $l, evalue $e, bitscore $b"
  fi
fi
if [[ $AUTO_COLS -eq 1 && $VF_COLS_SET -eq 0 && -n "$first_vf" ]]; then
  read -r p l e b g c < <(detect_layout "$first_vf") || true
  if [[ -n "${b:-}" && "$b" -gt 3 ]]; then
    V_PID_COL=$p; V_LEN_COL=$l; V_EVAL_COL=$e; V_BIT_COL=$b; V_GENE_COL=$g; V_CAT_COL=$c
    if [[ "$g" -gt 0 ]]; then
      log "VF layout          : ${b} cols -> Gene $g, Category $c, pident $p, length $l, evalue $e, bitscore $b"
    else
      log "VF layout          : ${b} cols, plain outfmt 6 - no Gene/Category columns."
      log "                     VF_gene from sseqid; VF_category needs --vf-map or it will be '$UNKNOWN'."
    fi
  fi
fi

log "samples            : ${#SAMPLES[@]}"
log "with MGE / VF file : ${n_mge} / ${n_vf}"
log "hit filters        : identity >= ${MIN_ID}%, length >= ${MIN_LEN}, e-value <= ${MAX_EVAL}"
log "best hit per read  : highest bitscore; ties -> lowest e-value -> named gene -> lowest sseqid"
log "VF gene source     : ${VF_GENE_SOURCE}"
log "output             : ${OUT_DIR}"
log ""

# ----------------------------------------------------------------- headers
printf 'Sample\tRead_ID\tARG_gene\tARG_class\tMGE_gene\tMGE_category\n' > "$ARG_MGE"
if [[ $INCLUDE_VF_GENE -eq 1 ]]; then
  printf 'Sample\tRead_ID\tARG_gene\tARG_class\tVF_gene\tVF_category\n' > "$ARG_VF"
else
  printf 'Sample\tRead_ID\tARG_gene\tARG_class\tVF_category\n' > "$ARG_VF"
fi
printf 'Sample\tRead_ID\tARG_gene\tARG_class\tMGE_gene\tMGE_category\tVF_gene\tVF_category\n' > "$TRIPLE"
printf 'Sample\tARG_reads\tMGE_reads\tVF_reads\tARG_MGE\tARG_VF\tARG_MGE_VF\tMGE_unclassified\tVF_unclassified\tMGE_tied_categories\n' > "$SUMMARY"

# ----------------------------------------------------------------- parsers
# DeepARG -> read_id \t ARG_gene \t ARG_class
parse_arg() {
  awk -F'\t' -v g="$A_GENE_COL" -v r="$A_READ_COL" -v c="$A_CLASS_COL" \
             -v pc="$A_PROB_COL" -v minp="$ARG_MIN_PROB" '
    /^#/ { next }
    tolower($r) ~ /^read[_-]?id$/ { next }
    NF >= 5 && $r != "" {
      if (minp > 0 && pc > 0 && pc <= NF && ($pc+0) < minp) next
      print $r "\t" $g "\t" $c
    }
  ' "$1"
}

# DIAMOND -> read_id \t gene \t category, one best hit per read.
# Shared by both databases; which columns and which annotation source to use is
# passed in, so the two very different file layouts go through the same logic.
best_hits() {
  local file="$1" db="$2"
  local q s pid len ev bit gcol ccol genefields catfield mapfile
  if [[ "$db" == "mge" ]]; then
    q=$M_Q_COL; s=$M_S_COL; pid=$M_PID_COL; len=$M_LEN_COL; ev=$M_EVAL_COL; bit=$M_BIT_COL
    gcol=0; ccol=0; genefields="$MGE_GENE_FIELDS"; catfield="$MGE_CAT_FIELD"; mapfile="$MGE_MAP"
  else
    q=$V_Q_COL; s=$V_S_COL; pid=$V_PID_COL; len=$V_LEN_COL; ev=$V_EVAL_COL; bit=$V_BIT_COL
    gcol=$V_GENE_COL; ccol=$V_CAT_COL; genefields=""; catfield=0; mapfile="$VF_MAP"
  fi

  awk -F'\t' \
    -v q="$q" -v s="$s" -v pid="$pid" -v len="$len" -v ev="$ev" -v bit="$bit" \
    -v gcol="$gcol" -v ccol="$ccol" -v genefields="$genefields" -v catfield="$catfield" \
    -v minid="$MIN_ID" -v minlen="$MIN_LEN" -v maxev="$MAX_EVAL" \
    -v mapfile="$mapfile" -v unknown="$UNKNOWN" -v napat="$NA_PATTERN" \
    -v vfsrc="$VF_GENE_SOURCE" -v db="$db" -v statsfile="$STATS" -v sample="$3" '
    function isna(x) { return (x ~ napat) }
    function accession(sid) {                       # VFG048797(gb|WP_004175261) -> gb|WP_004175261
      if (match(sid, /\(([^)]*)\)/)) return substr(sid, RSTART+1, RLENGTH-2)
      return sid
    }
    BEGIN {
      ngf = split(genefields, gf, ",")
      if (mapfile != "") {
        while ((getline line < mapfile) > 0) {
          if (line ~ /^[[:space:]]*$/ || line ~ /^#/) continue
          if (split(line, f, "\t") < 2) continue
          k = f[1]; v = f[2]
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", k); gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
          if (tolower(k) == "gene" || tolower(k) == "sseqid") continue
          if (!(k in cat)) cat[k] = v
        }
        close(mapfile)
      }
    }
    /^#/ { next }
    tolower($q) ~ /^qseqid$/ || tolower($q) == "query" { next }     # header line
    {
      total++
      if (($pid+0) < minid)  { drop_id++;  next }
      if (($len+0) < minlen) { drop_len++; next }
      if (($ev+0)  > maxev)  { drop_ev++;  next }
      kept++

      sid = $s
      # ---- gene
      gene = ""
      if (db == "vf") {
        if (vfsrc == "sseqid")          gene = sid
        else if (vfsrc == "accession")  gene = accession(sid)
        else {                                       # "gene", accession as fallback
          if (gcol > 0 && gcol <= NF && !isna($gcol)) gene = $gcol
          else gene = accession(sid)
        }
      } else {
        n = split(sid, p, "|")
        for (i = 1; i <= ngf; i++) {
          fi = gf[i] + 0
          if (fi > 0 && fi <= n && !isna(p[fi])) { gene = p[fi]; break }
        }
        if (gene == "") gene = sid
      }

      # ---- category
      category = ""
      if (ccol > 0 && ccol <= NF && !isna($ccol)) category = $ccol
      else if (catfield > 0) {
        n = split(sid, p, "|")
        if (catfield <= n && !isna(p[catfield])) category = p[catfield]
      }
      if (category == "") {
        if (sid in cat) category = cat[sid]
        else if (gene in cat) category = cat[gene]
        else {
          n = split(sid, p, "|")
          for (i = 1; i <= n; i++) if (p[i] in cat) { category = cat[p[i]]; break }
        }
      }
      if (category == "") category = unknown

      r = $q; b = $bit+0; e = $ev+0
      named = (gene == sid || isna(gene)) ? 0 : 1
      if (!(r in bestbit)) { take = 1 }
      else if (b > bestbit[r]) { take = 1 }
      else if (b < bestbit[r]) { take = 0 }
      else if (e < bestev[r]) { take = 1 }
      else if (e > bestev[r]) { take = 0 }
      else {
        if (category != bestcat[r]) tied[r] = 1        # tied hits disagree on category
        if (named > bestnamed[r]) take = 1
        else if (named < bestnamed[r]) take = 0
        else take = (sid < bestsid[r]) ? 1 : 0         # deterministic last resort
      }
      if (take) {
        bestbit[r] = b; bestev[r] = e; bestsid[r] = sid
        bestgene[r] = gene; bestcat[r] = category; bestnamed[r] = named
      }
    }
    END {
      nreads = 0; nun = 0; namb = 0
      for (r in bestgene) {
        print r "\t" bestgene[r] "\t" bestcat[r]
        nreads++
        if (bestcat[r] == unknown) nun++
        if (r in tied) namb++
      }
      printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", sample, db, total+0, kept+0, \
             drop_id+0, drop_len+0, drop_ev+0, nun, namb >> statsfile
    }
  ' "$file"
}

count_lines() { if [[ -s "$1" ]]; then wc -l < "$1" | tr -d ' '; else echo 0; fi; }

# ----------------------------------------------------------------- main loop
n_done=0
for sample in "${SAMPLES[@]}"; do
  arg_file="$ARG_DIR/${sample}${ARG_SUFFIX}"
  mge_file="$MGE_DIR/${sample}${MGE_SUFFIX}"
  vf_file="$VF_DIR/${sample}${VF_SUFFIX}"

  if [[ ! -f "$arg_file" ]]; then log "  [skip] $sample: no ARG file"; continue; fi
  parse_arg "$arg_file" | sort -u -k1,1 > "$TMP/arg.tsv"

  : > "$TMP/mge.tsv"; : > "$TMP/vf.tsv"
  if [[ -n "$MGE_DIR" ]]; then
    if [[ -f "$mge_file" ]]; then best_hits "$mge_file" mge "$sample" | sort -k1,1 > "$TMP/mge.tsv"
    else log "  [warn] $sample: no MGE file"; fi
  fi
  if [[ -n "$VF_DIR" ]]; then
    if [[ -f "$vf_file" ]]; then best_hits "$vf_file" vf "$sample" | sort -k1,1 > "$TMP/vf.tsv"
    else log "  [warn] $sample: no VF file"; fi
  fi

  awk -F'\t' -v S="$sample" '
    NR==FNR { g[$1]=$2; c[$1]=$3; next }
    ($1 in g) { print S "\t" $1 "\t" g[$1] "\t" c[$1] "\t" $2 "\t" $3 }
  ' "$TMP/arg.tsv" "$TMP/mge.tsv" > "$TMP/arg_mge.tsv"

  awk -F'\t' -v S="$sample" -v wg="$INCLUDE_VF_GENE" '
    NR==FNR { g[$1]=$2; c[$1]=$3; next }
    ($1 in g) {
      if (wg) print S "\t" $1 "\t" g[$1] "\t" c[$1] "\t" $2 "\t" $3
      else    print S "\t" $1 "\t" g[$1] "\t" c[$1] "\t" $3
    }
  ' "$TMP/arg.tsv" "$TMP/vf.tsv" > "$TMP/arg_vf.tsv"

  awk -F'\t' -v S="$sample" -v mgef="$TMP/mge.tsv" '
    FNR==NR { ag[$1]=$2; ac[$1]=$3; next }
    FILENAME==mgef { if ($1 in ag) { mg[$1]=$2; mc[$1]=$3 } next }
    ($1 in ag) && ($1 in mg) {
      print S "\t" $1 "\t" ag[$1] "\t" ac[$1] "\t" mg[$1] "\t" mc[$1] "\t" $2 "\t" $3
    }
  ' "$TMP/arg.tsv" "$TMP/mge.tsv" "$TMP/vf.tsv" > "$TMP/triple.tsv"

  sort -k2,2 "$TMP/arg_mge.tsv" >> "$ARG_MGE"
  sort -k2,2 "$TMP/arg_vf.tsv"  >> "$ARG_VF"
  sort -k2,2 "$TMP/triple.tsv"  >> "$TRIPLE"

  mge_un=$(awk -F'\t' -v s="$sample" '$1==s && $2=="mge" {print $8+0}' "$STATS"); mge_un=${mge_un:-0}
  vf_un=$(awk  -F'\t' -v s="$sample" '$1==s && $2=="vf"  {print $8+0}' "$STATS"); vf_un=${vf_un:-0}
  mge_amb=$(awk -F'\t' -v s="$sample" '$1==s && $2=="mge" {print $9+0}' "$STATS"); mge_amb=${mge_amb:-0}

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" \
    "$(count_lines "$TMP/arg.tsv")" "$(count_lines "$TMP/mge.tsv")" "$(count_lines "$TMP/vf.tsv")" \
    "$(count_lines "$TMP/arg_mge.tsv")" "$(count_lines "$TMP/arg_vf.tsv")" \
    "$(count_lines "$TMP/triple.tsv")" "$mge_un" "$vf_un" "$mge_amb" >> "$SUMMARY"

  n_done=$((n_done + 1))
  printf '\r  processed %d/%d' "$n_done" "${#SAMPLES[@]}" >&2
done
log ""; log ""

# ----------------------------------------------------------------- report
tot() { awk 'NR>1' "$1" | wc -l | tr -d ' '; }
log "ARG-MGE.tsv    : $(tot "$ARG_MGE") rows"
log "ARG-VF.tsv     : $(tot "$ARG_VF") rows"
log "ARG-MGE-VF.tsv : $(tot "$TRIPLE") rows"
log "summary.tsv    : per-sample counts"
log ""
log "hits removed by each filter (alignments, not reads):"
awk -F'\t' '
  { tot[$2]+=$3; kept[$2]+=$4; di[$2]+=$5; dl[$2]+=$6; de[$2]+=$7; un[$2]+=$8; amb[$2]+=$9 }
  END {
    for (d in tot)
      printf "  %-4s %d alignments -> %d kept | dropped: identity %d, length %d, e-value %d\n",
             toupper(d), tot[d], kept[d], di[d], dl[d], de[d]
  }
' "$STATS" >&2

un_mge=$(awk -F'\t' '$2=="mge" {s+=$8} END {print s+0}' "$STATS")
un_vf=$(awk  -F'\t' '$2=="vf"  {s+=$8} END {print s+0}' "$STATS")
amb=$(awk    -F'\t' '$2=="mge" {s+=$9} END {print s+0}' "$STATS")
[[ "$un_mge" -gt 0 ]] && log "" && log "NOTE: $un_mge MGE reads have category '$UNKNOWN' - check --mge-cat-field / --mge-map."
[[ "$un_vf"  -gt 0 ]] && log "NOTE: $un_vf VF reads have category '$UNKNOWN' - check --vf-cols / --vf-map."
[[ "$amb"    -gt 0 ]] && log "NOTE: $amb MGE reads had equally-scoring hits disagreeing on category;" \
  && log "      one was chosen deterministically. Column MGE_tied_categories in summary.tsv."
exit 0
