#!/usr/bin/env bash
#
# clean_colocation_tables.sh
#
# Drops rows carrying uninformative annotations from ARG-MGE.tsv, ARG-VF.tsv and
# ARG-MGE-VF.tsv. summary.tsv is left alone.
#
# Defaults (matched case-insensitively, surrounding whitespace ignored):
#   ARG_class    : unclassified
#   MGE_gene     : NA, NA:Keyword
#   VF_category  : Others
#
# Runs in place and keeps the originals in <indir>/pre_clean_backup/.
# Use --dry-run first to see what would go.
#
#   ./clean_colocation_tables.sh -i colocation/ --dry-run
#   ./clean_colocation_tables.sh -i colocation/
#   ./clean_colocation_tables.sh -i colocation/ -o colocation_clean/
#
# Add to a rule with a comma separated list, e.g.
#   --drop-mge-category unclassified --drop-vf-gene NA,unclassified
#
set -euo pipefail

IN_DIR=""; OUT_DIR=""; DRY=0
DROP_ARG_CLASS="unclassified"
DROP_ARG_GENE=""
DROP_MGE_GENE="NA,NA:Keyword"
DROP_MGE_CATEGORY=""
DROP_VF_GENE=""
DROP_VF_CATEGORY="Others"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log() { printf '%s\n' "$*" >&2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--indir)  IN_DIR="$2"; shift 2 ;;
    -o|--outdir) OUT_DIR="$2"; shift 2 ;;
    --dry-run)   DRY=1; shift ;;
    --drop-arg-class)    DROP_ARG_CLASS="$2"; shift 2 ;;
    --drop-arg-gene)     DROP_ARG_GENE="$2"; shift 2 ;;
    --drop-mge-gene)     DROP_MGE_GENE="$2"; shift 2 ;;
    --drop-mge-category) DROP_MGE_CATEGORY="$2"; shift 2 ;;
    --drop-vf-gene)      DROP_VF_GENE="$2"; shift 2 ;;
    --drop-vf-category)  DROP_VF_CATEGORY="$2"; shift 2 ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed 's/^#//; s/^ //' | head -n -1; exit 0 ;;
    *) die "unknown option: $1 (try --help)" ;;
  esac
done

[[ -n "$IN_DIR" ]] || die "-i/--indir is required"
[[ -d "$IN_DIR" ]] || die "not a directory: $IN_DIR"

if [[ -n "$OUT_DIR" ]]; then
  mkdir -p "$OUT_DIR"; BACKUP=""
else
  BACKUP="$IN_DIR/pre_clean_backup"
  # never clobber an existing backup: a second run would otherwise replace the
  # pristine originals with already-cleaned copies
  if [[ -d "$BACKUP" ]]; then
    BACKUP="$IN_DIR/pre_clean_backup_$(date +%Y%m%d-%H%M%S)"
    log "note: $IN_DIR/pre_clean_backup already exists, using $BACKUP instead"
  fi
  [[ $DRY -eq 1 ]] || mkdir -p "$BACKUP"
fi

log "rules (case-insensitive):"
[[ -n "$DROP_ARG_CLASS"    ]] && log "  ARG_class    drop: $DROP_ARG_CLASS"
[[ -n "$DROP_ARG_GENE"     ]] && log "  ARG_gene     drop: $DROP_ARG_GENE"
[[ -n "$DROP_MGE_GENE"     ]] && log "  MGE_gene     drop: $DROP_MGE_GENE"
[[ -n "$DROP_MGE_CATEGORY" ]] && log "  MGE_category drop: $DROP_MGE_CATEGORY"
[[ -n "$DROP_VF_GENE"      ]] && log "  VF_gene      drop: $DROP_VF_GENE"
[[ -n "$DROP_VF_CATEGORY"  ]] && log "  VF_category  drop: $DROP_VF_CATEGORY"
log ""
[[ $DRY -eq 1 ]] && log "DRY RUN - nothing will be written" && log ""

clean_one() {
  local file="$1" out="$2"
  awk -F'\t' -v OFS='\t' \
      -v r_ac="$DROP_ARG_CLASS" -v r_ag="$DROP_ARG_GENE" \
      -v r_mg="$DROP_MGE_GENE"  -v r_mc="$DROP_MGE_CATEGORY" \
      -v r_vg="$DROP_VF_GENE"   -v r_vc="$DROP_VF_CATEGORY" \
      -v out="$out" -v dry="$DRY" '
    function trim(x) { gsub(/^[[:space:]]+|[[:space:]]+$/, "", x); return x }
    function load(list, name,    n, a, i) {
      n = split(list, a, ",")
      for (i = 1; i <= n; i++) {
        v = tolower(trim(a[i]))
        if (v != "") rules[name, v] = 1
      }
    }
    BEGIN {
      load(r_ac, "ARG_class"); load(r_ag, "ARG_gene")
      load(r_mg, "MGE_gene");  load(r_mc, "MGE_category")
      load(r_vg, "VF_gene");   load(r_vc, "VF_category")
    }
    NR == 1 {
      for (i = 1; i <= NF; i++) { idx[$i] = i; hdr[i] = $i }
      if (!dry) print > out
      next
    }
    {
      total++
      drop = 0; why = ""
      for (h in idx) {
        c = idx[h]
        v = tolower(trim($c))
        if ((h, v) in rules) { drop = 1; why = h; dropped_value[h SUBSEP trim($c)]++; break }
      }
      if (drop) { removed++; by_col[why]++; next }
      kept++
      if (!dry) print > out
    }
    END {
      printf "%d\t%d\t%d\n", total+0, kept+0, removed+0 > "/dev/stdout"
      for (k in by_col) printf "      %-13s removed %d row(s)\n", k, by_col[k] > "/dev/stderr"
      for (k in dropped_value) {
        split(k, kk, SUBSEP)
        printf "        %s = \"%s\"  (%d)\n", kk[1], kk[2], dropped_value[k] > "/dev/stderr"
      }
    }
  ' "$file"
}

found=0
for base in ARG-MGE ARG-VF ARG-MGE-VF; do
  file="$IN_DIR/$base.tsv"
  [[ -f "$file" ]] || { log "  $base.tsv: not found, skipped"; continue; }
  found=$((found + 1))

  if [[ -n "$OUT_DIR" ]]; then target="$OUT_DIR/$base.tsv"; else target="$IN_DIR/.$base.cleaned.tmp"; fi

  log "  $base.tsv"
  read -r total kept removed < <(clean_one "$file" "$target")

  if [[ $DRY -eq 0 ]]; then
    if [[ -z "$OUT_DIR" ]]; then
      cp -p "$file" "$BACKUP/$base.tsv"
      mv -f "$target" "$file"
    fi
  else
    rm -f "$target" 2>/dev/null || true
  fi

  pct=$(awk -v r="$removed" -v t="$total" 'BEGIN { printf "%.1f", (t ? 100*r/t : 0) }')
  log "      $total rows -> $kept kept, $removed removed (${pct}%)"
  log ""
done

[[ $found -gt 0 ]] || die "none of ARG-MGE.tsv, ARG-VF.tsv, ARG-MGE-VF.tsv found in $IN_DIR"

if [[ $DRY -eq 1 ]]; then
  log "DRY RUN complete. Re-run without --dry-run to apply."
elif [[ -n "$OUT_DIR" ]]; then
  log "cleaned tables written to $OUT_DIR (originals untouched)"
else
  log "cleaned in place; originals in $BACKUP"
fi
log ""
log "Re-run make_gephi_networks.sh on the cleaned tables to rebuild the networks."
