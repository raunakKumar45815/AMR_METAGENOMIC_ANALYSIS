#!/usr/bin/env bash
#
# make_gephi_networks.sh
#
# Turns the colocation tables into Gephi edge lists. Every row in a colocation
# table is one read, so the weight of an edge is the number of reads on which
# the two annotations were found together.
#
#   Source	Target	Type	Weight
#   MLS	integration/excision	Undirected	108
#
# Inputs: the output directory of build_colocation_tables.sh. summary.tsv is
# ignored. Columns are located by header name, so a table written with
# --no-vf-gene works too.
#
# Outputs, per input table and per level:
#   edges_<table>_<pair>_<level>.tsv   Source Target Type Weight
#   nodes_<table>_<pair>_<level>.tsv   Id Label Class
#
# ARG-MGE-VF produces three pairings (ARG-MGE, ARG-VF, MGE-VF) from the reads
# that carry all three annotations, which is a different and much smaller set
# than the pairwise tables. Do not mix the two in one network.
#
# Levels:
#   class  ARG_class  x  MGE_category / VF_category   (your example)
#   gene   ARG_gene   x  MGE_gene     / VF_gene
#
# Usage:
#   ./make_gephi_networks.sh -i colocation/ -o gephi/
#   ./make_gephi_networks.sh -i colocation/ -o gephi/ --level class --min-weight 5
#   ./make_gephi_networks.sh -i colocation/ -o gephi/ --by-sample
#
set -euo pipefail

IN_DIR=""; OUT_DIR="gephi"
LEVEL="both"          # class | gene | both
MIN_WEIGHT=1
SEP=$'\t'; EXT="tsv"
WITH_NSAMPLES=0
PREFIX_LABELS=0
BY_SAMPLE=0

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log() { printf '%s\n' "$*" >&2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--indir)  IN_DIR="$2"; shift 2 ;;
    -o|--outdir) OUT_DIR="$2"; shift 2 ;;
    --level)     LEVEL="$2"; shift 2 ;;
    --min-weight) MIN_WEIGHT="$2"; shift 2 ;;
    --csv)       SEP=","; EXT="csv"; shift ;;
    --with-nsamples) WITH_NSAMPLES=1; shift ;;
    --prefix-labels) PREFIX_LABELS=1; shift ;;
    --by-sample) BY_SAMPLE=1; shift ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed 's/^#//; s/^ //' | head -n -1; exit 0 ;;
    *) die "unknown option: $1 (try --help)" ;;
  esac
done

[[ -n "$IN_DIR" ]] || die "-i/--indir is required (the colocation output directory)"
[[ -d "$IN_DIR" ]] || die "not a directory: $IN_DIR"
case "$LEVEL" in class|gene|both) ;; *) die "--level must be class, gene or both" ;; esac
mkdir -p "$OUT_DIR"

# ----------------------------------------------------------------- helpers
# index of a named column in a header line, or 0
col_of() {
  awk -F'\t' -v want="$2" 'NR==1 { for (i=1;i<=NF;i++) if ($i==want) { print i; exit } print 0; exit }' "$1"
}

# aggregate one pair into an edge list + node list
# args: file src_col tgt_col sample_col src_class tgt_class out_stem
aggregate() {
  local file="$1" sc="$2" tc="$3" smp="$4" scls="$5" tcls="$6" stem="$7"
  stem="${stem/ARG-MGE_ARG-MGE_/ARG-MGE_}"
  stem="${stem/ARG-VF_ARG-VF_/ARG-VF_}"
  local edges="$OUT_DIR/edges_${stem}.${EXT}" nodes="$OUT_DIR/nodes_${stem}.${EXT}"

  awk -F'\t' -v sc="$sc" -v tc="$tc" -v smp="$smp" \
             -v scls="$scls" -v tcls="$tcls" \
             -v minw="$MIN_WEIGHT" -v withn="$WITH_NSAMPLES" -v pfx="$PREFIX_LABELS" \
             -v sep="$SEP" -v edges="$edges" -v nodes="$nodes" '
    function lbl(v, cls) { return pfx ? cls ":" v : v }
    NR == 1 { next }
    $sc != "" && $tc != "" {
      s = lbl($sc, scls); t = lbl($tc, tcls)
      k = s SUBSEP t
      w[k]++
      if (!((k SUBSEP $smp) in seen)) { seen[k SUBSEP $smp] = 1; ns[k]++ }
      src[k] = s; tgt[k] = t
      if (!(s in ncls)) ncls[s] = scls; else if (ncls[s] != scls) clash[s] = 1
      if (!(t in ncls)) ncls[t] = tcls; else if (ncls[t] != tcls) clash[t] = 1
    }
    END {
      printf "Source%sTarget%sType%sWeight", sep, sep, sep > edges
      if (withn) printf "%sn_samples", sep > edges
      printf "\n" > edges

      n = 0
      for (k in w) {
        if (w[k] < minw) continue
        n++
        # sort by weight, descending, without relying on gawk extensions
        order[n] = k
      }
      for (i = 1; i < n; i++)
        for (j = i + 1; j <= n; j++)
          if (w[order[j]] > w[order[i]] ||
              (w[order[j]] == w[order[i]] && order[j] < order[i])) {
            tmpk = order[i]; order[i] = order[j]; order[j] = tmpk
          }

      for (i = 1; i <= n; i++) {
        k = order[i]
        printf "%s%s%s%sUndirected%s%d", src[k], sep, tgt[k], sep, sep, w[k] > edges
        if (withn) printf "%s%d", sep, ns[k] > edges
        printf "\n" > edges
      }

      printf "Id%sLabel%sClass\n", sep, sep > nodes
      for (v in ncls) printf "%s%s%s%s%s\n", v, sep, v, sep, ncls[v] > nodes

      nclash = 0
      for (v in clash) { nclash++; if (nclash <= 3) printf "  [!] label \"%s\" is used by more than one annotation type; Gephi will merge these into one node (use --prefix-labels to keep them apart)\n", v > "/dev/stderr" }
      printf "%d\t%d\n", n, nclash > "/dev/stdout"
    }
  ' "$file"
}

process_table() {
  local file="$1" table="$2"
  local c_arg_g c_arg_c c_mge_g c_mge_c c_vf_g c_vf_c c_smp
  c_smp=$(col_of "$file" "Sample")
  c_arg_g=$(col_of "$file" "ARG_gene");   c_arg_c=$(col_of "$file" "ARG_class")
  c_mge_g=$(col_of "$file" "MGE_gene");   c_mge_c=$(col_of "$file" "MGE_category")
  c_vf_g=$(col_of "$file" "VF_gene");     c_vf_c=$(col_of "$file" "VF_category")

  local rows; rows=$(awk 'NR>1' "$file" | wc -l | tr -d ' ')
  if [[ "$rows" -eq 0 ]]; then
    log "  $table: 0 data rows, skipped"
    return
  fi
  log "  $table ($rows reads)"

  local levels=()
  [[ "$LEVEL" == "class" || "$LEVEL" == "both" ]] && levels+=("class")
  [[ "$LEVEL" == "gene"  || "$LEVEL" == "both" ]] && levels+=("gene")

  for lv in "${levels[@]}"; do
    local a b acls bcls
    if [[ "$lv" == "class" ]]; then
      a=$c_arg_c; acls="ARG_class"
    else
      a=$c_arg_g; acls="ARG_gene"
    fi

    # ARG x MGE
    if [[ "$c_mge_c" -gt 0 || "$c_mge_g" -gt 0 ]]; then
      if [[ "$lv" == "class" ]]; then b=$c_mge_c; bcls="MGE_category"; else b=$c_mge_g; bcls="MGE_gene"; fi
      if [[ "$a" -gt 0 && "$b" -gt 0 ]]; then
        read -r ne nc < <(aggregate "$file" "$a" "$b" "$c_smp" "$acls" "$bcls" "${table}_ARG-MGE_${lv}")
        log "    ARG-MGE  $lv: $ne edges"
      fi
    fi

    # ARG x VF
    if [[ "$c_vf_c" -gt 0 || "$c_vf_g" -gt 0 ]]; then
      if [[ "$lv" == "class" ]]; then b=$c_vf_c; bcls="VF_category"; else b=$c_vf_g; bcls="VF_gene"; fi
      if [[ "$a" -gt 0 && "$b" -gt 0 ]]; then
        read -r ne nc < <(aggregate "$file" "$a" "$b" "$c_smp" "$acls" "$bcls" "${table}_ARG-VF_${lv}")
        log "    ARG-VF   $lv: $ne edges"
      fi
    fi

    # MGE x VF, only where both are present (the triple table)
    if [[ ( "$c_mge_c" -gt 0 || "$c_mge_g" -gt 0 ) && ( "$c_vf_c" -gt 0 || "$c_vf_g" -gt 0 ) ]]; then
      local m v mcls vcls
      if [[ "$lv" == "class" ]]; then
        m=$c_mge_c; v=$c_vf_c; mcls="MGE_category"; vcls="VF_category"
      else
        m=$c_mge_g; v=$c_vf_g; mcls="MGE_gene"; vcls="VF_gene"
      fi
      if [[ "$m" -gt 0 && "$v" -gt 0 ]]; then
        read -r ne nc < <(aggregate "$file" "$m" "$v" "$c_smp" "$mcls" "$vcls" "${table}_MGE-VF_${lv}")
        log "    MGE-VF   $lv: $ne edges"
      fi
    fi
  done
}

# ----------------------------------------------------------------- main
log "input      : $IN_DIR"
log "output     : $OUT_DIR"
log "level(s)   : $LEVEL"
log "min weight : $MIN_WEIGHT"
log "weight     : number of reads carrying both annotations"
log ""

found=0
for f in "$IN_DIR"/*.tsv "$IN_DIR"/*.csv; do
  [[ -f "$f" ]] || continue
  base="$(basename "$f")"
  name="${base%.*}"
  [[ "$name" == "summary" ]] && continue
  case "$name" in
    ARG-MGE|ARG-VF|ARG-MGE-VF) ;;
    *) log "  $base: not a recognised colocation table, skipped"; continue ;;
  esac
  process_table "$f" "$name"
  found=$((found + 1))
done
[[ $found -gt 0 ]] || die "no colocation tables found in $IN_DIR (expected ARG-MGE.tsv, ARG-VF.tsv, ARG-MGE-VF.tsv)"

# ----------------------------------------------------------------- per sample
if [[ $BY_SAMPLE -eq 1 ]]; then
  log ""
  log "writing per-sample edge lists to $OUT_DIR/by_sample/"
  mkdir -p "$OUT_DIR/by_sample"
  for f in "$IN_DIR"/ARG-MGE.tsv "$IN_DIR"/ARG-VF.tsv "$IN_DIR"/ARG-MGE-VF.tsv; do
    [[ -f "$f" ]] || continue
    name="$(basename "${f%.*}")"
    c_smp=$(col_of "$f" "Sample")
    mapfile -t samples < <(awk -F'\t' -v s="$c_smp" 'NR>1 {print $s}' "$f" | sort -u)
    for smp in "${samples[@]}"; do
      sub="$OUT_DIR/by_sample/.tmp_${name}_${smp}.tsv"
      { head -1 "$f"; awk -F'\t' -v s="$c_smp" -v want="$smp" 'NR>1 && $s==want' "$f"; } > "$sub"
      IN_SAVE="$OUT_DIR"; OUT_DIR="$OUT_DIR/by_sample"
      process_table "$sub" "${name}_${smp}" 2>/dev/null
      OUT_DIR="$IN_SAVE"
      rm -f "$sub"
    done
  done
fi

log ""
log "done. In Gephi: File > Import spreadsheet, choose the edges file as an"
log "Edges table, then the matching nodes file as a Nodes table. Weight is"
log "picked up automatically; colour nodes by the Class attribute."
