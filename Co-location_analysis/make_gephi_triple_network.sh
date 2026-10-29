#!/usr/bin/env bash
#
# make_gephi_triple_network.sh
#
# One combined Gephi network from ARG-MGE-VF.tsv, at class/category level.
# Every row of that table is a single read carrying an ARG, an MGE and a VF
# annotation, so each read contributes three undirected edges:
#
#   ARG_class  -- MGE_category
#   ARG_class  -- VF_category
#   MGE_category -- VF_category
#
# All three go into ONE edge file, which is what you want for localising MGE and
# VF on ARG-carrying reads: the ARG classes sit between the MGE and VF sides and
# you can see which combinations travel together on the same read.
#
# Outputs (2 files):
#   edges_ARG-MGE-VF_class.tsv   Source Target Type Weight Pair
#   nodes_ARG-MGE-VF_class.tsv   Id Label Class
#
# Pair says which of the three relationships an edge is, so you can colour or
# filter edges by it in Gephi. Drop it with --minimal for a strict 4-column file.
#
# Usage:
#   ./make_gephi_triple_network.sh -i colocation/ARG-MGE-VF.tsv -o gephi/
#   ./make_gephi_triple_network.sh -i colocation/ARG-MGE-VF.tsv -o gephi/ --min-weight 5
#
set -euo pipefail

IN=""; OUT_DIR="gephi"; MIN_WEIGHT=5; MINIMAL=0; SEP=$'\t'; EXT="tsv"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log() { printf '%s\n' "$*" >&2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)  IN="$2"; shift 2 ;;
    -o|--outdir) OUT_DIR="$2"; shift 2 ;;
    --min-weight) MIN_WEIGHT="$2"; shift 2 ;;
    --minimal)   MINIMAL=1; shift ;;
    --csv)       SEP=","; EXT="csv"; shift ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed 's/^#//; s/^ //' | head -n -1; exit 0 ;;
    *) die "unknown option: $1 (try --help)" ;;
  esac
done

[[ -n "$IN" ]] || die "-i/--input is required (path to ARG-MGE-VF.tsv)"
[[ -f "$IN" ]] || die "no such file: $IN"
mkdir -p "$OUT_DIR"

EDGES="$OUT_DIR/edges_ARG-MGE-VF_class.${EXT}"
NODES="$OUT_DIR/nodes_ARG-MGE-VF_class.${EXT}"

log "input      : $IN"
log "min weight : $MIN_WEIGHT (edges below this are dropped)"
log "weight     : reads carrying both annotations"
log ""

awk -F'\t' -v sep="$SEP" -v minw="$MIN_WEIGHT" -v minimal="$MINIMAL" \
           -v edges="$EDGES" -v nodes="$NODES" '
  function trim(x) { gsub(/^[[:space:]]+|[[:space:]]+$/, "", x); return x }
  function add(a, acls, b, bcls, pair,    k) {
    if (a == "" || b == "") return
    k = a SUBSEP b SUBSEP pair
    w[k]++; src[k] = a; tgt[k] = b; prs[k] = pair
    if (!(a in ncls)) ncls[a] = acls; else if (ncls[a] != acls) clash[a] = 1
    if (!(b in ncls)) ncls[b] = bcls; else if (ncls[b] != bcls) clash[b] = 1
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) idx[$i] = i
    for (h in idx) present[h] = 1
    if (!("ARG_class" in idx) || !("MGE_category" in idx) || !("VF_category" in idx)) {
      print "ERROR: input needs ARG_class, MGE_category and VF_category columns" > "/dev/stderr"
      exit 1
    }
    a = idx["ARG_class"]; m = idx["MGE_category"]; v = idx["VF_category"]
    next
  }
  {
    reads++
    A = trim($a); M = trim($m); V = trim($v)
    add(A, "ARG_class",    M, "MGE_category", "ARG-MGE")
    add(A, "ARG_class",    V, "VF_category",  "ARG-VF")
    add(M, "MGE_category", V, "VF_category",  "MGE-VF")
  }
  END {
    if (reads == 0) { print "ERROR: no data rows in input" > "/dev/stderr"; exit 1 }

    printf "Source%sTarget%sType%sWeight", sep, sep, sep > edges
    if (!minimal) printf "%sPair", sep > edges
    printf "\n" > edges

    n = 0
    for (k in w) if (w[k] >= minw) order[++n] = k
    for (i = 1; i < n; i++)
      for (j = i + 1; j <= n; j++)
        if (w[order[j]] > w[order[i]] ||
            (w[order[j]] == w[order[i]] && order[j] < order[i])) {
          t = order[i]; order[i] = order[j]; order[j] = t
        }

    for (i = 1; i <= n; i++) {
      k = order[i]
      printf "%s%s%s%sUndirected%s%d", src[k], sep, tgt[k], sep, sep, w[k] > edges
      if (!minimal) printf "%s%s", sep, prs[k] > edges
      printf "\n" > edges
      kept[prs[k]]++; keptw[prs[k]] += w[k]
      used[src[k]] = 1; used[tgt[k]] = 1
      deg[src[k]]++; deg[tgt[k]]++
    }

    printf "Id%sLabel%sClass\n", sep, sep > nodes
    nn = 0
    for (x in used) { printf "%s%s%s%s%s\n", x, sep, x, sep, ncls[x] > nodes; nn++; cls_n[ncls[x]]++ }

    total = 0; for (k in w) total++
    printf "reads in table   : %d\n", reads > "/dev/stderr"
    printf "edges before cut : %d\n", total > "/dev/stderr"
    printf "edges kept       : %d\n", n > "/dev/stderr"
    printf "nodes            : %d (ARG %d, MGE %d, VF %d)\n", nn, \
           cls_n["ARG_class"]+0, cls_n["MGE_category"]+0, cls_n["VF_category"]+0 > "/dev/stderr"
    printf "\nby relationship:\n" > "/dev/stderr"
    for (p in kept)
      printf "  %-8s %3d edges, %d read-pairs (%.0f%% of reads retained)\n", \
             p, kept[p], keptw[p], 100*keptw[p]/reads > "/dev/stderr"
    for (x in clash)
      printf "\n  [!] \"%s\" appears as more than one annotation type; Gephi will merge it into one node\n", x > "/dev/stderr"
  }
' "$IN"

log ""
log "wrote:"
log "  $EDGES"
log "  $NODES"
log ""
log "In Gephi: import the edges file as an Edges table, then the nodes file as a"
log "Nodes table. Colour nodes by Class (ARG / MGE / VF), size by degree, and"
log "colour or filter edges by Pair. Weight is picked up automatically."
