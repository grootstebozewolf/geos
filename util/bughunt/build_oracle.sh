#!/usr/bin/env bash
# Build bin/oracle_arc: the NetTopologySuite.Proofs reference oracle (Coq-extracted
# kernel + oracle/driver.ml), trimmed to its circular-arc modes.
#
# The proofs repo's checked-in extracted.ml is older than driver.ml, so the full
# oracle_bin does not link without re-running Rocq extraction.  Everything up to
# the I_CIRCULAR section (driver.ml lines 1-1644) plus the ARC_SEGMENT_XY /
# ARC_ARC_DISTANCE / ARC_SEGMENT_DISTANCE section (lines 2157-2436) compiles
# against the checked-in extraction unchanged; we splice exactly those line
# ranges and add a small dispatch loop.  No oracle code is modified.
#
# Needs: git, ocamlfind, ocamlopt, zarith  (Debian: ocaml-nox ocaml-findlib libzarith-ocaml-dev)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PROOFS="${NTS_PROOFS:-$HERE/.nts-proofs}"
PIN="${NTS_PROOFS_PIN:-57f82369}"
if [ ! -d "$PROOFS/.git" ]; then
  git clone -q https://github.com/grootstebozewolf/NetTopologySuite.Proofs "$PROOFS"
fi
git -C "$PROOFS" cat-file -e "$PIN^{commit}" 2>/dev/null || git -C "$PROOFS" fetch -q origin
W="$HERE/.oracle-build"; rm -rf "$W"; mkdir -p "$W" "$HERE/bin"
git -C "$PROOFS" show "$PIN:oracle/driver.ml" > "$W/driver.ml"
# extracted.ml/.mli are gitignored in the proofs repo: they are produced by
# `make -f Makefile.gen` (Rocq extraction of Validate_binary64_extract.v).
# Use the copy in the proofs working tree (or NTS_EXTRACTED_DIR).
EXT="${NTS_EXTRACTED_DIR:-$PROOFS/oracle}"
if [ ! -f "$EXT/extracted.ml" ]; then
  echo "no extracted.ml in $EXT; run 'make -f Makefile.gen' in the proofs repo first" >&2; exit 1
fi
cp "$EXT/extracted.ml" "$EXT/extracted.mli" "$W/"
{
  sed -n 1,1644p "$W/driver.ml"
  sed -n 2157,2436p "$W/driver.ml"
  cat <<'OCAML'
(* ----- trimmed dispatch (GEOS curve bug hunter): arc modes only ----- *)
let () =
  let rec loop () =
    match (try Some (String.trim (input_line stdin)) with End_of_file -> None) with
    | None -> ()
    | Some "" -> loop ()
    | Some mode ->
      (match mode with
       | "ARC_LENGTH_INVARIANTS_EXACT" -> run_arc_length_invariants_exact ()
       | "ARC_LENGTH" -> run_arc_length ()
       | "ARC_AREA_INVARIANTS_EXACT" -> run_arc_area_invariants_exact ()
       | "ARC_AREA" -> run_arc_area ()
       | "ARC_DISTANCE" -> run_arc_distance ()
       | "ARC_ARC_XY" -> run_arc_arc_xy ()
       | "ARC_SEGMENT_XY" -> run_arc_segment_xy ()
       | "ARC_ARC_DISTANCE" -> run_arc_arc_distance ()
       | "ARC_SEGMENT_DISTANCE" -> run_arc_segment_distance ()
       | "ORIENT_EXACT" -> run_orient_exact ()
       | m -> Printf.printf "UNKNOWN_MODE %s\n" m);
      flush stdout; loop ()
  in loop ()
OCAML
} > "$W/driver_arc.ml"
cd "$W"
O="ocamlfind ocamlopt -package str,zarith"
$O -c extracted.mli && $O -c extracted.ml && $O -c driver_arc.ml
$O -linkpkg extracted.cmx driver_arc.cmx -o "$HERE/bin/oracle_arc"
echo "built $HERE/bin/oracle_arc from NetTopologySuite.Proofs@$PIN"
