#!/usr/bin/env bash
# Static cell packing preview of one vessel geometry (no fluid steps).
# Usage (from this folder): ./run_preview.sh <bifurcation|stenosis|curved> [ranks] [packCells options]
# Builds the STL if missing, packs its bounding box, runs vessel_geometries,
# and leaves the result in output_<geometry>*/ (open *.xmf in ParaView).
set -euo pipefail
geom=${1:?geometry: bifurcation, stenosis or curved}
ranks=${2:-4}
shift $(( $# < 2 ? $# : 2 ))
# Most cells packed in the box touch the wall and are not placed; a WBC fits
# only within a few um of the axis, so many are packed for 1-2 to be placed
# (see README). The stenosis box holds fewer (packCells exits 3 with more).
case $geom in
  stenosis) default_opts="--hematocrit 0.3 --wbc 10 --seed 1" ;;
  *) default_opts="--hematocrit 0.3 --wbc 30 --seed 1" ;;
esac
pack_opts=${*:-$default_opts}
root=../..

[ -f "$geom.stl" ] || uv run -q "$root/scripts/geometry/make_vessels.py" "$geom" -o "$geom.stl"
# Pack the STL's bounding box (it starts at 0); positions are shifted by one
# dx, where the voxeliser puts the STL origin (see make_vessels.py).
dx_um=$(sed -n 's:.*<dx> *\([^ <]*\) *</dx>.*:\1:p' "config_$geom.xml" | awk '{print $1 * 1e6}')
read -r -a box < <(uv run -q --with trimesh python -c \
  "import sys,trimesh; print(*('%.4f' % v for v in trimesh.load(sys.argv[1]).bounds[1]))" "$geom.stl")
box+=(--walls xyz --offset "$dx_um" "$dx_um" "$dx_um")
# shellcheck disable=SC2086
"$root/tools/packCells/packCells" "${box[@]}" $pack_opts > "pack_$geom.log"
tr "\r" "\n" < "pack_$geom.log" | tail -n 3
mpirun -n "$ranks" ./vessel_geometries "config_$geom.xml" | tee "run_$geom.log" | grep -E "placed|VesselGeometries"
