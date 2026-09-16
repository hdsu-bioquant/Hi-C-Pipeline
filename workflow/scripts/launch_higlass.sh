#!/usr/bin/env bash
# launch_higlass.sh
#   Spin up a local HiGlass server (higlass/higlass-docker:latest), bind-mount
#   the pipeline's higlass/data/ directory, and auto-ingest every track listed
#   in higlass/manifest.yaml.
#
# Usage:
#   bash launch_higlass.sh --data /abs/.../results/higlass/data \
#                          --manifest /abs/.../results/higlass/manifest.yaml \
#                          [--port 8989] [--container higlass-server]
#
# Requires Docker OR Podman on the host.  On the BioQuant cluster this is
# typically only available on certain interactive nodes — not on `single`.

set -euo pipefail

DATA=""
MANIFEST=""
PORT=8989
CONTAINER_NAME="higlass-server"
IMAGE="higlass/higlass-docker:latest"

while [ $# -gt 0 ]; do
    case "$1" in
        --data)      DATA="$2"; shift 2;;
        --manifest)  MANIFEST="$2"; shift 2;;
        --port)      PORT="$2"; shift 2;;
        --container) CONTAINER_NAME="$2"; shift 2;;
        *) echo "unknown flag: $1" >&2; exit 2;;
    esac
done

if [ -z "$DATA" ] || [ -z "$MANIFEST" ]; then
    echo "Usage: bash $0 --data <data_dir> --manifest <manifest.yaml> [--port 8989]" >&2
    exit 2
fi

# Pick docker or podman
if command -v docker >/dev/null 2>&1; then RUNTIME=docker
elif command -v podman >/dev/null 2>&1; then RUNTIME=podman
else
    echo "ERROR: neither docker nor podman found on PATH" >&2
    exit 1
fi

echo "[higlass] using $RUNTIME · image=$IMAGE · port=$PORT"

# Pull image (no-op if already cached)
$RUNTIME pull "$IMAGE"

# Stop / remove any previous container of the same name
$RUNTIME rm -f "$CONTAINER_NAME" 2>/dev/null || true

# Run HiGlass server.
# The official image binds-volume mounts /data for the tileset directory.
$RUNTIME run -d \
    --name "$CONTAINER_NAME" \
    --publish "$PORT:80" \
    --volume "$DATA":/data:ro \
    --volume "$DATA":/tmp/hg-data \
    "$IMAGE"

echo "[higlass] container started; waiting 10s for the server to come up..."
sleep 10

# Ingest each track from the manifest
echo "[higlass] ingesting tracks from $MANIFEST ..."
# Crude YAML parser — just pulls 'uid' / 'filename' / 'filetype' / 'datatype' lines
python3 - "$MANIFEST" <<'PY' > /tmp/_higlass_ingest_cmds.sh
import sys, yaml, shlex
with open(sys.argv[1]) as fh:
    manifest = yaml.safe_load(fh)
for tr in manifest.get("tracks", []):
    uid      = tr.get("uid", "")
    filename = tr.get("filename", "")
    filetype = tr.get("filetype", "cooler")
    datatype = tr.get("datatype", "matrix")
    if not uid or not filename:
        continue
    print(f"docker exec {shlex.quote('higlass-server')} python higlass-server/manage.py "
          f"ingest_tileset --filename /data/{shlex.quote(filename)} "
          f"--filetype {filetype} --datatype {datatype} --uid {shlex.quote(uid)} || true")
PY

bash /tmp/_higlass_ingest_cmds.sh
rm -f /tmp/_higlass_ingest_cmds.sh

echo
echo "==========================================================================="
echo "HiGlass server is running at:    http://localhost:$PORT"
echo "Container name:                  $CONTAINER_NAME"
echo "To stop:                         $RUNTIME stop $CONTAINER_NAME"
echo "To remove:                       $RUNTIME rm $CONTAINER_NAME"
echo "==========================================================================="
