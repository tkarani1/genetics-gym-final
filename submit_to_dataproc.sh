#!/bin/bash
# General helper script to zip resources and submit any Python script to Dataproc
# Usage: ./submit_to_dataproc.sh <cluster_name> <script_path> [additional args]

set -e  # Exit on error

if [ $# -lt 2 ]; then
    echo "Usage: $0 <cluster_name> <script_path> [additional_args...]"
    echo ""
    echo "Example:"
    echo "  ./submit_to_dataproc.sh tk data_preprocessing/VSMs/3_coalesce_each_VSM/coalesce_each_VSM.py"
    echo ""
    exit 1
fi

CLUSTER_NAME=$1
SCRIPT_PATH=$2
shift 2  # Remove first two args, remaining are passed to hailctl

# Get script directory
SCRIPT_DIR=$(dirname "$0")
PROJECT_ROOT=$(cd "$SCRIPT_DIR" && pwd)

echo "Project root: $PROJECT_ROOT"
echo "Cluster: $CLUSTER_NAME"
echo "Script: $SCRIPT_PATH"
echo ""

# Navigate to project root
cd "$PROJECT_ROOT"

# Create fresh resources.zip
echo "Creating fresh resources.zip..."
rm -f resources.zip
cd resources && zip -q -r ../resources.zip . && cd ..
echo "✓ resources.zip created ($(du -h resources.zip | cut -f1))"

echo ""
echo "Submitting to Dataproc cluster '$CLUSTER_NAME'..."
hailctl dataproc submit "$CLUSTER_NAME" \
  "$SCRIPT_PATH" \
  --pyfiles=resources.zip \
  "$@"

echo ""
echo "✓ Job submitted successfully!"
