# Dataproc Submit Guide

Helper scripts to automatically zip resources and submit jobs to Dataproc.

## Quick Start

### Submit the coalesce job

```bash
cd /Users/tk508/Work/GenGym/genetics-gym-final
./submit_to_dataproc.sh tk data_preprocessing/VSMs/3_coalesce_each_VSM/coalesce_each_VSM.py
```

### Submit any Python script

```bash
./submit_to_dataproc.sh <cluster_name> <script_path> [additional_hailctl_args]
```

**Examples:**

```bash
# Submit coalesce script
./submit_to_dataproc.sh tk data_preprocessing/VSMs/3_coalesce_each_VSM/coalesce_each_VSM.py

# Submit import_AM script with extra memory
./submit_to_dataproc.sh tk data_preprocessing/VSMs/1_import_each_VSM_to_tsv/import_AM.py \
  --properties spark.driver.memory=16g

# Submit join_AM script
./submit_to_dataproc.sh tk data_preprocessing/VSMs/2_join_each_VSM_to_linker/join_AM.py
```

## What it does

The `submit_to_dataproc.sh` script:

1. Creates a fresh `resources.zip` from your `resources/` directory
2. Submits your Python script to the specified Dataproc cluster
3. Automatically includes `resources.zip` via `--pyfiles`
4. Passes any additional arguments to `hailctl dataproc submit`

## Why use this?

- ✓ Always uses latest code from `resources/`
- ✓ No need to manually zip resources
- ✓ Consistent submission command
- ✓ Easier to track what was submitted

## Files

- **`submit_to_dataproc.sh`** - General purpose submission script (recommended)
- **`data_preprocessing/VSMs/3_coalesce_each_VSM/submit_coalesce.sh`** - Specific script for coalesce job

## Alternative: Manual submission

If you prefer to manually control the zip:

```bash
cd /Users/tk508/Work/GenGym/genetics-gym-final
cd resources && zip -r ../resources.zip . && cd ..
hailctl dataproc submit tk <your_script.py> --pyfiles=resources.zip
```

## Troubleshooting

### "Permission denied"
Make the script executable:
```bash
chmod +x submit_to_dataproc.sh
```

### "resources.zip: No such file or directory"
The script creates it automatically, but if you see this error, run from project root:
```bash
cd /Users/tk508/Work/GenGym/genetics-gym-final
./submit_to_dataproc.sh tk <script>
```

### "ModuleNotFoundError: No module named 'resources'"
Make sure:
1. `resources/` directory has `__init__.py`
2. You're using the submit script (not raw `hailctl` command)
3. Script is run from project root
