"""
Convert TSV.bgz files to optimized Parquet format
Run this once to get 5-20x speedup on all subsequent operations
"""
import pandas as pd
import sys
from pathlib import Path
import time
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from resources.paths import WRITE_VSM_LINKER_TABLES_BASE, VSM_TABLE_NAMES


def convert_to_parquet(tsv_path, parquet_path, key_cols=['chrom', 'pos', 'ref', 'alt'], 
                      chunksize=5_000_000):
    """
    Convert TSV.bgz to optimized Parquet with:
    - Sorted by key columns
    - Categorical types for low-cardinality columns
    - Snappy compression (fast)
    - Optimal row group size
    """
    print(f"\nConverting: {tsv_path}")
    print(f"       to: {parquet_path}")
    
    start = time.time()
    
    # Read in chunks to handle large files
    print("Reading TSV...")
    chunks = []
    total_rows = 0
    
    for i, chunk in enumerate(pd.read_csv(tsv_path, sep='\t', compression='gzip', chunksize=chunksize)):
        total_rows += len(chunk)
        print(f"  Chunk {i+1}: {len(chunk):,} rows (total: {total_rows:,})", end='\r')
        chunks.append(chunk)
    
    print(f"\n  Total rows read: {total_rows:,}")
    
    # Combine chunks
    print("Combining chunks...")
    df = pd.concat(chunks, ignore_index=True)
    del chunks  # Free memory
    
    # Optimize dtypes
    print("Optimizing data types...")
    
    # Categorical for low-cardinality string columns
    categorical_cols = ['chrom', 'ref', 'alt']
    for col in categorical_cols:
        if col in df.columns:
            df[col] = df[col].astype('category')
            print(f"  {col}: {df[col].nunique()} unique values → category")
    
    # Float32 for score columns (saves 50% memory)
    float_cols = [col for col in df.columns if df[col].dtype == 'float64']
    for col in float_cols:
        df[col] = df[col].astype('float32')
    print(f"  Converted {len(float_cols)} float64 columns → float32")
    
    # Sort by key columns (huge speedup for subsequent groupby)
    print(f"Sorting by {key_cols}...")
    df = df.sort_values(key_cols)
    
    # Write to Parquet
    print("Writing Parquet...")
    df.to_parquet(
        parquet_path,
        engine='pyarrow',
        compression='snappy',  # Fast compression, good ratio
        index=False,
        row_group_size=1_000_000  # 1M rows per group (optimal for queries)
    )
    
    elapsed = time.time() - start
    
    # Get file sizes
    import os
    tsv_size = os.path.getsize(tsv_path) / (1024**3)  # GB
    parquet_size = os.path.getsize(parquet_path) / (1024**3)  # GB
    
    print(f"\n✓ Conversion complete!")
    print(f"  Time: {elapsed:.1f}s ({elapsed/60:.1f}min)")
    print(f"  TSV size: {tsv_size:.2f} GB")
    print(f"  Parquet size: {parquet_size:.2f} GB")
    print(f"  Compression ratio: {tsv_size/parquet_size:.1f}x")
    print(f"  Rows: {len(df):,}")


def convert_all_methods(methods=None):
    """
    Convert all VSM methods from TSV.bgz to Parquet
    """
    if methods is None:
        methods = ['AM', 'PRIMATEAI3D', 'POPEVE', 'MISFIT', 'MPC', 'POLYPHEN', 
                   'ESM1B', 'CPT', 'CADD', 'GPN_MSA', 'REVEL', 'RASP', 'PROTEINMPNN']
    
    print(f"{'='*70}")
    print(f"Converting {len(methods)} methods to Parquet")
    print(f"{'='*70}")
    
    total_start = time.time()
    success_count = 0
    
    for method in methods:
        try:
            tsv_file = f"{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[method]}.tsv.bgz"
            parquet_file = f"{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[method]}.parquet"
            
            convert_to_parquet(tsv_file, parquet_file)
            success_count += 1
            
        except FileNotFoundError:
            print(f"✗ File not found: {tsv_file}")
        except KeyError:
            print(f"✗ Method '{method}' not found in VSM_TABLE_NAMES")
        except Exception as e:
            print(f"✗ Error converting {method}: {e}")
            import traceback
            traceback.print_exc()
    
    total_time = time.time() - total_start
    
    print(f"\n{'='*70}")
    print(f"Conversion complete!")
    print(f"  Successful: {success_count}/{len(methods)}")
    print(f"  Total time: {total_time:.1f}s ({total_time/60:.1f}min)")
    print(f"{'='*70}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert TSV.bgz to optimized Parquet')
    parser.add_argument('--methods', nargs='+', help='Specific methods to convert')
    parser.add_argument('--all', action='store_true', help='Convert all methods')
    parser.add_argument('--file', help='Convert a single file (TSV path)')
    parser.add_argument('--output', help='Output parquet path (if using --file)')
    
    args = parser.parse_args()
    
    if args.file:
        # Convert single file
        output = args.output or args.file.replace('.tsv.bgz', '.parquet')
        convert_to_parquet(args.file, output)
    
    elif args.all or args.methods:
        # Convert multiple methods
        convert_all_methods(args.methods)
    
    else:
        # Default: convert common methods
        print("Converting common methods (use --all for all methods)")
        convert_all_methods(['AM', 'PRIMATEAI3D', 'POPEVE', 'MISFIT', 'MPC', 'POLYPHEN'])
