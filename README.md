# BOSS jr. Documentation
Author: Gurpinder 

Last update: 18-12-2024
## Overview
**B**rassica **O**rthologue **S**earch **S**cript
BOSS jr. is a Python script designed for gene mapping between different genomic assemblies, specifically for limited sets of genes (a few hundred at most). It was originally created to find gene mappings between *Brassica napus* assemblies and orthologues between *Arabidopsis* and Brassicas. The script performs reciprocal BLAST searches to establish gene correspondences. The script runs blastn, so is suitable for search using CDS sequences.

## Dependencies
It has been tested with
- python=3.11
- ncbi-blast=2.12
- Biopython=1.80

## Command-Line Interface (CLI):
```bash
Boss_jr.py [-h] [--evalue_threshold EVALUE_THRESHOLD] input_file query_database target_database {oto,otm}
```
### Positional Arguments:

**input_file**: Path to the file containing the list of query gene IDs (one per line).

**query_database**: Path to the query database (Nucleotide database created using `makeblastdb').

**target_database**: Path to the target database (Nucleotide database created using `makeblastdb').

**flag**: Specifies the operation mode:

 oto (One-to-One reciprocal BLAST)
 
 otm (One-to-Many reciprocal BLAST)

### Optional Arguments:

--evalue_threshold (default: 1e-50): The E-value threshold for filtering BLAST results.

## Example:
```bash
python boss_jr.py genes.txt query_db.fasta target_db.fasta oto --evalue_threshold 1e-30
```

## Notes
- The script prints the mapping results to the console in a structured format, indicating whether each query gene has a valid mapping in the target database.
 
- The script is inefficient for large-scale analyses due to the use of Python dictionaries and loops.

- For large-scale gene mappings, consider using pandas for improved performance.
