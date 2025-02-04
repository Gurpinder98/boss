# Overview of functions in Boss_jr.py

```shell build_sequence_dict(genes: list, database_path: str) -> dict```

Fetches sequences for specified genes from a given FASTA database.

Input: List of gene IDs and path to the FASTA file.

Output: Dictionary of gene IDs mapped to their sequences.

```shell blast_gene(Query_dict: dict, gene: str, database_to_search_in: str, evalue_threshold: float) -> list ```

Performs a BLAST search for a given gene against a target database.

Input: Query dictionary, gene ID, database path, and E-value threshold.

Output: List of BLAST alignments.

```shell pretty_printing(putative_orthologue_dict: dict) -> None ```

Prints orthologue mappings in a readable format.

```shell one_to_one_reciprocal_blast(query_genes: list, Query_database: str, Target_database: str, evalue_threshold: float) -> None ```

Performs a one-to-one reciprocal BLAST search.

Steps:

For each query gene, identifies the best hit in the target genome.

Performs a BLAST search for the best hit in the query genome.

If the best hit maps back to the original query gene, it is assigned as an orthologue.

```shell one_to_many_reciprocal_blast(query_genes: list, Query_database: str, Target_database: str, evalue_threshold: float) -> None ```

Performs a one-to-many reciprocal BLAST search.

Steps:

Identifies all hits above the E-value threshold for each query gene.

For each hit, performs a BLAST search against the query database.

If the best hit matches the original query gene, it is considered a valid mapping.

```shell read_query_genes(gene_list_file: str) -> list```

Reads a gene list from a file (one ID per line).
