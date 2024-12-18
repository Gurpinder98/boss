#!/usr/bin/env python
"""
BOSS jr. (More generalised(?) version of the bigger BOSS.py)

Originally created to find gene mappings between brassica napus assemblies, 
and orthologues between Arabidosis and Brassicas.

This version is more generalised, useful for a limited number of genes, 
a few hundered at max perhaps.

If you need to map more than that, the inefficient use of python 
dictionaries and for loops could be a bottleneck. Use atleast pandas
for sanity.

Gurpinder
gurpinder-singh.sidhu@jic.ac.uk
last big update: 18-12-2024
"""

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from tempfile import TemporaryDirectory
from pathlib import Path
import argparse


def build_sequence_dict(genes: list, database_path: str) -> dict:
    """
    A function to fetch sequence strings for a set of genes.
    Input: list of gene ids, path to the file containing those sequences.
    """
    Query_dict = {}
    with open(database_path) as in_f:
        for fasta in SeqIO.parse(in_f, 'fasta'):
            name, seq = fasta.id, str(fasta.seq)
            if name in genes:
                Query_dict[name] = seq
    # quick check to see if there are any gene ids not found in the database.
    for gene in genes:
        if gene not in list(Query_dict.keys()):
            # print a quick heads up!
            print(f"WARNING: {gene} not found in {database_path}.")
    return Query_dict


def blast_gene(Query_dict: dict, gene: str, database_to_seach_in: str, evalue_threshold: float) -> list:
    """
    A function that runs a blast search for a gene against a database. 
    Returns a list of blast alignments.
    """
    with TemporaryDirectory(prefix="blastn") as tmpdir:
        temp_path = Path(tmpdir)
        query_path = temp_path/"query.fa"
        result_path = temp_path/"blast.xml"

        #write sequence for the gene from query_dict into a file
        query_path.write_text(Query_dict[gene]) 

        blast_command_line = NcbiblastnCommandline(task = "blastn", query = str(query_path), db = str(database_to_seach_in), evalue = evalue_threshold, outfmt=5, out=str(result_path))

        try:
            blast_command_line()
        except:
            print(f"Problem occured at gene {gene}")
        
        with open(result_path) as in_f:
            blast_record = NCBIXML.read(in_f)
    return blast_record.alignments


def pretty_printing(putative_orthologue_dict: dict) -> None:
    """
    A function to print things off on the console.
    """
    print("\nOrthologues:")
    for gene in putative_orthologue_dict:
        if len(putative_orthologue_dict[gene]) > 0: print(f"{gene}: {', '.join(putative_orthologue_dict[gene])}")
        else: print(f"{gene}: NA")


def one_to_one_reciprocal_blast(query_genes: list, Query_database: str, Target_database: str, evalue_threshold: float) -> None:
    """
    A function to perform one to one reciprocal blast. This assumes that your target database has only one mapping for a query gene.
    How it works:
        1. For every gene in query genes list, it gets the best hit for it in the target genome.
        2. For every best hit gene, it runs a blast search for that gene within the query genome. 
        3. If the best hit in step 2 matches the query gene, it assigns it as a mapping.
    """

    dict_of_seqs_from_query_database = build_sequence_dict(query_genes, Query_database)

    # for every gene in query genes seq dict, get the best hit in the target genome
    dict_of_best_hits_for_query_genes = {}
    for gene in dict_of_seqs_from_query_database:
        result = blast_gene(dict_of_seqs_from_query_database, gene, Target_database, evalue_threshold)
        dict_of_best_hits_for_query_genes[gene] = [[res.hit_id for res in result][0]] #only take the best hit

    # ensure that blast result for each gene was correctly processed. 
    assert dict_of_seqs_from_query_database.keys() == dict_of_best_hits_for_query_genes.keys()
    
    putative_mappings = {}
    for gene in dict_of_best_hits_for_query_genes: #for every query gene
        putative_mappings[gene] = [] 
        dict_of_seqs_from_target_database = build_sequence_dict(dict_of_best_hits_for_query_genes[gene], Target_database) #get sequences for all the hits that it had
        for hit in dict_of_seqs_from_target_database: # for each hit
            result = blast_gene(dict_of_seqs_from_target_database, hit, Query_database, evalue_threshold) #blast it against query database
            best_hit = [res.hit_id for res in result][0] #only take the best hit
            if best_hit == gene:
                putative_mappings[gene].append(hit)

    # finally do the printing
    print("\nType of mapping: One-to-one")
    pretty_printing(putative_mappings)


def one_to_many_reciprocal_blast(query_genes: list, Query_database: str, Target_database: str, evalue_threshold: float):
    """
    A function to perform one to many reciprocal blast. This assumes that your target database has many possible mappings for a query gene.
    How it works:
        1. For every gene in query genes list, it gets all the hits above the e_value threshold.
        2. For every above e_value threshold hit gene, it runs a blast search for that gene within the query genome. 
        3. If the best hit for the e_value threshold gene in step 2 matches the query gene, it assigns it as a mapping.
        4. Step 3 is repeated for all e_value threshold genes.
    """
    dict_of_seqs_from_query_database = build_sequence_dict(query_genes, Query_database)

    # for every gene in query genes seq dict, get the best hit in the target genome
    dict_of_best_hits_for_query_genes = {}
    for gene in dict_of_seqs_from_query_database:
        result = blast_gene(dict_of_seqs_from_query_database, gene, Target_database, evalue_threshold)
        dict_of_best_hits_for_query_genes[gene] = [res.hit_id for res in result] #take all the hits in target database

    # ensure that blast result for each gene was correctly processed. 
    assert dict_of_seqs_from_query_database.keys() == dict_of_best_hits_for_query_genes.keys()
    
    putative_mappings = {}
    for gene in dict_of_best_hits_for_query_genes: #for every query gene
        putative_mappings[gene] = [] 
        dict_of_seqs_from_target_database = build_sequence_dict(dict_of_best_hits_for_query_genes[gene], Target_database) #get sequences for all the hits that it had
        for hit in dict_of_seqs_from_target_database: # for each hit
            result = blast_gene(dict_of_seqs_from_target_database, hit, Query_database, evalue_threshold) #blast it against query database
            best_hit = [res.hit_id for res in result][0] #only take the best hit
            if best_hit == gene:
                putative_mappings[gene].append(hit)

    # finally do the printing
    print("\nType of mapping: One-to-many")
    pretty_printing(putative_mappings)


def read_query_genes(gene_list_file):
    """
    A helper function to read the gene list file. It assumes one id per line.
    """
    genes = []
    try:
        with open(gene_list_file, "r") as in_f:
            genes = [line.strip().replace("\n", "") for line in in_f.readlines()]
    except FileNotFoundError:
        print(f"{gene_list_file} doesn't exist.")
    return genes


def cli_interface():
    parser = argparse.ArgumentParser(
        description="A script to perform one to one, or one to many reciprocal blast searches between query and target blast databases."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input file containing the list of query genes (one gene per line)."
    )
    parser.add_argument(
        "query_database",
        type=str,
        help="Path to the query database."
    )
    parser.add_argument(
        "target_database",
        type=str,
        help="Path to the target database."
    )
    parser.add_argument(
        "flag",
        type=str,
        choices=["oto", "otm"],
        help="Flag specifying the operation mode. Options: 'oto' (One to One), 'otm' (One to Many)."
    )
    parser.add_argument(
        "--evalue_threshold",
        type=float,
        default=1e-50,
        help="E-value threshold for processing (default: 1e-50)."
    )

    return parser.parse_args()


if __name__ == "__main__":

    args = cli_interface()
    input_file = args.input_file
    Query_database = args.query_database
    Target_database = args.target_database
    evalue_threshold = args.evalue_threshold
    flag = args.flag

    Query_genes = read_query_genes(input_file)
    print(f"BOSS jr.\n-----------------------\n\nTarget_database: {Target_database}\nQuery_database:{Query_database}\ne-value:{evalue_threshold}")
    if flag == "oto":
        print("Mode: One to One reciprocal blast")
        one_to_one_reciprocal_blast(Query_genes, Query_database, Target_database, evalue_threshold)
    elif flag == "otm":
        print("Mode: One to Many reciprocal blast")
        one_to_many_reciprocal_blast(Query_genes, Query_database, Target_database, evalue_threshold)
    