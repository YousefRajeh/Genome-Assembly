#!/usr/bin/env python3
"""
De Bruijn Graph (DBG) Genome Assembler with QUAST Integration

Usage:
    python DeBruijnGraph.py --input <input.fastq> --kmer_size <k> --output <output.fasta> 
                            [--reference <reference.fasta>] [--min_coverage <min_cov>] [--export_gfa <graph.gfa>]

Arguments:
    --input: Path to the input FASTQ file
    --kmer_size: Size of k-mers to use for graph construction
    --output: Path to the output FASTA file for contigs
    --reference: Path to the reference genome for evaluation (optional)
    --min_coverage: Minimum k-mer coverage (default: 1)
    --export_gfa: Export assembly graph in GFA format (optional)
"""

import argparse
import os
import sys
from collections import defaultdict, Counter
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import tempfile
import re
import shutil

class DeBruijnGraphAssembler:
    def __init__(self, k, min_coverage=1):
        """
        Initialize the De Bruijn Graph assembler.
        
        Args:
            k (int): k-mer size
            min_coverage (int): Minimum k-mer coverage
        """
        self.k = k
        self.min_coverage = min_coverage
        self.kmers = Counter()
        self.graph = nx.DiGraph()
        self.contigs = []
        
    def read_fastq(self, fastq_file):
        """
        Read sequences from a FASTQ file.
        
        Args:
            fastq_file (str): Path to the FASTQ file
        """
        reads = []
        for record in SeqIO.parse(fastq_file, "fastq"):
            reads.append(str(record.seq))
        return reads
    
    def extract_kmers(self, reads):
        """
        Extract k-mers from reads and count their occurrences.
        
        Args:
            reads (list): List of read sequences
        """
        for read in reads:
            # Skip reads shorter than k
            if len(read) < self.k:
                continue
                
            # Extract k-mers from the read
            for i in range(len(read) - self.k + 1):
                kmer = read[i:i + self.k]
                # Avoid k-mers with N (ambiguous bases)
                if 'N' not in kmer:
                    self.kmers[kmer] += 1
        
        # Filter k-mers by coverage
        self.kmers = Counter({kmer: count for kmer, count in self.kmers.items() 
                              if count >= self.min_coverage})
        
        print(f"Extracted {len(self.kmers)} k-mers with coverage >= {self.min_coverage}")
    
    def build_graph(self):
        """
        Build the De Bruijn graph from k-mers.
        """
        # Create nodes and edges
        for kmer, count in self.kmers.items():
            # k-1 prefixes and suffixes become nodes
            prefix = kmer[:-1]
            suffix = kmer[1:]
            
            # Add nodes if they don't exist
            if prefix not in self.graph:
                self.graph.add_node(prefix)
            if suffix not in self.graph:
                self.graph.add_node(suffix)
            
            # Add edge with weight = k-mer count
            if self.graph.has_edge(prefix, suffix):
                self.graph[prefix][suffix]['weight'] += count
                self.graph[prefix][suffix]['kmers'].append(kmer)
            else:
                self.graph.add_edge(prefix, suffix, weight=count, kmers=[kmer])
        
        print(f"Built graph with {self.graph.number_of_nodes()} nodes and {self.graph.number_of_edges()} edges")
    
    def find_contigs(self):
        """
        Find contigs by identifying Eulerian paths in the graph.
        """
        # Find non-branching paths (simple paths)
        def get_in_degree(node):
            return self.graph.in_degree(node)
        
        def get_out_degree(node):
            return self.graph.out_degree(node)
        
        # Start from nodes that have out_degree != in_degree or out_degree > 1
        start_nodes = [node for node in self.graph.nodes() 
                       if get_in_degree(node) != get_out_degree(node) or get_out_degree(node) > 1]
        
        # If no suitable start nodes, use any node (circular genome case)
        if not start_nodes and self.graph.nodes():
            start_nodes = [list(self.graph.nodes())[0]]
        
        contigs = []
        visited_edges = set()
        
        # Process each start node
        for start_node in start_nodes:
            if get_out_degree(start_node) == 0:
                continue
                
            for _, end_node in self.graph.out_edges(start_node):
                edge = (start_node, end_node)
                if edge in visited_edges:
                    continue
                
                # Start a new path
                path = [start_node]
                current_node = end_node
                
                while True:
                    path.append(current_node)
                    visited_edges.add((path[-2], current_node))
                    
                    # If we reached a junction or end of path
                    if get_out_degree(current_node) != 1 or get_in_degree(current_node) != 1:
                        break
                    
                    # Continue path
                    next_nodes = list(self.graph.successors(current_node))
                    if not next_nodes:
                        break
                    next_node = next_nodes[0]
                    if (current_node, next_node) in visited_edges:
                        break
                    current_node = next_node
                
                # Convert path to contig sequence
                if len(path) > 1:
                    contig = path[0]
                    for i in range(1, len(path)):
                        contig += path[i][-1]
                    contigs.append(contig)
        
        # Add isolated cycles (completely circular paths)
        components = list(nx.weakly_connected_components(self.graph))
        for component in components:
            if all(get_in_degree(node) == get_out_degree(node) == 1 for node in component):
                cycle = list(nx.simple_cycles(self.graph.subgraph(component)))[0]
                if cycle:
                    contig = cycle[0]
                    for i in range(1, len(cycle)):
                        contig += cycle[i][-1]
                    contig += cycle[0][-1]  # Close the cycle
                    contigs.append(contig)
        
        # Sort contigs by length (descending)
        self.contigs = sorted(contigs, key=len, reverse=True)
        print(f"Generated {len(self.contigs)} contigs")
    
    def export_gfa(self, gfa_file):
        """
        Export the assembly graph in GFA format.
        
        Args:
            gfa_file (str): Path to the output GFA file
        """
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(gfa_file), exist_ok=True)
        
        with open(gfa_file, 'w') as f:
            # Write header
            f.write("H\tVN:Z:1.0\n")
            
            # Write nodes (segments)
            for node in self.graph.nodes():
                f.write(f"S\t{node}\t{node}\n")
            
            # Write edges (links)
            for u, v, data in self.graph.edges(data=True):
                # In GFA, we connect from the end of u to the start of v
                f.write(f"L\t{u}\t+\t{v}\t+\t{self.k-2}M\tRC:i:{data['weight']}\n")
    
    def write_contigs(self, output_file):
        """
        Write contigs to a FASTA file.
        
        Args:
            output_file (str): Path to the output FASTA file
        """
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with open(output_file, 'w') as f:
            for i, contig in enumerate(self.contigs):
                f.write(f">contig_{i+1} length={len(contig)}\n")
                # Write sequence with 60 characters per line
                for j in range(0, len(contig), 60):
                    f.write(f"{contig[j:j+60]}\n")
    
    def assemble(self, fastq_file, output_file, export_gfa=None, reference_file=None):
        """
        Run the complete assembly pipeline.
        
        Args:
            fastq_file (str): Path to the input FASTQ file
            output_file (str): Path to the output FASTA file
            export_gfa (str, optional): Path to export GFA file
            reference_file (str, optional): Path to reference genome for evaluation
        """
        print(f"Starting assembly with k={self.k}, min_coverage={self.min_coverage}")
        
        # Step 1: Read FASTQ file
        print(f"Reading reads from {fastq_file}")
        reads = self.read_fastq(fastq_file)
        print(f"Read {len(reads)} reads")
        
        # Step 2: Extract k-mers
        print("Extracting k-mers...")
        self.extract_kmers(reads)
        
        # Step 3: Build De Bruijn graph
        print("Building De Bruijn graph...")
        self.build_graph()
        
        # Step 4: Find contigs
        print("Finding contigs...")
        self.find_contigs()
        
        # Step 5: Write contigs to output file
        print(f"Writing {len(self.contigs)} contigs to {output_file}")
        self.write_contigs(output_file)
        
        # Export GFA if requested
        if export_gfa:
            print(f"Exporting graph to GFA format: {export_gfa}")
            self.export_gfa(export_gfa)
        
        # Calculate basic assembly statistics
        assembly_stats = self.calculate_assembly_stats()
        
        # Print assembly statistics
        print("\nAssembly Statistics:")
        print(f"Total assembly length: {assembly_stats['total_length']} bp")
        print(f"Number of contigs: {assembly_stats['num_contigs']}")
        print(f"Largest contig: {assembly_stats['largest_contig']} bp")
        print(f"GC content: {assembly_stats['gc_content']:.2f}%")
        print(f"N50: {assembly_stats['n50']} bp")
        print(f"L50: {assembly_stats['l50']}")
        print(f"N90: {assembly_stats['n90']} bp")
        
        # Try running external QUAST if reference is provided
        if reference_file and os.path.exists(reference_file):
            try:
                self.run_quast(output_file, reference_file)
            except Exception as e:
                print(f"\nNote: Could not run external QUAST: {e}")
        
        return self.contigs
    
    def calculate_assembly_stats(self):
        """
        Calculate basic assembly statistics.
        
        Returns:
            dict: Dictionary with assembly statistics
        """
        contig_lengths = [len(contig) for contig in self.contigs]
        total_length = sum(contig_lengths)
        largest_contig = max(contig_lengths) if contig_lengths else 0
        
        stats = {
            'total_length': total_length,
            'num_contigs': len(self.contigs),
            'largest_contig': largest_contig,
            'gc_content': self.calculate_gc_content(),
            'n50': self.calculate_n50(),
            'l50': self.calculate_l50(),
            'n90': self.calculate_nx(0.9)
        }
        
        return stats
    
    def calculate_gc_content(self):
        """
        Calculate GC content of the assembly.
        
        Returns:
            float: GC content percentage
        """
        if not self.contigs:
            return 0.0
            
        total_gc = sum(contig.count('G') + contig.count('C') for contig in self.contigs)
        total_length = sum(len(contig) for contig in self.contigs)
        
        if total_length == 0:
            return 0.0
            
        return (total_gc / total_length) * 100
    
    def calculate_n50(self):
        """
        Calculate N50 of the assembly.
        
        Returns:
            int: N50 value
        """
        contig_lengths = sorted([len(contig) for contig in self.contigs], reverse=True)
        total_length = sum(contig_lengths)
        cumulative_length = 0
        
        for length in contig_lengths:
            cumulative_length += length
            if cumulative_length >= total_length / 2:
                return length
        
        return 0
    
    def calculate_nx(self, x):
        """
        Calculate Nx (e.g., N90) for the assembly.
        
        Args:
            x (float): Fraction between 0 and 1
            
        Returns:
            int: Nx value
        """
        contig_lengths = sorted([len(contig) for contig in self.contigs], reverse=True)
        total_length = sum(contig_lengths)
        cumulative_length = 0
        
        for length in contig_lengths:
            cumulative_length += length
            if cumulative_length >= total_length * x:
                return length
        
        return 0
    
    def calculate_l50(self):
        """
        Calculate L50 of the assembly.
        
        Returns:
            int: L50 value
        """
        contig_lengths = sorted([len(contig) for contig in self.contigs], reverse=True)
        total_length = sum(contig_lengths)
        cumulative_length = 0
        
        for i, length in enumerate(contig_lengths):
            cumulative_length += length
            if cumulative_length >= total_length / 2:
                return i + 1
        
        return len(contig_lengths)
    
    def run_quast(self, assembly_file, reference_file):
        """
        Run QUAST for comprehensive assembly evaluation.
        
        Args:
            assembly_file (str): Path to assembly FASTA file
            reference_file (str): Path to reference genome FASTA file
        """
        if not reference_file:
            print("Reference file not provided, skipping QUAST analysis.")
            return
        
        # Create output directory
        output_dir = os.path.join(os.path.dirname(assembly_file), "quast_results")
        os.makedirs(output_dir, exist_ok=True)
        
        # Build QUAST command
        quast_cmd = ["quast.py", "-o", output_dir, assembly_file]
        if reference_file:
            quast_cmd.extend(["-r", reference_file])
        
        # Try to run QUAST
        try:
            print(f"Running QUAST: {' '.join(quast_cmd)}")
            result = subprocess.run(quast_cmd, 
                                   capture_output=True, 
                                   text=True, 
                                   check=True)
            
            print("\nQUAST Analysis Completed Successfully")
            print(f"Results available in: {output_dir}")
            
            # Parse and display key QUAST results
            report_file = os.path.join(output_dir, "report.txt")
            if os.path.exists(report_file):
                with open(report_file, 'r') as f:
                    report_content = f.read()
                    print("\nQUAST Report Highlights:")
                    for line in report_content.split('\n'):
                        if any(metric in line for metric in [
                            "Total length", "Genome fraction", "Duplication ratio",
                            "# misassemblies", "# mismatches per 100 kbp",
                            "# indels per 100 kbp", "Largest contig", "GC (%)",
                            "N50", "L50", "N90"
                        ]):
                            print(line)
        except FileNotFoundError:
            print("QUAST not found. Please install QUAST or check if it's in your PATH.")
            raise
        except subprocess.CalledProcessError as e:
            print(f"QUAST execution failed: {e}")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            raise
        except Exception as e:
            print(f"Error running QUAST: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(description='De Bruijn Graph Assembly with QUAST Integration')
    parser.add_argument('--input', required=True, help='Input FASTQ file')
    parser.add_argument('--kmer_size', type=int, required=True, help='Size of k-mers')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    parser.add_argument('--reference', help='Reference genome for evaluation')
    parser.add_argument('--min_coverage', type=int, default=1, help='Minimum k-mer coverage (default: 1)')
    parser.add_argument('--export_gfa', help='Export assembly graph in GFA format')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' does not exist")
        sys.exit(1)
    
    # Check if reference file exists if provided
    if args.reference and not os.path.isfile(args.reference):
        print(f"Error: Reference file '{args.reference}' does not exist")
        sys.exit(1)
    
    # Check if k-mer size is valid
    if args.kmer_size < 2:
        print("Error: k-mer size must be at least 2")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Run assembly
    assembler = DeBruijnGraphAssembler(args.kmer_size, args.min_coverage)
    assembler.assemble(args.input, args.output, args.export_gfa, args.reference)

if __name__ == "__main__":
    main()