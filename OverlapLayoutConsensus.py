#!/usr/bin/env python3
"""
Overlap-Layout-Consensus (OLC) Genome Assembler with QUAST Integration

Usage:
    python olc_assembler.py --input <input.fastq> --min_overlap <n> --output <output.fasta> 
                           [--reference <reference.fasta>] [--min_coverage <min_cov>] [--export_gfa <graph.gfa>]

Arguments:
    --input: Path to the input FASTQ file
    --min_overlap: Minimum overlap length between reads
    --output: Path to the output FASTA file for contigs
    --reference: Path to the reference genome for evaluation (optional)
    --min_coverage: Minimum read coverage for consensus (default: 1)
    --export_gfa: Export assembly graph in GFA format (optional)
"""

import argparse
import os
import sys
from collections import defaultdict
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from itertools import combinations
import subprocess
import tempfile
import re

class OLCAssembler:
    def __init__(self, min_overlap, min_coverage=1):
        """
        Initialize the OLC assembler.
        
        Args:
            min_overlap (int): Minimum overlap length between reads
            min_coverage (int): Minimum read coverage for consensus
        """
        self.min_overlap = min_overlap
        self.min_coverage = min_coverage
        self.reads = []
        self.read_ids = {}
        self.overlap_graph = nx.DiGraph()
        self.contigs = []
    
    def read_fastq(self, fastq_file):
        """
        Read sequences from a FASTQ file.
        
        Args:
            fastq_file (str): Path to the FASTQ file
        """
        reads = []
        read_ids = {}
        
        for i, record in enumerate(SeqIO.parse(fastq_file, "fastq")):
            read_seq = str(record.seq)
            read_id = f"read_{i}"
            reads.append(read_seq)
            read_ids[read_id] = i
        
        self.reads = reads
        self.read_ids = read_ids
        return reads
    
    def find_overlaps(self):
        """
        Find all-vs-all read overlaps to build the overlap graph.
        """
        print("Finding overlaps between reads...")
        
        # For efficiency with large datasets, we could implement a k-mer based pre-filter
        # to only compare reads that are likely to overlap
        
        total_reads = len(self.reads)
        processed = 0
        
        # For each pair of reads
        for i in range(len(self.reads)):
            # Print progress every 100 reads
            processed += 1
            if processed % 100 == 0 or processed == total_reads:
                print(f"Progress: {processed}/{total_reads} reads processed ({processed/total_reads:.1%})")
                
            for j in range(len(self.reads)):
                if i == j:  # Skip self-comparison
                    continue
                
                read1 = self.reads[i]
                read2 = self.reads[j]
                
                # Find suffix-prefix overlap
                max_overlap = 0
                for overlap_len in range(min(len(read1), len(read2)), self.min_overlap - 1, -1):
                    # Check if suffix of read1 matches prefix of read2
                    if read1[-overlap_len:] == read2[:overlap_len]:
                        max_overlap = overlap_len
                        break
                
                # If significant overlap found, add edge to graph
                if max_overlap >= self.min_overlap:
                    read1_id = f"read_{i}"
                    read2_id = f"read_{j}"
                    
                    # Add nodes if they don't exist
                    if read1_id not in self.overlap_graph:
                        self.overlap_graph.add_node(read1_id, sequence=read1)
                    if read2_id not in self.overlap_graph:
                        self.overlap_graph.add_node(read2_id, sequence=read2)
                    
                    # Add directed edge with overlap length
                    self.overlap_graph.add_edge(read1_id, read2_id, overlap=max_overlap)
        
        print(f"Built overlap graph with {self.overlap_graph.number_of_nodes()} nodes and {self.overlap_graph.number_of_edges()} edges")
    
    def find_paths(self):
        """
        Find non-branching paths in the overlap graph to form contigs.
        """
        print("Finding non-branching paths...")
        
        # Identify start nodes (in-degree = 0 or junction nodes)
        start_nodes = [node for node in self.overlap_graph.nodes() 
                      if self.overlap_graph.in_degree(node) == 0 or 
                      self.overlap_graph.in_degree(node) > 1 or 
                      self.overlap_graph.out_degree(node) > 1]
        
        # If no start nodes (circular graph), use any node
        if not start_nodes and self.overlap_graph.nodes():
            start_nodes = [list(self.overlap_graph.nodes())[0]]
        
        paths = []
        visited_edges = set()
        
        # Process each start node
        for start_node in start_nodes:
            if self.overlap_graph.out_degree(start_node) == 0:
                continue
                
            # Start path from each outgoing edge
            for _, end_node in self.overlap_graph.out_edges(start_node):
                edge = (start_node, end_node)
                if edge in visited_edges:
                    continue
                
                # Start a new path
                path = [start_node]
                current_node = end_node
                
                while True:
                    path.append(current_node)
                    visited_edges.add((path[-2], current_node))
                    
                    # If we reached a junction or end
                    if self.overlap_graph.out_degree(current_node) != 1 or self.overlap_graph.in_degree(current_node) != 1:
                        break
                    
                    # Continue path
                    next_nodes = list(self.overlap_graph.successors(current_node))
                    if not next_nodes:
                        break
                    next_node = next_nodes[0]
                    
                    # Check if edge already visited
                    if (current_node, next_node) in visited_edges:
                        break
                    
                    current_node = next_node
                
                # If valid path found, save it
                if len(path) > 1:
                    paths.append(path)
        
        # Check for isolated cycles
        components = list(nx.weakly_connected_components(self.overlap_graph))
        for component in components:
            subgraph = self.overlap_graph.subgraph(component)
            # If all nodes have in-degree = out-degree = 1, it's a cycle
            if all(subgraph.in_degree(node) == subgraph.out_degree(node) == 1 for node in component):
                # Find a cycle
                cycle = find_cycle(subgraph)
                if cycle and len(cycle) > 0:
                    paths.append(cycle)
        
        print(f"Found {len(paths)} non-branching paths")
        return paths
    
    def generate_consensus(self, path):
        """
        Generate consensus sequence for a path of reads.
        
        Args:
            path (list): List of read IDs forming a contig
            
        Returns:
            str: Consensus sequence
        """
        if not path:
            return ""
        
        # Start with the first read
        first_node = path[0]
        consensus = self.overlap_graph.nodes[first_node]['sequence']
        
        # For each subsequent read in the path
        for i in range(1, len(path)):
            prev_node = path[i-1]
            curr_node = path[i]
            
            # Get overlap length between reads
            overlap_len = self.overlap_graph[prev_node][curr_node]['overlap']
            
            # Get the next read sequence
            next_seq = self.overlap_graph.nodes[curr_node]['sequence']
            
            # Append the non-overlapping part to the consensus
            consensus += next_seq[overlap_len:]
        
        return consensus
    
    def assemble(self, fastq_file, output_file, export_gfa=None, reference_file=None):
        """
        Run the complete assembly pipeline.
        
        Args:
            fastq_file (str): Path to the input FASTQ file
            output_file (str): Path to the output FASTA file
            export_gfa (str, optional): Path to export GFA file
            reference_file (str, optional): Path to reference genome for evaluation
        """
        print(f"Starting assembly with min_overlap={self.min_overlap}, min_coverage={self.min_coverage}")
        
        # Step 1: Read FASTQ file
        print(f"Reading reads from {fastq_file}")
        self.read_fastq(fastq_file)
        print(f"Read {len(self.reads)} reads")
        
        # Step 2: Find overlaps and build graph
        self.find_overlaps()
        
        # Step 3: Find non-branching paths
        paths = self.find_paths()
        
        # Step 4: Generate consensus for each path
        print("Generating consensus sequences...")
        contigs = []
        for i, path in enumerate(paths):
            contig = self.generate_consensus(path)
            contigs.append(contig)
            print(f"Generated contig {i+1}: {len(contig)} bp from {len(path)} reads")
        
        # Sort contigs by length (descending)
        contigs = sorted(contigs, key=len, reverse=True)
        self.contigs = contigs
        
        # Step 5: Write contigs to output file
        self.write_contigs(output_file)
        
        # Export GFA if requested
        if export_gfa:
            print(f"Exporting graph to GFA format: {export_gfa}")
            self.export_gfa(export_gfa)
        
        # Calculate basic assembly statistics
        assembly_stats = self.calculate_assembly_stats()
        
        # Print assembly statistics (moved to quest instead)
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
        
        return contigs
    
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
        
        print(f"Wrote {len(self.contigs)} contigs to {output_file}")
    
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
            for node in self.overlap_graph.nodes():
                sequence = self.overlap_graph.nodes[node]['sequence']
                f.write(f"S\t{node}\t{sequence}\n")
            
            # Write edges (links) with more detailed information
            for u, v, data in self.overlap_graph.edges(data=True):
                overlap_len = data['overlap']
                # Add more information to the edge representation
                f.write(f"L\t{u}\t+\t{v}\t+\t{overlap_len}M\tOL:i:{overlap_len}\n")
        
        print(f"Exported overlap graph to {gfa_file}")
    
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

def find_cycle(graph):
    """Find a cycle in the graph."""
    try:
        return list(nx.simple_cycles(graph))[0]
    except:
        return None

def main():
    parser = argparse.ArgumentParser(description='Overlap-Layout-Consensus Assembly with QUAST Integration')
    parser.add_argument('--input', required=True, help='Input FASTQ file')
    parser.add_argument('--min_overlap', type=int, required=True, help='Minimum overlap length')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    parser.add_argument('--reference', help='Reference genome for evaluation')
    parser.add_argument('--min_coverage', type=int, default=1, help='Minimum read coverage (default: 1)')
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
    
    # Check if minimum overlap is valid
    if args.min_overlap < 1:
        print("Error: Minimum overlap must be at least 1")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Run assembly
    assembler = OLCAssembler(args.min_overlap, args.min_coverage)
    assembler.assemble(args.input, args.output, args.export_gfa, args.reference)

if __name__ == "__main__":
    main()