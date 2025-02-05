from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from difflib import SequenceMatcher
import numpy as np

print("The program will predict possible amplicons generated using the input primers, against a supplied genome sequence.")

def find_primer_matches_with_mismatches(reference_seq_str, primer_seq, max_mismatches=2):
    matches = []
    for i in range(len(reference_seq_str) - len(primer_seq) + 1):
        window = reference_seq_str[i:i + len(primer_seq)]
        s = SequenceMatcher(None, primer_seq, window)
        ratio = s.ratio()
        if ratio >= (1 - (max_mismatches / len(primer_seq))):
            matches.append((i, i + len(primer_seq)))
    return matches
    
def predict_amplicons(reference_fasta, primer_csv, min_amplicon_size=100, max_amplicon_size=500, max_mismatches=3):
    reference_seq = SeqIO.read(reference_fasta, "fasta")
    reference_seq_str = str(reference_seq.seq)

    primers_df = pd.read_csv(primer_csv)
    primers = [(row['primer_label'], row['primer_sequence']) for _, row in primers_df.iterrows()]

    primer_matches = {}
    for primer_label, primer_seq in primers:
        forward_matches = find_primer_matches_with_mismatches(reference_seq_str, primer_seq, max_mismatches)
        reverse_matches = find_primer_matches_with_mismatches(reference_seq_str, str(Seq(primer_seq).reverse_complement()), max_mismatches)
        primer_matches[primer_label] = {'forward': forward_matches, 'reverse': reverse_matches}

    amplicons = []
    for primer_label1, matches1 in primer_matches.items():
        for primer_label2, matches2 in primer_matches.items():
            if primer_label1 == primer_label2:
                continue
            for forward_match1 in matches1['forward']:
                for reverse_match2 in matches2['reverse']:
                    if forward_match1[1] < reverse_match2[0]:
                        amplicon_size = reverse_match2[0] - forward_match1[1]
                        if min_amplicon_size <= amplicon_size <= max_amplicon_size:
                            amplicon = {
                                'name': f'{primer_label1}-{primer_label2}',
                                'start': forward_match1[1],
                                'end': reverse_match2[0],
                                'primers': [primer_label1, primer_label2]
                            }
                            amplicons.append(amplicon)

    G = nx.Graph()

    for primer_label, _ in primers:
        G.add_node(primer_label, type='primer')
    for amplicon in amplicons:
        G.add_node(amplicon['name'], type='amplicon', start=amplicon['start'], end=amplicon['end'])

    for amplicon in amplicons:
        for primer_label in amplicon['primers']:
            G.add_edge(primer_label, amplicon['name'])

    return G, reference_seq_str, amplicons

def calculate_coverage_without_overlap(amplicon_positions):
    amplicon_positions.sort(key=lambda x: x[1])

    covered_regions = []
    for _, start, end in amplicon_positions:
        new_region = (start, end)
        merged = False
        for i in range(len(covered_regions)):
            existing_start, existing_end = covered_regions[i]
            overlap_start = max(start, existing_start)
            overlap_end = min(end, existing_end)

            if overlap_start < overlap_end:
                merged = True
                covered_regions[i] = (min(start, existing_start), max(end, existing_end))
                break
        if not merged:
            covered_regions.append(new_region)

    total_covered_length = sum(end - start for start, end in covered_regions)
    max_position = max(end for _, _, end in amplicon_positions) if amplicon_positions else 0
    coverage_percentage = (total_covered_length / max_position) * 100 if max_position>0 else 0

    return coverage_percentage


reference_fasta = input("Path to Reference Genome file.fasta: ")
primer_csv = input("Path to primer sequence file.csv (Must contain two headers- primer_label,primer_sequence): ")

amplicon_graph, reference_seq_str, amplicons = predict_amplicons(reference_fasta, primer_csv, max_mismatches=3)

# Output a Visualisation Graph
plt.figure(figsize=(12, 6))

amplicon_positions = [(n, d['start'], d['end']) for n, d in amplicon_graph.nodes(data=True) if d['type'] == 'amplicon']
amplicon_positions.sort(key=lambda x: x[1])
y_coords = range(len(amplicon_positions))

for i, (name, start, end) in enumerate(amplicon_positions):
    plt.plot([start, end], [i, i], color='blue', linewidth=2)

coverage = calculate_coverage_without_overlap(amplicon_positions)

plt.xlabel('Reference Sequence Position')
plt.ylabel('Amplicon')
plt.yticks(y_coords, [amp[0] for amp in amplicon_positions])
plt.title('Predicted Amplicons across the genome')

min_coord = min(amp[1] for amp in amplicon_positions)
max_coord = max(amp[2] for amp in amplicon_positions)
plt.xlim(min_coord - 10, max_coord + 10)

plt.show()
print("Predicted Coverage Percentage for input genome sequence is approximately ",round(coverage,2),"%.")
