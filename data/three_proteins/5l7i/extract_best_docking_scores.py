#!/usr/bin/env python3
"""
Script to extract the highest docking score for each ligand from PDBQT files.
Each ligand may have multiple conformers (e.g., 9_prepared-1_out.pdbqt, 9_prepared-2_out.pdbqt).
This script finds the best score for each unique ligand.
"""

import os
import re
import glob
from collections import defaultdict

def extract_ligand_number(filename):
    """Extract the ligand number from filename like '9_prepared-1_out.pdbqt' -> '9'"""
    match = re.match(r'ligand(\d+)_prepared-\d+_out\.pdbqt', filename)
    if match:
        return match.group(1)
    return None

def extract_docking_score(filepath):
    """Extract all docking scores from a PDBQT file and return the best one"""
    scores = []
    vina_result_re = re.compile(r'^REMARK VINA RESULT:\s*([-\d\.]+)')
    try:
        with open(filepath, 'r') as f:
            for line in f:
                match = vina_result_re.match(line)
                if match:
                    try:
                        score = float(match.group(1))
                        scores.append(score)
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None
    if scores:
        return min(scores)  # Lower (more negative) scores are better
    return None

def main():
    # Directory containing the PDBQT files
    docking_dir = "/home/yang2531/Documents/Bo_toolbox/PatWalters/Benchmarking_gene_model/data/three_proteins/7wf5/7wf5_cpu_result_top500"
    
    # Dictionary to store the best score for each ligand
    ligand_scores = defaultdict(list)
    
    # Get all PDBQT files
    pdbqt_files = glob.glob(os.path.join(docking_dir, "*_out.pdbqt"))
    
    print(f"Found {len(pdbqt_files)} PDBQT files")
    
    # Process each file
    for filepath in pdbqt_files:
        filename = os.path.basename(filepath)
        ligand_number = extract_ligand_number(filename)
        
        if ligand_number is None:
            print(f"Warning: Could not extract ligand number from {filename}")
            continue
        
        score = extract_docking_score(filepath)
        if score is not None:
            ligand_scores[ligand_number].append((score, filename))
            print(f"Ligand {ligand_number}: {filename} -> Score: {score:.3f}")
        else:
            print(f"Warning: No valid score found in {filename}")
    
    # Find the best score for each ligand
    best_scores = {}
    for ligand_number, scores_list in ligand_scores.items():
        if scores_list:
            # Sort by score (lower is better) and take the first
            best_score, best_file = min(scores_list, key=lambda x: x[0])
            best_scores[ligand_number] = (best_score, best_file)
    
    # Sort ligands by number and print results
    print("\n" + "="*60)
    print("BEST DOCKING SCORES FOR EACH LIGAND")
    print("="*60)
    print(f"{'Ligand':<8} {'Best Score':<12} {'File':<25} {'Conformers':<12}")
    print("-" * 60)
    
    for ligand_number in sorted(best_scores.keys(), key=int):
        best_score, best_file = best_scores[ligand_number]
        num_conformers = len(ligand_scores[ligand_number])
        print(f"{ligand_number:<8} {best_score:<12.3f} {best_file:<25} {num_conformers:<12}")
    
    # Save results to a CSV file
    output_file = "best_vina_docking_scores_top500.csv"
    if best_scores:
        scores_only = [score for score, _ in best_scores.values()]
        avg_score = sum(scores_only) / len(scores_only)
        best_score_val = min(scores_only)
        worst_score_val = max(scores_only)
    with open(output_file, 'w') as f:
        f.write("Ligand,Best_Score,Best_File,Total_Conformers,Average_Score,Best_Overall_Score,Worst_Overall_Score\n")
        first_row = True
        for ligand_number in sorted(best_scores.keys(), key=int):
            best_score, best_file = best_scores[ligand_number]
            num_conformers = len(ligand_scores[ligand_number])
            if first_row:
                f.write(f"{ligand_number},{best_score:.3f},{best_file},{num_conformers},{avg_score:.3f},{best_score_val:.3f},{worst_score_val:.3f}\n")
                first_row = False
            else:
                f.write(f"{ligand_number},{best_score:.3f},{best_file},{num_conformers},,,\n")
    print(f"\nResults and summary statistics saved to {output_file}")
    
    # Save summary statistics to a text file
    # summary_file = "docking_score_summary.txt"
    # if best_scores:
    #     scores_only = [score for score, _ in best_scores.values()]
    #     avg_score = sum(scores_only) / len(scores_only)
    #     best_score_val = min(scores_only)
    #     worst_score_val = max(scores_only)
    #     with open(summary_file, 'w') as f:
    #         f.write(f"Total ligands processed: {len(best_scores)}\n")
    #         f.write(f"Best overall score: {best_score_val:.3f}\n")
    #         f.write(f"Worst overall score: {worst_score_val:.3f}\n")
    #         f.write(f"Average score: {avg_score:.3f}\n")
    #     print(f"Summary statistics saved to {summary_file}")

if __name__ == "__main__":
    main() 