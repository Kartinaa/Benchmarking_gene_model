import os
import gzip
import pickle
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from joblib import Parallel, delayed
from typing import Set

class BMScaffoldChecker:
    
    def __init__(self, bm_scaffold_file: str, n_jobs: int = 8):
        self.bm_scaffold_smiles_set = self.load_scaffold_smiles_parallel(bm_scaffold_file, n_jobs=n_jobs)

    @staticmethod
    def load_scaffold_smiles_parallel(file_path: str, n_jobs: int = 8) -> Set[str]:
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        def parse_block(block_lines):
            return [line.split()[0] for line in block_lines if line]

        # Split into N chunks
        chunk_size = len(lines) // n_jobs + 1
        chunks = [lines[i:i+chunk_size] for i in range(0, len(lines), chunk_size)]

        parsed_chunks = Parallel(n_jobs=n_jobs)(
            delayed(parse_block)(chunk) for chunk in chunks
        )

        # Flatten and deduplicate
        return set(smiles for chunk in parsed_chunks for smiles in chunk)

    def _extract_bm_scaffold_smiles(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if not scaffold or scaffold.GetNumAtoms() == 0:
            return None
        return Chem.MolToSmiles(scaffold, isomericSmiles=True)

    def check_smiles(self, smiles: str) -> bool:
        """
        Return True if the BM scaffold of the given SMILES exists in the source scaffold set.
        """
        bm_smi = self._extract_bm_scaffold_smiles(smiles)
        return bm_smi in self.bm_scaffold_smiles_set if bm_smi else False

    def check_file(self, smi_file: str) -> float:
        """
        Return the % of BM scaffolds in a .smi file that occur in the source scaffold set.
        Only valid molecules with extractable scaffolds are considered.
        """
        total = 0
        matched = 0
        with open(smi_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                smiles = parts[0]
                bm_smi = self._extract_bm_scaffold_smiles(smiles)
                if bm_smi:
                    total += 1
                    if bm_smi in self.bm_scaffold_smiles_set:
                        matched += 1
        return (matched / total * 100) if total > 0 else 0.0, matched, total
