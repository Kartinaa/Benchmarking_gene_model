import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolToSmiles, MolToMolBlock
from rdkit.Chem.Scaffolds import MurckoScaffold
from typing import List, Optional, Union

class BMScaffoldMatcher:
    """
    Class to check if the BM scaffold extracted from a molecule is present in a reference set of BM scaffolds.
    """
    def __init__(self, bm_scaffold_file: str):
        """
        Initialize with a BM scaffold .smi file (one SMILES per line, optional title column).
        """
        self.bm_scaffolds = self._load_scaffolds(bm_scaffold_file)
        self.bm_scaffold_smiles = set([Chem.MolToSmiles(mol, isomericSmiles=True) for mol in self.bm_scaffolds if mol])

    def _load_scaffolds(self, file_path: str) -> List[Chem.Mol]:
        scaffolds = []
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                mol = Chem.MolFromSmiles(parts[0])
                if mol:
                    scaffolds.append(mol)
        return scaffolds

    def _read_smi_file(self, smi_file: str) -> List[dict]:
        data = []
        with open(smi_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                entry = {'smiles': parts[0]}
                if len(parts) > 1:
                    entry['title'] = parts[1]
                else:
                    entry['title'] = None
                data.append(entry)
        return data

    def _extract_bm_scaffold(self, mol: Chem.Mol) -> Optional[str]:
        """
        Extract the BM scaffold (Murcko scaffold) as a SMILES string from a molecule.
        Returns None if extraction fails.
        """
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            if scaffold and scaffold.GetNumAtoms() > 0:
                return Chem.MolToSmiles(scaffold, isomericSmiles=True)
        except Exception:
            return None
        return None

    def _match_scaffold(self, mol: Chem.Mol) -> (Optional[str], bool):
        """
        Extract the BM scaffold from the molecule and check if it is in the reference BM scaffolds.
        Returns (bm_scaffold_smiles, occurs)
        """
        bm_scaffold_smiles = self._extract_bm_scaffold(mol) if mol else None
        occurs = bm_scaffold_smiles in self.bm_scaffold_smiles if bm_scaffold_smiles else False
        return bm_scaffold_smiles, occurs

    def process(self, input_data: Union[str, List[Union[str, tuple]]], is_file: bool = True) -> pd.DataFrame:
        """
        Process a .smi file, a list of SMILES/(SMILES, title) tuples, or a single SMILES string.
        Returns a DataFrame with columns: ['smiles', 'mol', 'bm_scaffold', 'occurs', 'title']
        """
        if is_file:
            entries = self._read_smi_file(input_data)
        else:
            # Accept a single SMILES string, a list of SMILES, or a list of (SMILES, title) tuples
            entries = []
            if isinstance(input_data, str):
                entries.append({'smiles': input_data, 'title': None})
            elif isinstance(input_data, (list, tuple)):
                for item in input_data:
                    if isinstance(item, str):
                        entries.append({'smiles': item, 'title': None})
                    elif isinstance(item, (list, tuple)):
                        entries.append({'smiles': item[0], 'title': item[1] if len(item) > 1 else None})

        results = []
        for entry in entries:
            mol = Chem.MolFromSmiles(entry['smiles'])
            bm_scaffold, occurs = self._match_scaffold(mol) if mol else (None, False)
            results.append({
                'smiles': entry['smiles'],
                'mol': mol,
                'bm_scaffold': bm_scaffold,
                'occurs': int(occurs),
                'title': entry['title']
            })
        df = pd.DataFrame(results)
        return df

    @staticmethod
    def visualize(df: pd.DataFrame, subset: Optional[list] = None, **kwargs):
        """
        Visualize a DataFrame of molecules using mols2grid.
        Parameters:
            df: DataFrame containing a 'mol' column (RDKit Mol objects).
            subset: list of columns to display (default ['img', 'Count'] if present, else ['img']).
            kwargs: additional arguments passed to mols2grid.display.
        """
        try:
            import mols2grid
        except ImportError:
            raise ImportError("mols2grid is not installed. Please install it with 'pip install mols2grid'.")
        # Generate 'img' column if not present
        if 'img' not in df.columns:
            from rdkit.Chem import Draw
            df = df.copy()
            df['img'] = df['mol'].apply(lambda m: Draw.MolToImage(m) if m else None)
        # Default subset
        if subset is None:
            subset = ['img', 'occurs'] if 'occurs' in df.columns else ['img']
        return mols2grid.display(df, subset=subset, **kwargs)
