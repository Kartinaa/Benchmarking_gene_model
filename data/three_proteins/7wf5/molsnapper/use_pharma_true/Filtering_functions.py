from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def atom_in_multiple_small_rings(mol):
    """
    Check if any atom is in two or more small (3- or 4-membered) rings.
    
    Args:
        mol: RDKit molecule object.
    
    Returns:
        bool: True if any atom is in two or more small rings, False otherwise.
    """
    # Get all ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Get the atom indices for each ring

    # Create a dictionary to count how many rings each atom is part of
    atom_ring_count = {}
    for ring in atom_rings:
        # Only count small rings (3- or 4-membered)
        if len(ring) <= 4:
            for atom_idx in ring:
                if atom_idx in atom_ring_count:
                    atom_ring_count[atom_idx] += 1
                else:
                    atom_ring_count[atom_idx] = 1

    # Check if any atom is in more than one small ring
    for count in atom_ring_count.values():
        if count > 1:
            return True  # Atom is in multiple small rings
    
    return False  # No atom is in multiple small rings

def double_bond_in_small_ring(mol):
    """
    Check if any small (3- or 4-membered) ring contains a double bond.
    
    Args:
        mol: RDKit molecule object.
    
    Returns:
        bool: True if a double bond exists inside a small ring, False otherwise.
    """
    # Get all ring information
    ring_info = mol.GetRingInfo()
    bond_rings = ring_info.BondRings()  # Get the bond indices for each ring

    # Iterate through the rings
    for ring in bond_rings:
        if len(ring) <= 4:  # Only small rings (3- or 4-membered)
            for bond_idx in ring:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    return True  # Double bond found in a small ring
    
    return False  # No double bonds in small rings

def filter_molecules(mol):
    """
    Filter molecules based on the two rules:
    1. No atom is part of two small rings.
    2. No double bond inside small rings.
    
    Args:
        molecules: List of RDKit molecule objects.
    
    Returns:
        List of filtered molecule objects.
    """
    if not atom_in_multiple_small_rings(mol) and not double_bond_in_small_ring(mol):
        return True
    else:
        return False
    

