# %%
from rdkit import Chem
import sys
import os
import useful_rdkit_utils as uru
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem.Scaffolds import MurckoScaffold
import itertools
import pandas as pd
from collections import Counter
import mols2grid
from tdc.single_pred import ADME
from Filtering_functions import filter_molecules
import glob
from rdkit.Chem.Scaffolds import MurckoScaffold
import molvs as mv
from tqdm import tqdm
import concurrent.futures

# %%
# Get the current working directory
current_dir = os.getcwd()
print(current_dir)
# Get the parent directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
print(parent_dir)
# Add the parent directory to sys.path
sys.path.append(parent_dir)

# %% [markdown]
# Set the default image size

# %%
uru.rd_set_image_size(300,300)

# %% [markdown]
# ### Reading the data from Targetdiff

# %% [markdown]
# See if they are valid or not

# %%
smi_list = []
with open('/home/yang2531/Documents/Bo_toolbox/PatWalters/Benchmarking_gene_model/Targetdiff/combined_smiles_target_diff.txt') as f:
    for smi in f:
        smi_list.append(smi.strip())
standardized_smiles_list = []
for smi in smi_list:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        standardized_smiles_list.append(Chem.MolToSmiles(mol))

print(len(smi_list))
len(standardized_smiles_list)

# %% [markdown]
# Convert SMILES and name of them to a dataframe

# %%
df = pd.DataFrame(standardized_smiles_list, columns=['SMILES'])
df

# %% [markdown]
# Add an RDKit molecule to the dataframe

# %%
df['mol'] = df.SMILES.apply(Chem.MolFromSmiles)

# %% [markdown]
# Strip salts from the molecule

# %%
df.mol = df.mol.apply(uru.get_largest_fragment)
df.mol[0]

# %% [markdown]
# Instantiate a RingSystemFinder object and find ring systems for the molecules in df

# %%
ring_system_finder = uru.RingSystemFinder()
df['ring_sytems'] = df.mol.apply(ring_system_finder.find_ring_systems)
df.ring_sytems

# %%
df.ring_sytems.values

# %% [markdown]
# See how many times each ring system occurs

# %%
ring_system_list = list(itertools.chain.from_iterable(df.ring_sytems.values))
ring_count_df = pd.DataFrame(Counter(ring_system_list).items(),columns=["SMILES","Count"]) ### Convert a dictionary to a DataFrame.
ring_count_df.sort_values("Count",ascending=False,inplace=True)
ring_count_df

# %% [markdown]
# View the ring system frequencies

# %%
mols2grid.display(ring_count_df,subset=["img","Count"])

# %% [markdown]
# The RingSystemLookup object has a dictionary of how many times each ring system occurs in the ChEMBL database.  We can use this object to evaluate the molecules in df.

# %%
ring_system_lookup = uru.RingSystemLookup.default()
res = df.mol.apply(ring_system_lookup.process_mol)

# %%
res

# %%
df[['min_ring','min_freq']] = res.apply(uru.get_min_ring_frequency).tolist()
df

# %%
mols2grid.display(df.sort_values("min_freq"),mol_col="mol",subset=["img","min_freq"])

# %% [markdown]
# ### Check how many molecules are unique the ring frequency of them

# %% [markdown]
# How many molecules are unique using inChI?

# %%
df['inchi'] = df.mol.apply(Chem.MolToInchi)
df = df.drop_duplicates("inchi", keep="first", ignore_index=True)
df

# %% [markdown]
# How many molecules contains ring structure?

# %%
df_ring = df[df.min_freq != -1]
df_ring

# %%
filtered_df_ring_freq = df[(df.min_freq > 100) & (df.min_freq != -1)]
filtered_df_ring_freq

# %% [markdown]
# ### How many of them could pass PAINS filter?

# %% [markdown]
# Get a list of rules

# %%
reos = uru.REOS()
reos.get_available_rule_sets()

# %% [markdown]
# Get the currently active rule sets

# %%
reos.get_active_rule_sets()

# %% [markdown]
# Set active rule set to PAINS

# %%
reos.set_active_rule_sets(['PAINS'])
reos.get_active_rule_sets()

# %% [markdown]
# Apply PAINS filter to unique mols

# %%
reos.pandas_mols(df.mol)

# %%
df_PAINS_filter = pd.concat([df, reos.pandas_mols(df.mol)], axis=1)
df_PAINS_filter = df_PAINS_filter[df_PAINS_filter.description == 'ok']
df_PAINS_filter

# %% [markdown]
# ### How many of them could pass filters suggested by Dr.Reymond?

# %%
df['Reymond'] = df.mol.apply(filter_molecules)
df

# %%
df_Reymond_filter = df[df.Reymond == True]
df_Reymond_filter

# %% [markdown]
# ### How many of them has BM scaffold appears in ZINC20 druglike molecules?

# %%
with open('/home/yang2531/Documents/Bo_toolbox/PatWalters/Benchmarking_gene_model/combined_BM_scaffolds_ZINC20_unique.smi', 'r') as f:
    bm_scaffolds = {line.strip() for line in f}

largest_frag_chooser = mv.fragment.LargestFragmentChooser()

def contains_bm_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        bm_scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        largest_frag = largest_frag_chooser.choose(bm_scaffold)
    for scaffold in bm_scaffolds:
        scaffold_mol = Chem.MolFromSmiles(scaffold)
        if scaffold_mol is None:
            continue
        if largest_frag.HasSubstructMatch(scaffold_mol) and scaffold_mol.HasSubstructMatch(largest_frag):
            return True
    return False

# %%
with concurrent.futures.ProcessPoolExecutor(max_workers=28) as executor:
    results = list(tqdm(executor.map(contains_bm_scaffold, df['SMILES']), total=len(df)))

# Add results to DataFrame
df['contains_bm_scaffold'] = results

# Filter the DataFrame
df_bm_zinc20 = df[df['contains_bm_scaffold']]


