import os
import shutil

def get_unique_filename(dst_folder, filename):
    """
    Generates a unique filename by appending a counter if the file already exists in the destination folder.
    
    Args:
        dst_folder (str): The destination folder.
        filename (str): The original filename.
    
    Returns:
        str: A unique filename.
    """
    base, extension = os.path.splitext(filename)
    counter = 1
    new_filename = filename
    
    # Loop until we find a filename that does not exist
    while os.path.exists(os.path.join(dst_folder, new_filename)):
        new_filename = f"{base}_{counter}{extension}"
        counter += 1
        
    return new_filename

def gather_shifted_sdfs(src_root, dst_folder):
    """
    Gather all *_shifted.sdf files from subfolders and copy them to a single destination folder.
    Renames files with a counter if they already exist.
    
    Args:
        src_root (str): The root folder containing subfolders.
        dst_folder (str): The folder where all the *_shifted.sdf files will be copied.
    """
    if not os.path.exists(dst_folder):
        os.makedirs(dst_folder)

    # Walk through all directories and subdirectories
    for root, dirs, files in os.walk(src_root):
        for file in files:
            if file.endswith('_shifted.sdf'):
                # Get the unique filename to avoid overwriting
                unique_filename = get_unique_filename(dst_folder, file)
                
                # Build the full path to the source file
                full_file_path = os.path.join(root, file)
                
                # Copy the file to the destination folder with the unique filename
                shutil.copy(full_file_path, os.path.join(dst_folder, unique_filename))
                print(f"Copied: {full_file_path} to {os.path.join(dst_folder, unique_filename)}")

if __name__ == "__main__":
    src_root = '/home/yang2531/Documents/Bo_toolbox/PatWalters/Benchmarking_gene_model/Molsnapper/sample_MolDiff_20241002_214910_clash_rate_0.1_SDF'  # Change this to your parent folder
    dst_folder = '/home/yang2531/Documents/Bo_toolbox/PatWalters/Benchmarking_gene_model/Molsnapper/sample_MolDiff_20241002_214910_clash_rate_0.1_SDF/combined_shifted_sdf'  # Change this to your destination folder
    
    gather_shifted_sdfs(src_root, dst_folder)
