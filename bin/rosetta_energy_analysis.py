#!/usr/bin/env python3
"""
Rosetta Energy Analysis for AlphaFold3 CIF Files

This script performs comprehensive energy analysis on protein structures using PyRosetta.
It processes CIF files from AlphaFold3 and calculates various energy terms and structural metrics.

Usage:
    python rosetta_energy_analysis.py --input_dir /path/to/cif/files --output results.csv

Requirements:
    - PyRosetta installed and properly configured
    - Python 3.6+
"""

import os
import sys
import argparse
import re
import polars as pl
from pathlib import Path
import logging
import tempfile
from Bio.PDB import MMCIF2Dict, PDBParser, PDBIO
from Bio.PDB.mmcifio import MMCIFIO

try:
    import pyrosetta
    from pyrosetta import pose_from_pdb
    from pyrosetta.rosetta.core.scoring import CA_rmsd, all_atom_rmsd
    from pyrosetta.rosetta.core.scoring import calc_total_sasa
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
except ImportError as e:
    print(f"Error importing PyRosetta: {e}")
    print("Please ensure PyRosetta is properly installed and configured.")
    sys.exit(1)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def initialize_pyrosetta():
    """Initialize PyRosetta with standard options."""
    try:
        pyrosetta.init("-ignore_unrecognized_res -ignore_zero_occupancy false -obey_ENDMDL false")
        logger.info("PyRosetta initialized successfully")
    except Exception as e:
        logger.error(f"Failed to initialize PyRosetta: {e}")
        sys.exit(1)


def convert_cif_to_pdb(cif_path, temp_dir=None):
    """
    Convert a CIF file to PDB format.
    
    Args:
        cif_path (str): Path to the CIF file
        temp_dir (str): Directory for temporary PDB file (optional)
        
    Returns:
        str: Path to the converted PDB file
    """
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        
        # Parse CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('structure', str(cif_path))
        
        # Create temporary PDB file
        if temp_dir is None:
            temp_dir = tempfile.gettempdir()
        
        pdb_filename = Path(cif_path).stem + '_converted.pdb'
        pdb_path = Path(temp_dir) / pdb_filename
        
        # Write to PDB format
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_path))
        
        logger.info(f"Converted {cif_path} to {pdb_path}")
        return str(pdb_path)
        
    except Exception as e:
        logger.error(f"Failed to convert CIF to PDB {cif_path}: {e}")
        return None


def load_pose_from_cif(cif_path, temp_dir=None):
    """
    Load a protein structure from a CIF file by first converting to PDB format.
    
    Args:
        cif_path (str): Path to the CIF file
        temp_dir (str): Directory for temporary files (optional)
        
    Returns:
        pyrosetta.Pose: Loaded pose object
    """
    try:
        # First convert CIF to PDB
        pdb_path = convert_cif_to_pdb(cif_path, temp_dir)
        if pdb_path is None:
            logger.error(f"Failed to convert CIF to PDB: {cif_path}")
            return None
        
        # Load pose from the converted PDB file
        pose = pose_from_pdb(pdb_path)
        logger.info(f"Successfully loaded pose from converted PDB: {pdb_path}")
        
        # Clean up temporary PDB file
        try:
            os.remove(pdb_path)
            logger.debug(f"Cleaned up temporary file: {pdb_path}")
        except Exception as cleanup_e:
            logger.warning(f"Could not clean up temporary file {pdb_path}: {cleanup_e}")
        
        return pose
        
    except Exception as e:
        logger.error(f"Failed to load pose from {cif_path}: {e}")
        return None


def interface_ddg_IAM(pose, partners_str, scorefxn=None,
                      pack_separated=True, pack_input=True, verbose=False):
    """
    Compute interface ΔΔG using InterfaceAnalyzerMover.

    Args:
        pose: PyRosetta Pose
        partners_str (str): Docking-style partners string, e.g. "A_BC"
        scorefxn: optional ScoreFunction (defaults to ref2015 if None)
        pack_separated (bool): repack each partner in separated state
        pack_input (bool): repack the complex/input state
        verbose (bool): IAM verbose output

    Returns:
        float: interface ΔΔG (REU)
    """
    if scorefxn is None:
        scorefxn = pyrosetta.create_score_function("ref2015")

    iam = InterfaceAnalyzerMover(partners_str)
    iam.set_scorefunction(scorefxn)

    # Common, robust settings for meaningful ΔΔG:
    # - Repack in the bound and separated states to let sidechains optimize
    # - Keep defaults for minimization unless you want to enable relax
    iam.set_pack_separated(pack_separated)
    iam.set_pack_input(pack_input)

    iam.apply(pose)
    return iam.get_interface_dG()


def collect_energy_info(pose, reference_pose=None, interface_mode=None):
    """
    Collects comprehensive energy terms and structural metrics from a pose.
    
    Args:
        pose: PyRosetta pose object
        reference_pose: Reference pose for RMSD calculation (optional)
        interface_mode: Interface analysis mode - "peptide_antibody" for 3-chain system (A vs B+C)
                       or tuple of chain IDs for 2-chain analysis (e.g., ("A", "B"))
        
    Returns:
        dict: Dictionary containing energy metrics and structural information
    """
    try:
        # Create score function using ref2015 scoring scheme
        scorefxn = pyrosetta.create_score_function("ref2015")
        
        # Compute the total energy for the whole pose
        total_energy = scorefxn(pose)
        
        # Collect key energy terms for the whole structure
        energy_data = pose.energies().total_energies()
        fa_atr = energy_data[pyrosetta.rosetta.core.scoring.fa_atr]
        fa_rep = energy_data[pyrosetta.rosetta.core.scoring.fa_rep]
        fa_sol = energy_data[pyrosetta.rosetta.core.scoring.fa_sol]
        fa_elec = energy_data[pyrosetta.rosetta.core.scoring.fa_elec]
        hbond_sr_bb = energy_data[pyrosetta.rosetta.core.scoring.hbond_sr_bb]
        hbond_lr_bb = energy_data[pyrosetta.rosetta.core.scoring.hbond_lr_bb]
        hbond_sc = energy_data[pyrosetta.rosetta.core.scoring.hbond_sc]

        # Calculate per-chain energy
        per_chain_energies = {}
        num_chains = pose.num_chains()
        for chain in range(1, num_chains + 1):
            try:
                chain_pose = pose.split_by_chain(chain)
                chain_energy = scorefxn(chain_pose)
                per_chain_energies[f"chain_{chain}"] = chain_energy
            except Exception as e:
                logger.warning(f"Could not calculate energy for chain {chain}: {e}")
                per_chain_energies[f"chain_{chain}"] = None
        
        # Calculate RMSD if a reference pose is provided
        rmsd_ca = None
        rmsd_all_atom = None
        rmsd_per_chain = {}
        if reference_pose:
            try:
                rmsd_ca = CA_rmsd(reference_pose, pose)
                rmsd_all_atom = all_atom_rmsd(reference_pose, pose)
                
                # Calculate per-chain RMSD
                for chain in range(1, min(num_chains, reference_pose.num_chains()) + 1):
                    try:
                        ref_chain = reference_pose.split_by_chain(chain)
                        pose_chain = pose.split_by_chain(chain)
                        chain_rmsd_ca = CA_rmsd(ref_chain, pose_chain)
                        rmsd_per_chain[f"rmsd_ca_chain_{chain}"] = chain_rmsd_ca
                    except Exception as e:
                        logger.warning(f"Could not calculate RMSD for chain {chain}: {e}")
                        rmsd_per_chain[f"rmsd_ca_chain_{chain}"] = None
            except Exception as e:
                logger.warning(f"Could not calculate RMSD: {e}")

        # Calculate SASA (Solvent Accessible Surface Area)
        total_sasa = None
        try:
            probe_radius = 1.4  # Standard probe radius for water
            total_sasa = calc_total_sasa(pose, probe_radius)
        except Exception as e:
            logger.warning(f"Could not calculate SASA: {e}")
        
        # Calculate interface binding energy
        binding_energy = None
        peptide_antibody_binding = None
        
        if interface_mode == "peptide_antibody" and num_chains >= 3:
            # OLD simple subtraction kept if you still want it:
            try:
                peptide_energy = per_chain_energies.get("chain_1", 0) or 0
                heavy_energy = per_chain_energies.get("chain_2", 0) or 0
                light_energy = per_chain_energies.get("chain_3", 0) or 0
                peptide_antibody_binding = total_energy - (peptide_energy + heavy_energy + light_energy)
            except Exception as e:
                logger.warning(f"Could not calculate simple peptide-antibody binding energy: {e}")

            # NEW: IAM-based ΔΔG for A vs (B+C)
            try:
                partners = "A_BC"   # peptide (A) vs antibody (B+C)
                binding_energy = interface_ddg_IAM(
                    pose,
                    partners_str=partners,
                    scorefxn=scorefxn,
                    pack_separated=True,  # repack each partner when separated
                    pack_input=True,      # repack the complex too
                    verbose=False
                )
                logger.info(f"IAM interface ΔΔG (A vs B+C): {binding_energy:.2f} REU")
            except Exception as e:
                logger.warning(f"InterfaceAnalyzerMover failed for A vs B+C: {e}")
                
        elif isinstance(interface_mode, tuple) and len(interface_mode) == 2 and num_chains >= 2:
            # Existing 2-chain case (A vs B, etc.)
            try:
                partners = f"{interface_mode[0]}_{interface_mode[1]}"
                binding_energy = interface_ddg_IAM(
                    pose, partners_str=partners, scorefxn=scorefxn,
                    pack_separated=True, pack_input=True, verbose=False
                )
            except Exception as e:
                logger.warning(f"Could not calculate 2-chain IAM binding energy: {e}")

        # Return comprehensive energy and structural information
        results = {
            "total_energy": total_energy,
            "fa_atr": fa_atr,                 # van der Waals attraction
            "fa_rep": fa_rep,                 # van der Waals repulsion
            "fa_sol": fa_sol,                 # solvation energy
            "fa_elec": fa_elec,               # electrostatic energy
            "hbond_sr_bb": hbond_sr_bb,       # short-range backbone H-bonds
            "hbond_lr_bb": hbond_lr_bb,       # long-range backbone H-bonds
            "hbond_sc": hbond_sc,             # sidechain H-bonds
            "rmsd_ca": rmsd_ca,               # C-alpha RMSD to reference
            "rmsd_all_atom": rmsd_all_atom,   # All-atom RMSD to reference
            "total_sasa": total_sasa,         # Solvent accessible surface area
            "binding_energy": binding_energy, # Interface binding energy (2-chain)
            "peptide_antibody_binding": peptide_antibody_binding, # Peptide-Antibody binding energy
            "num_chains": num_chains,         # Number of chains
            "num_residues": pose.total_residue() # Total number of residues
        }
        
        # Add per-chain energies and per-chain RMSD
        results.update(per_chain_energies)
        results.update(rmsd_per_chain)
        
        return results
        
    except Exception as e:
        logger.error(f"Error in collect_energy_info: {e}")
        return None


def process_simple_cif_directory(input_dir, reference_file=None, interface_mode=None):
    """
    Simple mode: Process all CIF files found recursively in a directory.

    Args:
        input_dir (str): Directory containing CIF files
        reference_file (str): Path to reference structure for RMSD calculation (optional)
        interface_mode: Interface analysis mode - "peptide_antibody" for 3-chain system or tuple for 2-chain

    Returns:
        dict: Dictionary with filenames as keys and energy results as values
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        logger.error(f"Input directory {input_dir} does not exist")
        return {}

    # Find all CIF files recursively
    cif_files = list(input_path.rglob("*.cif"))

    if not cif_files:
        logger.warning(f"No CIF files found in {input_dir}")
        return {}

    logger.info(f"Found {len(cif_files)} CIF files to process")

    # Load reference pose if provided
    reference_pose = None
    if reference_file:
        reference_pose = load_pose_from_cif(reference_file)
        if reference_pose is None:
            logger.warning("Could not load reference file, proceeding without RMSD calculations")

    # Process each CIF file
    results = {}
    for cif_file in cif_files:
        logger.info(f"Processing {cif_file.name}")

        pose = load_pose_from_cif(cif_file)
        if pose is None:
            continue

        energy_info = collect_energy_info(pose, reference_pose, interface_mode)
        if energy_info:
            # Create a descriptive filename from the directory structure
            # For nested paths like ATP9B_fragment37_mAb_CSF2_full/atp9b_fragment37_mab_csf2_full/atp9b_fragment37_mab_csf2_full_model.cif
            # We want to use: ATP9B_fragment37_mAb_CSF2_full_model
            relative_path = cif_file.relative_to(input_path)
            if len(relative_path.parts) > 1:
                # Use the top-level directory name + file stem
                result_key = f"{relative_path.parts[0]}_{cif_file.stem}"
            else:
                # Direct file in input directory
                result_key = cif_file.stem
            results[result_key] = energy_info
        else:
            logger.warning(f"Failed to collect energy info for {cif_file.name}")

    return results


def process_cif_directory(input_dir, reference_file=None, interface_mode=None, pattern='cd320_fab320', expected_kmer_size=None):
    """
    Process all CIF files in a directory and calculate energy metrics.
    
    Args:
        input_dir (str): Directory containing CIF files
        reference_file (str): Path to reference structure for RMSD calculation (optional)
        interface_mode: Interface analysis mode - "peptide_antibody" for 3-chain system or tuple for 2-chain
        pattern (str): Directory naming pattern prefix (e.g., "cd320_fab320" or "bcr")

    Returns:
        dict: Dictionary with filenames as keys and energy results as values
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        logger.error(f"Input directory {input_dir} does not exist")
        return {}
    
    # Find all final/best CIF files (exclude individual seed files)
    # Look for *_model.cif files in subdirectories (AF3 final outputs)
    cif_files = []
    
    # First check for direct CIF files in the input directory
    direct_cifs = list(input_path.glob("*.cif"))
    cif_files.extend(direct_cifs)
    
    # Then check for AF3-style nested structure: subdir/subdir/*_model.cif
    for subdir in input_path.iterdir():
        if subdir.is_dir():
            # Filter directories with EXACT BCR match and extract kmer size
            # Pattern format: {pattern}_{KMER}mer_window_{N}
            # This ensures BCR name is followed immediately by kmer size, not another BCR variant
            pattern_regex = re.escape(f'{pattern}_') + r'(\d+)mer_window_\d+'
            match = re.match(pattern_regex, subdir.name)
            if not match:
                continue  # Skip directories that don't match the exact pattern

            # Extract kmer size from directory name (captured group 1 from regex)
            kmer_size = match.group(1) + 'mer'

            # Filter by expected kmer size if specified
            # Normalize comparison: strip trailing 's' to handle both "15mer" and "15mers"
            if expected_kmer_size and kmer_size.rstrip('s') != expected_kmer_size.rstrip('s'):
                continue  # Skip this directory - wrong kmer size

            # Look for nested subdirectory with same/similar name
            # Handle window directories - nested dirs drop leading zeros from window numbers
            if 'window_' in subdir.name:
                # Extract kmer size and window number using regex
                # e.g., "BCAS3_epitope_mAb1_HA1C_15mer_window_00"
                kmer_match = re.search(r'_(\d+mer)_window_(\d+)', subdir.name)
                if not kmer_match:
                    continue  # Skip if no valid kmer and window found
                kmer_size = kmer_match.group(1)  # e.g., "15mer"

                # Filter by expected kmer size if specified
                # Normalize comparison: strip trailing 's' to handle both "15mer" and "15mers"
                if expected_kmer_size and kmer_size.rstrip('s') != expected_kmer_size.rstrip('s'):
                    continue  # Skip this directory - wrong kmer size

                window_num_int = int(kmer_match.group(2))  # Extracts window number as int
                # Nested directory uses lowercase pattern and drops leading zeros
                nested_subdir = subdir / f"{pattern.lower()}_{kmer_size}_window_{window_num_int}"
            else:
                # For non-window patterns, use simple lowercase transformation
                nested_subdir = subdir / subdir.name.lower().replace('-', '_')
            if nested_subdir.exists() and nested_subdir.is_dir():
                # Find the final model CIF file (not individual seeds)
                model_cifs = list(nested_subdir.glob("*_model.cif"))
                if model_cifs:
                    cif_files.append(model_cifs[0])  # Take only first CIF from canonical directory
                else:
                    # Canonical dir exists but has no CIF, search recursively in window dir
                    subdir_cifs = list(subdir.rglob("*_model.cif"))
                    if subdir_cifs:
                        cif_files.append(subdir_cifs[0])
            else:
                # Canonical dir doesn't exist, search recursively for any CIF
                subdir_cifs = list(subdir.rglob("*_model.cif"))
                if subdir_cifs:
                    cif_files.append(subdir_cifs[0])
    
    if not cif_files:
        logger.warning(f"No CIF files found in {input_dir} or its subdirectories")
        return {}
    
    logger.info(f"Found {len(cif_files)} CIF files to process")
    
    # Load reference pose if provided
    reference_pose = None
    if reference_file:
        reference_pose = load_pose_from_cif(reference_file)
        if reference_pose is None:
            logger.warning("Could not load reference file, proceeding without RMSD calculations")
    
    # Process each CIF file
    results = {}
    for cif_file in cif_files:
        logger.info(f"Processing {cif_file.name}")
        
        pose = load_pose_from_cif(cif_file)
        if pose is None:
            continue
        
        energy_info = collect_energy_info(pose, reference_pose, interface_mode)
        if energy_info:
            # Use parent directory name + file stem for unique naming
            if cif_file.parent.name != input_path.name:
                result_key = f"{cif_file.parent.parent.name}_{cif_file.stem}"
            else:
                result_key = cif_file.stem
            results[result_key] = energy_info
        else:
            logger.warning(f"Failed to collect energy info for {cif_file.name}")
    
    return results


def save_results_to_csv(results, output_file):
    """
    Save energy analysis results to a CSV file.
    
    Args:
        results (dict): Dictionary containing energy analysis results
        output_file (str): Path to output CSV file
    """
    if not results:
        logger.warning("No results to save")
        return
    
    # Convert results to DataFrame
    df = pl.DataFrame(list(results.values()), schema=list(next(iter(results.values())).keys()))
    # Add the index (filename) as the first column
    df = df.with_columns(pl.Series("filename", list(results.keys())))
    df = df.select(["filename"] + [col for col in df.columns if col != "filename"])
    
    # Reorder columns for better readability
    basic_cols = ['total_energy', 'num_chains', 'num_residues']
    energy_cols = ['fa_atr', 'fa_rep', 'fa_sol', 'fa_elec', 'hbond_sr_bb', 'hbond_lr_bb', 'hbond_sc']
    structural_cols = ['rmsd_ca', 'rmsd_all_atom', 'total_sasa', 'binding_energy', 'peptide_antibody_binding']
    
    # Get chain energy and RMSD columns
    chain_energy_cols = [col for col in df.columns if col.startswith('chain_') and not col.startswith('rmsd')]
    chain_rmsd_cols = [col for col in df.columns if col.startswith('rmsd_ca_chain_')]
    chain_energy_cols.sort()
    chain_rmsd_cols.sort()
    chain_cols = chain_energy_cols + chain_rmsd_cols
    
    # Reorder columns
    column_order = ['filename'] + basic_cols + energy_cols + structural_cols + chain_cols
    available_cols = [col for col in column_order if col in df.columns]
    remaining_cols = [col for col in df.columns if col not in available_cols]
    final_cols = available_cols + remaining_cols
    
    df = df.select(final_cols)
    df = df.sort("filename")
    
    # Save to CSV
    df.write_csv(output_file)
    logger.info(f"Results saved to {output_file}")
    
    # Print summary
    print(f"\nEnergy Analysis Summary:")
    print(f"Processed {len(results)} structures")
    print(f"Results saved to: {output_file}")
    
    if 'total_energy' in df.columns:
        total_energy_col = df.select('total_energy')
        min_energy = total_energy_col.min().item(0, 0)
        max_energy = total_energy_col.max().item(0, 0)
        print(f"Total energy range: {min_energy:.2f} to {max_energy:.2f}")
    
    print(f"\nFirst few rows:")
    print(df.head())


def main():
    """Main function to run the energy analysis."""
    parser = argparse.ArgumentParser(description='Analyze energy terms of AlphaFold3 CIF files using Rosetta')
    parser.add_argument('--input_dir', required=True, help='Directory containing CIF files')
    parser.add_argument('--output', default='rosetta_energy_results.csv', help='Output CSV file')
    parser.add_argument('--reference', help='Reference CIF file for RMSD calculations')
    parser.add_argument('--interface', nargs='*', metavar='CHAIN', 
                      help='Interface analysis: use "peptide_antibody" for 3-chain system (A vs B+C) or two chain IDs for 2-chain analysis (e.g., A B)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')
    parser.add_argument('--pattern', default='cd320_fab320',
                      help='Directory naming pattern prefix (e.g., "cd320_fab320" or "bcr")')
    parser.add_argument('--no-windows', action='store_true',
                      help='Simple mode: find all CIF files in input directory and process them directly')
    parser.add_argument('--kmer-size', type=str,
                      help='Filter to only process structures with this kmer size (e.g., "15mer" or "25mer")')

    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize PyRosetta
    initialize_pyrosetta()
    
    # Process interface argument
    interface_mode = None
    if args.interface:
        if len(args.interface) == 1 and args.interface[0] == "peptide_antibody":
            interface_mode = "peptide_antibody"
        elif len(args.interface) == 2:
            interface_mode = tuple(args.interface)
        else:
            logger.error("Interface argument must be either 'peptide_antibody' or two chain IDs")
            sys.exit(1)
    
    # Process CIF files
    if args.no_windows:
        # Simple mode: find all CIF files recursively
        results = process_simple_cif_directory(
            args.input_dir,
            args.reference,
            interface_mode
        )
    else:
        # Window-based mode with pattern matching
        results = process_cif_directory(
            args.input_dir,
            args.reference,
            interface_mode,
            args.pattern,
            args.kmer_size
        )
    
    if results:
        # Save results
        save_results_to_csv(results, args.output)
    else:
        logger.error("No results to save. Check your input directory and files.")
        sys.exit(1)


if __name__ == "__main__":
    main()
