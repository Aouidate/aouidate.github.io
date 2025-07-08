---
title: "A Beginner’s Guide to Molecular Docking with Smina"
collection: tutorials
type: "Undergraduate course"
permalink: /tutorials/2023-Molecular-Docking
venue: "School of Applied Sciences, Ait-Melloul, Chemistry Department"
date: 2023-07-18
location: "Agadir, Morocco"
author: Aouidate
tags: [molecular docking, smina, computational chemistry, tutorial, undergraduate]
---

**A Beginner’s Guide to Molecular Docking with Smina**  
*by Aouidate*  

*Originally on Medium : [Medium Article](https://medium.com/loops-strands/a-beginners-guide-to-molecular-docking-with-smina-e1e4360950c2)*

---

Molecular docking is an essential technique in computational drug discovery, allowing us to predict how small molecules (ligands) interact with biological targets (proteins). With the rise of open-source tools, this field is more accessible than ever, even for undergraduate students. One such tool is **Smina**—a fork of AutoDock Vina with enhanced scoring and user features.

In this tutorial, we’ll walk through a beginner-friendly workflow for running molecular docking experiments with Smina.

---

## Table of Contents

1. [What is Molecular Docking?](#what-is-molecular-docking)
2. [Why Use Smina?](#why-use-smina)
3. [Setting Up Your Environment](#setting-up-your-environment)
4. [Preparing the Protein and Ligand](#preparing-the-protein-and-ligand)
5. [Running Docking with Smina](#running-docking-with-smina)
6. [Analyzing the Results](#analyzing-the-results)
7. [Tips and Best Practices](#tips-and-best-practices)
8. [Further Reading](#further-reading)

---

## What is Molecular Docking?

Molecular docking simulates the interaction between a receptor (usually a protein) and a ligand (a small molecule or drug candidate). The goal is to predict the preferred orientation of the ligand when bound to the protein and to estimate the binding affinity.

---

## Why Use Smina?

- **Open Source:** Smina is free and actively maintained.
- **Improved Scoring:** Better scoring functions than vanilla AutoDock Vina.
- **Customizability:** Allows user-defined scoring functions.
- **Cross-Platform:** Runs on Windows, macOS, and Linux.

---

## Setting Up Your Environment

### 1. Install Smina

You can download pre-built binaries or build from source.  
**Example (Linux):**
```bash
wget https://github.com/mwojcikowski/smina/releases/download/2020-12-17/smina.static -O smina
chmod +x smina
```
Place the executable in a directory included in your `$PATH`.

### 2. Install Support Tools

- **Open Babel:** For format conversion (`sudo apt install openbabel`)
- **PDBFixer or PyMOL:** For protein preparation

---

## Preparing the Protein and Ligand

### 1. Protein Preparation

- Download your protein structure (PDB format) from the [Protein Data Bank](https://www.rcsb.org/).
- Remove water molecules and unwanted chains.
- Add hydrogens and assign correct protonation states (using PDBFixer, PyMOL, or Chimera).
- Save as PDBQT format (required by Smina):
  ```bash
  obabel protein.pdb -O protein.pdbqt
  ```

### 2. Ligand Preparation

- Draw or download your ligand structure.
- Convert to 3D, add hydrogens.
- Save as PDBQT:
  ```bash
  obabel ligand.sdf -O ligand.pdbqt --gen3d
  ```

---

## Running Docking with Smina

Set the binding site (center and size) based on your protein structure (e.g., using PyMOL or Chimera to find the coordinates):

```bash
./smina -r protein.pdbqt -l ligand.pdbqt --center_x 10 --center_y 15 --center_z 20 --size_x 20 --size_y 20 --size_z 20 --out docked_ligand.pdbqt --log docking.log
```

- `-r`: receptor (protein)
- `-l`: ligand
- `--center_*` and `--size_*`: define the docking box
- `--out`: output file
- `--log`: log file

---

## Analyzing the Results

- Use PyMOL, Chimera, or UCSF ChimeraX to visualize the docked poses.
- The log file provides binding affinity scores (in kcal/mol).
- Look for low (more negative) scores and reasonable binding modes.

---

## Tips and Best Practices

- Always verify the preparation of your structures.
- Try multiple docking runs for reliability.
- Visual inspection is key—scores are not everything!
- Compare results with literature or experimental data when possible.

---

## Further Reading

- [Smina GitHub Repository](https://github.com/mwojcikowski/smina)
- [AutoDock Vina Manual](http://vina.scripps.edu/manual.html)
- [PyMOL Documentation](https://pymol.org/)
- [Open Babel Documentation](http://openbabel.org/wiki/Main_Page)

---

*Happy Docking! If you have questions or want to share your results, feel free to comment below or reach out.*
