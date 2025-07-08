---
title: "Incorporating Stereochemical Information into ML-Driven QSAR and QSPR Models"
date: 2025-07-08
permalink: /posts/2025/08/Encoding_Sterochemistry/
tags:
  - Cheminformatics
  - Stereochemistry
  - Drug Design
---

<div style="text-align: justify;">

In this blog, we will explore how to encode molecular stereochemistry for QSAR models. We’ll begin by discussing the importance of QSAR in drug discovery, the critical role stereochemistry plays in drug effectiveness and safety, and then dive into practical ways to represent stereochemistry using accessible tools like molecular fingerprints. Finally, we’ll compare how different fingerprint types capture stereochemical differences and which ones are best suited for this purpose.

</div>

# What is QSAR and Why Is It Important?

<div style="text-align: justify;">

QSAR/QSPR (Quantitative Structure-Activity/Property Relationship) modeling is a computational technique used to predict the biological activity or properties of chemical compounds based on their molecular structures. By building mathematical relationships between chemical features (molecular descriptors or fingerprints) and experimental activity data, QSAR helps identify promising drug candidates, optimize their properties, and assess chemical safety without extensive lab testing. This accelerates drug discovery and reduces costs and risks associated with chemical exposure, as recommended by OECD guidelines.

</div>

# The Crucial Role of Stereochemistry in Real Life and QSAR

<div style="text-align: justify;">

Stereochemistry - the 3D spatial arrangement of atoms in a molecule - often plays a vital role in determining biological activity. Different stereoisomers, such as enantiomers or diastereomers, can interact very differently with biological targets, leading to significant variations in efficacy, potency, or toxicity.

A well-known example is <b> thalidomide** </b>, which was withdrawn from the market in 1961. Its R-enantiomer acted as a sedative, while the S-enantiomer caused severe birth defects.

Imagine a classroom of 30 students and a key that opens the door to their study room. I, as the professor, make 30 copies of the key: 15 are exact copies (R-enantiomers) and 15 are mirror images (S-enantiomers). I randomly hand out these keys to the students. Only the 15 students with the exact copies can open the door and enter, while the other 15 cannot. However, these mirror-image keys might accidentally open other doors - unknown to us - and if so, may cause damage.

This metaphor illustrates how the S-enantiomer of thalidomide interacts with proteins unrelated to its intended target, leading to harmful effects.

<div style="text-align: center;">

<img src="/images/Encoding_Stereochemistry/thalidomide.png" alt="Thalidomide enantiomers" width="600" height="400" class="img-fluid rounded mx-auto d-block mb-4" loading="lazy" />

</div>

Another example is <b> limonene </b>: the R-enantiomer smells like orange, whereas its mirror image, the S-enantiomer, has a lemon-like scent. These stark differences highlight why accurately capturing stereochemical information in QSAR models is essential for reliably predicting biological activity and designing safer, more effective drugs.

One of the main challenges in QSAR modeling is how to accurately represent the 3D stereochemical features of molecules. Capturing these subtle 3D differences is essential but not straightforward.

3D-QSAR methods, such as CoMFA and CoMSIA, inherently consider stereochemistry and can be effective solutions. However, these approaches are generally limited to series of analogs sharing the same scaffold, restricting their broader applicability.

Traditional (or classical) QSAR methods rely mainly on 2D descriptors or fingerprints, which often ignore or oversimplify stereochemical details. This can result in models unable to distinguish between stereoisomers, leading to poor predictions of activity, properties, or toxicity. On the other hand, fully 3D approaches require reliable 3D structures or conformers, which are computationally expensive to generate and highly dependent on the quality of conformer sampling and alignment.

Furthermore, integrating stereochemistry into descriptors or fingerprints in a way that machine learning algorithms can effectively use remains an active area of research. Striking the right balance between computational efficiency and stereochemical accuracy is a key challenge in building robust, predictive QSAR models.

</div>

# Comparing Fingerprints of R- and S- Thalidomide to Detect Stereochemistry Differences

To illustrate how different fingerprint types capture stereochemical differences, we will compare the R- and S-enantiomers of thalidomide using various molecular fingerprints. We will use RDKit, a powerful cheminformatics toolkit, to generate and visualize these fingerprints.

Here we will import the necessary libraries and define the SMILES strings for the R- and S-enantiomers of thalidomide. and then generate the Rdkit molecular objects from these SMILES strings.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors, DataStructs, rdFingerprintGenerator, Draw
from mapchiral.mapchiral import encode, jaccard_similarity


# Define the SMILES strings for the R- and S-enantiomers of thalidomide
smiles_R = "O=C(N1)CC[C@@H](N2C(C3=CC=CC=C3C2=O)=O)C1=O"
smiles_S = "O=C(N1)CC[C@H](N2C(C3=CC=CC=C3C2=O)=O)C1=O"


mol_R = Chem.MolFromSmiles(smiles_R)
mol_S = Chem.MolFromSmiles(smiles_S)
```
Th next we generate the fingerprints for both enantiomers using various fingerprint generators from RDKit. We will use the Morgan fingerprint (It is worth to note that it is important to use includeChirality=True), RDKit fingerprint, Topological Torsion fingerprint, MACCS keys, and MapFps (a chiral fingerprint generator). These fingerprints will help us understand how well each type captures the stereochemical differences between the two enantiomers.

```python 

# Function to calculate different fingerprints

mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048, includeChirality=True)
ttgen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)

def get_fingerprints(mol):
    fps = {}
    # Morgan fingerprint (radius 2), chirality enabled
    fps['Morgan'] = mfpgen.GetFingerprint(mol)
    # RDKit fingerprint (default)
    fps['RDKit'] = Chem.RDKFingerprint(mol)
    # Topological torsion fingerprint
    fps['TopologicalTorsion'] = ttgen.GetFingerprint(mol)
    # MACCS keys
    fps['MACCS'] = MACCSkeys.GenMACCSKeys(mol)
    #Map Chiral fingerprints
    fps['MapFps'] = encode(mol, max_radius=2, n_permutations=2048, mapping=False)

    return fps

fps_R = get_fingerprints(mol_R)
fps_S = get_fingerprints(mol_S)
```
Once we have the fingerprints, we can check the similarity between the R- and S-enantiomers, which will help us understand how well each fingerprint type captures the stereochemical differences. If they capture it, it means that the fingerprint can distinguish between the two enantiomers, which is crucial for accurate QSAR modeling. if not, it means that the fingerprint is not sensitive to stereochemistry and may not be suitable for applications where stereochemical differences are important.
```python

# Calculate Tanimoto similarity for bit vector fingerprints
def tanimoto(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def dice(fp1, fp2):
    # Convert ExplicitBitVect to set of bits
    return DataStructs.DiceSimilarity(fp1, fp2)

print("Similarity between R- and S-limonene:")

for key in ['Morgan', 'RDKit', 'MACCS']:
    sim = tanimoto(fps_R[key], fps_S[key])
    print(f"{key} fingerprint similarity: {sim:.3f}")

sim_tt = dice(fps_R['TopologicalTorsion'], fps_S['TopologicalTorsion'])
print(f"Topological Torsion fingerprint similarity: {sim_tt:.3f}")

sim_jac = jaccard_similarity(fps_R["MapFps"], fps_S["MapFps"])
print(f"MinHashed Atom-Pair Fingerprint Chiral similarity: {sim_jac:.3f}")
```
So the output will be something like this:

```
Similarity between R- and S-limonene:

Morgan fingerprint similarity: 0.714

RDKit fingerprint similarity: 1.000

MACCS fingerprint similarity: 1.000

Topological Torsion fingerprint similarity: 1.000

MinHashed Atom-Pair Fingerprint Chiral similarity: 0.879
```

As we can see from the similarity scores, the Morgan fingerprint is able to capture some of the stereochemical differences between the two enantiomers, while the RDKit, MACCS and Topological Torsion fingerprints are not. However, the MapFps fingerprint, which is specifically designed to capture stereochemical information, shows a high similarity score of 0.879, indicating that it can effectively distinguish between the two enantiomers. This demonstrates the importance of using appropriate fingerprint types that can capture stereochemical differences for accurate QSAR modeling. and even without the need to generate 3D conformers, which can be computationally expensive and time-consuming.