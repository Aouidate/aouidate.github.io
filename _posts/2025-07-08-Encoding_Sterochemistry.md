---
title: 'Incorporating Stereochemical Information into ML-Driven QSAR and QSPR Models'
date: 2025-07-08
permalink: /posts/2025/08/Encoding_Sterochemistry/
tags:
  - Cheminformatics
  - Sterochemistry
  - Drug Design
---

In this blog, we will explore how to encode molecular stereochemistry for QSAR models. We’ll begin by discussing the importance of QSAR in drug discovery, the critical role stereochemistry plays in drug effectiveness and safety, and then dive into practical ways to represent stereochemistry using accessible tools like molecular fingerprints. Finally, we’ll compare how different fingerprint types capture stereochemical differences and which ones are best suited for this purpose.

# What is QSAR and Why Is It Important?

QSAR/QSPR Quantitative Structure-Activity/Property Relationship modeling is a computational technique used to predict the biological activity or properties of chemical compounds based on their molecular structures. By building mathematical relationships between chemical features (Moleuclar descriptors or fingerprints) and experimental activity data, QSAR can help to identify promising drug candidates, optimise their properties, and assess chemical safety without extensive lab testing. This accelerates drug discovery and reduces costs and risks associated with chemical exposure, which is recommended by OECD

# The Crucial Role of Stereochemistry in Real Life and QSAR

Stereochemistry - the 3D spatial arrangement of atoms in a molecule - often plays a vital role in determining biological activity. Different stereoisomers, such as enantiomers or diastereomers, can interact very differently with biological targets, leading to significant variations in efficacy, potency, or even toxicity.

A well-known example is **thalidomide**, which was withdrawn from the market in 1961. Its R-enantiomer acted as a sedative, while the S-enantiomer caused severe birth defects. Imagine a classroom of 30 students and a key that opens the door to their study room. I, as the professor, make 30 copies of the key: 15 are exact copies (R-enantiomers) and 15 are mirror images (S-enantiomers). I randomly hand out these keys to the students. Only the 15 students with the exact copies can open the door and enter, while the other 15 cannot. However, these mirror-image keys might accidentally open other doors — unknown to us — and if so, may cause damage. This metaphor illustrates how the S-enantiomer of thalidomide interacts with proteins unrelated to its intended target, leading to harmful effects.

<img src="/images/Encoding_Stereochemistry/thalidomide.png" alt="Thalidomide enantiomers" width="600" height="400" class="img-fluid rounded mx-auto d-block mb-4" loading="lazy" />

Another example is **limonene**: the R-enantiomer smells like orange, whereas its mirror image, the S-enantiomer, has a lemon-like scent. These stark differences highlight why accurately capturing stereochemical information in QSAR models is essential for reliably predicting biological activity and designing safer, more effective drugs.

One of the main challenges in QSAR modeling is how to accurately represent the 3D stereochemical features of molecules. Capturing these subtle 3D differences is essential but not straightforward.

3D-QSAR methods, such as CoMFA and CoMSIA, inherently consider stereochemistry and can be effective solutions. However, these approaches are generally limited to series of analogs sharing the same scaffold, which restricts their broader applicability.

Traditional (or classical) QSAR methods rely mainly on 2D descriptors or fingerprints, which often ignore or oversimplify stereochemical details. This can result in models unable to distinguish between stereoisomers, leading to poor predictions of activity, properties, or toxicity. On the other hand, fully 3D approaches require reliable 3D structures or conformers, which are computationally expensive to generate and highly dependent on the quality of conformer sampling and alignment.

Furthermore, integrating stereochemistry into descriptors or fingerprints in a way that machine learning algorithms can effectively use remains an active area of research. Striking the right balance between computational efficiency and stereochemical accuracy is a key challenge in building robust, predictive QSAR models.

# Comparing Fingerprints of R- and S- Thalidomide to Detect Stereochemistry Differences
Coming...

------