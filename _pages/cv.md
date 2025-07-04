---
layout: archive
title: "CV"
permalink: /cv/
author_profile: true
redirect_from:
  - /resume
---

{% include base_path %}

# Profiles
* [Google Scholar](https://scholar.google.com/citations?user=nzQ-UzAAAAAJ&hl=en&oi=ao)
<!-- * [ORCID](http://orcid.org/yourorcidurl) -->
* [GitHub](https://github.com/Aouidate)
* [ResearchGate](https://www.researchgate.net/profile/Adnane-Aouidate?ev=hdr_xprf)

# Experience
## Assistant Professor (Computational Chemist)  
*School of Applied Sciences of Ait Melloul, Ibn Zohr University, Agadir, Morocco*  
01/2023 – Present  
- Developed web applications for chemoinformatics tools using Streamlit & Django  
- Explored chemical space of BRAF inhibitors through cheminformatics and AI  
- Supervised PhD students working on ADMET predictions using ML and DL-based QSAR models  
- Delivered lectures and tutorials to students  
- Led junior researcher supervision and feedback  

## Postdoctoral Researcher (Chemoinformatics)  
*Structural Bioinformatic & Chemoinformatic Team, ICOA, Orleans, France*  
01/2021 – 12/2022  
- Encoded ACE2/SARS-CoV-2 interactions as fingerprints using RDKit and ODDT  
- Developed methods to dock and score protein-protein interaction inhibitors using Tanimoto coefficient and SPLIF-score  
- Conducted molecular dynamics simulations of protein-ligand complexes  
- Applied Fragment-Based Drug Discovery to generate LIMKs protein kinase inhibitors  

## Postdoctoral Researcher (Computational Chemist)  
*CADD Center, SIAT, Shenzhen, China*  
2019 – 2020  
- Automated AI-based drug discovery workflows with KNIME to accelerate R&D  
- Generated large chemical libraries and performed scaffold hopping  
- Predicted molecular properties and biological activities of small molecules  
- Curated chemical databases and performed in silico ADMET/toxicity assessment  

# Skills
- Structure-based drug design (SBDD): Molecular Docking, Virtual Screening, Molecular Dynamics, Homology Modeling  
- Ligand-based drug design (LBDD): Pharmacophore modeling, 2D & 3D QSAR (CoMFA, CoMSIA)  
- Chemical database management and retrieval  
- Automation of cheminformatics and ML workflows with KNIME  
- Machine learning models for drug discovery, risk assessment, toxicology  
- Database preparation for virtual screening  
- In silico ADMET assessment  
- Software: Schrödinger Suite, Gaussian09, MOE, Material Studio, Autodock, Gromacs, KNIME, DataWarrior, RDKit, DeepChem  
- Programming: Python (Data Science), Streamlit, Git, Jupyter Notebook, Docker, Django, Singularity  
- Operating Systems: Linux, Windows  

# Education
- **Ph.D in Chemistry**, UMI University, 2019  
  Focus: New bioactive organic molecules related to inhibition of protein kinases using 3D-QSAR, Molecular Docking, and ADMET  
- **M.S. in Chemistry**, UMI University, 2012–2014  
  Focus: Molecular Chemistry and Natural Products  
- **B.S. in Chemistry**, UMI University, 2009–2012  
  Focus: Organic Chemistry  

# Languages
- Arabic: Fluent  
- French: Fluent  
- English: Fluent  
- Spanish: Elementary  

# Interests
- Running  
- Traveling  
- Data science  

# Selected Publications
<ul>{% for post in site.publications reversed %}
  {% include archive-single-cv.html %}
{% endfor %}</ul>

# Presentations
<ul>{% for post in site.talks reversed %}
  {% include archive-single-talk-cv.html %}
{% endfor %}</ul>

# Teaching
<ul>{% for post in site.teaching reversed %}
  {% include archive-single-cv.html %}
{% endfor %}</ul>
