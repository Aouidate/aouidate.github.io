---
title: "An AI Agent for Virtual Screening with Smina and Ollama"
date: 2026-03-16
author: "Adnane Aouidate"
permalink: /posts/2026/03/AI_Agent_Virtual_Screening/
tags:
  - Cheminformatics
  - AI Agents
  - Drug Design
  - Molecular Docking
  - Ollama
---

# AI Agent for Virtual Screening with Smina

## Overview

This tutorial demonstrates how to build an **AI Agent for Virtual Screening** using Ollama with local LLMs. The agent acts as a molecular docking expert that can autonomously perform virtual screening tasks using natural language commands.

## Tutorial: ADN_T018_AI_Agent_Virtual_Screening_with_Smina.ipynb

### What You'll Learn

1. **AI Agent Architecture**: Understand how to build an AI agent with tools
2. **Ollama Integration**: Use Ollama with `gpt-oss:latest` for local LLM inference
3. **Molecular Docking Expert**: Create a domain-specific AI persona
4. **Tool Integration**: Connect computational tools (Smina, OpenBabel) to the agent
5. **Function Calling**: Implement autonomous tool usage by the AI
6. **Web Deployment**: Convert to a Streamlit web application

### Key Components

#### 🤖 AI Brain

- **LLM**: Ollama with `gpt-oss:latest` model
- **Framework**: Function calling for tool orchestration
- **Persona**: Molecular docking expert with deep domain knowledge

#### 🛠️ Tools

1. **OpenBabel**: Protein and ligand preparation (hydrogen addition, format conversion)
2. **Smina**: Molecular docking engine
3. **Pandas**: Data analysis and pose ranking
4. **py3Dmol**: 3D visualization of protein-ligand complexes

#### 💻 Web Interface

- **Streamlit**: Interactive web app for the AI agent
- **Features**: File upload, chat interface, 3D visualization, results download

## Prerequisites

### Software Requirements

1. **Ollama** - Local LLM runtime

   ```bash
   # Download from https://ollama.ai/download
   # For macOS/Linux, install and pull the model:
   ollama pull gpt-oss:latest
   ```

2. **Conda** - For package management (recommended)

### Python Packages

```bash
# Create environment
conda create -n vs_agent python=3.9
conda activate vs_agent

# Install core packages
conda install -c conda-forge rdkit openbabel smina py3dmol
pip install pandas ollama streamlit stmol streamlit-ketcher
```

## Getting Started

### 1. Run the Jupyter Notebook

```bash
cd notebooks
jupyter notebook ADN_T018_AI_Agent_Virtual_Screening_with_Smina.ipynb
```

Work through the cells sequentially:

- **Section 1**: Setup Ollama and verify model
- **Section 2**: Define molecular docking expert persona
- **Section 3**: Create tool functions
- **Section 4**: Build the AI agent
- **Section 5**: Demo the agent with virtual screening tasks
- **Section 6**: Deploy as Streamlit web app

### 2. Example Agent Interactions

The agent can understand natural language commands like:

```python
# Example 1: Protein preparation
agent.chat("""
I have a PDB file called 5R82.pdb. Please:
1. Extract the protein and ligand (residue name RZS) separately
2. Add hydrogens to the protein at pH 7.4
3. Tell me about the structure
""")

# Example 2: Virtual screening
agent.chat("""
Please perform virtual screening:
1. Prepare ligands from '../databases/data-test.sdf'
2. Dock them against 5R82_proteinH.pdb using reference ligand 5R82_ligand.pdb
3. Analyze the top 5 compounds and tell me about their binding affinities
""")
```

### 3. Deploy as Web App

After completing the notebook, save the required files:

```python
# Save the tool classes to vs_tools.py
# Save the Streamlit code to streamlit_vs_agent.py

# Run the app
streamlit run streamlit_vs_agent.py
```

Access at `http://localhost:8501`

## Notebook Structure

```
ADN_T018_AI_Agent_Virtual_Screening_with_Smina.ipynb
├── 1. Introduction & Overview
├── 2. Setup Ollama and Verification
├── 3. Define Molecular Docking Expert Persona
├── 4. Create Tool Functions
│   ├── prepare_protein()
│   ├── prepare_ligands()
│   ├── perform_docking()
│   ├── analyze_results()
│   └── extract_protein_ligand()
├── 5. Build AI Agent with Function Calling
│   ├── MolecularDockingAgent class
│   ├── Tool specifications
│   └── Chat interface
├── 6. Demo: Virtual Screening Workflow
│   ├── Protein preparation
│   ├── Ligand preparation
│   ├── Molecular docking
│   ├── Results analysis
│   └── 3D visualization
├── 7. Convert to Streamlit Web App
│   ├── Complete app code
│   ├── Deployment instructions
│   └── Usage guide
└── 8. Summary and Extensions
```

## Features

### ✅ What the AI Agent Can Do

- **Autonomous Operation**: Break down complex tasks into steps
- **Smart Tool Selection**: Choose appropriate tools based on context
- **Expert Insights**: Provide domain knowledge about binding modes and interactions
- **Interactive Chat**: Natural language interface
- **Error Handling**: Gracefully handle issues and provide helpful feedback
- **Result Visualization**: Generate 3D views of docking results

### 🔄 Example Workflow

1. User: "I want to screen a library of compounds against SARS-CoV-2 main protease"
2. Agent:
   - Extracts protein and reference ligand from PDB
   - Prepares protein (adds hydrogens)
   - Prepares ligand library (protonation, format conversion)
   - Sets up and runs docking with optimal parameters
   - Analyzes results and ranks compounds
   - Visualizes top hits
   - Provides expert commentary on binding modes

## Advantages Over Traditional Workflows

| Traditional Approach         | AI Agent Approach               |
| ---------------------------- | ------------------------------- |
| Manual command execution     | Natural language requests       |
| Need to remember parameters  | Agent suggests optimal settings |
| Sequential scripting         | Autonomous task orchestration   |
| Technical expertise required | Accessible to non-experts       |
| Static workflows             | Adaptive and flexible           |

## Extensions and Next Steps

### Additional Tools to Add

- AutoDock Vina integration
- ADMET prediction (SwissADME, pkCSM)
- Protein-ligand interaction fingerprints (PLIF)
- Molecular dynamics preparation
- Multi-target docking

### Enhanced Analysis

- Binding site identification
- SAR (Structure-Activity Relationship) analysis
- Pharmacophore modeling
- Fragment-based drug design

### Advanced Features

- Batch processing of large libraries (100K+ compounds)
- Cloud deployment (AWS, Google Cloud, Azure)
- Database integration (ChEMBL, PubChem)
- Automated report generation
- Email notifications for long jobs

### Model Improvements

- Fine-tune on chemistry/docking data
- Use larger models (GPT-4, Claude)
- Implement RAG (Retrieval-Augmented Generation) with scientific literature

## Troubleshooting

### Common Issues

1. **Ollama not found**

   ```bash
   # Verify installation
   which ollama
   ollama --version
   ```

2. **Model not available**

   ```bash
   # Pull the model
   ollama pull gpt-oss:latest

   # List available models
   ollama list
   ```

3. **Smina not found**

   ```bash
   # Install via conda
   conda install -c conda-forge smina

   # Or download static binary from sourceforge
   ```

4. **OpenBabel issues**
   ```bash
   # Reinstall
   conda install -c conda-forge openbabel --force-reinstall
   ```

## Resources

### Documentation

- **Ollama**: https://ollama.ai/
- **Smina**: https://sourceforge.net/projects/smina/
- **OpenBabel**: http://openbabel.org/
- **RDKit**: https://www.rdkit.org/
- **Streamlit**: https://streamlit.io/

### Related Tutorials

- ADN_T013: Molecular Docking with Smina (basic)
- ADN_T017: Virtual Screening with Smina (intermediate)

### Papers and References

- Koes et al. (2013). "Lessons Learned in Empirical Scoring with smina"
- Trott & Olson (2010). "AutoDock Vina: improving the speed and accuracy of docking"

## Citation

If you use this tutorial in your research, please cite:

```
Aouidate, A. (2023-2026). AI Agent for Virtual Screening with Smina and Ollama.
Chemoinformatics Tutorials. Ibn Zohr University, Agadir, Morocco.
```

## Contact

For questions or feedback:

- Open an issue on GitHub
- Contact: [Your contact info]

## License

This tutorial is released under the [LICENSE] specified in the repository.

---

**Happy Docking! 🧬🚀**
