# asi_blender

## Table of Contents
- [Introduction](#introduction)
- [Set-up](#set-up)
- [How to Run](#how-to-run)
- [Cite this Work](#cite-this-work)

## Introduction 

This repository contains code to automatically identify the axon–spine interface (ASI) using 3D reconstructions of presynaptic axons and postsynaptic dendritic spines, and to quantify the degree of astroglial apposition at the ASI.

Read more in our [preprint](https://doi.org/10.1101/2025.05.13.653827).

## Set-up 

1. **Blender version**: This script currently supports **Blender 3.6.5**.

2. **Enable LoopTools Add-on**:
   - In Blender, go to `Edit` → `Preferences` → `Add-ons`.
   - Search for "LoopTools" and check the box to enable it.
   - 📄 [LoopTools Documentation](https://docs.blender.org/manual/en/3.5/addons/mesh/looptools.html)

3. **Object naming conventions** (required for the script to identify structures):
   - `d##c##` → **c-object**
   - `d##sp##` → **spine object**
   - `d##ax##` → **axon object**
   - `astroAll` → **astroglia object**
   - `1um_rad_ref` → **1 µm reference circle**

## How to Run 

1. Edit the path to your `.blend` file in the script:
   ```bash
   blender_asi_analysis.sh
2. Set the correct path to your Blender executable within the same script. For example:
   ```bash
   /path/to/blender --background your_file.blend --python run_asi_analysis.py
3. Make the shell script executable:
   ```bash
   chmod u+x blender_asi_analysis.sh
4. Run the script:
   ```bash
   ./blender_asi_analysis.sh

## Cite this Work 

If you use this code or data, please cite our preprint:
Nam, A. J., Kuwajima, M., Parker, P. H., Bowden, J. B., Abraham, W. C., & Harris, K. M. (2025). Perisynaptic astroglial response to in vivo long-term potentiation and concurrent long-term depression in the hippocampal dentate gyrus. bioRxiv. https://doi.org/10.1101/2025.05.13.653827

```bibtex
@article{nam2025perisynaptic,
  author    = {Andrea J. Nam and Masaaki Kuwajima and Patrick H. Parker and Jared B. Bowden and Wickliffe C. Abraham and Kristen M. Harris},
  title     = {Perisynaptic astroglial response to in vivo long-term potentiation and concurrent long-term depression in the hippocampal dentate gyrus},
  journal   = {bioRxiv},
  year      = {2025},
  doi       = {10.1101/2025.05.13.653827},
  url       = {https://doi.org/10.1101/2025.05.13.653827}
}
