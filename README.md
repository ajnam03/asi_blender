### asi_blender  
- [Introduction](#Introduction)
- [Set-up](#Set-up)
- [Cite](#Cite)
#### Introduction
This repository contains code to automatically identify the axon–spine interface (ASI) using 3D reconstructions of presynaptic axons and postsynaptic dendritic spines, and to quantify the degree of astroglial apposition at the ASI. Read more in our [preprint]([https://doi.org/10.1101/2025.05.13.653827])

#### 🔧 Set-up
1. **Blender version**: This script currently works with **Blender 3.6.5**.  
2. **LoopTools Add-on**: Ensure the LoopTools add-on is activated in Blender.  
   ➤ [LoopTools Documentation](https://docs.blender.org/manual/en/3.5/addons/mesh/looptools.html)  
3. **Naming conventions**:
   - `d##c##`: **c-object**
   - `d##sp##`: **spine object**
   - `d##ax##`: **axon object**
   - `astroAll`: **astroglia object**
   - `1um_rad_ref`: **1 µm reference circle**

#### ▶️ How to Run
1. Edit the path to your `.blend` file in `blender_asi_analysis.sh`.
2. Set the path to the Blender executable.
3. Make the script executable:
   ```bash
   chmod u+x blender_asi_analysis.sh

### Cite
## 📖 Cite this work

If you use this code or data, please cite our preprint:

> Nam, A. J., Kuwajima, M., Parker, P. H., Bowden, J. B., Abraham, W. C., & Harris, K. M. (2025). Perisynaptic astroglial response to in vivo long-term potentiation and concurrent long-term depression in the hippocampal dentate gyrus. bioRxiv. https://doi.org/10.1101/2025.05.13.653827

<details>
<summary>BibTeX</summary>

```bibtex
@article{nam2025perisynaptic,
  author    = {Andrea J. Nam and Masaaki Kuwajima and Patrick H. Parker and Jared B. Bowden and Wickliffe C. Abraham and Kristen M. Harris},
  title     = {Perisynaptic astroglial response to in vivo long-term potentiation and concurrent long-term depression in the hippocampal dentate gyrus},
  journal   = {bioRxiv},
  year      = {2025},
  doi       = {10.1101/2025.05.13.653827},
  url       = {https://doi.org/10.1101/2025.05.13.653827}
}
