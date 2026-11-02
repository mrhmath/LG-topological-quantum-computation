# Links--Gould-topological-quantum-computation
**Authors:** Matthew Harper and Ben-Michael Kohli

This repository provides the official Python implementation of the paper \
[An example of Links--Gould topological quantum computation](https://arxiv.org/abs/).

## Repository Description

We provide scripts for computing the matrix units of the braiding $$\sigma_1$$ and the right Jucys-Murphy element $$\sigma_1^2$$ on the right-associated morphism space $$\mathrm{Hom}_{U_q}(\mathbf{V}, \mathbf{W}^{\otimes 2}\otimes \mathbf{V})$$ where $$\mathbf{V}$$ and $$\mathbf{W}$$ are the Kac modules $$V(0,\alpha)$$ and $$V(0,-\frac{1}{2})$$, respectively.
   
## Scripts

The 'LG_TQC.py' Python program file computes the actions of the generators of the bosonization of the quantum group for any Kac module, as well as the braiding on any tensor product. For a given $$n$$, it computes explicit highest weight vectors in the tensor product decomposition of $$V(0,\beta)\otimes V(n,\alpha)$$ and the associated inclusion and projection maps. These are then used to compute matrix units for the actions of braid group elements.

## Citation
If you use this script in your research, please cite our preprint.

```bibtex
@misc{LGTQC,
  title = {An example of universal braiding in {L}inks--{G}ould topological quantum computation},
  author = {Harper, Matthew and Kohli, Ben-Michael},
  note = {\href{https://arxiv.org/abs/}{arXiv:}},
}
```
## Acknowledgments
MH thanks the organizers of the 2026 Moab Topology Conference and to Sung Kim whose talk made him aware of the works on braiding universality in an ADO category, which inspired this work. MH was partially supported through the NSF-RTG grant \#DMS-2135960.
