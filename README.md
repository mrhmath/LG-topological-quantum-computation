# Links--Gould-topological-quantum-computation
**Authors:** Matthew Harper and Ben-Michael Kohli

This repository provides the official Python implementation of the paper \
[An example of universal braiding in Links--Gould topological quantum computation](https://arxiv.org/abs/2611.04290).

## Repository Description

We provide a script for computing the matrix units of the braiding $$\sigma_1$$ and the right Jucys-Murphy element $$\sigma_2^2$$ on the right-associated morphism space $$\mathrm{Hom}_{U_q}(\mathbf{V}, \mathbf{W}^{\otimes 2}\otimes \mathbf{V})$$ where $$\mathbf{V}$$ and $$\mathbf{W}$$ are the Kac modules $$V(0,\alpha)$$ and $$V(0,-\frac{1}{2})$$, respectively.
   
## Scripts

The 'LG_TQC.py' Python program file computes the actions of the generators of the bosonization of the quantum group for any Kac module, as well as the braiding on any tensor product of Kac modules. It outputs the braiding applied to $$V(0,\beta)\otimes V(0,\alpha)$$. For a given $$n$$, it computes explicit highest weight vectors in the tensor product decomposition of $$V(0,\beta)\otimes V(n,\alpha)$$ and the associated inclusion and projection maps. These are then used to compute matrix units for the actions of braid group elements.

## Citation
If you use this script in your research, please cite our preprint.

```bibtex
@misc{LGTQC,
  title = {An example of universal braiding in {L}inks--{G}ould topological quantum computation},
  author = {Harper, Matthew and Kohli, Ben-Michael},
  note = {\href{https://arxiv.org/abs/2611.04290}{arXiv:2611.04290}},
}
```
## Acknowledgments
MH was partially supported through the NSF-RTG grant \#DMS-2135960.
