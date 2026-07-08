# MultiCMPSKit.jl

Research code for the numerical study of boson mixtures with multi-component continuous matrix product states.

This repository contains code associated with the paper:

> Wei Tang, Benoît Tuybens, and Jutho Haegeman,
> *Numerical study of boson mixtures with multi-component continuous matrix product states*,
> arXiv:2512.24998.

The code implements routines for variational calculations with multi-component continuous matrix product states (cMPS), with a focus on bosonic mixtures and the regularity constraints required for finite kinetic energy.

## Relation to CMPSKit.jl

This repository is derived from [`CMPSKit.jl`](https://github.com/Jutho/CMPSKit.jl), a Julia package for variational simulations with continuous matrix product states.

At the time of publication, the multi-component functionality contained here has not yet been merged into the main `CMPSKit.jl` repository. This repository is therefore provided separately to make the code used for the paper visible, linkable, and reproducible.

The intention is that the relevant functionality may later be integrated into `CMPSKit.jl`. This repository should therefore be viewed as a research-code snapshot associated with the paper, rather than as a replacement for the upstream package.

## Main features

* Multi-component cMPS ansatz for one-dimensional bosonic quantum field theories.
* Regularity-preserving parametrisations for bosonic mixtures.
* Variational optimisation routines for multi-component cMPS.
* Numerical tools used for the two-component Lieb-Liniger calculations.
* Routines for extracting observables such as densities, energies, correlation data, and low-energy quantities.
* Scripts associated with the numerical results reported in the paper.

## Installation

This repository is written in Julia.

Clone the repository:

```bash
git clone https://github.com/BenoitTuybens/MultiCMPSKit.jl.git
cd MultiCMPSKit.jl
```

Start Julia in the repository directory and instantiate the environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will install the dependencies specified in `Project.toml` and `Manifest.toml`.

## Data availability

The numerical data associated with the paper are archived separately on Zenodo.

Data DOI:

```text
https://doi.org/10.5281/zenodo.18109141
```

The repository contains the code used to generate and analyse the data. The Zenodo archive contains the corresponding data files used for the results reported in the paper.

## Repository structure

```text
src/            Core implementation, built on the devnew branch of CMPSKit.jl
examples/       Examples and usage demonstrations
test/           Tests and consistency checks
Project.toml    Julia environment
Manifest.toml   Exact dependency versions
README.md       This file
```

The exact structure may evolve as parts of the implementation are prepared for possible integration into `CMPSKit.jl`.

## Citation

If you use this code, please cite the accompanying paper:

```bibtex
@misc{TangTuybensHaegeman2025,
  author        = {Tang, Wei and Tuybens, Beno\^{i}t and Haegeman, Jutho},
  title         = {Numerical study of boson mixtures with multi-component continuous matrix product states},
  howpublished  = {arXiv:2512.24998},
  year          = {2025},
  eprint        = {2512.24998},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.quant-gas},
  doi           = {10.48550/arXiv.2512.24998}
}
```

The article is currently available as an arXiv preprint. Once the article is published in a journal, the recommended citation will be updated accordingly.

Please also cite the original `CMPSKit.jl` / cMPS optimisation work where appropriate.

## Requirements

The code requires Julia and the package dependencies listed in `Project.toml`.

The calculations may depend on packages from the Julia tensor-network ecosystem, including packages developed in the group of Jutho Haegeman and collaborators. Exact dependency versions are specified in `Manifest.toml`.

## License

This repository follows the license of CMPSKit.jl, unless stated otherwise.

See LICENSE for details.

## Acknowledgements

This code builds on `CMPSKit.jl` and related Julia tensor-network tools developed by Jutho Haegeman.

The scientific work associated with this repository was carried out by Wei Tang, Benoît Tuybens, and Jutho Haegeman.

## Contact

For questions about the paper or this code, please contact:

* Wei Tang: wei.tang.phys@gmail.com
* Benoît Tuybens: benoit.tuybens@gmail.com
* Jutho Haegeman: jutho.haegeman@ugent.be
