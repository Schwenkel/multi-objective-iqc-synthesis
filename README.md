# Multi-objective robust controller synthesis with integral quadratic constraints in discrete-time

## Description

This repository contains the MATLAB code that accompanies the paper 
> Lukas Schwenkel, Johannes Köhler, Matthias A. Müller, Carsten W. Scherer, Frank Allgöwer "Multi-objective robust controller synthesis with integral quadratic constraints in discrete-time", 2025, arXiv: .


## Prerequisites

- MATLAB (tested with version R2023b)
- YALMIP (tested with version 22-June-2023)
- MOSEK (tested with version 10.1.21)

## Usage

Run `example1_p2p.m`, `example1_e2p.m`, `example2.m`, and `example3.m` to get the results of Example 1, 2 and 3 in the article. These examples call the functions in `iqc_analysis.m` and `iqc_synthesis.m` that contain the analysis and synthesis algorithms proposed in the article. Integral quadratic constraints are implemented using the class described in `IQC.m`. The synthesis algorithm iterates between the analysis step `ana_step.m` and the synthesis step `syn_step.m`. The factorization from Theorem 4 that is needed in the synthesis step is written in `factorization.m`. To warm start the synthesis algorithm a nominal synthesis using `nom_synthesis.m` is performed. The computation of the lower bounds for Example 2 are computed for each of the three controllers using the function `lower_bound_exmp2.m`. Finally, the results of the examples are also available in `example1_p2p_result.m`, `example1_e2p_result.m`, `example2_result.m`, and `example3_result.m`.

## License

This project is licensed under the MIT License.

## Citation

If you use this code in your research, please cite our work:

```text
@article{Schwenkel2025,
  title={Multi-objective robust controller synthesis with integral quadratic constraints in discrete-time},
  author={L. Schwenkel and J. K{\"o}hler and M. A. M{\"u}ller and C. W. Scherer and F. Allg{\"o}wer}},
  year={2025}
}
```
  
## Support and Contact

For any questions or issues related to this code, please contact the author:

- Lukas Schwenkel schwenkel(at)ist.uni-stuttgart.de
