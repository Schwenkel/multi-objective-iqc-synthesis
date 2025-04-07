# Multi-objective robust controller synthesis with integral quadratic constraints in discrete-time

## Description

This repository contains the MATLAB code that accompanies the paper 
> Lukas Schwenkel, Johannes Köhler, Matthias A. Müller, Carsten W. Scherer, Frank Allgöwer "Multi-objective robust controller synthesis with integral quadratic constraints in discrete-time", 2025, [arXiv:2503.22429](https://arxiv.org/abs/2503.22429).


## Prerequisites

- [MATLAB](https://de.mathworks.com/products/matlab.html) (tested with version R2023b) 
- [YALMIP](https://yalmip.github.io/) (tested with version 22-June-2023)
- [Mosek](https://www.mosek.com/) (tested with version 10.1.21)

## Usage

* Run `example1_p2p.m`, `example1_e2p.m`, `example2.m`, and `example3.m` to get the results of Examples 1, 2, and 3 in the article.
* These examples call the functions in `iqc_analysis.m` and `iqc_synthesis.m` that contain the analysis and synthesis algorithms proposed in the article.
* Integral quadratic constraints are implemented using the class described in `IQC.m`.
* The synthesis algorithm iterates between the analysis step `ana_step.m` and the synthesis step `syn_step.m`.
* The factorization from Theorem 4 that is needed in the synthesis step is written in `factorization.m`.
* A nominal synthesis using `nom_synthesis.m` is performed to warm start the synthesis algorithm.
* The computation of the lower bounds for Example 2 are computed for each of the three controllers using the function `lower_bound_exmp2.m`.
* Finally, the results of the examples are also available in `results/example1_p2p_result.m`, `results/example1_e2p_result.m`, `results/example2_result.m`, and `results/example3_result.m`.

## License

This project is licensed under the MIT License.

## Citation

If you use this code in your research, please cite our work:

```text
@article{Schwenkel2025,
  title={Multi-objective robust controller synthesis with integral quadratic constraints in discrete-time},
  author={L. Schwenkel and J. K{\"o}hler and M. A. M{\"u}ller and C. W. Scherer and F. Allg{\"o}wer}},
  year={2025},
  journal={arxiv:2503.22429},
  doi={10.48550/arXiv.2503.22429},
}
```
  
## Contact

For any questions or issues related to this code, don't hesitate to get in touch with the author:

- Lukas Schwenkel schwenkel(at)ist.uni-stuttgart.de

## References
Some parts of the code implement methods from the following articles. These parts include a reference in the comments.

- V. Ionescu and M. Weiss, On computing the stabilizing solution of the discrete-time Riccati equation, Linear Algebra and its Applications 174 (1992), 229–238
- C. W. Scherer, Dissipativity and integral quadratic constraints: Tailored computational robustness tests for complex interconnections, IEEE Control Systems 42 (2022), no. 3, 115–139
- C. W. Scherer, J. Veenman, Stability analysis by dynamic dissipation inequalities: On merging frequency-domain techniques with time-domain conditions. Systems & Control Letters 121 (2018), 7–15
- C. W. Scherer, S. Weiland, Linear Matrix Inequalities in Control - Lecture Notes, Delft University of Technology (2005), Delft, Netherlands
- L. Schwenkel, J. Köhler, M. A. Müller, and F. Allgöwer, Robust peak-to-peak gain analysis using integral quadratic constraints, Proc. 22nd IFAC World Congress (2023), 11564–11569
- J. Veenman, C. W. Scherer, and H. Köroglu, Robust stability and performance analysis based on integral quadratic constraints, European J. Control 31 (2016), 1–32

