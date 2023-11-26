## Fairness-driven Private Collaborative Machine Learning

### About
This repository supports the following paper:
> Pessach, D., Tassa, T., & Shmueli, E. (2021). Fairness-driven private collaborative machine learning. arXiv preprint arXiv:2109.14376. ([https://doi.org/10.48550/arXiv.2109.14376](https://doi.org/10.48550/arXiv.2109.14376)).

This repository implements a privacy-preserving pre-process mechanism for enhancing fairness of collaborative Machine Learning algorithms. In particular, the technique improves fairness by decreasing distances between the distributions of attributes of the privileged and unprivileged groups. It uses a binning approach that enables the implementation of privacy-preserving enhancements, by means of Multi-Party Computation.

### Usage
To run the experiment run the notebook file [Main.ipynb](https://github.com/Danapessach/MultiDistExperiment/blob/main/Main.ipynb).

### Datasets
The datasets utilized for the experiments can be found in the [datasets](https://github.com/Danapessach/MultiDistExperiment/tree/main/MyDistExperiment/datasets) folder.
(The ProPublica Recidivism dataset, the ProPublica Violent Recidivism dataset, and the Bank Marketing dataset).

The experiment is assisted by the pre-processed dataset files introduced by:
> Friedler, S. A., Scheidegger, C., Venkatasubramanian, S., Choudhary, S., Hamilton, E. P., & Roth, D. (2019, January). A comparative study of fairness-enhancing interventions in machine learning. In Proceedings of the conference on fairness, accountability, and transparency (pp. 329-338).‏
([https://dl.acm.org/doi/pdf/10.1145/3287560.3287589](https://dl.acm.org/doi/pdf/10.1145/3287560.3287589)).

### Packages
The experiment is implemented in Python assisted with the Virtual Ideal Functionality Framework (VIFF) library for secure multi-party computations, based on:
> Damgård, I., Geisler, M., Krøigaard, M., & Nielsen, J. B. (2009, March). Asynchronous multiparty computation: Theory and implementation. In International workshop on public key cryptography (pp. 160-179). Berlin, Heidelberg: Springer Berlin Heidelberg.‏
([https://link.springer.com/chapter/10.1007/978-3-642-00468-1_10](https://link.springer.com/chapter/10.1007/978-3-642-00468-1_10)).

### Reference
If you find the code useful, please cite our paper.

<pre>@article{pessach2021fairness,
  title={Fairness-driven private collaborative machine learning},
  author={Pessach, Dana and Tassa, Tamir and Shmueli, Erez},  
  journal={arXiv preprint arXiv:2109.14376},
  year={2021}  
}
</pre>

2023/11/26
