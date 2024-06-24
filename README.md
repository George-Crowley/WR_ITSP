<!-- Useful github stuff 
https://hackmd.io/@vivek-blog/github_article
https://viveikjha.github.io
https://jekyllrb.com
https://github.com/ManufacturingInformatics/marl_fixture_planner/blob/main/README.md?plain=1
-->

# Information-Theoretic Sensor Placement for Large-Scale Sewer Networks

This is the repository to go along with the paper ["Information-Theoretic Sensor Placement for Large-Scale Sewer Networks"](link-here). This paper provides an information-theoretic approach to sensor selection in wastewater collection networks and uses estimation techniques to validate the approach.

The repository contains the code to reproduce the results from Case Study 2 (CS2) in the aforementioned paper. We note that for the estimation simulations for CS2, the general regression neural network (GRNN) results were simulated on the Sheffield University high-performance computer due to the RAM required for the amount of training data we used.  

If you want to replicate CS2 results, the code and relevant data sets will be provided. We will also provide step-by-step documentation in another file which will walk you through our MATLAB code.

## Contents

- [Getting Started - datasets](#1)
- [Description of source code](#8)
- [Workflow](#2)
- [Dependencies](#3)
- [Installing](#4)
- [Executing program](#5)
- [Mutual Information calculations](#6)
- [Reproducing figures without running the simulations and configuration files](#9)
- [Authors](#7)
  
<!-- - [Brief Synopsis](#2)
- [Manual Installation](#3)
  - [MATLAB Runtime Installation](#3a)
  - [Package Installation](#3b)
  - [Training](#3c)
  - [Inference](#3d)
- [Docker Container](#4)
  - [Setup](#4a)
  - [Usage](#4b)
- [Troubleshooting](#5)
  - [MATLAB Runtime Issues](#5a)
  - [Docker Execution Issues](#5b)
- [Citing This Work](#8) -->

<!-- <a id='1'></a> -->

## Description

* This repo can serve several purposes. The first is to replicate the results as in the paper. The second is to be able to use these coded simulations to run on your own networks.
<a id='1'></a>
## Getting Started - datasets

* To get started, you will need to download the datasets. Note that for both volumetric flow and water level, there are two datasets each 1gb approximately in size saved in a `.mat` file.
* To download the flow datasets, visit: `https://zenodo.org/doi/10.5281/zenodo.11442174`.
* To download the water level dataset, visit: `https://zenodo.org/doi/10.5281/zenodo.11636092`.

* Once the files are downloaded, download the source code files from this repo and move the files into the Data folder or any relevant folder that is added to the Matlab path.
<a id='8'></a>
## Description of source code

* There are three folders to download. The first is the Code file, which contains the scripts to be run. The Data folder currently contains the random sensor selections the paper uses in its estimation, the Rule based sensor selection vector, and the shape files for Bellinge (CS2) saved in a structured array in a `.mat` file. To download the datasets, see [Getting Started - datasets](#1). The last folder is the Step-by-step code instructions, which include extra steps (not commented out) on what is going on in the code. 

<a id='2'></a>
## Workflow  

* The workflow of this code is as follows: there are 7 `.m` files.
### Sensor Placement scripts
* `Run_code.m` This houses the script that runs and combines all the other matlab scripts.
* `MI_comparison_Bellinge.m` This file runs mutual information calculations for the heuristic sensor selections that aren't the one-step modified greedy algorithm.
* `sensor_selection.m` This file houses the one-step modified greedy algorithm proposed in the paper.
### Estimation scripts 
* `GLM_estimation1.m` This file houses the general linear model estimation method where different sensor selections are inputted into the file etc.
* `GRNNET_estimation1.m` This file houses the general regression neural network estimation method where different sensor selections are inputted into the file etc.
* `error_calculations.m` This file contains the script for calculating the error matrices and normalized mean square error and is used at the end of both the GLM and GRNN scripts.
* `Run_estimation.m` For different inputted numbers of sensors, this script runs the estimation for each of the different sensor placements mentioned in the paper.

### Running order
* In `run_code.m`, the relevant datasets will be loaded (choose which you would like to use i.e. flow or water level). Then, the relevant heuristic sensor selections will be loaded - to reproduce the results as in the paper for CS2, the rule based sensor selection for Bellinge will also be loaded. Then, the proposed Algorithm will be run in `sensor_selection.m`. A plot of the results will follow. Note that depending on the size of the network and how many sensors you want to implement in the network this could take some time, so we would suggest starting with a small number and building up.

* Then, the plot of mutual information comparing the different sensor selections is plotted.

* Following that, we then look at the estimation results for all of the sensor selections for different numbers of sensors. We call this vector for the number of sensors we are going to use `k_list_vector`. In the paper, this vector has the values `25,50,...,200,225,250`, and will be pre-loaded.

* Following that, a loop is created which runs the file `Run_estimation.m`. In this file, all of the relevant sensor selections are included alongside the relevant estimation techniques we are going to use. The loop runs for all values in the vector `k_list_vector`. For each sensor selection heuristic and our own proposed, the GLM (general linear model) in the file `GLM_estimation1.m` and GRNN (general regression neural network) in the file `GRNNET_estimation1.m` are called, and input the relevant datasets and modify them according to the sensor selection inputted into the function. In both the estimation files, once the estimation is done, the `error_calculations.m` file is then called, which calculates the error matrix between true and estimated values, and the NMSE (normalized mean square error).

* Once this is done, all of the results are saved in each loop and saved in the folders created in the first few lines in `run_code.m`. The final remains of the code show the results obtained in the remainder of the paper.

### NOTE: if you try to run the `run_estimation.m` in the `Run_code.m` file with the current data files, the GRNN will not run and the file will error due to the amount of RAM needed to run the algorithm. We used The University of Sheffield high-performance computer to run the estimation simulations. To run smaller simulations, before the run_estimation.m code is run, change the number of data points contained in the time series for the training dataset. This is added into the dataset (we take the first 4000 data points of both the training and validation, but you can amend as necessary.)
<a id='3'></a>
### Dependencies

* All simulations were run on MATLAB 2023a.
* Packages needed:
  - Statistics and Machine Learning Toolbox
  - Deep Learning Toolbox
  - Mapping Toolbox
<a id='4'></a>
### Installing

* To download MATLAB: `https://uk.mathworks.com/`
* To download/clone the repository, please follow these steps outlined in ["Cloning a Github repo"](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).
* Once the folder has been downloaded, add it to the path in Matlab and make sure the downloaded data sets are also added to the path.
<a id='5'></a>
### Executing program

* Download relevant datasets in the Getting Started section.
* See step-by-step folder for explanations of running of code. Otherwise, download the relevant functions and files, add everything to the MATLAB path, and run the `run_code.m` file.

  

<a id='6'></a>
## Mutual Information calculations

In Theorem 3 of the paper, we state that

$$
\begin{align} 
        & I(X^n; {\bf{H}} X^n +  Z^k) =  \dfrac{1}{2} \log \left( \dfrac{1}{\sigma^{2k}} \hbox{det} \left( {\bf{H}} \hbox{$\mathbf{\Sigma}$} {\bf{H}}^{{\sf T}} + \sigma^2 {{\bf I}}_k \right) \right).  
\end{align}
$$

By using the usual log laws, we obtain

$$
\begin{align*}
I(X^n; {\bf{H}} X^n +  Z^k) = \dfrac{1}{2} \log \hbox{det} \left( {\bf{H}} \hbox{$\mathbf{\Sigma}$} {\bf{H}}^{{\sf T}} + \sigma^2 {{\bf I}}_k \right) - \dfrac{k}{2} \log \sigma^2. 
\end{align*}
$$

Now as the matrix inside the determinant operator grows as we select more sensor locations, dependent on the state measurements we are interested in, the eigenvalues of the matrix will tend towards 0 and give an undefined value for our calculation. To mitigate this issue in the numerical simulation, we will take advantage of some determinant properties. Recall that for a matrix ${\bf{A}} \in \mathbb{R}^{k \times k}$ and scalar $\gamma \in \mathbb{R}$, we have that

$$
\begin{align}
\hbox{det}\left(\gamma {\bf{A}} \right) = \gamma^k \hbox{det}({\bf{A}}) \implies \log \hbox{det}\left({\bf{A}}\right) = \log \hbox{det}\left(\gamma {\bf{A}} \right) - k \log \gamma.
\end{align}
$$

This substitution will allow for larger matrix determinant calculations for appropriately chosen $\gamma$. This is implemented in the `MI_comparion_Bellinge.m` and `sensor_selection.m` code with 

$$
{\bf{A}} = {\bf{H}} \hbox{$\mathbf{\Sigma}$} {\bf{H}}^{{\sf T}} + \sigma^2 {{\bf I}}_k.
$$

<a id='9'></a>
## Reproducing figures without running the simulations

<a id='7'></a>
## Authors
<!-- Contributors names and contact info -->

ex. George Crowley  
<!-- ex. [@DomPizzie](https://twitter.com/dompizzie) -->

<!--
## Citing This Work

To cite this work, please refer to [Paper name](Paper_link) and cite us using the format below:

```bibtex
 @article{, 
    type={Journal}, 
    title={}, 
    DOI={.....}, 
    publisher={.....}, 
    author={...}, 
    year={2024}, 
    month={....}, 
    language={en} 
 }
```
-->

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the [BSD 3-Clause] License - see the LICENSE.md file for details
<!--
## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46) -->
