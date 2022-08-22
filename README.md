# recetox-aplcms
![GitHub R package version](https://img.shields.io/github/r-package/v/RECETOX/recetox-aplcms)
[![Conda](https://img.shields.io/conda/v/bioconda/r-recetox-aplcms)](https://anaconda.org/bioconda/r-recetox-aplcms)

This is a fork of the official [aplcms repo](https://github.com/tianwei-yu/apLCMS) that takes the project towards large-scale MS analyses.

## Usage
Install through conda: https://anaconda.org/bioconda/r-recetox-aplcms

Use as a Galaxy Tool: https://github.com/RECETOX/galaxytools/tree/master/tools/recetox_aplcms

## Testing
Before being able to run the tests, it is necessary to fetch the required data using the following commands:

```
wget -P tests/testdata/adjusted -i tests/remote-files/adjusted.txt
wget -P tests/testdata/aligned -i tests/remote-files/aligned.txt
wget -P tests/testdata/extracted -i tests/remote-files/extracted.txt
wget -P tests/testdata/input -i tests/remote-files/input.txt
wget -P tests/testdata/recovered -i tests/remote-files/recovered.txt
wget -P tests/testdata/recovered/recovered-extracted -i tests/remote-files/recovered-extracted.txt
wget -P tests/testdata/recovered/recovered-corrected -i tests/remote-files/recovered-corrected.txt
wget -P tests/testdata/filtered -i tests/remote-files/filtered.txt
wget -P tests/testdata/features -i tests/remote-files/features.txt
wget -P tests/testdata/clusters -i tests/remote-files/clusters.txt
wget -P tests/testdata/hybrid -i tests/remote-files/hybrid.txt
```

The `hybrid` and `unsupervised` tests of recetox-aplcms are [reported](https://github.com/RECETOX/recetox-aplcms/issues/24) to be OS specific and may fail depending on the platrform they are run on. To ensure reproducibility during development process you can run the tests in a designated Docker container as follows:
```
# from the repository root run
$ docker build -t recetox-aplcms .
```
After `docker-build` has built the image run:
```
$ docker run --rm -t -v $(pwd):/usr/src/recetox-aplcms recetox-aplcms
```
This will create a container and automatically run all the tests from the **tests** folder.

# Documentation for developers

## Setting up a develeopment environment
The development environment can be set up in two ways, either via **VSCode's devcontainer** extension or a **docker container**.

### Devcontainer
To use a devcontainer you need VSCode with [Remote - Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension and docker installed on your machine:
- Clone your fork of the repository and open the folder in VSCode;
- From VSCode's command pallete run `Remote-Containers: Open Folder in Container`. VSCode may take a few minutes building a container;
- After container is ready, open a **new** terminal and type `conda activate recetox-aplcms` to activate Conda environment;
- Run `R` or `radian` to enter R terminal (we recommend `radian` due to its ease of use);
- A good starting point would be fetching the test data as described above, running `devtools::test()` and waiting until all tests pass to ensure the environment is set correctly.

### Docker container
To use a docker development environment you need **Docker** installed on your machine. If you don't have **Docker** you can follow installation instructions on Docker's [web](https://docs.docker.com/engine/install/).
- Clone your fork of the repository;
- From the package root folder run `docker build -t recetox-aplcms .` to build an image. This may take a few minutes.
- After the image is build start the container:
    ```bash
    $ docker run -it \
    -v $(pwd):/usr/src/recetox-aplcms \
    --entrypoint '/bin/bash' \
    recetox-aplcms
    ```
- Once in container, finish setting up the environment by running:
    ```bash
    $ apt update && apt upgrade
    ```
    ```shell
    $ apt install git && git config --global --add safe.directory /usr/src/recetox-aplcms
    ```
- Enter a Conda environment by running `conda activate recetox-aplcms`
- Run `R` or `radian` to enter R terminal (we recommend `radian` due to its ease of use);
- A good starting point would be fetching the test data as described above, running `devtools::test()` and waiting until all tests pass to ensure the environment is set correctly.


## Reference
Yu, T., Park, Y., Johnson, J. M. & Jones, D. P. apLCMS—adaptive processing of high-resolution LC/MS data. Bioinformatics 25, 1930–1936 (2009). DOI: [10.1093/bioinformatics/btp291](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp291)
