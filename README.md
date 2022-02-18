# recetox-aplcms
![GitHub R package version](https://img.shields.io/github/r-package/v/RECETOX/recetox-aplcms)
[![Conda](https://img.shields.io/conda/v/bioconda/r-recetox-aplcms)](https://anaconda.org/bioconda/r-recetox-aplcms)

This is a fork of the official [aplcms repo](https://github.com/tianwei-yu/apLCMS) that takes the project towards large-scale MS analyses.

## Usage
Install through conda: https://anaconda.org/bioconda/r-recetox-aplcms

Use as a Galaxy Tool: https://github.com/RECETOX/galaxytools/tree/master/tools/recetox_aplcms

## Testing
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



## Reference
Yu, T., Park, Y., Johnson, J. M. & Jones, D. P. apLCMS—adaptive processing of high-resolution LC/MS data. Bioinformatics 25, 1930–1936 (2009). DOI: [10.1093/bioinformatics/btp291](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp291)
