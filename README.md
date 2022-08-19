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



## Reference
Yu, T., Park, Y., Johnson, J. M. & Jones, D. P. apLCMS—adaptive processing of high-resolution LC/MS data. Bioinformatics 25, 1930–1936 (2009). DOI: [10.1093/bioinformatics/btp291](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp291)
