<div style="text-align:center"><img style="width:60%; height: auto" src="https://github.com/yanailab/enhance/raw/master/images/splash.jpg"/></div>

## ENHANCE: Accurate denoising of single-cell RNA-Seq data

This repository contains a Python implementation of the ENHANCE algorithm for denoising single-cell RNA-Seq data ([Wagner et al., 2019](https://www.biorxiv.org/content/10.1101/655365v1)).

The **R implementation** can be found in a [separate repository](https://github.com/yanailab/enhance-R).

### Running ENHANCE from the command-line

Follow these instructions to run the Python implementation of ENHANCE from the command-line.

1. Install dependencies

   Make sure you have Python 3 and the Python packages `scikit-learn`, `pandas`, and `click` installed. The easiest way to install Python 3 as well as these packages is to download and install [Anaconda](https://www.anaconda.com/distribution/) (select "Python 3.7 version").

2. Download the GitHub repository

   [Download ENHANCE](https://github.com/yanailab/enhance/archive/master.zip), and extract the contents into a folder.

3. Test running the script

   To run the script, change into the folder where you extracted the files, and run (on Linux/Mac):
    
   ``` bash
   python3 enhance.py --help
   ```

   You should see the following output:

    ```
    Usage: enhance.py [OPTIONS]

    Options:
    -f, --fpath TEXT            The input UMI-count matrix.
    -o, --saveto TEXT           The output matrix.
    --transcript-count INTEGER  The target median transcript count for
                                determining thenumber of neighbors to use for
                                aggregation.(Ignored if "--num-neighbors" is
                                specified.)
    --max-neighbor-frac FLOAT   The maximum number of neighbors to use for
                                aggregation, relative to the total number of
                                cells in the dataset. (Ignored if "--num-
                                neighbors" is specified.)
    --pc-var-fold-thresh FLOAT  The fold difference in variance required for
                                relevant PCs, relative to the variance of the
                                first PC of a simulated dataset containing only
                                noise.
    --max-components INTEGER    The maximum number of principal components to
                                use.
    --num-neighbors INTEGER     The number of neighbors to use for aggregation.
    --sep TEXT                  Separator used in input file. The output file
                                will use this separator as well.  [default: \t]
    --use-double-precision      Whether to use double-precision floating point
                                format. (This doubles the amount of memory
                                required.)
    --seed INTEGER              Seed for pseudo-random number generator.
                                [default: 0]
    --test                      Test if results for test data are correct.
    --help                      Show this message and exit.
    ```

4. Make sure your expression matrix file is formatted correctly

   By default, the script expects your expression matrix to be stored as a tab-separated plain-text file, with gene labels contained in the first column, and cell labels contained in the first row (the top-left "cell" in the matrix can either be empty or contain the first cell label). You can compress the text file using gzip. A properly formatted example dataset (`data/pbmc-4k_expression.tsv.gz`) is included in this repository.

   If your file uses a separator other than the tab character, you must specify it by passing the `--sep` argument to the script. For example, if you're using comma-separated values (csv), pass `--sep ,`.  This will also affect the separator used in the output file.

5. Run ENHANCE!

   Let's say your (tab-separated) expression matrix file is called `expression.tsv.gz`, and you saved it in the same directory as the `enhance.py` script. Then, to ENHANCE, you would use:

   ``` bash
   python3 enhance.py -f expression.tsv.gz -o denoised_expression.tsv
   ```

   This will produce a denoised matrix called `denoised_expression.tsv`.


### Example

  Running ENHANCE from the command-line, on the test dataset included
  in this repository (`data/pbmc-4k_expression.tsv.gz`):

  ``` bash
  $ python3 enhance.py -f data/pbmc-4k_expression.tsv.gz -o denoised_pbmc-4k_expression.tsv
  ```

  Output:
  ```
	Loading the data... done. (Took 0.1 s.)
	The expression matrix contains 7145 genes and 100 cells.

	Performing kNN-smoothing with k=32, d=2, and dither=0.030...
	Step 1/5: Smooth using k=2
		PCA took 0.1 s.
		The fraction of variance explained by the top 2 PCs is 4.6 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 0.0 s.
	Step 2/5: Smooth using k=4
		PCA took 0.0 s.
		The fraction of variance explained by the top 2 PCs is 8.0 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 0.0 s.
	Step 3/5: Smooth using k=8
		PCA took 0.0 s.
		The fraction of variance explained by the top 2 PCs is 14.3 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 0.0 s.
	Step 4/5: Smooth using k=16
		PCA took 0.0 s.
		The fraction of variance explained by the top 2 PCs is 25.4 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 0.0 s.
	Step 5/5: Smooth using k=32
		PCA took 0.0 s.
		The fraction of variance explained by the top 2 PCs is 52.4 %.
		Calculating pair-wise distance matrix took 0.0 s.
		Calculating the smoothed expression matrix took 0.1 s.
	kNN-smoothing finished in 0.3 s.

	Writing results to "test_smoothing_smoothed.tsv"... done. (Took 0.5 s.)
  ```

The results are shown using UMAP below. (The UMAP result is included in this repository, under `data/umap_result.tsv`).

<div style="text-align:center"><img style="width:60%; height: auto" src="https://github.com/yanailab/enhance/raw/master/images/github_example.png"/></div>

### Changelog

We will note all changes to the code here.

#### 6/3/2019 - Python implementation of ENHANCE released
This is the first release of ENHANCE (algorithm version "0.1"), which is the version we used to generate the results presented in our bioRxiv paper, including all benchmark analyses.
