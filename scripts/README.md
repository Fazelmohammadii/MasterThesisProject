
# Scripts Directory

This directory contains the R scripts used throughout the thesis project on investigating unique mutations in DNMT3A and TET2 genes within clonal hematopoiesis (CH) and acute myeloid leukemia (AML) populations. Each script serves a specific purpose in the data processing, analysis, and modeling pipeline.

## Structure

```plaintext
scripts/
├── PreProcessing.R
├── MutationalSpectrum.R
├── MutationalSignatures.R
├── FeatureEngineering.R
├── ModelTraining.R
└── README.md
```

### Scripts

- **PreProcessing.R**: Preprocesses the raw mutation data, which are already mixed and harmonized.
  - **Usage**: Prepares the dataset for subsequent analysis by cleaning and organizing the mutation data.

- **MutationalSpectrum.R**: Creates the mutational spectrum and performs statistical analysis of the mutations in CH and AML conditions.
  - **Usage**: Generates mutational spectrum plots and statistical summaries to compare CH and AML mutations.

- **MutationalSignatures.R**: Analyzes the mutational signatures based on the COSMIC catalogues.
  - **Usage**: Identifies and compares mutational signatures in the dataset using established COSMIC signatures.

- **FeatureEngineering.R**: Creates new features and prepares the dataset for machine learning model training.
  - **Usage**: Engineers features from the mutation data to be used as input for machine learning models.

- **ModelTraining.R**: Trains machine learning models and tests them for mutation classification.
  - **Usage**: Implements and evaluates machine learning models to classify mutations as CH-like or AML-like.

## Usage

Each script is designed to be run independently. Below are general usage instructions for running the scripts:

### PreProcessing

To preprocess the raw mutation data:

```r
Rscript scripts/PreProcessing.R
```

### Mutational Spectrum Analysis

To create the mutational spectrum and perform statistical analysis:

```r
Rscript scripts/MutationalSpectrum.R
```

### Mutational Signatures Analysis

To analyze mutational signatures:

```r
Rscript scripts/MutationalSignatures.R
```

### Feature Engineering

To create new features and prepare the dataset for machine learning:

```r
Rscript scripts/FeatureEngineering.R
```

### Model Training

To train and test the machine learning models for mutation classification:

```r
Rscript scripts/ModelTraining.R
```

## Notes

- Ensure that all dependencies are installed before running the scripts. Refer to `requirements.txt` for a list of required R packages.
- Review each script's header comments for detailed usage instructions and parameters.
- The scripts are designed to be run in sequence, starting with `PreProcessing.R` and ending with `ModelTraining.R`.

This directory is an integral part of the thesis project, providing the necessary scripts to preprocess data, analyze mutations, engineer features, and develop machine learning models.
