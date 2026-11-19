# deepFPlearn prediction of AR, ER, and ED endpoints

This repository contains a Jupyter Notebook used to predict androgen receptor (AR), estrogen receptor (ER), and endocrine disruption (ED)-related activities using pretrained deepFPlearn models.

## Workflow

1. Extract SMILES from the input Excel file.
2. Predict AR, ER, and ED activities using five pretrained models per endpoint.
3. Calculate the mean and standard deviation of predictions across the five models.
4. Match the predictions to the original identifiers using SMILES.
5. Export the combined results as CSV files.

## Input

The Notebook expects an Excel file named `input.xlsx`, containing a sheet named `general` with the following columns:

- `group`: Original identifier
- `SMILES.sirius`: SMILES assigned by SIRIUS
- `score.sirius`: SIRIUS score

The input dataset is not included in this repository.

## Dependencies

This workflow uses pretrained models from [deepFPlearn](https://github.com/yigbt/deepFPlearn). The original software, pretrained models, and required configuration files must be obtained separately.

## Output

The combined prediction results are saved to:

`results_predictgnn/ensemble/ALL_endpoints_ensemble_with_id.csv`