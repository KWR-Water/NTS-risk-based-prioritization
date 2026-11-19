# Scripts

This repository contains R-scripts and Python-scripts, which can be used to execute the workflow as described in the manuscript. Required R- and Python versions, packages and libraries are described in the scripts. 

## Workflow

1. 1_patRoon_feature_finding.R -> perform feature parameter optimization and find and group features.
2. 2_patRoon_QC_check.R -> QC check of spiked samples
3. 3_patRoon_NTS_and_prioritization.R -> feature finding using patRoon and prioritization using MS2Tox, MS2Quant and structural alerts
4. 4_MS2bioactivity_nrahr.iypnb -> script to run MS2Tox nr.ahr activity model
5. 5_rearrange_SIRIUS_output_for_MLinvitroTox.R -> re-arrange SIRIUS output generated through patRoon, so that it is compatible with MLinvitroTox

## Dependencies
This workflow depends on external software, including:
- [patRoon](https://rickhelmus.github.io/patRoon/) 
- [MLinvitroTox](https://github.com/LiliGasser/MLinvitroTox-tutorial/tree/main)
- [deepFPlearn](https://github.com/yigbt/deepFPlearn)
- [SIRIUS](https://bright-giant.com/sirius-features/)
- [MS2Tox nr.ahr](https://github.com/RaHub4AI/MS2bioactivity/tree/SIRIUS5)
