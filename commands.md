# MRMkit pipeline
## DDA quantitation and library generation
```
cd {path_to_folder_with_mzML_files,_param.txt,_assay_info_and_batch_info}

# internal standards (ISTDs) peak extraction
python3 {path_to_MRMkit}/MRMistd.py

# peak detection
python3 {path_to_MRMkit}/MRMgetpeak.py

# draw ion chromatograms
Rscript {path_to_MRMkit}/MRMionc.r

# batch correction
python3 {path_to_MRMkit}/MRMcorrect.py

# RT selection (optional)
# users to input RT of true peaks in user_RT.txt
python3 {path_to_MRMkit}/MRM_RT.py
```



