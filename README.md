# DICOM to NIfTI Conversion for Breast MRI

A Python pipeline that converts raw DICOM medical image slices into organized 3D and 4D NIfTI volumes ready for research analysis. Built and tested on the QIN-BREAST-02 collection from The Cancer Imaging Archive (TCIA), a publicly available Dynamic Contrast-Enhanced (DCE) breast MRI dataset.

This is the second project in a learning portfolio focused on medical imaging pipeline development. It builds on a previous metadata extraction pipeline by moving from inspection to transformation: converting clinical-format data into the research-format expected by tools like 3D Slicer, ITK-SNAP, LIFEx, MONAI, and PyRadiomics.

---

## Why This Project Exists

Medical scanners output images in DICOM format. Each slice of a scan is saved as a separate `.dcm` file with rich clinical metadata attached. While this format is excellent for hospital systems and clinical workflows, it is impractical for research.

A single DCE breast MRI study can produce 2,000 to 3,000 individual DICOM files, scattered across nested folders by patient, study, and series. Before any meaningful analysis can begin, those slices must be:

1. Grouped together by acquisition sequence
2. Sorted by physical position in space
3. Stacked into proper 3D or 4D volumes
4. Saved in a format the analysis tools understand

NIfTI is the standard format in medical imaging research. One NIfTI file represents one full 3D volume (or 4D in the case of dynamic studies), with spatial geometry preserved in a compact, self-contained file.

This pipeline automates that conversion.

---

## Dataset

**QIN-BREAST-02**, from The Cancer Imaging Archive

| Property | Value |
|---|---|
| Modality | Dynamic Contrast-Enhanced MRI |
| Subjects | 13 patients |
| Studies | 34 imaging visits |
| Series | 235 acquisition sequences |
| Total DICOM files | 31,790 |
| Scanner | Philips Achieva 3T |
| License | CC BY 4.0 |
| Source | [TCIA Collection Page](https://www.cancerimagingarchive.net/collection/qin-breast-02/) |

The dataset contains a mix of sequence types: DCE-MRI with multiple temporal positions per series, diffusion weighted imaging, T1 mapping, quantitative magnetization transfer, and field maps. Each requires either 3D or 4D representation depending on whether temporal information is present.

---

## What This Pipeline Does

For every series folder in the dataset, the pipeline performs the following steps:

1. Reads all DICOM files in the folder using `pydicom`
2. Groups slices by `TemporalPositionIdentifier`, which is relevant for DCE series with multiple time points
3. Sorts each group by physical Z position using `ImagePositionPatient`
4. Stacks slices into 3D volumes with correct voxel spacing, origin, and direction cosines
5. For multi-timepoint series, joins individual 3D volumes into a single 4D volume
6. Saves the result as a compressed NIfTI file (`.nii.gz`)

The output structure mirrors the patient hierarchy and uses meaningful filenames based on series number and sequence description.

---

## Input and Output Structure

**Input** (DICOM hierarchy from TCIA):

```
qin_breast_02/
└── PatientID/
    └── StudyInstanceUID/
        └── SeriesInstanceUID/
            ├── slice_001.dcm
            ├── slice_002.dcm
            └── ...
```

**Output** (NIfTI volumes):

```
nifti_output/
├── QIN-BREAST-02-0001/
│   ├── 8101_THRIVE_SENSE.nii.gz
│   ├── 7401_DWI_EPI_b0200800.nii.gz
│   ├── 8401_qMT_3D_Pulsed_Trig.nii.gz
│   └── ...
├── QIN-BREAST-02-0002/
│   └── ...
...
└── QIN-BREAST-02-0013/
```

For this dataset, the conversion produces:

- 31,790 DICOM slices, reduced to
- 170 NIfTI volumes (99 are 3D, 71 are 4D)
- approximately 13 sequences per patient

---

## Installation

```bash
git clone https://github.com/<your-username>/MRI_breast_dicom_to_nifti.git
cd MRI_breast_dicom_to_nifti
pip install -r requirements.txt
```

Requirements:

```
pydicom>=2.4.0
SimpleITK>=2.3.0
numpy>=1.24.0
```

---

## Usage

The default data and output paths are hardcoded at the top of the script. Most users only need to run:

```bash
python dicom_to_nifti.py
```

To override paths or enable pixel statistics logging:

```bash
# Override paths
python dicom_to_nifti.py -d /path/to/dicoms -o /path/to/output

# Log intensity statistics for each volume after conversion
python dicom_to_nifti.py --pixel-stats
```

The script skips files that have already been converted, so it is safe to re-run.

---

## Key Implementation Details

### Sorting by Physical Position, Not by Filename

DICOM filenames are arbitrary hashes. The only reliable way to order slices is by the `ImagePositionPatient` tag, which gives the physical (x, y, z) coordinate of the slice in scanner space. Within each temporal group, slices are sorted by their Z coordinate.

### Voxel Spacing from Actual Z Differences

The `SliceThickness` DICOM tag is not always reliable, especially when slices have gaps between them. The pipeline computes slice spacing directly from the difference between consecutive Z positions of the sorted slices. This produces accurate voxel geometry in the output NIfTI.

### Handling DCE Temporal Positions

A DCE-MRI series interleaves slices from multiple time points within the same folder. The pipeline groups them by `TemporalPositionIdentifier` first, builds one 3D volume per time point, then joins these into a 4D volume.

### Solving the Physical Space Mismatch

During DCE acquisition, patients breathe between scans. This causes each time point to be at a marginally different physical position, often submillimeter differences. SimpleITK's `JoinSeries` strictly requires identical physical space across all inputs, and rejects the join with a "physical space mismatch" error even for tiny discrepancies.

The pipeline resolves this by copying spatial metadata from the first time point to all subsequent time points before joining. The pixel data itself is not modified. If precise motion correction is required for downstream analysis, that is handled as a separate registration step.

### Rescale Slope and Intercept

DICOM pixel values are often stored as integers with a slope and intercept that must be applied to recover the true intensity values:

```
true_intensity = stored_value * RescaleSlope + RescaleIntercept
```

The pipeline applies this conversion automatically for every slice before stacking.

---

## Example Output

```
14:23:11 | INFO | DICOM to NIfTI Conversion Pipeline (4D-aware)
14:23:11 | INFO | Source : /home/usama/dicom_learning/data/qin_breast_02
14:23:11 | INFO | Output : /home/usama/dicom_learning/nifti_output
14:23:11 | INFO | Found 235 series folders to convert.
14:23:12 | INFO | [1/235] 1.3.6.1.4.1.14519.5.2.1.8162.5966.347764622180583097
14:23:12 | INFO |   multi-flip_T1-map | 10 timepoint(s) x 10 slices
14:23:13 | INFO |   4D volume shape: (10, 10, 192, 192)
14:23:13 | INFO |   Saved: 8501_multi-flip_T1-map.nii.gz
...
14:47:11 | INFO | Conversion complete.
14:47:11 | INFO |   Successful : 235
14:47:11 | INFO |   Failed     : 0
```

---

## Verifying Output

A quick check that 4D volumes are correctly structured:

```python
import SimpleITK as sitk

img = sitk.ReadImage("nifti_output/QIN-BREAST-02-0001/8501_multi-flip_T1-map.nii.gz")
print("Size       :", img.GetSize())     # (192, 192, 10, 10)
print("Spacing    :", img.GetSpacing())  # (1.33, 1.33, 5.0, 1.0)
print("Dimensions :", img.GetDimension())
```

A 4D size of (192, 192, 10, 10) reads as: 192 by 192 in-plane resolution, 10 slices through the breast, 10 temporal positions across the DCE acquisition.

---

## What This Project Demonstrates

For someone reviewing this as part of a portfolio, this project shows:

1. Understanding of DICOM internal structure and the DICOM hierarchy (Patient, Study, Series, Image)
2. Knowledge of the practical differences between clinical (DICOM) and research (NIfTI) formats
3. Familiarity with `pydicom` and `SimpleITK`, the two foundational libraries for medical imaging in Python
4. Awareness of real-world data issues such as breathing motion, sorting by physical position, and metadata inconsistencies
5. Ability to debug iteratively: the first version of this pipeline produced incorrect 3D output by collapsing temporal positions into the Z axis, which was identified and corrected
6. Clean code structure with hardcoded defaults, optional CLI overrides, error handling, and informative logging

---

## Possible Extensions

This pipeline is a foundation. Logical next steps include:

- Adding intensity normalization (z-score, min-max, or percentile-based)
- Adding resampling to isotropic voxel spacing using `SimpleITK.Resample`
- Adding bias field correction with N4 for MRI volumes
- Integrating with PyRadiomics to extract quantitative imaging features from each volume
- Adding a CLI flag to extract only specific sequence types (for example, only DCE)
- Building a 3D visualization viewer using matplotlib or napari

---

## Author

**Usama Khan**, PhD in Industrial and Information Engineering
Research focus: medical image analysis, breast ultrasound radiomics, deep learning pipelines, PyTorch

---

## License

This project is released under the MIT License.
The QIN-BREAST-02 dataset is distributed by The Cancer Imaging Archive under CC BY 4.0.
