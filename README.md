# DICOM to NIfTI Conversion Pipeline for Breast MRI

A conversion and quality-control pipeline that transforms raw clinical DICOM image series into research-ready NIfTI volumes. The pipeline uses `dcm2niix` as its conversion engine and adds patient-level organization, automated QC validation, and summary reporting around it.

Built and tested on the QIN-BREAST-02 collection from The Cancer Imaging Archive (TCIA), a publicly available Dynamic Contrast-Enhanced (DCE) breast MRI dataset of 13 subjects. It follows a DICOM metadata extraction pipeline and moves from inspecting data to transforming it into the format expected by research tools such as 3D Slicer, ITK-SNAP, LIFEx, MONAI, and PyRadiomics.

---

## Table of Contents

- [Motivation](#motivation)
- [Design Decision: Why dcm2niix](#design-decision-why-dcm2niix)
- [Dataset](#dataset)
- [Pipeline Overview](#pipeline-overview)
- [Input and Output Structure](#input-and-output-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Implementation Notes](#implementation-notes)
- [What This Project Demonstrates](#what-this-project-demonstrates)
- [Possible Extensions](#possible-extensions)

---

## Motivation

Medical scanners store images in DICOM format, where every individual slice is a separate `.dcm` file carrying a full header of clinical metadata. This format suits hospital systems but is impractical for research. A single DCE breast MRI study can produce several thousand DICOM files spread across nested folders.

Before analysis can begin, those slices must be grouped by acquisition sequence, ordered by physical position, separated by acquisition sub-dimension (time point, b-value, echo), and stacked into coherent 3D or 4D volumes. NIfTI is the research standard: one file holds one complete volume with spatial geometry preserved.

This pipeline automates that transformation across an entire dataset.

---

## Design Decision: Why dcm2niix

This project initially used a custom converter written with `pydicom` and `SimpleITK`. That version grouped slices only by `TemporalPositionIdentifier`. It worked for DCE and multi-flip series but failed silently on others.

The failure surfaced during visual inspection in ITK-SNAP and LIFEx: diffusion-weighted (DWI) and field-map volumes showed a distinctive horizontal striping artifact. The cause was identified as multi-volume interleaving. A DWI series contains several b-values stored in one folder; the custom converter did not separate them, so slices from different b-values were stacked together in the wrong order.

Correctly handling every acquisition sub-dimension (DCE time points, DWI b-values, multi-echo field maps, multi-flip series, qMT offsets) along with coordinate-system conversion is a substantial engineering problem. `dcm2niix` is the established open-source tool built specifically for this, refined over years of edge-case handling.

The pipeline was therefore redesigned around a clear separation of responsibilities:

| Responsibility | Handled by |
|---|---|
| DICOM to NIfTI conversion | `dcm2niix` |
| Sub-dimension separation, orientation | `dcm2niix` |
| Batch orchestration over the dataset | This pipeline |
| Per-patient output organization | This pipeline |
| QC validation of every output volume | This pipeline |
| Summary reporting | This pipeline |

Choosing an established tool over a custom reimplementation is a deliberate engineering decision, not a shortcut. The value of this project lies in robust orchestration and validation around a trusted conversion core.

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

The dataset mixes several sequence types: DCE-MRI, diffusion-weighted imaging, T1 mapping, quantitative magnetization transfer, and field maps. Each produces 3D or 4D volumes depending on its acquisition structure.

---

## Pipeline Overview

```
Step 1   For each series folder, read PatientID for output organization
Step 2   Run dcm2niix on the series folder (the conversion itself)
Step 3   QC check every produced NIfTI file (dimensions, size, spacing)
Step 4   Repeat across all 235 series folders
Step 5   Write a QC summary report
```

---

## Input and Output Structure

**Input** (DICOM hierarchy as downloaded from TCIA):

```
qin_breast_02/
└── PatientID/
    └── StudyInstanceUID/
        └── SeriesInstanceUID/
            ├── slice_001.dcm
            ├── slice_002.dcm
            └── ...
```

**Output** (NIfTI volumes plus metadata sidecars):

```
nifti_output/
├── QIN-BREAST-02-0001/
│   ├── 8101_THRIVE_SENSE.nii.gz      image volume
│   ├── 8101_THRIVE_SENSE.json        metadata sidecar
│   ├── 7401_DWI_EPI.nii.gz
│   ├── 7401_DWI_EPI.json
│   └── ...
├── QIN-BREAST-02-0002/
│   └── ...
├── ...
└── conversion_report.txt             QC summary
```

For each series, `dcm2niix` produces a compressed `.nii.gz` volume and a `.json` sidecar containing acquisition metadata in BIDS format. When a series contains multiple sub-volumes (for example several b-values), `dcm2niix` writes each as a separate correctly ordered file.

---

## Installation

```bash
git clone https://github.com/Usamakhan843/MRI_breast_dicom_to_nifti.git
cd MRI_breast_dicom_to_nifti
pip install -r requirements.txt
```

Requirements:

```
dcm2niix>=1.0.20240000
SimpleITK>=2.3.0
pydicom>=2.4.0
```

The pip `dcm2niix` package bundles the conversion binary, so no separate system installation is required.

---

## Usage

Default data and output paths are set at the top of the script. Routine use requires only:

```bash
python dicom_to_nifti.py
```

To override the paths:

```bash
python dicom_to_nifti.py -d /path/to/dicoms -o /path/to/output
```

The script reports progress per series and writes a QC summary on completion.

---

## Implementation Notes

### Locating the dcm2niix binary

The pip `dcm2niix` package does not place a command on the system PATH. Instead it bundles the binary inside the package directory and exposes its location through `dcm2niix.bin_path`. The pipeline imports the package and resolves the binary from there, with a fallback to any system-installed `dcm2niix` on the PATH.

### Success is judged by output, not by exit code

`dcm2niix` returns a non-zero exit status in situations that are not failures. Its `--version` call exits non-zero by design, and some successful conversions also return non-zero. Relying on the exit code therefore produces false failures.

The pipeline judges success by inspecting the actual result instead:

- For the availability check, it confirms the binary runs and emits recognizable version output.
- For each conversion, it snapshots the output folder before and after the call and treats any newly created `.nii.gz` files as the success signal.

This is a general principle worth noting: an exit code is a claim, not a guarantee. Verifying the real artifact is more reliable.

### Quality control on every volume

After each conversion, every produced NIfTI file is opened with `SimpleITK` and checked for readability, dimensionality (3D or 4D), size, and voxel spacing. These results feed the final summary report, giving an at-a-glance audit of the whole conversion.

---

## What This Project Demonstrates

For a reviewer assessing this as part of a portfolio, the project shows:

1. Understanding of the DICOM hierarchy and the practical differences between clinical and research image formats.
2. Knowledge of real acquisition structure: temporal positions, b-values, echoes, and why naive slice stacking corrupts multi-volume series.
3. The judgment to identify a flawed approach, diagnose it through visual inspection, and replace it with an established tool rather than persisting with a custom reimplementation.
4. Practical engineering habits: separating responsibilities cleanly, validating results rather than trusting return codes, and building QC into the pipeline.
5. Familiarity with the standard medical imaging toolchain: `dcm2niix`, `SimpleITK`, `pydicom`, and viewers such as ITK-SNAP and LIFEx.

---

## Possible Extensions

- Intensity normalization (z-score, percentile-based) for cross-subject comparability
- Resampling to isotropic voxel spacing with `SimpleITK.Resample`
- N4 bias field correction for MRI intensity uniformity
- Parsing the `dcm2niix` JSON sidecars into a dataset-wide metadata table
- A CLI flag to convert only selected sequence types, for example DCE only
- Integration with PyRadiomics to extract quantitative features from each volume

---


## Author

**Usama Khan**, PhD in Industrial and Information Engineering

Research focus: medical image analysis, breast ultrasound radiomics, deep learning pipelines, PyTorch

---

## License

This project is released under the MIT License.
The QIN-BREAST-02 dataset is distributed by The Cancer Imaging Archive under CC BY 4.0.
