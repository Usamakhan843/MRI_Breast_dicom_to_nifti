"""
DICOM to NIfTI Conversion Pipeline — QIN-BREAST-02 (DCE Breast MRI)
=====================================================================
Converts DICOM series into correctly structured NIfTI volumes:
  - Non-DCE series   : 3D volume  (192, 192, N_slices)
  - DCE series       : 4D volume  (192, 192, N_slices, N_timepoints)

This study grouped by TemporalPositionIdentifier (192, 192, 10, 10) first, sorts each
group by ImagePositionPatient Z, then stacks into proper 4D.

Folder structure expected:
    data_dir/
    └── PatientID/
        └── StudyInstanceUID/
            └── SeriesInstanceUID/
                └── *.dcm

Output structure:
    output_dir/
    └── PatientID/
        └── SeriesNumber_SeriesDescription.nii.gz

Usage:
    python dicom_to_nifti.py                   # uses default paths
    python dicom_to_nifti.py --pixel-stats     # also log volume statistics
    python dicom_to_nifti.py -d /path -o /path # override paths

Author: Usama Khan
"""

import argparse
import logging
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pydicom
import SimpleITK as sitk

# ---------------------------------------------------------------------------
# DEFAULT PATHS
# ---------------------------------------------------------------------------
DEFAULT_DATA_DIR   = "/home/usama/dicom_learning/dicom_reading/data/qin_breast_02"
DEFAULT_OUTPUT_DIR = "/home/usama/dicom_learning/dicom_to_Nifit/nifti_output"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def slugify(text: str) -> str:
    """Convert a string to a safe filename."""
    text = text.strip().replace(" ", "_")
    text = re.sub(r"[^\w\-]", "", text)
    return text[:60]


def get_spacing_origin_direction(ds: pydicom.Dataset):
    """
    Extract spatial geometry from a DICOM slice.
    Returns voxel spacing (x, y), image origin, and direction cosines.
    """
    pixel_spacing = [float(v) for v in ds.PixelSpacing]       # [row, col] spacing in mm
    origin        = [float(v) for v in ds.ImagePositionPatient]
    direction     = [float(v) for v in ds.ImageOrientationPatient]
    return pixel_spacing, origin, direction


def log_volume_stats(arr: np.ndarray, label: str) -> None:
    """Log basic intensity statistics for a volume."""
    logger.info(
        f"  {label} | shape: {arr.shape} | "
        f"mean: {arr.mean():.2f} | std: {arr.std():.2f} | "
        f"min: {arr.min():.2f} | max: {arr.max():.2f}"
    )


# ---------------------------------------------------------------------------
# Core: read and group DICOM slices in one series folder
# ---------------------------------------------------------------------------

def read_series_folder(series_dir: Path) -> tuple[list, dict]:
    """
    Read all DICOM files in a series folder.
    Returns:
        datasets : list of pydicom Datasets (with pixel data)
        info     : dict with patient_id, series_num, series_desc
    """
    dcm_files = list(series_dir.glob("*.dcm"))
    if not dcm_files:
        return [], {}

    datasets = []
    for f in dcm_files:
        try:
            ds = pydicom.dcmread(str(f))
            datasets.append(ds)
        except Exception as e:
            logger.warning(f"  Could not read {f.name}: {e}")

    if not datasets:
        return [], {}

    # Extract naming info from first file
    ds0  = datasets[0]
    info = {
        "patient_id":  getattr(ds0, "PatientID",        "UNKNOWN"),
        "series_num":  getattr(ds0, "SeriesNumber",      "000"),
        "series_desc": getattr(ds0, "SeriesDescription", "unknown"),
    }
    return datasets, info


def group_by_temporal_position(datasets: list) -> dict:
    """
    Group DICOM slices by TemporalPositionIdentifier.
    If tag is absent (non-DCE series), all slices go into group '1'.
    Within each group, sort slices by ImagePositionPatient Z coordinate.
    Returns: {temporal_position: [sorted datasets]}
    """
    groups = defaultdict(list)
    for ds in datasets:
        tp = str(getattr(ds, "TemporalPositionIdentifier", "1"))
        groups[tp].append(ds)

    # Sort each group by Z position (physical slice location)
    for tp in groups:
        groups[tp].sort(
            key=lambda ds: float(ds.ImagePositionPatient[2])
        )

    return dict(sorted(groups.items(), key=lambda x: int(x[0])))


def build_sitk_volume(sorted_slices: list) -> sitk.Image:
    """
    Stack a sorted list of DICOM slices into a 3D SimpleITK image.
    Applies RescaleSlope and RescaleIntercept per slice.
    Sets correct voxel spacing, origin, and direction.
    """
    pixel_arrays = []
    for ds in sorted_slices:
        pixels    = ds.pixel_array.astype(np.float32)
        slope     = float(getattr(ds, "RescaleSlope",     1))
        intercept = float(getattr(ds, "RescaleIntercept", 0))
        pixel_arrays.append(pixels * slope + intercept)

    # Stack slices along Z axis: shape (N_slices, rows, cols)
    volume = np.stack(pixel_arrays, axis=0)

    # Convert to SimpleITK image
    # Note: SimpleITK uses (x, y, z) but numpy uses (z, y, x) — flip needed
    image = sitk.GetImageFromArray(volume)

    # Set spatial geometry from first slice
    ds0           = sorted_slices[0]
    ds1           = sorted_slices[1] if len(sorted_slices) > 1 else sorted_slices[0]
    pixel_spacing = [float(v) for v in ds0.PixelSpacing]

    # Compute slice thickness from actual Z positions (more reliable than tag)
    z0             = float(ds0.ImagePositionPatient[2])
    z1             = float(ds1.ImagePositionPatient[2])
    slice_spacing  = abs(z1 - z0) if z0 != z1 else float(getattr(ds0, "SliceThickness", 5.0))

    image.SetSpacing([pixel_spacing[1], pixel_spacing[0], slice_spacing])
    image.SetOrigin([float(v) for v in ds0.ImagePositionPatient])

    # Set direction cosines from ImageOrientationPatient
    iop    = [float(v) for v in ds0.ImageOrientationPatient]
    row    = iop[:3]   # direction of rows
    col    = iop[3:]   # direction of columns
    normal = list(np.cross(row, col))   # normal to the slice plane
    image.SetDirection(row + col + normal)

    return image


# ---------------------------------------------------------------------------
# Core: convert one series folder to NIfTI
# ---------------------------------------------------------------------------

def convert_series(series_dir: Path, output_dir: Path,
                   pixel_stats: bool = False) -> bool:
    """
    Convert one DICOM series folder to a 3D or 4D NIfTI volume.
    DCE series with multiple temporal positions produce 4D output.
    """
    datasets, info = read_series_folder(series_dir)
    if not datasets:
        logger.warning(f"  Skipping empty series: {series_dir.name}")
        return False

    patient_id  = info["patient_id"]
    series_num  = info["series_num"]
    series_desc = slugify(info["series_desc"])

    patient_out = output_dir / patient_id
    patient_out.mkdir(parents=True, exist_ok=True)

    out_filename = f"{series_num}_{series_desc}.nii.gz"
    out_path     = patient_out / out_filename

    if out_path.exists():
        logger.info(f"  Already exists, skipping: {out_path.name}")
        return True

    try:
        # Group slices by temporal position
        groups = group_by_temporal_position(datasets)
        n_timepoints = len(groups)
        n_slices     = len(list(groups.values())[0])

        logger.info(
            f"  {series_desc} | "
            f"{n_timepoints} timepoint(s) x {n_slices} slices"
        )

        # Build one 3D volume per temporal position
        volumes = []
        for tp, slices in groups.items():
            vol = build_sitk_volume(slices)
            volumes.append(vol)

        if n_timepoints == 1:
            # Non-DCE: save as 3D
            final_image = volumes[0]
            arr         = sitk.GetArrayFromImage(final_image)
            logger.info(f"  3D volume shape: {arr.shape}")
        else:
            # DCE: force all volumes to share the same physical space as
            # the first temporal position before joining into 4D.
            # Patient motion between time points causes slightly different
            # origins/directions which makes JoinSeries fail.
            ref_spacing   = volumes[0].GetSpacing()
            ref_origin    = volumes[0].GetOrigin()
            ref_direction = volumes[0].GetDirection()
            for vol in volumes:
                vol.SetSpacing(ref_spacing)
                vol.SetOrigin(ref_origin)
                vol.SetDirection(ref_direction)
            final_image = sitk.JoinSeries(volumes)
            arr         = sitk.GetArrayFromImage(final_image)
            logger.info(f"  4D volume shape: {arr.shape}  (slices x timepoints x rows x cols)")

        if pixel_stats:
            log_volume_stats(sitk.GetArrayFromImage(final_image), out_filename)

        sitk.WriteImage(final_image, str(out_path))
        logger.info(f"  Saved: {out_path.name}")
        return True

    except Exception as e:
        logger.error(f"  Conversion failed for {series_dir.name}: {e}")
        return False


# ---------------------------------------------------------------------------
# Batch: walk entire dataset
# ---------------------------------------------------------------------------

def convert_dataset(data_dir: str, output_dir: str,
                    pixel_stats: bool = False) -> None:
    data_path   = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    series_dirs = [
        p for p in data_path.rglob("*")
        if p.is_dir() and any(p.glob("*.dcm"))
    ]

    total   = len(series_dirs)
    success = 0
    failed  = 0

    logger.info(f"Found {total} series folders to convert.")

    for i, series_dir in enumerate(series_dirs, 1):
        logger.info(f"[{i}/{total}] {series_dir.name[:50]}")
        ok = convert_series(series_dir, output_path, pixel_stats=pixel_stats)
        if ok:
            success += 1
        else:
            failed += 1

    logger.info("")
    logger.info("=" * 50)
    logger.info("Conversion complete.")
    logger.info(f"  Successful : {success}")
    logger.info(f"  Failed     : {failed}")
    logger.info(f"  Output dir : {output_dir}")
    logger.info("=" * 50)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert QIN-BREAST-02 DICOM series to NIfTI volumes."
    )
    parser.add_argument("--data-dir",    "-d", default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir",  "-o", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--pixel-stats", "-p", action="store_true",
                        help="Log intensity statistics per volume")
    return parser.parse_args()


def main():
    args = parse_args()
    logger.info("DICOM to NIfTI Conversion Pipeline (4D-aware)")
    logger.info(f"Source : {args.data_dir}")
    logger.info(f"Output : {args.output_dir}")
    convert_dataset(args.data_dir, args.output_dir, pixel_stats=args.pixel_stats)


if __name__ == "__main__":
    main()
