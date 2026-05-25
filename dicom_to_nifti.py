"""
DICOM to NIfTI Conversion Pipeline — QIN-BREAST-02 (DCE Breast MRI)
=====================================================================
A hybrid pipeline that uses dcm2niix as the conversion engine and
adds batch orchestration, per-patient organization, QC validation,
and summary reporting on top.

Why dcm2niix:
    A hand-written converter must correctly separate every acquisition
    sub-dimension — DCE temporal positions, DWI b-values, multi-echo
    field maps, multi-flip series, qMT offsets. dcm2niix is the
    industry-standard tool built specifically to handle all of these
    plus coordinate-system (LPS/RAS) conversion. This pipeline delegates
    the conversion to dcm2niix and focuses on orchestration and QC.

Folder structure expected:
    data_dir/
    └── PatientID/
        └── StudyInstanceUID/
            └── SeriesInstanceUID/
                └── *.dcm

Output structure:
    output_dir/
    └── PatientID/
        ├── <series>.nii.gz      (image volume)
        ├── <series>.json        (metadata sidecar from dcm2niix)
        └── ...

Usage:
    python dicom_to_nifti.py                    # uses default paths
    python dicom_to_nifti.py -d /path -o /path  # override paths

Requirements:
    dcm2niix   (pip install dcm2niix)
    SimpleITK  (for QC dimension checks)
    pydicom    (to read PatientID for output organization)

Author: Usama
"""

import argparse
import json
import logging
import shutil
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path

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
# Step 0 — Verify dcm2niix is available
# ---------------------------------------------------------------------------

def check_dcm2niix() -> list:
    """
    Confirm dcm2niix is available and return the command list to invoke it.

    Resolution order:
      1. The pip 'dcm2niix' package, which bundles the binary and exposes
         its path via dcm2niix.bin_path
      2. A dcm2niix binary already on the system PATH

    Note: 'dcm2niix --version' exits with a non-zero status code by
    design, so check=True must NOT be used for the probe. We only
    confirm the binary is executable and produces version output.

    Returns the invoker as a list, e.g. ['/path/to/dcm2niix'].
    Raises RuntimeError if neither is found.
    """
    candidates = []

    # Option 1: pip package that bundles the binary
    try:
        import dcm2niix
        candidates.append(str(dcm2niix.bin_path))
    except Exception:
        pass

    # Option 2: dcm2niix already on PATH (e.g. apt install)
    on_path = shutil.which("dcm2niix")
    if on_path:
        candidates.append(on_path)

    for binary in candidates:
        try:
            result = subprocess.run(
                [binary, "--version"], capture_output=True, text=True
            )
            # dcm2niix prints its version then exits non-zero; we only
            # need the binary to run and emit recognizable output.
            if "dcm2nii" in (result.stdout + result.stderr).lower():
                logger.info(f"dcm2niix found: {binary}")
                return [binary]
        except (OSError, FileNotFoundError):
            continue

    raise RuntimeError(
        "dcm2niix not found. Install it with:  pip install dcm2niix"
    )


# ---------------------------------------------------------------------------
# Step 1 — Identify the PatientID for a series folder
# ---------------------------------------------------------------------------

def get_patient_id(series_dir: Path) -> str:
    """
    Read one DICOM file in the series folder to extract PatientID.
    Used to organize outputs into per-patient folders.
    Returns 'UNKNOWN' if no readable file is found.
    """
    for dcm_file in series_dir.glob("*.dcm"):
        try:
            ds = pydicom.dcmread(str(dcm_file), stop_before_pixels=True)
            return str(getattr(ds, "PatientID", "UNKNOWN"))
        except Exception:
            continue
    return "UNKNOWN"


# ---------------------------------------------------------------------------
# Step 2 — Convert one series folder with dcm2niix
# ---------------------------------------------------------------------------

def convert_series(invoker: list, series_dir: Path, patient_out: Path) -> list[Path]:
    """
    Run dcm2niix on one series folder.

    dcm2niix automatically:
      - separates DWI b-values, multi-echo, multi-flip, DCE dynamics
      - corrects LPS to RAS coordinate orientation
      - writes a .json metadata sidecar next to each .nii.gz

    Filename format flags:
      %s = series number     %d = series description
      %e = echo number       (appended automatically when relevant)

    Returns a list of generated .nii.gz file paths.
    """
    patient_out.mkdir(parents=True, exist_ok=True)

    cmd = invoker + [
        "-o", str(patient_out),     # output directory
        "-f", "%s_%d",              # filename: seriesNumber_seriesDescription
        "-z", "y",                  # compress output to .nii.gz
        "-b", "y",                  # write .json BIDS sidecar
        str(series_dir),            # input series folder
    ]

    # Snapshot existing files so we only return newly created ones
    before = set(patient_out.glob("*.nii.gz"))

    try:
        # dcm2niix may return a non-zero exit code even on a successful
        # conversion, so success is judged by whether new .nii.gz files
        # were actually produced, not by the return code.
        subprocess.run(cmd, capture_output=True, text=True)
        after    = set(patient_out.glob("*.nii.gz"))
        produced = list(after - before)
        return produced
    except (OSError, FileNotFoundError) as e:
        logger.error(f"  dcm2niix could not run for {series_dir.name}: {e}")
        return []


# ---------------------------------------------------------------------------
# Step 3 — QC check on a converted NIfTI file
# ---------------------------------------------------------------------------

def qc_check(nifti_path: Path) -> dict:
    """
    Read a converted NIfTI file and return basic QC information:
    dimensionality, size, voxel spacing.
    """
    try:
        img = sitk.ReadImage(str(nifti_path))
        return {
            "file":      nifti_path.name,
            "dimension": img.GetDimension(),
            "size":      img.GetSize(),
            "spacing":   tuple(round(s, 3) for s in img.GetSpacing()),
            "ok":        True,
        }
    except Exception as e:
        return {"file": nifti_path.name, "ok": False, "error": str(e)}


# ---------------------------------------------------------------------------
# Step 4 — Batch process the whole dataset
# ---------------------------------------------------------------------------

def convert_dataset(data_dir: str, output_dir: str) -> list[dict]:
    """
    Walk the dataset, convert every series folder, and run QC.
    Returns a list of QC records, one per generated NIfTI file.
    """
    invoker     = check_dcm2niix()
    data_path   = Path(data_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Collect all series-level folders (those directly containing .dcm files)
    series_dirs = [
        p for p in data_path.rglob("*")
        if p.is_dir() and any(p.glob("*.dcm"))
    ]
    total = len(series_dirs)
    logger.info(f"Found {total} series folders to convert.")

    qc_records   = []
    converted    = 0
    failed       = 0

    for i, series_dir in enumerate(series_dirs, 1):
        patient_id  = get_patient_id(series_dir)
        patient_out = output_path / patient_id

        logger.info(f"[{i}/{total}] {patient_id} | {series_dir.name[:40]}")

        produced = convert_series(invoker, series_dir, patient_out)

        if not produced:
            failed += 1
            continue

        # QC each produced file
        for nifti_path in produced:
            qc = qc_check(nifti_path)
            qc["patient_id"] = patient_id
            qc_records.append(qc)
            if qc["ok"]:
                logger.info(f"  {qc['file']} | dim {qc['dimension']} | size {qc['size']}")
            else:
                logger.warning(f"  QC FAILED: {qc['file']} | {qc.get('error')}")

        converted += 1

    logger.info("")
    logger.info("=" * 55)
    logger.info("Conversion complete.")
    logger.info(f"  Series converted : {converted}/{total}")
    logger.info(f"  Series failed    : {failed}")
    logger.info(f"  NIfTI files made : {len(qc_records)}")
    logger.info("=" * 55)

    return qc_records


# ---------------------------------------------------------------------------
# Step 5 — Summary report
# ---------------------------------------------------------------------------

def write_summary(qc_records: list[dict], data_dir: str,
                  output_dir: str, report_path: str) -> None:
    """Generate a plain-text QC summary report."""

    per_patient = defaultdict(int)
    dim_counts  = defaultdict(int)
    qc_failures = []

    for r in qc_records:
        per_patient[r["patient_id"]] += 1
        if r.get("ok"):
            dim_counts[r["dimension"]] += 1
        else:
            qc_failures.append(r)

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "=" * 60,
        "  QIN-BREAST-02 — DICOM TO NIfTI CONVERSION REPORT",
        "=" * 60,
        f"  Generated        : {now}",
        f"  Source dir       : {data_dir}",
        f"  Output dir       : {output_dir}",
        f"  Conversion engine: dcm2niix",
        f"  Total NIfTI files: {len(qc_records)}",
        "",
        "  FILES PER PATIENT",
        "-" * 60,
    ]
    for pid, count in sorted(per_patient.items()):
        lines.append(f"  {pid:<35} {count} volumes")

    lines += [
        "",
        "  VOLUME DIMENSIONALITY",
        "-" * 60,
        f"  3D volumes       : {dim_counts.get(3, 0)}",
        f"  4D volumes       : {dim_counts.get(4, 0)}",
        "",
        "  QC STATUS",
        "-" * 60,
        f"  Passed QC        : {len(qc_records) - len(qc_failures)}",
        f"  Failed QC        : {len(qc_failures)}",
    ]
    for r in qc_failures:
        lines.append(f"    - {r['file']}: {r.get('error', 'unknown error')}")

    lines += ["", "=" * 60]

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    logger.info(f"Summary report saved → {report_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert QIN-BREAST-02 DICOM series to NIfTI using dcm2niix."
    )
    parser.add_argument("--data-dir",   "-d", default=DEFAULT_DATA_DIR,
                        help=f"Root DICOM directory (default: {DEFAULT_DATA_DIR})")
    parser.add_argument("--output-dir", "-o", default=DEFAULT_OUTPUT_DIR,
                        help=f"NIfTI output directory (default: {DEFAULT_OUTPUT_DIR})")
    return parser.parse_args()


def main():
    args = parse_args()
    logger.info("DICOM to NIfTI Conversion Pipeline (dcm2niix engine)")
    logger.info(f"Source : {args.data_dir}")
    logger.info(f"Output : {args.output_dir}")

    qc_records  = convert_dataset(args.data_dir, args.output_dir)

    report_path = str(Path(args.output_dir) / "conversion_report.txt")
    write_summary(qc_records, args.data_dir, args.output_dir, report_path)

    print("\n✓ Pipeline complete.")
    print(f"  NIfTI files → {args.output_dir}")
    print(f"  QC report   → {report_path}")


if __name__ == "__main__":
    main()
