#!/usr/bin/env python3
"""Synthetic gantry-tilt CT series for In/Synthetic, and a geometric check: the
raw (non-_Tilt) NIfTI sform must place every slice at its ImagePositionPatient,
and pixdim x*y*z must equal the voxel volume (|det| of the sform). Classic slices, Enhanced CT without
SpacingBetweenSlices (the parser stores the IPP distance) and Enhanced CT with a
perpendicular SpacingBetweenSlices, each at 0 and +/-20 degrees; classic slices
whose instance order is descending (z flip) or permuted, and feet-first (slice
column pointing down). Cases from external and internal audits (2026-11-18).

Regenerate the fixtures (deterministic UIDs, so the files do not change), then
convert them with dcm2niix and check the geometry; batch.sh compares Ref as usual:
  uv run --no-project --with pydicom --with numpy python3 make_synthetic.py [path/to/dcm2niix]
"""
import hashlib, json, math, struct, subprocess, sys, tempfile
from pathlib import Path

import numpy as np
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian

IN = Path(__file__).resolve().parent / "In" / "Synthetic"


def uid(*parts):
    """A stable 2.25 UID from the case name and role (fixtures regenerate byte-identical)."""
    return "2.25." + str(int(hashlib.sha256("/".join(map(str, parts)).encode()).hexdigest()[:30], 16))

EXE = sys.argv[1] if len(sys.argv) > 1 else "dcm2niix"
STEP = 2.0  # bed step (mm) between slice origins


def item(**kw):
    ds = Dataset()
    for k, v in kw.items():
        setattr(ds, k, v)
    return ds


def write_series(folder, number, enhanced, spacing, tilt, inst=(1, 2, 3), row_x=1, position="HFS", maker="SIEMENS", zs=None):
    name = folder.name
    c, s = math.cos(math.radians(tilt)), math.sin(math.radians(tilt))
    orient = [row_x, 0, 0, 0, round(c, 10), round(s, 10)]
    ds = Dataset()
    ds.file_meta = FileMetaDataset()
    ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds.SOPClassUID = ds.file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2" + (".1" if enhanced else "")
    ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID = uid(name, "sop")
    ds.StudyInstanceUID, ds.SeriesInstanceUID, ds.FrameOfReferenceUID = uid("study"), uid(name, "series"), uid(name, "for")
    ds.Modality, ds.Manufacturer, ds.PatientPosition = "CT", maker, position
    ds.SeriesNumber, ds.InstanceNumber, ds.SeriesDescription, ds.ProtocolName = number, 1, name, name
    ds.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
    ds.StudyDate, ds.StudyTime, ds.SeriesTime = "20260101", "120000", "120000"
    ds.Rows = ds.Columns = 4
    ds.SamplesPerPixel, ds.PhotometricInterpretation = 1, "MONOCHROME2"
    ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 16, 16, 15, 1
    ds.RescaleSlope, ds.RescaleIntercept = 1, 0
    ds.GantryDetectorTilt = tilt
    measures = item(PixelSpacing=[1, 1], SliceThickness=round(STEP * c, 10))
    if spacing:
        measures.SpacingBetweenSlices = round(STEP * c, 10)  # perpendicular
    n = len(inst)
    pixels = np.stack([np.full((4, 4), 100 + k, dtype="<i2") for k in range(n)])
    if enhanced:
        ds.NumberOfFrames = 3
        ds.SharedFunctionalGroupsSequence = [
            item(PlaneOrientationSequence=[item(ImageOrientationPatient=orient)], PixelMeasuresSequence=[measures])
        ]
        org = uid(name, "org")
        ds.DimensionOrganizationSequence = [item(DimensionOrganizationUID=org)]
        ds.DimensionIndexSequence = [
            item(DimensionOrganizationUID=org, DimensionIndexPointer=0x00209057, FunctionalGroupPointer=0x00209111)
        ]
        ds.PerFrameFunctionalGroupsSequence = [
            item(
                FrameContentSequence=[item(DimensionIndexValues=k + 1, InStackPositionNumber=k + 1, StackID="1")],
                PlanePositionSequence=[item(ImagePositionPatient=[0, 0, k * STEP])],
            )
            for k in range(3)
        ]
        ds.PixelData = pixels.tobytes()
        ds.save_as(folder / "enh.dcm", enforce_file_format=True)
        return
    ds.ImageOrientationPatient, ds.PixelSpacing, ds.SliceThickness = orient, [1, 1], round(STEP * c, 10)
    for k in range(n):
        ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID = uid(name, "sop", k)
        ds.ImagePositionPatient, ds.InstanceNumber = [0, 0, zs[k] if zs else k * STEP], inst[k]
        ds.PixelData = pixels[k].tobytes()
        ds.save_as(folder / f"{k}.dcm", enforce_file_format=True)


def ipp_error(nii, n, zs=None):
    """Worst distance from slice k's ImagePositionPatient (pixel value 100+k) to
    the sform position of its nearest corner voxel (rows and columns may flip);
    with unequal spacing (zs) also slice 0's error and the worst error of the
    sidecar's SlicePositions placement."""
    b = nii.read_bytes()
    m = np.array([struct.unpack("<4f", b[280 + 16 * r : 296 + 16 * r]) for r in range(3)] + [(0, 0, 0, 1)])
    off = int(struct.unpack("<f", b[108:112])[0])
    img = np.frombuffer(b[off : off + 2 * 16 * n], "<i2").reshape(n, 4, 4)
    pos = json.loads(nii.with_suffix(".json").read_text()).get("SlicePositions") if zs else None
    unit = m[:3, 2] / np.linalg.norm(m[:3, 2])
    err = first = sp = 0.0
    for k in range(n):
        kk = int(np.nonzero(img[:, 0, 0] == 100 + k)[0][0])
        want = [0, 0, zs[k] if zs else k * STEP]
        best = min(np.linalg.norm((m @ [i, j, kk, 1])[:3] - want) for i in (0, 3) for j in (0, 3))
        err = max(err, best)
        if kk == 0:
            first = best
        if pos:
            at = lambda i, j: (m @ [i, j, 0, 1])[:3] + pos[kk] * unit
            sp = max(sp, min(np.linalg.norm(at(i, j) - want) for i in (0, 3) for j in (0, 3)))
    if zs and not pos:
        sp = float("inf")  # unequal spacing must write SlicePositions
    # pixdim x*y*z must be the voxel volume, |det| of the sform (sheared or not).
    px = struct.unpack("<8f", b[76:108])
    det = abs(np.linalg.det(m[:3, :3]))
    if abs(px[1] * px[2] * px[3] - det) > 1e-4 * det:
        err = float("inf")
    return m[2, 2], err, first, sp


CASES = [
    (f"tilt_{'enhanced' if e else 'classic'}_{'spacing' if sp else 'nospacing'}_{t if t >= 0 else 'm' + str(-t)}", dict(enhanced=e, spacing=sp, tilt=t))
    for e, sp in [(False, False), (True, False), (True, True)]
    for t in [0, 20, -20]
] + [
    ("tilt_classic_descending_20", dict(enhanced=False, spacing=False, tilt=20, inst=(4, 3, 2, 1))),
    ("tilt_classic_permuted_20", dict(enhanced=False, spacing=False, tilt=20, inst=(1, 2, 4, 3), maker="TOSHIBA")),
    ("tilt_classic_ffs_20", dict(enhanced=False, spacing=False, tilt=20, inst=(1, 2, 3, 4), row_x=-1, position="FFS")),
    ("tilt_classic_ffsdescending_20", dict(enhanced=False, spacing=False, tilt=20, inst=(4, 3, 2, 1), row_x=-1, position="FFS")),
    # Unequal spacing without tilt: slice 0 exact, minimax step (worst 2 mm for
    # these positions; the mean would give 3.5), SlicePositions exact.
    ("unequal_classic_0", dict(enhanced=False, spacing=False, tilt=0, inst=(1, 2, 3, 4, 5), zs=[0, 2, 3, 4, 10])),
    ("unequal_classic_descending_0", dict(enhanced=False, spacing=False, tilt=0, inst=(5, 4, 3, 2, 1), zs=[0, 2, 3, 4, 10])),
]
MINIMAX_WORST = 2.0


failures = 0
with tempfile.TemporaryDirectory() as tmp:
    for number, (name, kw) in enumerate(CASES, start=1):
        src, out = IN / name, Path(tmp, name)
        src.mkdir(parents=True, exist_ok=True)
        out.mkdir()
        write_series(src, number, **kw)
        subprocess.run([EXE, "-z", "n", "-b", "y", "-f", "out", "-o", str(out), str(src)], capture_output=True)
        nii = out / "out.nii"
        if not nii.exists():
            print(f"FAIL {name}: no output")
            failures += 1
            continue
        zs = kw.get("zs")
        col_z, err, first, sp = ipp_error(nii, len(kw.get("inst", (1, 2, 3))), zs)
        ok = (first < 1e-3 and err <= MINIMAX_WORST + 1e-3 and sp < 1e-3) if zs else err < 1e-3
        failures += not ok
        extra = f", slice 0 {first:.6f}, SlicePositions {sp:.6f} mm" if zs else ""
        print(f"{'ok  ' if ok else 'FAIL'} {name}: sform slice column z {col_z:.4f}, max IPP error {err:.6f} mm{extra}")
sys.exit(1 if failures else 0)
