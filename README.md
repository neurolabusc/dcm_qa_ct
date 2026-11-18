
## About

dcm_qa_ct is a simple DICOM to NIfTI validator script and dataset to test conversion of Computerized Axial Tomography (CT/CAT) scans. For example, it can be used to identify regressions in the performance of [dcm2niix](https://github.com/rordenlab/dcm2niix). This repository includes data from two CT sessions and synthetic series.

## DataSets

* Philips
  * Description: Phantom scan
  * Model: Ingenuity
  * Comments: Variable gantry tilt courtesy of Paul Morgan (Nottingham Medical School). 

* GE
  * Description: Human scan
  * Model: HiSpeed Dual
  * Comments: Variable gantry tilt and variable inter-slice distances (1.08..6.99mm between slices), provided by [Aliza Medical Imaging](https://www.aliza-dicom-viewer.com/download/datasets). The raw output's slice step minimizes the worst slice position error (15 mm), and its sidecar lists each slice's distance from slice 0 in `SlicePositions` (dcm2niix, 2026-11-18).

* Synthetic
  * Description: 4x4-pixel CT series written by `make_synthetic.py` (deterministic UIDs: rerunning it reproduces `In/Synthetic` byte for byte).
  * Comments: gantry tilt 0 and +/-20 degrees with a 2 mm bed step: classic slices, Enhanced CT without SpacingBetweenSlices (dcm2niix stores the IPP distance) and with a perpendicular SpacingBetweenSlices; classic slices whose instance order is descending (z flip) or permuted, and feet-first (slice column pointing down). `make_synthetic.py [dcm2niix]` also converts them and checks that the raw (non-`_Tilt`) sform puts every slice at its ImagePositionPatient (exact since dcm2niix's bed-step fix, 2026-11-18). Known: the `_Tilt` output of an Enhanced series without SpacingBetweenSlices is resampled 1/cos too far (pre-existing, recorded as is). Also unequal slice spacing without tilt (`unequal_classic_0`, `unequal_classic_descending_0`, positions [0,2,3,4,10]): slice 0 exact, the worst error within the minimax bound (2 mm), and every slice placed by the sidecar's `SlicePositions`. Also overlay planes on a 4D stack (`overlay_4d`, two frames at one position, as a surview): `_ROI1` holds one volume per frame (it was an empty 3D plane before 2026-11-18).

## Links

More CT DICOM datasets are available from these pages:

  - [dcm2nii NITRC web page](https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage#Computed_Tomography_.28CT.2C_CAT.29).
  - [List of DICOM datasets](https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage#Sample_DICOM_Images).

