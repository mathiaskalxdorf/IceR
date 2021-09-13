# IceR 1.2.0

## New Features

* Added capability to process input data which is supplied in a specific input format enabling use of IceR for data preprocessed by any pipeline besides MaxQuant

## Fixes

* Fixed bug in requantify_features() during gam fitting of decoy feature intensities which can result in a crash if too few/no intensities are available (issue #4 on Github)
* Fixed bug in requantify_features() during preparation of some alignment performance QC plots (issue #7 on Github)
* Fixed bug in align_features() during feature alignment which could result in an unhandled error.

---

# IceR 1.1.0

## New Features

* Added capability to handle SILAC data acquired on Thermo devices and analyzed by MaxQuant 

---

# IceR 0.9.9

Initial release.
