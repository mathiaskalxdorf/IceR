# IceR 0.9.13

## Fixes

* Added additional console prints for better monitoring progress of IceR

---

# IceR 0.9.12

## Fixes

* Added up-front checking that samples were acquired with same gradient lengths. If major differences in gradient lengths are detected, IceR returns an error
* Added report of observed m/z and RT deviations between samples to console allowing detection if high inconsistencies are present
* Added warning in case of major deviations in m/z and RT between samples which will hamper reliable results from IceR

---

# IceR 0.9.11

## Fixes

* Check that selected quantification column in load_MaxQ_data() by parameter `intensity_used` is available and otherwise report an error

---

IceR 0.9.10

## Fixes

* Fixed bug in requantify_features() during gam fitting of decoy feature intensities which can result in a crash if too few/no intensities are available (issue #4 on Github)
* Fixed bug in requantify_features() during preparation of some alignment performance QC plots (issue #7 on Github)
* Fixed bug in align_features() during feature alignment which could result in an unhandled error.

---

# IceR 0.9.9

Initial release.
