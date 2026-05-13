# Third-Party Dependencies and Licensing

The source code in this repository is released under the **MIT License**
(see [`LICENSE`](./LICENSE)).

However, this software depends on external tools that are **not distributed
with this repository** and which carry their own licensing restrictions.
Users are responsible for obtaining these tools directly from their
respective providers and for complying with all applicable terms.

The MIT License applies **only to the original code in this repository**.
It does not, and cannot, grant any rights to the third-party software
described below.

---

## 1. AlphaFold 3

This project uses [AlphaFold 3](https://github.com/google-deepmind/alphafold3),
developed by Google DeepMind and Isomorphic Labs.

### Licensing
The AlphaFold 3 source code, model parameters, and outputs are made available
by Google **free of charge for certain non-commercial uses only**, subject to:

- [AlphaFold 3 Model Parameters Terms of Use](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)
- [AlphaFold 3 Model Parameters Prohibited Use Policy](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_PROHIBITED_USE_POLICY.md)
- [AlphaFold 3 Output Terms of Use](https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md)

### Key restrictions
- **Non-commercial use only.** Use on behalf of, or in connection with,
  any commercial organization (including research conducted on behalf of
  commercial organizations) is prohibited.
- Only non-commercial organizations (universities, non-profit research
  institutes, educational and government bodies, and journalism) may use
  the AlphaFold 3 assets.
- You may not share the AlphaFold 3 assets with any commercial organization
  or use them in a manner that grants any rights to such an organization.
- Google retains the right to revoke access at any time.

### Obtaining model parameters
The AlphaFold 3 model parameters are **not** redistributed with this
repository. They must be requested directly from Google DeepMind by
submitting the official request form using an institutional email address.

### Required citation
If you use this software with AlphaFold 3, you must cite:

> Abramson, J. et al. Accurate structure prediction of biomolecular
> interactions with AlphaFold 3. *Nature* (2024).

---

## 2. Rosetta

This project uses [Rosetta](https://rosettacommons.org/), developed by the
Rosetta Commons and licensed through the University of Washington (UW).

### Licensing
Rosetta is governed by the
[Rosetta Software License Agreement](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md).

- **Free** for academic, non-profit, and government use under the
  Rosetta Software Academic License Agreement.
- **Commercial use requires a separate paid license** obtained from
  UW CoMotion (license@uw.edu). Fee-for-service work by academic users
  is also considered commercial under this license.

### Key restrictions
- Rosetta and any modifications to Rosetta **may not be redistributed**
  by licensees, except as forks of official Rosetta Commons repositories
  on the same version control platform.
- All forks of the Rosetta code must maintain the current licensing
  restrictions.
- Copyright and license notices in the Rosetta source must be retained.

---

## 3. Effective use restrictions

Although the original code in this repository is MIT-licensed, the
**practical ability to run this software end-to-end is limited to
non-commercial use** unless you have separately obtained:

1. Authorization from Google DeepMind to use AlphaFold 3 (or a future
   commercial arrangement, if one becomes available), **and**
2. A commercial license for Rosetta from UW CoMotion.

You are solely responsible for ensuring that your use of AlphaFold 3,
Rosetta, and any other third-party dependencies complies with their
respective licenses. The authors of this repository make no
representations about the suitability of any external tool for any
particular use case and disclaim all liability for license compliance
relating to third-party software.

---

## 4. Other dependencies

Standard open-source Python and system dependencies (e.g., those installed
via `pip` or your package manager) are governed by their own licenses,
which can be found in their respective package metadata or repositories.
