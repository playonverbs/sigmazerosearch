# HyperonProduction NTuple Specification

:::{card}
:link: https://github.com/playonverbs/HyperonProduction

This table describes the format of output NTuples produced by the LArSoft
analyser at <github:playonverbs/HyperonProduction>.
:::

Subtypes listed in the 'Type' column are described in separate tables. However
in the ntuple they are unpacked into primitive component types.

:::{important}
Currently this spec applies to a single-slice analysis.

All `std::*` types are implied unless specified.
:::

| Variable                      | Type                    | Notes                     | Bounds     | Units  |
| ----------------------------- | ----------------------- | ------------------------- | --------   | -----  |
| run                           | `unsigned int`          | Event run number          | n/a        | n/a    |
| subrun                        | `unsigned int`          | Event subrun number       | n/a        | n/a    |
| event                         | `unsigned int`          | Event number              | n/a        | n/a    |
| **mc**_nu_pdg                 | `unsigned int`          | MC nu PDG code            | n/a        | n/a    |
| **mc**_nu_q2                  | `double`                | MC nu interaction q2      | $[0,...]$  | n/a    |
| **mc**_ccnc                   | `string`                | MC event CC/NC            | CC/NC/None | n/a    |
| **mc**_mode                   | `string`                | MC event mode             | n/a        | n/a    |
| **mc**\_nu_pos_{x,y,z}        | `double`                | MC nu pos coord           | n/a        | cm     |
| **mc**_lepton_pdg             | `int`                   | MC lepton PDG code        | n/a        | n/a    |
| **mc**_lepton_mom             | `double`                | MC lepton momentum mag    | n/a        | n/a    |
| true_nu_slice_ID              | `int`                   |                           | $[0,...]$  | n/a    |
| true_nu_slice_completeness    | `double`                |                           | $[0,1]$    | n/a    |
| true_nu_slice_purity          | `double`                |                           | $[0,1]$    | n/a    |
| n_slices                      | `int`                   |                           | $[0,...]$  | n/a    |
| flash_match_nu_slice_ID       | `int`                   |                           | $[0,...]$  | n/a    |
| pandora_nu_slice_ID           | `int`                   |                           | $[0,...]$  | n/a    |
| **pfp**_purity                | `vector<double>`        |                           | $[0,1]$    | n/a    |
| **pfp**_completeness          | `vector<double>`        |                           | $[0,1]$    | n/a    |
| **pfp**_has_truth             | `vector<bool>`          |                           | n/a        | n/a    |
| **pfp**_trackID               | `vector<int>`           |                           | $[0,...]$  | n/a    |
| **pfp**_true_energy           | `vector<double>`        |                           | n/a        | n/a    |
| **pfp**_true_ke               | `vector<double>`        |                           | $[0,...]$  | n/a    |
| **pfp**_true_px               | `vector<double>`        |                           | $[0,...]$  | n/a    |
| **pfp**_true_py               | `vector<double>`        |                           | $[0,...]$  | n/a    |
| **pfp**_true_pz               | `vector<double>`        |                           | $[0,...]$  | n/a    |
| **pfp**_true_length           | `vector<double>`        |                           | $[0,...]$  | n/a    |
| **pfp**_true_origin           | `vector<int>`           |                           | $[0,5]$    | n/a    |
| **pfp**_true_pdg              | `vector<int>`           | Backtracked PDG value     | $[0,5]$    | n/a    |
| **pfp**_trk_shr_score         | `vector<double>`        |                           | $[0,1]$    | n/a    |
| **pfp**_x                     | `vector<double>`        |                           | $[0,1]$    | n/a    |
| **pfp**_y                     | `vector<double>`        |                           | $[0,1]$    | n/a    |
| **pfp**_z                     | `vector<double>`        |                           | $[0,1]$    | n/a    |
| **trk**_llrpid                | `vector<double>`        | Log-Likelihood ratio PID  | $[-1,1]$   | n/a    |
| **trk**_three_plane_mean_dedx | `vector<double>`        | Three-plane mean dE/dx    | $[0,...]$  | MeV/cm |

