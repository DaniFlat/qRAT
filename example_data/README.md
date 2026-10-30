# qRAT example data

Example qPCR result files to try out qRAT. Download all files at once as
[qRAT_example_data.zip](https://github.com/DaniFlat/qRAT/raw/main/example_data/qRAT_example_data.zip).

| File | Use in qRAT | Content |
|---|---|---|
| `single_plate.csv` | *Single Plate* | One plate, samples `C` (control) and `T` (treatment), reference gene `ref G` and genes of interest `GOI 1`–`GOI 7`, technical replicates |
| `plate_1.csv`, `plate_2.csv` | *Multiple Plates* (upload both) | Two plates with samples `WT`/`M` (control `C`, time points `T1`–`T9`), `ref G`, `GOI 1`, `GOI 2`. Both plates contain the inter-plate calibrator `IPC` and the control `WT C` |

All files are comma separated with the columns `Well`, `Sample Name`, `Target Name` and `Cq`.
For the analysis steps see the [qRAT User Guide](https://github.com/DaniFlat/qRAT/blob/main/www/qRAT_User_Guide.pdf).
