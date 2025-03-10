# ==============================================================================
# This script is used to merge all JKC synthetic photometry data file generated by gaiaxpy_phot.py into a single CSV file
# ==============================================================================

import pathlib
import pandas as pd
import tqdm

base_path = pathlib.Path("./jkc_photometry/tables")

file_ls = list(base_path.glob("XpContinuousMeanSpectrum*.csv"))

_table = pd.read_csv(file_ls[0])
base_table = pd.DataFrame(data={"source_id": _table["source_id"], "Jkc_mag_B": _table["Jkc_mag_B"], "Jkc_mag_V": _table["Jkc_mag_V"]})

for i in tqdm.tqdm(file_ls[1:]):  # append to base_table
    _table = pd.read_csv(i)
    base_table = pd.concat([base_table, pd.DataFrame(data={"source_id": _table["source_id"], "Jkc_mag_B": _table["Jkc_mag_B"], "Jkc_mag_V": _table["Jkc_mag_V"]})], ignore_index=True)

# save to file
base_table.to_csv(base_path / "Gaia_XP_JKC.csv", index=False)
