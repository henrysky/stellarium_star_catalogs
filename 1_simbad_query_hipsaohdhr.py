import logging
import pathlib

import tqdm
from astropy.table import Table, vstack

from py.utils import custom_simbad

if __name__ == "__main__":
    base_path = pathlib.Path("simbad_query_results")
    base_path.mkdir(parents=True, exist_ok=True)
    hip_subdir = base_path / "hip"
    hip_subdir.mkdir(parents=True, exist_ok=True)
    sao_subdir = base_path / "sao"
    sao_subdir.mkdir(parents=True, exist_ok=True)
    hd_subdir = base_path / "hd"
    hd_subdir.mkdir(parents=True, exist_ok=True)
    hr_subdir = base_path / "hr"
    hr_subdir.mkdir(parents=True, exist_ok=True)
    hip_combined_path = base_path / "hip_combined.dat"
    sao_combined_path = base_path / "sao_combined.dat"
    hd_combined_path = base_path / "hd_combined.dat"
    hr_combined_path = base_path / "hr_combined.dat"
    max_hip_id = 120416
    max_sao_id = 258997
    max_hd_id = 272150
    max_hr_id = 9110
    query_batch_size = 2000

    catalog_configs = [
        (hip_subdir, max_hip_id, hip_combined_path),
        (sao_subdir, max_sao_id, sao_combined_path),
        (hd_subdir, max_hd_id, hd_combined_path),
        (hr_subdir, max_hr_id, hr_combined_path),
    ]

    for subdir, max_id, combined_path in catalog_configs:
        if combined_path.exists():  # if the combined file already exists, skip
            logging.warning(
                f"Skipping {subdir.name} as the combined file already exists. If you want to re-query, delete the combined file."
            )
            continue
        for batch in tqdm.tqdm(
            range(max_id // query_batch_size + 1), desc=f"Querying {subdir.name}"
        ):
            max_id_clipped = min(max_id, batch * query_batch_size + query_batch_size)
            ids = [
                f"{subdir.name.upper()} {str(i)}".strip()
                for i in range(1 + batch * query_batch_size, 1 + max_id_clipped)
            ]
            result = custom_simbad.query_objects(ids)
            result.write(
                subdir / f"simbad_{subdir.name}_{str(batch)}.dat",
                format="ascii",
                overwrite=True,
            )
        # merge all the tables
        table_list = []
        counter = 0
        curr_path = subdir / f"simbad_{subdir.name}_{str(counter)}.dat"
        while (
            curr_path.exists()
        ):  # need to loop through all the files in the exact order. can't use glob("*")
            table_list.append(Table.read(curr_path, format="ascii"))
            counter += 1
            curr_path = subdir / f"simbad_{subdir.name}_{str(counter)}.dat"
        simbad_table = vstack(table_list)
        simbad_table.write(combined_path, format="ascii", overwrite=True)
