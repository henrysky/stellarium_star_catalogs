# Many of these HIP stars are in a binary system, unresolved by Hipparcos but resolved by Gaia which will cause cross-matching issue.
# Usually by the fact that they are missing `FLUX_G` but not `FLUX_V`, `FLUX_J` because Gaia did resolve them.
# The underlying issue is these binary stars are closed enough to be observed by 2MASS/Hipparcos as a single source hence having
# V, JHK band magnitude but resolved by Gaia hence not having a G-band magnitude because such single object does not exist in Gaia.

import hashlib
import pathlib
import re
import struct
import textwrap
import time
import warnings

import astropy.units as u
import polars as pl
from astropy.coordinates.name_resolve import NameResolveError
from astropy.table import Table
from astroquery.gaia import Gaia
from tqdm import tqdm

from py.utils import custom_simbad
from astropy.coordinates import SkyCoord

base_path = pathlib.Path("simbad_query_results")
hip_combined_path = base_path / "hip_combined.ecsv"
sao_combined_path = base_path / "sao_combined.ecsv"
hd_combined_path = base_path / "hd_combined.ecsv"
hr_combined_path = base_path / "hr_combined.ecsv"
# check if the combined file already exists, if not then raise
for combined_path in [
    hip_combined_path,
    sao_combined_path,
    hd_combined_path,
    hr_combined_path,
]:
    if not combined_path.exists():
        raise FileNotFoundError(
            f"{combined_path} does not exist. Please run simbad_query_hipsaohdhr.py first."
        )
hip_combined_df = Table.read(hip_combined_path, format="ecsv").to_df("polars")
sao_combined_table = Table.read(sao_combined_path, format="ecsv")
hd_combined_table = Table.read(hd_combined_path, format="ecsv")
hr_combined_table = Table.read(hr_combined_path, format="ecsv")

gaia_id_path = base_path / "gaia_source_id_lookup.csv"
gaia_id_schema = {"dr2_sourceid": pl.Int64, "dr3_sourceid": pl.Int64}
try:
    gaia_id_df = pl.read_csv(gaia_id_path, schema=gaia_id_schema)
except FileNotFoundError:
    gaia_id_df = pl.DataFrame(schema=gaia_id_schema)

simbad_heriarchy_cache = base_path / "cache_star_hierarchy"
simbad_heriarchy_cache.mkdir(parents=True, exist_ok=True)

# regex to parse ID of catalogs
re_gaiadr2 = re.compile(r"Gaia DR2\s+(\d+)")
re_gaiadr3 = re.compile(r"Gaia DR3\s+(\d+)")
re_hip = re.compile(r"HIP\s+(\d+)([A-Z]?)")
re_sao = re.compile(r"SAO\s+(\d+)")
re_hd = re.compile(r"HD\s+(\d+)")
re_hr = re.compile(r"HR\s+(\d+)")


def parse_cross_id(df: pl.DataFrame) -> pl.DataFrame:
    """"
    Parse cross-identification catalog numbers from a concatenated ID string.
    Extracts Hipparcos (HIP), SAO, HD, and HR catalog identifiers from a 
    combined ID column, along with the component designation. Converts 
    component letters (A, B, C, etc.) to numeric values (1, 2, 3, etc.),
    with 0 representing missing or null values.

    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame containing an 'ids' column with concatenated 
        catalog identifiers in a standard format.
    
    Returns
    -------
    pl.DataFrame
        DataFrame with four new integer columns appended:
        - hip : Hipparcos catalog number (Int64)
        - sao : SAO catalog number (Int64)
        - hd : HD catalog number (Int64)
        - hr : HR catalog number (Int64)
        - component : Component designation as integer (Int64),
                      where A=1, B=2, C=3, etc., and null=0
    """
    return df.with_columns(
    [
        pl.col("ids")
        .str.extract(re_hip.pattern, group_index=1)
        .fill_null(0)
        .cast(pl.Int64)
        .alias("hip")
        .fill_null(0),
        pl.col("ids").str.extract(re_sao.pattern).fill_null(0).cast(pl.Int64).alias("sao"),
        pl.col("ids").str.extract(re_hd.pattern).fill_null(0).cast(pl.Int64).alias("hd"),
        pl.col("ids").str.extract(re_hr.pattern).fill_null(0).cast(pl.Int64).alias("hr"),
        # replace the NaN with 0, A to 1, B to 2, etc
        pl.col("ids")
        .str.extract(re_hip.pattern, group_index=2)
        .alias("component")
        .map_elements(
            lambda x: 0 if x is None or x == "" else ord(x) - 64, return_dtype=pl.Int64
        ),
    ]
)

hip_combined_df = parse_cross_id(hip_combined_df)

# remove the known bad HIP stars
hip_to_remove = [
    1902,  # HIP 1902 is a globular cluster NGC 104
    24647,  # HIP 24647 is HIP 24648
    29116,  # HIP 29116 is HIP 29119
    35194,  # HIP 35194 is HIP 35195
    54948,  # HIP 54948 is a cluster of star (cant be sure which one they are referring, doubt Hipparcos has resolved them)
    88759,  # HIP 88759 is HIP 88762
    91906,  # HIP 91906 is HIP 91924
    98623,  # HIP 98623 is HIP 98625
]  # indices of rows
hip_combined_df = (
    hip_combined_df.filter(~pl.col("hip").is_in(hip_to_remove))
)

def extract_gaia_ids(
    df: pl.DataFrame, strict_dr: bool = False, cone_search: bool = False
) -> pl.DataFrame:
    """
    This function will extract the Gaia DR IDs from `ids` of a dataframe and append the result in a new column `source_id`.
    Try Gaia DR3 first, if not found try Gaia DR2.

    Parameters
    ----------
    df : pl.DataFrame
        DataFrame containing an ``ids`` column from which Gaia DR source IDs will be extracted.
    strict_dr : bool
        If True, only extract Gaia DR3 number. Otherwise, try to extract other Gaia DR number (except DR1) if Gaia DR3 is not found.
    cone_search : bool
        If True, try to cone search the source_id if not found
    """
    # more strict on resolving the source_id in case of binary stars
    df = df.with_columns(
        [
            pl.col("ids")
            .str.extract(re_gaiadr3.pattern)
            .fill_null(0)
            .cast(pl.Int64)
            .alias("source_id"),
            pl.col("ids")
            .str.extract(re_gaiadr2.pattern)
            .fill_null(0)
            .cast(pl.Int64)
            .alias("dr2_source_id"),
        ]
    )
    if (
        not cone_search and strict_dr
    ):  # getting the latest Gaia DR without additional resolving methods
        return df.drop("dr2_source_id")
    else:
        # in case we are upgrading to DR4, DR2 need to be upgraded to DR3
        # =========== the case where DR2 equals DR3, need to do this check before cone_search ===========
        if strict_dr:  # if any DR2 fit in this criteria
            # For rows where DR3 source_id is null but DR2 source_id exists, try to resolve via Gaia first with DR2 ID directly
            null_dr3_mask = (
                df["source_id"].is_null() & df["dr2_source_id"].is_not_null()
            )
            dr2_ids = df.filter(null_dr3_mask)["dr2_source_id"].drop_nulls().to_list()
            ids_str = ",".join(str(i) for i in dr2_ids)
            gaia_result = Gaia.launch_job(
                f"SELECT source_id FROM gaiadr3.gaia_source WHERE source_id IN ({ids_str})"
            ).results
            matched_dr3_ids = set()
            if len(gaia_result) > 0:
                matched_dr3_ids = set(int(r["source_id"]) for r in gaia_result)
            df = df.with_columns(
                pl.when(
                    pl.col("source_id").is_null()
                    & pl.col("dr2_source_id").is_in(list(matched_dr3_ids))
                )
                .then(pl.col("dr2_source_id"))
                .otherwise(pl.col("source_id"))
                .alias("source_id")
            )
    # if cone_search:
    #     # try to cone search in case SIMBAD missed the Gaia DR number, stricter in radius to prevent false positive
    #     null_source_id_mask = df["source_id"] == 0
    #     for idx, row in enumerate(df.filter(null_source_id_mask).iter_rows(named=True)):
    #         try:
    #             coord = SkyCoord(ra=row["ra"], dec=row["dec"], unit="deg")
    #             gaia_result = Gaia.cone_search(coord, radius=7.5 * u.arcsecond).results
    #             if len(gaia_result) > 0:
    #                 found_source_id = int(gaia_result[0]["source_id"])
    #                 df = df.with_columns(
    #                 pl.when(pl.col("ra") == row["ra"])
    #                 .then(pl.lit(found_source_id))
    #                 .otherwise(pl.col("source_id"))
    #                 .alias("source_id")
    #                 )
    #         except Exception as e:
    #             print(f"Cone search failed for row {idx}: {e}")
    #         continue
    #     # if len(gaia_source) > 0:
    #     #     return int(gaia_source[0]["source_id"])
    return df


def parse_result(simbad_df: pl.DataFrame) -> pl.DataFrame:
    """
    This function will parse the result from the query:
    - Remove rows that are not stars
    - Add additional rows to the table for stars that are in a binary system

    To attempt to resolve potential binary system, the logic is as follows:
    1. If the stars have V and J magnitude dimmer than 3.0, then it should have Gaia source_id
    3. If children have no Gaia source_id, then assume this is not a binary system. Just for whatever reason Gaia did not observe the parent star
    3. If childen have no Gaia source_id, then assume this is not a binary system. Just for whatever reason Gaia did not observe the parent star
    4. If the star has children with source_id, then add the children to the table and remove the parent star

    Parameters
    ----------
    simbad_df : polars.Dataframe
        Dataframe from the query result

    Returns
    -------
    polars.Dataframe
        Dataframe with additional rows from the query result
    """
    simbad_df = simbad_df.with_columns(pl.col("otype").fill_null("")).filter(
        ~pl.col("otype").str.contains("err")
    )
    # more strict on resolving the source_id in case of binary stars
    simbad_df = extract_gaia_ids(simbad_df, strict_dr=True)

    # TODO: resolve Gaia ID and use the new source_id col
    potential_binary_df = simbad_df.filter(
        ((pl.col("source_id") == 0) | (pl.col("source_id").is_null()))
        # & (pl.col("ids").is_not_null())  # need to at least have some id
        & (pl.col("V") > 3.0)
        & (pl.col("J") > 2.0)
    )
    print("Potential Binary: ", potential_binary_df.height)

    result = []

    # missing Gaia ID might indicate Gaia resolved binary not resolved by HIP
    for i, irow in enumerate(
        tqdm(
            potential_binary_df.iter_rows(named=True),
            total=potential_binary_df.height,
            desc="Resolving potential binary",
        )
    ):
        # if children have no SAO/HD/HR id, then use the parent id
        parent_sao = irow["sao"]
        parent_hd = irow["hd"]
        parent_hr = irow["hr"]

        query_str = custom_simbad.query_hierarchy(
            irow["main_id"],
            hierarchy="children",
            get_query_payload=True,
            detailed_hierarchy=False,
        )["QUERY"]
        tmp_hash = hashlib.sha1(query_str.encode("utf-8")).hexdigest()
        tmp_file = simbad_heriarchy_cache / f"{tmp_hash}.dat"
        # if not exist then query and save to file for caching
        if tmp_file.exists():
            children_df = pl.read_csv(tmp_file, schema=simbad_df.schema)
        else:
            temp_table: Table = custom_simbad.query_tap(query_str)
            children_df = temp_table.to_df("polars").with_columns(
                pl.lit(irow["user_specified_id"].strip()).alias("user_specified_id")
            )
            if len(children_df) <= 1:
                # in case some stars only have one component which is itself or no result
                continue
            children_df = extract_gaia_ids(children_df)
            # add letter to the end of the user_specified_id, so first row in temp will end with "A", second row in temp will end with "B", etc
            df_list = []
            for idx, row in enumerate(children_df.iter_rows(named=True)):
                user_specified_id = f"{row['user_specified_id']}{chr(65 + idx)}"
                ids = f"{user_specified_id}|{row['ids']}"
                # add temp_source_id for each row to the ids column if it is not null
                if row["source_id"] is not None:
                    ids = f"{ids}|Gaia DR3 {row['source_id']}"
                row["user_specified_id"] = user_specified_id
                row["ids"] = ids
                df_list.append(row)
            children_df = pl.DataFrame(df_list, schema=simbad_df.schema)
            children_df.write_csv(tmp_file)
        source_id = children_df["source_id"]
        # the second row and beyond should have source_id, if not then assume it is not a binary system
        # also if all of them has the same source_id, then assume it is not a binary system
        if source_id[1:].is_null().any() or source_id.is_unique().all():
            continue
        # add the parent SAO/HD/HR id to the ids column if it is not NaN
        if parent_sao:
            children_df = children_df.with_columns(
                (pl.col("ids") + f"|SAO {parent_sao}").alias("ids")
            )
        if parent_hd:
            children_df = children_df.with_columns(
                (pl.col("ids") + f"|HD {parent_hd}").alias("ids")
            )
        if parent_hr:
            children_df = children_df.with_columns(
                (pl.col("ids") + f"|HR {parent_hr}").alias("ids")
            )
        children_df.with_columns(pl.lit(irow["user_specified_id"].strip()).alias("user_specified_id"))
        # children_df.drop("matched_id")
        children_df = extract_gaia_ids(parse_cross_id(children_df), strict_dr=True)
        result.append(children_df)
        simbad_df = (
            simbad_df.with_row_index("_idx").filter(pl.col("_idx") != i).drop("_idx")
        )

    if len(result) > 0:
        # append list of result to simbad_df
        # print(simbad_df.schema)
        print(pl.concat(result, how="vertical").schema)
        simbad_df = pl.concat([simbad_df, *result], how="diagonal_relaxed")

    # drop empty rows (if missing ra, no use of keeping the row)
    simbad_df = simbad_df.filter(pl.col("ra").is_not_nan())
    # # try to resolve the source_id for the remaining stars
    # wo_source_id_mask = (
    #     (simbad_df["source_id"].is_null() | (simbad_df["source_id"] == 0))
    #     & (simbad_df["V"].fill_null(99.99) > 3.0)
    #     & simbad_df["ids"].is_not_null()  # need to at least have some id
    # )
    # # apply extract_gaia_number with cone_search=True to the remaining stars
    # star_wo_source_id_df = simbad_df.filter(wo_source_id_mask)
    # for idx, row in enumerate(
    #     tqdm(
    #         star_wo_source_id_df.iter_rows(named=True),
    #         total=star_wo_source_id_df.height,
    #         desc="Resolving Gaia DR3 ID xmatch",
    #     )
    # ):
    #     time.sleep(1.0)  # sleep for 1.0 seconds to prevent rate limit
    #     try:
    #         new_ids = f"{row['ids']}|Gaia DR3 {extract_gaia_ids(pl.DataFrame(row), cone_search=True)}"
    #         simbad_df = simbad_df.with_columns(
    #             pl.when(pl.col("main_id") == row["main_id"])
    #             .then(pl.lit(new_ids))
    #             .otherwise(pl.col("ids"))
    #             .alias("ids")
    #         )
    #         print(f"Reolsved {row['main_id']} with cone_search")
    #     except NameResolveError as e:
    #         print(f"Failed to resolve {row['main_id']} due to {e}")
    #         continue
    return simbad_df

hip_with_binary_df = parse_result(hip_combined_df)  # with binary stars resolved

# deal with Sirius A and Sirius B, because SIMBAD said Sirius does not have children
hip_with_binary_df = hip_with_binary_df.with_columns(
    pl.when(pl.col("main_id") == "* alf CMa")
    .then(pl.lit("HIP 32349A|") + pl.col("ids"))
    .otherwise(pl.col("ids"))
    .alias("ids")
)
result = custom_simbad.query_object("* alf CMa B")  # Sirius B
result["ids"] = "HIP 32349B|" + result["ids"]
result["user_specified_id"] = "* alf CMa B"
result.remove_column("matched_id")
restult_df = extract_gaia_ids(parse_cross_id(result.to_df("polars")), strict_dr=True)
# cast to Dataframe and back to prevent potential dtype incompatibility issue
simbad_table = pl.concat([hip_with_binary_df, restult_df])
simbad_table = simbad_table.drop_nulls(subset=["ids"])
simbad_table = (
    extract_gaia_ids(simbad_table)
    .rename({"source_id": "gaia_dr3"})
)
simbad_table.write_csv(base_path / "hip_processed_with_binary.dat")


# # drop rows with NaN "ids" column
# cross_id_df = cross_id_df.dropna(how="all")

# # find which SAO id is missing, should be SAO 1 to SAO max_sao_id. Only query the missing SAO id
# existing_sao_ids = set(cross_id_df["sao"].dropna().astype(int))
# all_sao_ids = set(range(1, len(sao_combined_table) + 1))
# missing_sao_idx = list(all_sao_ids - existing_sao_ids)
# # get the missing SAO id with corresponding row
# missing_sao = sao_combined_table[[i - 1 for i in missing_sao_idx]]
# sao_cross_id_df = parse_cross_id(missing_sao.to_df("polars"))
# cross_id_df = pd.concat([cross_id_df, sao_cross_id_df])


# # find which HD id is missing, should be HD 1 to HD max_hd_id. Only query the missing HD id
# existing_hd_ids = set(cross_id_df["hd"].dropna().astype(int))
# all_hd_ids = set(range(1, len(hd_combined_table) + 1))
# missing_hd_idx = list(all_hd_ids - existing_hd_ids)
# # get the missing HD id with corresponding row
# missing_hd = hd_combined_table[[i - 1 for i in missing_hd_idx]]
# cross_id_df = pd.concat([cross_id_df, parse_cross_id(missing_hd.to_df("polars"))])
