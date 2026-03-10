from astroquery.simbad import Simbad

custom_simbad = Simbad(timeout=3600)

cols_to_add = [
    "ids",
    "pmra",
    "pmdec",
    "plx_value",
    "plx_err",
    "U",
    "B",
    "V",
    "R",
    "I",
    "G",
    "J",
    "H",
    "K",
    "otype",
    "sp_type",
    "rvz_radvel",
    "rvz_err",
]

# you can add multiple columns at once, but the order of
# the columns remaining the same is important for caching
for col in cols_to_add:
    custom_simbad.add_votable_fields(col)
