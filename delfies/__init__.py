from importlib import metadata

__version__ = metadata.version("delfies")

ID_DELIM = "__"
REGION_DELIM1 = ":"
REGION_DELIM2 = "-"
REGION_CLICK_HELP = (
    f"Region to focus on, format: 'contig{REGION_DELIM1}start{REGION_DELIM2}stop'."
)
