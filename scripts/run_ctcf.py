import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.ctcf_analysis import parse_fu2008_sites, ctcf_ms_profile

sites = parse_fu2008_sites("data/fu2008_ctcf_sites.txt", which="occupied")