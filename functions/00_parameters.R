
############################################################################
#  Main article selection
############################################################################

# Initial read parameters
TC_year_min <- 1 # Min number of citations
PY_min <- 2016 # Start year
PY_max <- 2021 # End year

# Select variables to keep
#vars <- c("AU", "Author.s..ID", "TI", "PY", "SO", "VL", "IS", "PP", "TC", "DI", "Affiliations", "C1", "AB", "DE", "ID", "FU", "FX",
#          "CR", "RP", "LA", "JI", "DT", "DB", "UT", "J9", "AU_UN", "AU1_UN", "AU_UN_NR", "SR_FULL", "SR")

###########################################################################
# Network Biblio
###########################################################################

# Initial Filter
cutof_edge_bib <- 2
cutof_node_bib <- 5

cutof_edge_pct_bib <- 0.05
cutof_node_pct_bib <- 0.25

# community detection
com_max_bib <- 12

############################################################################
#  Network Co-Citation
############################################################################

# Initial Filter
cutof_edge_cit <- 2
cutof_node_cit <- 5

cutof_edge_pct_cit <- 0.05
cutof_node_pct_cit <- 0.25

# community detection
com_max_cit <- 12

############################################################################
#  Datavix
############################################################################

pal_kb = 'Set3'
pal_ra = 'Paired'
pal_tp = NULL
