###########################################################################################
########################### Workflow description
###########################################################################################

# For the whole project, the following workflow is applied (-> we are here)
# 
# 1. get data from institutes and departments on their publications
# 2. Run 00_preproccess_institute_data script on the institution
# 3. Upload to Scival by matching DOIs (record nonmatching DOIs)
# 4. Generate standard scientific evaluation based in Scival
#
# 5. Download a.) Bibliography of whole institute from Scopus as CSV, b.) Scival records of departments as CSV
# 6. Filter Institute from Scopus by Department in SciVal by EID
# 7. Run 11_preprocessing_seed to identify the department seeds
# 8. Take seed EIDs and query the 2k most relevant articles to each
# 9. Run 12_preprocessing_all
# 10. Run 91_descriptives on it.

###########################################################################################
########################### Preamble
###########################################################################################

### Generic preamble
rm(list=ls())
set.seed(1337)

library(tidyverse)

###########################################################################################
########################### Preamble
###########################################################################################


###########################################################################################
########################### Create department 
###########################################################################################

# 1. get data from institutes and departments on their publications
# --> 2. Run 00_preproccess_institute_data script on the institution

source('R/00_preprocess_institute.R')

###########################################################################################
########################### 1. Select Seed
###########################################################################################

# 6. Filter Institute from Scopus by Department in SciVal by EID
# 7. Run 11_preprocessing_seed to identify the department seeds

source('R/11_preprocess_seed.R')


###########################################################################################
########################### Field Mapping
###########################################################################################
 
# 9. Run 12_preprocessing_all
# --> 10. Run 91_descriptives on it.

source('R/12_preprocess_all.R')

###########################################################################################
########################### Report creation
###########################################################################################

# 9. Run 12_preprocessing_all
# --> 11. Create all reports

##########
### Field mapping general
##########

rmarkdown::render("R/91_descriptives_general.Rmd", output_file = '../output/field_mapping/field_mapping_general.html')

##########
### Field mapping bibliometric categorization
##########

rm(list=ls())

rmarkdown::render("R/92_descriptives_mapping.Rmd", output_file = '../output/field_mapping/field_mapping')





