# Compile information on phenophases, intensities, and abundances
# Erin Zylstra
# erinz@usanpn.org
# 2024-11-22

require(rnpn)
require(dplyr)
require(stringr)

rm(list = ls())

#- Phenophases ----------------------------------------------------------------#
# Extract list of phenophases in NPN database
ph <- npn_phenophases() %>% 
  select(-color) %>%
  rename(class_id = pheno_class_id,
         category = phenophase_category) %>% 
  data.frame() 
# Returns dataframe with phenophase id, name, category, and class_id
# Note: phenophase_name does not always match up with phenophase_description
# in the status-intensity datasets because the descriptions often tack on a 
# parenthetical reference to a taxonomic group (eg, forbs, birds).
# Note: many trailing white spaces in phenophase_name entries

# Extract information about pheno classes (higher-level order of phenophases)
phc <- npn_pheno_classes() %>% 
  rename(class_id = id,
         class_name = name,
         class_description = description) %>%
  data.frame()
# Returns dataframe with pheno class id, name, description, sequence

# Merge the two
ph <- left_join(ph, phc, by = "class_id")

# count(ph, category)
# filter(ph, category == "Activity")     # Animal categories 
# filter(ph, category == "Development")  # Animal, by stage (includes Dead categories!)
# filter(ph, category == "Method")       # Animals captured or baited
# filter(ph, category == "Reproduction") # Any animal reproduction-related stage (including eggs)
# filter(ph, category == "Flowers")      # Any flower-related stage (buds to senescense)
# filter(ph, category == "Fruits")       # Any fruit-related stage (first to drop)
# filter(ph, category == "Leaves")       # Any leaf-related stage (breaking buds to drop)
# filter(ph, category == "Needles")      # Any needle-reated stage (budburst to drop)
# filter(ph, category == "Pollen cones") # Any cone-related stage (first to pollen release)
# filter(ph, category == "Seed cones")   # Any seed cone-related stage (first to drop)

# For summarizing animal data, 
# Dead observation = Any dead observations
# Live observation = Any non-dead phenophase in Activity, Development, Reproduction categories
#TODO: Decide whether to leave observations in Methods category separate?
#TODO: Decide whether it's okay to lump eggs, pupae, caterpillars with adults

# For summarizing plant data, the classes seem too broad, but phenophase_name
# seems to detailed. Creating an intermediate classification for LPP reports 
# (mostly based on categories that were used previously by Janet Prevey):

ph %>%
  filter(category %in% c("Flowers", "Pollen cones", "Fruits", "Seed cones")) %>%
  arrange(class_id)
# Flower buds: class = Flowers or pollen cones
# Flowering: class = Open flowers or pollen cones, Pollen release
# Flower senescence: class = End of flowering
# Unripe fruits: class = Fruits or seed cones, Unripe fruits or seed cones
# Ripe fruits: class = Ripe fruits or seed cones, Recent fruit, cone or seed drop

ph %>%
  filter(category %in% c("Leaves", "Needles")) %>%
  arrange(class_id)
# Leaf out: class = Initial shoot or leaf growth, Young leaves or needles, Leaves or needles
# Leaf senescence: class = Colored leaves or needles, Falling leaves or needles

# Create new "pheno_group"
ph <- ph %>%
  mutate(pheno_group = case_when(
    grepl("dead|Dead", phenophase_name) ~ "Dead observation",
    category %in% c("Activity", 
                    "Development", 
                    "Reproduction") ~ "Live observation",
    category == "Method" ~ "Trapped/baited",
    class_name == "Flowers or pollen cones" ~ "Flower buds",
    class_name %in% c("Open flowers or pollen cones", 
                      "Pollen release") ~ "Flowering",
    class_name == "End of flowering" ~ "Flower senescence",
    class_name %in% c("Fruits or seed cones", 
                      "Unripe fruits or seed cones") ~ "Unripe fruits",
    class_name %in% c("Ripe fruits or seed cones", 
                      "Recent fruit, cone or seed drop") ~ "Ripe fruits",
    class_name %in% c("Initial shoot or leaf growth",
                      "Young leaves or needles",
                      "Leaves or needles") ~ "Leaf out",
    class_name %in% c("Colored leaves or needles",
                      "Falling leaves or needles") ~ "Leaf senescence",
    .default = NA
  ))

# Check:
count(ph, pheno_group, category, class_name)

# Write phenophase table to file
write.csv(ph, "phenophases.csv", row.names = FALSE)

#- Intensities and Abundances -------------------------------------------------#
# Extract list of intensity and abundance categroes in NPN database
ia <- npn_abundance_categories() 
# Returns tibble with category ids, each with a dataframe that lists the value
# ids and value names (which is what appears in intensity_value column in a
# status-intensity dataset)

# Remove blank entry 
ia <- ia %>%
  filter(!(is.na(category_name) | category_name == ""))

# Create a 2-dim dataframe that contains all the data
cat_ids <- select(ia, category_id, category_name) %>% data.frame()
cat_values <- ia$category_values
cat_values <- mapply(cbind, cat_values, "category_id" = cat_ids$category_id, 
                     SIMPLIFY = FALSE)
cat_values <- mapply(cbind, cat_values, "category_name" = cat_ids$category_name, 
                     SIMPLIFY = FALSE)
ia2 <- bind_rows(cat_values) %>%
  select(category_id, category_name, value_id, value_name, value_description)

# Identify whether values are a number, percent, or qualitative and 
# extract if not qualitative, extract bounding values to calculate a midpoint
# or representative number
value_names <- sort(unique(ia2$value_name))
values <- data.frame(value_names = value_names)
values <- values %>%
  mutate(value1 = NA,
         value2 = NA,
         type = case_when(
           str_detect(value_names, "%") ~ "percent",
           str_detect(value_names, "[0-9]") ~ "numeric",
           .default = "qualitative"
         ))
for (i in 1:nrow(values)) {
  if (str_detect(value_names[i], " to ")) {
    values[i, 2:3] <- str_split_fixed(value_names[i], " to ", 2)
    values[i, 2:3] <- as.numeric(str_remove(values[i, 2:3], ","))
  } else if (str_detect(value_names[i], "-")) {
    values[i, 2:3] <- str_split_fixed(value_names[i], "-", 2)
    values[i, 3] <- str_remove(values[i, 3], "%")
  } else if (str_detect(value_names[i], "% or more")) {
    values[i, 2:3] <- str_remove(value_names[i], "% or more")
  } else if (str_detect(value_names[i], "Less than ")) {
    values[i, 2] <- 0
    values[i, 3] <- str_remove(value_names[i], "Less than ")
    values[i, 3] <- str_remove(values[i, 3], "%")
  } else if (str_detect(value_names[i], "More than ")) {
    values[i, 2:3] <- str_remove(value_names[i], "More than ")
    values[i, 2] <- as.numeric(str_remove(values[i, 2], ",")) + 1
    values[i, 3] <- as.numeric(str_remove(values[i, 3], ",")) + 1
  }
}

#- PICK UP HERE ---------------------------------------------------------------#
# Calculate a midpoint... (or easy middlish value)






