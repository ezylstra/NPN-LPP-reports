# Compile information on phenophases, intensities, and abundances
# Erin Zylstra
# erinz@usanpn.org
# 2024-11-22

require(rnpn)
require(dplyr)
require(stringr)

rm(list = ls())

# Logical indicating whether to replace csvs that list phenophase, intensity, 
# and abundance values if they already exist
replace <- FALSE

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
# Returns dataframe with pheno class id, name, description, sequence.
  # class_id = 1:13 are the phenophase classes associated with plants. Sequence
  # is just class_id * 10, and both generally reflect the order in which things
  # occur within a category (eg, Leaves). Based on exploration of 
  # status-intensity data for a particular plant species and class_id, there is
  # only one phenophase per phenophase class for a given species, which is
  # helpful. 
  # class_id's > 60 are associated with animals. Many more phenophase classes
  # for animals than plants and the relationship between class_id and sequence
  # is not obvious.

# Merge the two
ph <- left_join(ph, phc, by = "class_id") %>%
  arrange(category, class_id, phenophase_id)

# Look at phenophases, classes associated with each category
  count(ph, category)
  filter(ph, category == "Development")  # Animal, by stage (includes Dead categories!)
  filter(ph, category == "Method")       # Animals captured or baited
  filter(ph, category == "Reproduction") # Any animal reproduction-related stage (including eggs)
  filter(ph, category == "Activity")     # All other animal phenophases
  filter(ph, category == "Flowers")      # Any flower-related stage (buds to senescense)
  filter(ph, category == "Fruits")       # Any fruit-related stage (first to drop)
  filter(ph, category == "Leaves")       # Any leaf-related stage (breaking buds to drop)
  filter(ph, category == "Needles")      # Any needle-reated stage (budburst to drop)
  filter(ph, category == "Pollen cones") # Any cone-related stage (first to pollen release)
  filter(ph, category == "Seed cones")   # Any seed cone-related stage (first to drop)

# Proposed pheno "groups" for animals
  # Dead observation = Any dead observations
  # Live observation = Any non-dead phenophase in Activity, Development, Reproduction categories
  #TODO: Decide whether to leave observations in Methods category separate?
  #TODO: Decide whether it's okay to lump eggs, pupae, caterpillars with adults

# Prosed pheno "groups" for plants
  # For summarizing plant data, the classes seem too broad, but phenophase_name
  # seems too detailed. Creating an intermediate classification for LPP reports 
  # (largely based on categories that were used previously by Janet Prevey):
  
  ph %>%
    filter(category %in% c("Leaves", "Needles")) %>%
    arrange(class_id)
  ph %>%
    filter(category %in% c("Flowers", "Pollen cones", "Fruits", "Seed cones")) %>%
    arrange(class_id)
  
  # Leaves: classes 1:3 = Initial shoot or leaf growth, Young leaves or needles, Leaves or needles (Janet called this Leaf out)
  # Leaf senescence: class 4:5 = Colored leaves or needles, Falling leaves or needles
  # Flowers: class 6 = Flowers or pollen cones (Janet called this Flower buds)
  # Open flowers: classes 7:8 = Open flowers or pollen cones, Pollen release (Janet called this Flowering)
  # Flower senescence: class 9 = End of flowering
  # Fruits: classes 10:11 = Fruits or seed cones, Unripe fruits or seed cones (Janet called this Unripe fruit)
  # Ripe fruits: classes 12:13 = Ripe fruits or seed cones, Recent fruit, cone or seed drop

# Create new "pheno_group" 
ph <- ph %>%
  mutate(pheno_group = case_when(
    grepl("dead|Dead", phenophase_name) ~ "Dead observation",
    category %in% c("Activity", 
                    "Development", 
                    "Reproduction") ~ "Live observation",
    category == "Method" ~ "Trapped/baited",
    class_name == "Flowers or pollen cones" ~ "Flowers",
    class_name %in% c("Open flowers or pollen cones", 
                      "Pollen release") ~ "Open flowers",
    class_name == "End of flowering" ~ "Flower senescence",
    class_name %in% c("Fruits or seed cones", 
                      "Unripe fruits or seed cones") ~ "Fruits",
    class_name %in% c("Ripe fruits or seed cones", 
                      "Recent fruit, cone or seed drop") ~ "Ripe fruits",
    class_name %in% c("Initial shoot or leaf growth",
                      "Young leaves or needles",
                      "Leaves or needles") ~ "Leaves",
    class_name %in% c("Colored leaves or needles",
                      "Falling leaves or needles") ~ "Leaf senescence",
    .default = NA
  ))

# Check:
count(ph, pheno_group, category, class_name)

# Write phenophase table to file
pheno_file <- "phenophases.csv"
if (!file.exists(pheno_file) | replace) {
  write.csv(ph, pheno_file, row.names = FALSE)
}

#- Intensities ----------------------------------------------------------------#
# Extract list of intensity categories in NPN database
ia <- npn_abundance_categories() 
# Returns tibble with category ids, each with a dataframe that lists the value
# ids and value names (which is what appears in intensity_value column in a
# status-intensity dataset)

# Remove any blank entries 
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

# Identify whether values are a number/count, percent, or qualitative and 
# if not qualitative, extract bounding values to calculate a midpoint
# or representative number
value_name <- sort(unique(ia2$value_name))
values <- data.frame(value_name = value_name)
values <- values %>%
  mutate(value1 = NA,
         value2 = NA,
         type = case_when(
           str_detect(value_name, "%") ~ "percent",
           str_detect(value_name, "[0-9]") ~ "number",
           .default = "qualitative"
         ))
for (i in 1:nrow(values)) {
  if (str_detect(value_name[i], " to ")) {
    values[i, 2:3] <- str_split_fixed(value_name[i], " to ", 2)
    values[i, 2:3] <- as.numeric(str_remove(values[i, 2:3], ","))
  } else if (str_detect(value_name[i], "-")) {
    values[i, 2:3] <- str_split_fixed(value_name[i], "-", 2)
    values[i, 3] <- str_remove(values[i, 3], "%")
  } else if (str_detect(value_name[i], "% or more")) {
    values[i, 2:3] <- str_remove(value_name[i], "% or more")
  } else if (str_detect(value_name[i], "Less than ")) {
    values[i, 2] <- 0
    values[i, 3] <- str_remove(value_name[i], "Less than ")
    values[i, 3] <- str_remove(values[i, 3], "%")
  } else if (str_detect(value_name[i], "More than ")) {
    values[i, 2:3] <- str_remove(value_name[i], "More than ")
    values[i, 2] <- as.numeric(str_remove(values[i, 2], ",")) + 1
    values[i, 3] <- as.numeric(str_remove(values[i, 3], ",")) + 1
  }
}

values <- values %>%
  mutate_at(c("value1", "value2"), as.numeric) %>%
  arrange(type, value1)

# Assigning a middle-ish value for each range (For count ranges, keeping it to 
# nice numbers like 5, 50, 500, and 5000)
values <- values %>%
  mutate(mag = nchar(value1) - 1) %>%
  mutate(value = case_when(
    value1 == value2 ~ round(value1),
    type == "number" & value1 == 0 ~ 1,
    type == "number" & value1 != 0 ~ 
      plyr::round_any(rowMeans(across(value1:value2)), 5 * (10 ^ mag)),
    type == "percent" ~ round(rowMeans(across(value1:value2))),
    .default = NA
  )) %>%
  select(-mag)

# Append middle values to full list of intensities
ia2 <- ia2 %>%
  left_join(select(values, value_name, type, value), by = "value_name") %>%
  arrange(category_id, value_id)

# For qualitative categories, create integer values from 1 to n
#TODO: Decide whether we want to assign more meaningful values to qualitative categories
ia2 <- ia2 %>%
  group_by(category_id) %>%
  mutate(value = case_when(
    is.na(value) ~ value_id - min(value_id) + 1,
    .default = value
  ))
# There are a few categories with "(peak)" in the name that have number ranges
# along with a "peak" value_id without a number, which has been labeled as
# quantitative using the rules above. Not sure how often these categories are 
# being used, but for now, make value 1 greater than the largest number in 
# category
for (i in 2:nrow(ia2)) {
  ia2$value[i] <- ifelse(ia2$category_id[i] == ia2$category_id[i - 1] & 
                         ia2$type[i] != ia2$type[i - 1],
                         ia2$value[i - 1] + 1, ia2$value[i])
} 
ia2 <- data.frame(ia2)

ia_file <- "intensities.csv"
if (!file.exists(ia_file) | replace) {
  write.csv(ia2, ia_file, row.names = FALSE)
}
