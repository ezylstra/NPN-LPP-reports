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
# Returns dataframe with phenophase id, name, category, and class_id
ph <- npn_phenophases() %>% 
  select(-color) %>%
  rename(class_id = pheno_class_id,
         category = phenophase_category) %>% 
  data.frame() 
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

#- Intensities ----------------------------------------------------------------#




#- Abundances -----------------------------------------------------------------#



