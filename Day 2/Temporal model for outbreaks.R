library(INLA)
library(readr)
library(dplyr)
library(ggplot2)

#Load the data
#https://figshare.com/articles/dataset/A_global_dataset_of_pandemic-_and_epidemic-prone_disease_outbreaks/17207183
load("~/Downloads/17207183/data-monthly-updated-1996-2025/Last update/outbreaks_03102025.RData")
# --- Prepare dataset ---
# Filter out missing or incomplete entries
df <- outbreaks %>%
  filter(!is.na(Country), !is.na(Year)) %>%
  group_by(Country, Year) %>%
  summarise(outbreaks = n(), .groups = "drop")

# Encode as factors
df$country_id <- as.numeric(as.factor(df$Country))
df$year_id <- as.numeric(as.factor(df$Year))

# Check distribution
ggplot(df, aes(x = Year, y = outbreaks)) +
  geom_point(alpha = 0.6) +
  labs(title = "Outbreaks per Year (All Countries)", y = "Count")

# --- INLA model specification ---
formula <- outbreaks ~ 1 +
  f(country_id, model = "iid") +     # random intercept per country
  f(year_id, model = "rw1")          # smooth temporal trend

# --- Run INLA ---
res <- inla(formula,
            family = "poisson",
            data = df,
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))

# --- Results summary ---
summary(res)

# Fixed effects (intercept)
res$summary.fixed

# Random effects by country (top 10 countries with highest effects)
country_effects <- data.frame(
  Country = levels(as.factor(df$Country)),
  Mean = res$summary.random$country_id$mean
) %>%
  arrange(desc(Mean)) %>%
  head(10)

print(country_effects)

# --- Plot temporal random effect (trend over years) ---
year_trend <- data.frame(
  Year = sort(unique(df$Year)),
  Mean = res$summary.random$year_id$mean,
  Lower = res$summary.random$year_id$`0.025quant`,
  Upper = res$summary.random$year_id$`0.975quant`
)

ggplot(year_trend, aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.4) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Temporal Trend in Global Disease Outbreaks (INLA)",
       x = "Year", y = "Posterior Mean (RW1 effect)") +
  theme_minimal()

#Ethiopia specifically

eth_out <- outbreaks %>% filter(Country == "Ethiopia")

# Check how many outbreaks by disease type, by year
eth_summary <- eth_out %>%
  group_by(Year, Disease) %>%
  summarise(n_outbreaks = n(), .groups = "drop")

ggplot(eth_summary, aes(x = Year, y = n_outbreaks, color = Disease)) +
  geom_line() +
  labs(title = "Number of reported outbreaks in Ethiopia by disease & year",
       y = "Outbreak count", x = "Year") +
  theme_minimal()

# --- Prepare for Poisson model (example) ---
# We might count outbreaks per year in Ethiopia
eth_yearly <- eth_out %>%
  group_by(Year) %>%
  summarise(count = n(), .groups = "drop")

# Make sure Year is numeric
eth_yearly <- eth_yearly %>% mutate(Year = as.integer(Year))

# Plot
ggplot(eth_yearly, aes(x = Year, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title="Outbreaks per year in Ethiopia", x="Year", y="Count")

formula <- eth_yearly$count ~ 1 +   
  f(Year, model = "rw2")          # smooth temporal trend


res <- inla(formula,
            family = "poisson",
            data = eth_yearly,
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))


# --- Plot temporal random effect (trend over years) ---
year_trend <- data.frame(
  Year = sort(unique(eth_yearly$Year)),
  Mean = res$summary.random$Year$mean,
  Lower = res$summary.random$Year$`0.025quant`,
  Upper = res$summary.random$Year$`0.975quant`
)

ggplot(year_trend, aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.4) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Ethiopia: Trend in Global Disease Outbreaks (INLA)",
       x = "Year", y = "Posterior Mean (RW2 effect)") +
  theme_minimal()

# --- Posterior predictive check ---
plot(res$summary.fitted.values$mean, eth_yearly$count,
     xlab = "Fitted values", ylab = "Observed",
     main = "Observed vs Fitted Outbreak Counts (INLA)")
abline(0, 1, col = "red", lwd = 2)


#Adding space

df <- outbreaks %>%
  filter(!is.na(Country), !is.na(Year)) %>%
  group_by(Country, Year) %>%
  summarise(outbreaks = n(), .groups = "drop")

# --- Load world shapefile (for adjacency) ---
# You can use rnaturalearth data (included globally)
# install.packages("rnaturalearth")
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Match dataset countries to map
df <- df %>% filter(Country %in% world$name_long)
world <- world %>% filter(name_long %in% df$Country)

# Merge IDs
world <- world %>%
  mutate(country_id = as.numeric(as.factor(name_long)))

df <- df %>%
  mutate(country_id = as.numeric(as.factor(Country)))

# --- Build adjacency matrix ---
nb <- poly2nb(world, queen = TRUE)
nb2INLA("world_adj.graph", nb)
g <- inla.read.graph(filename = "world_adj.graph")

# --- Encode year as factor ---
df$year_id <- as.numeric(as.factor(df$Year))

# --- Model formula ---
formula <- outbreaks ~ 1 +
  f(country_id, model = "bym2", graph = g) +   # spatially structured + unstructured
  f(year_id, model = "rw1")                    # temporal trend

# --- Run INLA ---
res <- inla(formula,
            family = "poisson",
            data = df,
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))

# --- Summaries ---
summary(res)

# --- Plot spatial random effects ---
spatial_effects <- data.frame(
  country_id = res$summary.random$country_id$ID,
  mean = res$summary.random$country_id$mean
)

# Merge with shapefile
world_map <- merge(world, spatial_effects, by.x = "country_id", by.y = "country_id", all.x = TRUE)

ggplot(world_map) +
  geom_sf(aes(fill = mean)) +
  scale_fill_viridis_c(option = "viridis") +
  labs(title = "Spatial random effects (BYM2-INLA)",
       fill = "Posterior mean") +
  theme_minimal()

# --- Plot temporal random effect ---
year_trend <- data.frame(
  Year = sort(unique(df$Year)),
  Mean = res$summary.random$year_id$mean,
  Lower = res$summary.random$year_id$`0.025quant`,
  Upper = res$summary.random$year_id$`0.975quant`
)

ggplot(year_trend, aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.4) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Temporal trend in outbreak intensity (RW1-INLA)",
       x = "Year", y = "Posterior mean (v_t)") +
  theme_minimal()

