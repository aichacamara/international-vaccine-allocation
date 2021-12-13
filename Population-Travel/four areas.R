#title: "Population and Travel DF"
#author: "Abraham Holleran"
#We group the countries of the world by high/low vaccination rates and high/low travel rates.

library(tidyverse)
library(lubridate)

#passengers is the number of passengers carried by country per year.
#This data was collected via International Civil Aviation Organization,
#Civil Aviation Statistics of the World and ICAO staff estimates. details here:
#https://data.worldbank.org/indicator/IS.AIR.PSGR?end=2019&name_desc=false&start=2019&view=bar
#This includes domestic travel.

#vaccinations is a dataframe detailing covid cases and vaccinations per country per day. It comes from ourworldindata, link here:
#https://ourworldindata.org/covid-vaccinations?country=USA

#pop is simply the population by region. It comes from our world in data.
#https://ourworldindata.org/grapher/world-population-by-world-regions-post-1820



passengers <- read.csv("passengers.csv")
vaccinations <- read.csv("vaccinations.csv")
pop <- read.csv("pop.csv")

#passengers is sparse. Let's replace the 2019 passenger number with the maximum
#of the 2000-2019 passengers. If there's no data between 2000-2019, we drop the country.

fill <- passengers[is.na(passengers$X2019),45:64] #Years 2000 to 2019
fill_max <- apply(X = fill, FUN = max, na.rm = T, MARGIN = 1)
passengers$X2019[is.na(passengers$X2019)] <- fill_max
count2019 <- passengers %>%
  filter(X2019 > 0) %>%
  select(1:2, X2019)


#Let's take the maximum population in each area.
area_pop <- pop %>%
  group_by(Code) %>%
  slice(which.max(Year)) %>%
  select(-c(3))
#Let's take the people vaccinated, per hundred, for each country.
#We use the most recent day for our calcuations.
top_v <- vaccinations %>%
  mutate(date = ymd(date)- ymd(20200101)) %>%
  group_by(iso_code) %>%
  slice(which.max(date)) %>%
  arrange(desc(people_vaccinated_per_hundred)) %>%
  select(iso_code, location, people_vaccinated_per_hundred) %>%
  ungroup()
top_v
#Let's join our country, vaccination, population, and travel data together.
#Then, we drop NAs and calculate the flights per person.

df <- top_v %>%
  inner_join(count2019, by = c("iso_code" = "Country.Code")) %>%
  inner_join(area_pop, by = c("iso_code" = "Code")) %>%
  select(location, vax = people_vaccinated_per_hundred, X2019, pop = Population..historical.estimates.) %>%
  drop_na() %>%
  mutate(flights_per_pop = X2019/pop)
df

#Next, we hand-pick 6 countries from each of 4 categories: high/low vaccination rates and high/low travel rates.

six_of_each <- df %>%
  slice(c(1,13,19,28,38,14,84,96, 106,68,66,75,4,7,20,30,12,54,118,130,131,132,135,139)) %>%
  mutate(area = rep(1:4, each = 6)) %>%
  arrange(area, location) %>%
  select(area, location, vax, flights_per_pop, pop)
six_of_each

#The lines between these 4 areas are blurred.
#A better way to pick would be to split the dataframe into four categories and systematically pick from them.

high_v_high_f <- df %>%
  filter(vax >= mean(vax), flights_per_pop >= mean(flights_per_pop)) %>%
  mutate(area = 1)
low_v_high_f <- df %>%
  filter(vax < mean(vax), flights_per_pop >= mean(flights_per_pop)) %>%
  mutate(area = 2)
high_v_low_f <- df %>%
  filter(vax >= mean(vax), flights_per_pop < mean(flights_per_pop)) %>%
  mutate(area = 3)
low_v_low_f <- df %>%
  filter(vax < mean(vax), flights_per_pop < mean(flights_per_pop)) %>%
  mutate(area = 4)
four_areas <- union(high_v_high_f, low_v_high_f) %>%
  union(high_v_low_f) %>%
  union(low_v_low_f)
four_areas

#Unfortunately, the low vax, high flight category only has 3 observations.
#More work needs to be done to split these systematically.