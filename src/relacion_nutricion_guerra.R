library(dplyr)

whi <- read.csv("data/world_hunger_index.csv", sep = '|')
population <- read.csv("data/WB_population/population.csv")
refugees <- read.csv("data/WB_refugees/refugees.csv")

refugees <- refugees[c('Country.Name', 'X2019')]
population <- population[c('Country.Name', 'X2019')]
names(refugees)[1] <- 'Country'
names(population)[1] <- 'Country'

head(whi, 10)
head(refugees, 10)
head(population, 10)

refugees$Country <- toupper(refugees$Country)
population$Country <- toupper(population$Country)
whi$Country <- toupper(whi$Country)

df <- merge(x=whi, y=refugees, by='Country')
df <- merge(x=df, y=population, by='Country')

df$Perc <- df$X2019.x / df$X2019.y

plot(log(df$Perc), df$X2019..14..18)
