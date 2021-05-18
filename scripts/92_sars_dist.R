library(tidyverse)
library(ggplot2)
library(lubridate)


# ----- read the data & clean -> df
df = read.csv("data/92_sars/india_18thmay.txt", sep = "")
colnames(df) = c("c_date", "lineage", "")
df = df[,1:2]
df$c_date = as.Date(df$c_date)



# tally
df = df %>% 
  group_by(year_week = floor_date(x = c_date, unit = "1 week"), lineage) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n))

# force tobble -> data.frame
df = as.data.frame(df)
df$year_week = as.character.Date(df$year_week)

# handle variants with 0 at time points
for(w in unique(df$year_week)){
  for(l in unique(df$lineage)){
    if(nrow(df[df$lineage == l & df$year_week == w,]) == 0){ # if theres no data, impute 0
      df = rbind(df, data.frame(year_week = w,lineage = l,n = 0,percentage = 0))  
    }
  }
}

# as we floored by week, we can safely re-convert to date
df$year_week = as.Date(df$year_week)

# ----- plot
ggplot(df, aes(x = year_week, y = percentage, fill = as.factor(lineage))) +
  geom_area(alpha=0.6 , size=1, colour="black") +
  theme_classic()
