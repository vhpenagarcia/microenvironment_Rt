#### Finding relationship between indoor and outdoor temperatures

library(tidyverse)
library(lubridate)
library(broom)

### Getting the data
filename <- "cities_hobo_filtered.csv"
fullpath <- file.path("https://github.com/vhpenagarcia/microenvironment_Rt/blob/main/",filename)
data <- data.frame(read_csv(fullpath))
head(data) # To explore a bit the data


#### Exploring regression on months
data <- data %>% separate(time, c("date","time"), " ")
data$date <- as.Date(data$date)

monthly_data <- data %>% 
  mutate(month = floor_date(date, unit = "month")) 
head(monthly_data)

### monthly regressions: spreading data
indoor_t_month <- monthly_data %>% filter(location == "in") %>% 
  select(city, date, month, temp)
colnames(indoor_t_month)[colnames(indoor_t_month) == "temp"] <- "indoor"
head(indoor_t_month)
outdoor_t_month <- monthly_data %>% filter(location == "out") %>%
  select(city, date, month, temp)
colnames(outdoor_t_month)[colnames(outdoor_t_month) == "temp"] <- "outdoor"
head(outdoor_t_month)

tempers_month <- inner_join(indoor_t_month,outdoor_t_month, by = c("city","date","month"))
head(tempers_month)

### Obtaining monthly models from averaged daily data
monthly_lm <- tempers_month %>%  
  group_by(city,date,month) %>%
  summarize(indoor=mean(indoor), outdoor=mean(outdoor)) %>%
  group_by(city,month) %>%
  do(tidy(lm(indoor ~ outdoor, data = .), conf.int = TRUE)) %>%
  select(city, month, estimate, conf.low, conf.high)
head(monthly_lm)

monthly_lm$coefficient <- rep(c("intercept","slope"), 32)
monthly_lm_summ <- monthly_r2_2 %>% 
  select(city,month,estimate,coefficient) %>% 
  spread(., coefficient, estimate)
head(monthly_lm_summ)

## Obtaining R2
monthly_lm_r2 <- tempers_month %>%  
  group_by(city,date,month) %>%
  summarize(indoor=mean(indoor), outdoor=mean(outdoor)) %>%
  group_by(city,month) %>%
  do(tidy(summary(lm(indoor ~ outdoor, data = .))$r.squared))
head(monthly_lm_r2)
colnames(monthly_lm_r2)[colnames(monthly_lm_r2) == "x"] <- "R_squared"

monthly_lm_summary <- inner_join(monthly_lm_summ,monthly_lm_r2, by = c("city","month"))
head(monthly_lm_summary)


### Obtaining plot
## Filtering july
tempers_month %>% 
  filter(month > "2019-07-01") %>%
  mutate(monthyear = format(month,format="%Y - %m")) %>%
  group_by(city, date, monthyear) %>% summarize(indoor = mean(indoor),
                                                outdoor = mean(outdoor)) %>%
  ggplot(aes(outdoor,indoor, col = city)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  facet_wrap(.~monthyear)
