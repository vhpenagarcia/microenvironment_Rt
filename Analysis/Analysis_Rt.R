#### Analysis of data from Hobo

library(tidyverse)
library(lubridate)
library(matrixStats)

### Functions
source("Functions_Caminade.R") # To execute code with Caminade's method
#source("Functions_liu.R") # To execute code with Liu-Helmersson's method
#source("Functions_Mordecai.R") # To execute code with Mordecai's method

### Obtaining data
## File with data from the HOBOs
filename <- "cities_hobo_filtered.csv"
fullpath <- file.path("https://github.com/vhpenagarcia/microenvironment_Rt/blob/main/",filename)  # Can set path of your computer
data <- data.frame(read_csv(fullpath))
data <- data %>% separate(time, c("date","time"), " ")
data$date <- as.Date(data$date)

## Data from weather stations
filename2 <- "cities_meteorologic.csv"
fullpath2 <- file.path("https://github.com/vhpenagarcia/microenvironment_Rt/blob/main/",filename2)  # Can set path of your computer
data_ws <- data.frame(read_csv(fullpath2))
data_ws$Fecha <- as.Date(data_ws$Fecha)
data_ws <- data_ws %>% mutate(week = floor_date(Fecha, unit = "week", week_start = 7))

## Archivo epidemiol?gico
filename_epi <- "epid_data.csv"
fullpath_epi <- file.path("https://github.com/vhpenagarcia/microenvironment_Rt/blob/main/",filename_epi) # Can set path of your computer
data_epi <- data.frame(read_csv(fullpath_epi))
data_epi$date <- as.Date(data_epi$date)


### Convert HOBO data to rounded weekly data
week_hobo <- data %>% 
  mutate(week = floor_date(date, unit = "week", week_start = 7)) %>%
  group_by(week, location, city) %>%
  summarize(mean = mean(temp), sd = sd(temp))
head(week_hobo)

### Filtering HOBO data to each city
nei_data <- week_hobo %>% filter(city == "Neiva")
sin_data <- week_hobo %>% filter(city == "Sincelejo")
sol_data <- week_hobo %>% filter(city == "Soledad")

### filter and Convert data to rounded weekly data: A little wrangling
##Neiva
ws_neiva <- data_ws %>% 
  filter(DescripcionSerie == "Temperatura seca de las 700, 1300 y 1800" &
           Municipio == "Neiva" & week %in% nei_data$week) %>%
  group_by(week, Municipio) %>%
  summarize(mean = mean(Valor), sd = sd(Valor))
ws_neiva <- ws_neiva %>% rename(c("city" = "Municipio")) %>%
  select(week, city, mean, sd)
ws_neiva$location <- "ws"
ws_neiva <- ws_neiva[,c(1,5,2,3,4)]

##Sincelejo
ws_sin <- data_ws %>% 
  filter(DescripcionSerie == "Temperatura seca de las 700, 1300 y 1800" &
           Municipio == "Corozal" & week %in% sin_data$week) %>%
  group_by(week, Municipio) %>%
  summarize(mean = mean(Valor), sd = sd(Valor))
ws_sin <- ws_sin %>% rename(c("city" = "Municipio")) %>%
  select(week, city, mean, sd)
ws_sin$location <- "ws"
ws_sin$city <- ifelse(ws_sin$city == "Corozal", "Sincelejo", NA)
ws_sin <- ws_sin[, c(1,5,2,3,4)]

## Bind data bases into one for each city
nei_data <- rbind(nei_data,ws_neiva)
sin_data <- rbind(sin_data,ws_sin)


### Split by location to analyze them separately
##Neiva
nei_data_in <- nei_data %>% filter(location == "in")
nei_data_out <- nei_data %>% filter(location == "out")
nei_data_ws <- nei_data %>% filter(location == "ws")
##Sincelejo
sin_data_in <- sin_data %>% filter(location == "in")
sin_data_out <- sin_data %>% filter(location == "out")
sin_data_ws <- sin_data %>% filter(location == "ws")
##Soledad
sol_data_in <- sol_data %>% filter(location == "in")
sol_data_out <- sol_data %>% filter(location == "out")

#### Create matrices of cases accross 1000 simulations
### Obtaining initial cases
nei_index <- min(nei_data$week)
previous_week_nei <- which(data_epi$city == "Neiva" & data_epi$date == nei_index)-1
cases_initial_nei <- data_epi$arbovirus[previous_week_nei]

sin_index <- min(sin_data$week)
previous_week_sin <- which(data_epi$city == "Sincelejo" & data_epi$date == sin_index)-1
cases_initial_sin <- data_epi$arbovirus[previous_week_sin]

sol_index <- min(sol_data$week)
previous_week_sol <- which(data_epi$city == "Soledad" & data_epi$date == nei_index)-1
cases_initial_sol <- data_epi$arbovirus[previous_week_sol]

#### Running Monte Carlo simulations
set.seed(2, sample.kind = "Rounding") # For reproducibility
B <- 1000

# Depending if you are running code for Liu-Helmersson, Mordecai or Caminade 
# method, you should pick for running the appropriate one piece of code from the
# following sections:

### For the Caminade's method, run the following piece of code:
## Neiva
# In
nei_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_cam(Temp=nei_data_in$mean, sd_T = nei_data_in$sd, cases_initial = cases_initial_nei,
                 m1=731.721069, m2=668.480164)
}), ncol=B, nrow=nrow(nei_data_in))
# Out
nei_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_cam(Temp=nei_data_out$mean, sd_T = nei_data_out$sd, cases_initial = cases_initial_nei,
                 m1=731.721069, m2=668.480164)
}), ncol=B, nrow=nrow(nei_data_out))
# Ws
nei_matrix_cases_ws <- matrix(replicate(B, {
  temp_cases_cam(Temp=nei_data_ws$mean, sd_T = nei_data_ws$sd, cases_initial = cases_initial_nei,
                 m1=731.721069, m2=668.480164)
}), ncol=B, nrow=nrow(nei_data_ws))

## Sincelejo
# In
sin_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_cam(Temp=sin_data_in$mean, sd_T = sin_data_in$sd, cases_initial = cases_initial_sin,
                 m1=752.021729, m2=383.251526)
}), ncol=B, nrow=nrow(sin_data_in))
# Out
sin_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_cam(Temp=sin_data_out$mean, sd_T = sin_data_out$sd, cases_initial = cases_initial_sin,
                 m1=752.021729, m2=383.251526)
}), ncol=B, nrow=nrow(sin_data_out))
# Ws
sin_matrix_cases_ws <- matrix(replicate(B, {
  temp_cases_cam(Temp=sin_data_ws$mean, sd_T = sin_data_ws$sd, cases_initial = cases_initial_sin,
                 m1=752.021729, m2=383.251526)
}), ncol=B, nrow=nrow(sin_data_ws))

## Soledad
# In
sol_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_cam(Temp=sol_data_in$mean, sd_T = sol_data_in$sd, cases_initial = cases_initial_sol,
                 m1=894.557617, m2=177.526230)
}), ncol=B, nrow=nrow(sol_data_in))
# Out
sol_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_cam(Temp=sol_data_out$mean, sd_T = sol_data_out$sd, cases_initial = cases_initial_sol,
                 m1=894.557617, m2=177.526230)
}), ncol=B, nrow=nrow(sol_data_out))


### For the Liu-Helmersson´s method, run the following piece of code:
## Neiva
# In
nei_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_liu(Temp=nei_data_in$mean, sd_T = nei_data_in$sd, cases_initial = cases_initial_nei)
}), ncol=B, nrow=nrow(nei_data_in))
# Out
nei_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_liu(Temp=nei_data_out$mean, sd_T = nei_data_out$sd, cases_initial = cases_initial_nei)
}), ncol=B, nrow=nrow(nei_data_out))
# WS
nei_matrix_cases_ws <- matrix(replicate(B, {
  temp_cases_liu(Temp=nei_data_ws$mean, sd_T = nei_data_ws$sd, cases_initial = cases_initial_nei)
}), ncol=B, nrow=nrow(nei_data_ws))

## Sincelejo
#In
sin_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_liu(Temp=sin_data_in$mean, sd_T = sin_data_in$sd, cases_initial = cases_initial_sin)
}), ncol=B, nrow=nrow(sin_data_in))
# Out
sin_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_liu(Temp=sin_data_out$mean, sd_T = sin_data_out$sd, cases_initial = cases_initial_sin)
}), ncol=B, nrow=nrow(sin_data_out))
# WS
sin_matrix_cases_ws <- matrix(replicate(B, {
  temp_cases_liu(Temp=sin_data_ws$mean, sd_T = sin_data_ws$sd, cases_initial = cases_initial_sin)
}), ncol=B, nrow=nrow(sin_data_ws))

## Soledad
# In
sol_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_liu(Temp=sol_data_in$mean, sd_T = sol_data_in$sd, cases_initial = cases_initial_sol)
}), ncol=B, nrow=nrow(sol_data_in))
# Out
sol_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_liu(Temp=sol_data_out$mean, sd_T = sol_data_out$sd, cases_initial = cases_initial_sol)
}), ncol=B, nrow=nrow(sol_data_out))


### For the Mordecai's method, run the following piece of code:
## Neiva
# In
nei_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_mor(Temp=nei_data_in$mean, sd_T = nei_data_in$sd, cases_initial = cases_initial_nei)
}), ncol=B, nrow=nrow(nei_data_in))
# Out
nei_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_mor(Temp=nei_data_out$mean, sd_T = nei_data_out$sd, cases_initial = cases_initial_nei)
}), ncol=B, nrow=nrow(nei_data_out))
# WS
nei_matrix_cases_ws <- matrix(replicate(B, {
  temp_cases_mor(Temp=nei_data_ws$mean, sd_T = nei_data_ws$sd, cases_initial = cases_initial_nei)
}), ncol=B, nrow=nrow(nei_data_ws))

## Sincelejo
#In
sin_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_mor(Temp=sin_data_in$mean, sd_T = sin_data_in$sd, cases_initial = cases_initial_sin)
}), ncol=B, nrow=nrow(sin_data_in))
# Out
sin_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_mor(Temp=sin_data_out$mean, sd_T = sin_data_out$sd, cases_initial = cases_initial_sin)
}), ncol=B, nrow=nrow(sin_data_out))
# WS
sin_matrix_cases_ws <- matrix(replicate(B, {
  temp_cases_mor(Temp=sin_data_ws$mean, sd_T = sin_data_ws$sd, cases_initial = cases_initial_sin)
}), ncol=B, nrow=nrow(sin_data_ws))

## Soledad
# In
sol_matrix_cases_in <- matrix(replicate(B, {
  temp_cases_mor(Temp=sol_data_in$mean, sd_T = sol_data_in$sd, cases_initial = cases_initial_sol)
}), ncol=B, nrow=nrow(sol_data_in))
# Out
sol_matrix_cases_out <- matrix(replicate(B, {
  temp_cases_mor(Temp=sol_data_out$mean, sd_T = sol_data_out$sd, cases_initial = cases_initial_sol)
}), ncol=B, nrow=nrow(sol_data_out))


### Create data frames with matrix statistics and IC
##Neiva
# In
cases_in_nei <- data.frame(cases = round(rowMeans(nei_matrix_cases_in)),
                           sd_cases = rowSds(nei_matrix_cases_in)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))
# Out
cases_out_nei <- data.frame(cases = round(rowMeans(nei_matrix_cases_out)),
                            sd_cases = rowSds(nei_matrix_cases_out)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))
# Ws
cases_ws_nei <- data.frame(cases = round(rowMeans(nei_matrix_cases_ws)),
                           sd_cases = rowSds(nei_matrix_cases_ws)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))


## Sincelejo
# In
cases_in_sin <- data.frame(cases = round(rowMeans(sin_matrix_cases_in)),
                           sd_cases = rowSds(sin_matrix_cases_in)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))
# Out
cases_out_sin <- data.frame(cases = round(rowMeans(sin_matrix_cases_out)),
                            sd_cases = rowSds(sin_matrix_cases_out)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))
# WS
cases_ws_sin <- data.frame(cases = round(rowMeans(sin_matrix_cases_ws)),
                           sd_cases = rowSds(sin_matrix_cases_ws)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))

## Soledad
# In
cases_in_sol <- data.frame(cases = round(rowMeans(sol_matrix_cases_in)),
                           sd_cases = rowSds(sol_matrix_cases_in)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))
# Out
cases_out_sol <- data.frame(cases = round(rowMeans(sol_matrix_cases_out)),
                            sd_cases = rowSds(sol_matrix_cases_out)) %>%
  mutate(se_cases = sd_cases/sqrt(B), lower95_IC = round(cases - (1.96*se_cases)),
         upper95_IC = round(cases + (1.96*se_cases)))


### Create complete data frames
## Neiva
cases_nei_in <- cbind(nei_data_in,cases_in_nei)
cases_nei_out <- cbind(nei_data_out,cases_out_nei)
cases_nei_ws <- cbind(nei_data_ws, cases_ws_nei)
cases_nei <- rbind(cases_nei_in,cases_nei_out,cases_nei_ws)

## Sincelejo
cases_sin_in <- cbind(sin_data_in,cases_in_sin)
cases_sin_out <- cbind(sin_data_out,cases_out_sin)
cases_sin_ws <- cbind(sin_data_ws,cases_ws_sin)
cases_sin <- rbind(cases_sin_in,cases_sin_out,cases_sin_ws)

## Soledad
cases_sol_in <- cbind(sol_data_in,cases_in_sol)
cases_sol_out <- cbind(sol_data_out,cases_out_sol)
cases_sol <- rbind(cases_sol_in,cases_sol_out)

### Truncate negative values to 0
## Neiva
cases_nei$lower95_IC <- ifelse(cases_nei$lower95_IC < 0, 0, cases_nei$lower95_IC)
cases_nei$cases <- ifelse(cases_nei$cases < 0, 0, cases_nei$cases)
## Sincelejo
cases_sin$lower95_IC <- ifelse(cases_sin$lower95_IC < 0, 0, cases_sin$lower95_IC)
cases_sin$cases <- ifelse(cases_sin$cases < 0, 0, cases_sin$cases)
## Soledad
cases_sol$lower95_IC <- ifelse(cases_sol$lower95_IC < 0, 0, cases_sol$lower95_IC)
cases_sol$cases <- ifelse(cases_sol$cases < 0, 0, cases_sol$cases)


### plots for each city
## Neiva
cases_nei %>% ggplot(aes(week,cases, col = location)) +
  geom_ribbon(aes(ymin=lower95_IC, ymax=upper95_IC), linetype = 0, fill = "gray70") +
  geom_line() +
  geom_point()

## Sincelejo
cases_sin %>% ggplot(aes(week,cases, col = location)) +
  geom_ribbon(aes(ymin=lower95_IC, ymax=upper95_IC), linetype = 0, fill = "gray70") +
  geom_line() +
  geom_point()

## Soledad
cases_sol %>% ggplot(aes(week,cases, col = location)) +
  geom_ribbon(aes(ymin=lower95_IC, ymax=upper95_IC), linetype = 0, fill = "gray70") +
  geom_line() +
  geom_point()

### Summary statistics
## HOBOs
whole_df_cam <- rbind(cases_nei, cases_sin, cases_sol) # For Caminade's method
#whole_df_liu <- rbind(cases_nei, cases_sin, cases_sol) # For Liu-Helmersson's method
#whole_df_mor <- rbind(cases_nei, cases_sin, cases_sol) # For Mordecai's method
summary_cases <- whole_df_mor %>% group_by(city, location) %>%
  summarize(total_cases = sum(cases), lowerIC = sum(lower95_IC), 
            upperIC = sum(upper95_IC), avg_cases = mean(cases))
summary_cases

whole_df_cam$method <- "Caminade" # For Caminade's method
# whole_df_liu$method <- "Liu-Helmersson" # For Liu-Helmersson's method
# whole_df_mor$method <- "Mordecai" # For Mordecai's method


## Epidemics
nei_epid <- data_epi %>% filter(city == "Neiva" & date %in% nei_data$week)
sin_epid <- data_epi %>% filter(city == "Sincelejo" & date %in% sin_data$week)
sol_epid <- data_epi %>% filter(city == "Soledad" & date %in% sol_data$week)
cities_epid <- rbind(nei_epid, sin_epid, sol_epid)
epid_total <- cities_epid %>% group_by(city) %>% 
  summarize(total_cases = sum(arbovirus), avg_cases = mean(arbovirus))
epid_total
  


