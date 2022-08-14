
library(tidyverse)

## To work with this code it is necessary to  run before "Analysis_Rt" with the
# three methods and have generated a data_frame from each of them named 
# whole_df_method

## putting data frames together
whole_df <- rbind(whole_df_liu,whole_df_cam,whole_df_mor)
whole_df %>% ggplot(aes(week,cases, col=location)) +
  geom_ribbon(aes(ymin=lower95_IC, ymax=upper95_IC), fill="gray70") +
  geom_line()+
  geom_point()+
  facet_grid(method~city)+
  scale_y_continuous(trans = "log10")


### Plots for each city
whole_df %>% filter(city == "Neiva") %>%
  ggplot(aes(week, cases, col = location)) +
  geom_ribbon(aes(ymin = lower95_IC, ymax = upper95_IC), linetype = 0, 
              fill = "gray70") +
  geom_line() +
  facet_wrap(.~method, scales = "free_y", dir = "v", strip.position = "right") +
  ylab("Estimated cases") +
  xlab("Time")


whole_df %>% filter(city == "Sincelejo") %>%
  ggplot(aes(week, cases, col = location)) +
  geom_ribbon(aes(ymin = lower95_IC, ymax = upper95_IC), linetype = 0, 
              fill = "gray70") +
  geom_line() +
  facet_wrap(.~method, scales = "free_y", dir = "v", strip.position = "right") +
  ylab("Estimated cases") +
  xlab("Time")


whole_df %>% filter(city == "Soledad") %>%
  ggplot(aes(week, cases, col = location)) +
  geom_ribbon(aes(ymin = lower95_IC, ymax = upper95_IC), linetype = 0, 
              fill = "gray70") +
  geom_line() +
  facet_wrap(.~method, scales = "free_y", dir = "v", strip.position = "right") +
  ylab("Estimated cases") +
  xlab("Time")


### Epidemiological data
summary_cases <- whole_df %>% group_by(city,method,location) %>%
  summarize(total_cases = sum(cases), lower_IC = sum(lower95_IC), 
            upper_IC = sum(upper95_IC), avg_cases = mean(cases))
summary_cases


## Epidemics
nei_epid <- data_epi %>% filter(city == "Neiva" & date %in% nei_data$week)
sin_epid <- data_epi %>% filter(city == "Sincelejo" & date %in% sin_data$week)
sol_epid <- data_epi %>% filter(city == "Soledad" & date %in% sol_data$week)
cities_epid <- rbind(nei_epid, sin_epid, sol_epid)
epid_total <- cities_epid %>% group_by(city) %>% 
  summarize(total_cases = sum(arbovirus), avg_cases = mean(arbovirus))
epid_total


### Comparison of modelled cases
liu_nei_in <- whole_df_liu %>% filter(city == "Neiva" & location == "in")
liu_nei_out <- whole_df_liu %>% filter(city == "Neiva" & location == "out")
liu_nei_ws <- whole_df_liu %>% filter(city == "Neiva" & location == "ws")
liu_sin_in <- whole_df_liu %>% filter(city == "Sincelejo" & location == "in")
liu_sin_out <- whole_df_liu %>% filter(city == "Sincelejo" & location == "out")
liu_sin_ws <- whole_df_liu %>% filter(city == "Sincelejo" & location == "ws")
liu_sol_in <- whole_df_liu %>% filter(city == "Soledad" & location == "in")
liu_sol_out <- whole_df_liu %>% filter(city == "Soledad" & location == "out")

cam_nei_in <- whole_df_cam %>% filter(city == "Neiva" & location == "in")
cam_nei_out <- whole_df_cam %>% filter(city == "Neiva" & location == "out")
cam_nei_ws <- whole_df_cam %>% filter(city == "Neiva" & location == "ws")
cam_sin_in <- whole_df_cam %>% filter(city == "Sincelejo" & location == "in")
cam_sin_out <- whole_df_cam %>% filter(city == "Sincelejo" & location == "out")
cam_sin_ws <- whole_df_cam %>% filter(city == "Sincelejo" & location == "ws")
cam_sol_in <- whole_df_cam %>% filter(city == "Soledad" & location == "in")
cam_sol_out <- whole_df_cam %>% filter(city == "Soledad" & location == "out")

mor_nei_in <- whole_df_mor %>% filter(city == "Neiva" & location == "in")
mor_nei_out <- whole_df_mor %>% filter(city == "Neiva" & location == "out")
mor_nei_ws <- whole_df_mor %>% filter(city == "Neiva" & location == "ws")
mor_sin_in <- whole_df_mor %>% filter(city == "Sincelejo" & location == "in")
mor_sin_out <- whole_df_mor %>% filter(city == "Sincelejo" & location == "out")
mor_sin_ws <- whole_df_mor %>% filter(city == "Sincelejo" & location == "ws")
mor_sol_in <- whole_df_mor %>% filter(city == "Soledad" & location == "in")
mor_sol_out <- whole_df_mor %>% filter(city == "Soledad" & location == "out")

wx_liu_nei_inxout <-wilcox.test(liu_nei_in$cases,liu_nei_out$cases)
wx_liu_nei_inxws <-wilcox.test(liu_nei_in$cases,liu_nei_ws$cases)
wx_liu_nei_outxws <-wilcox.test(liu_nei_out$cases,liu_nei_ws$cases)
wx_liu_sin_inxout <-wilcox.test(liu_sin_in$cases,liu_sin_out$cases)
wx_liu_sin_inxws <-wilcox.test(liu_sin_in$cases,liu_sin_ws$cases)
wx_liu_sin_outxws <-wilcox.test(liu_sin_out$cases,liu_sin_ws$cases)
wx_liu_sol_inxout <-wilcox.test(liu_sol_in$cases,liu_sol_out$cases)

wx_cam_nei_inxout <-wilcox.test(cam_nei_in$cases,cam_nei_out$cases)
wx_cam_nei_inxws <-wilcox.test(cam_nei_in$cases,cam_nei_ws$cases)
wx_cam_nei_outxws <-wilcox.test(cam_nei_out$cases,cam_nei_ws$cases)
wx_cam_sin_inxout <-wilcox.test(cam_sin_in$cases,cam_sin_out$cases)
wx_cam_sin_inxws <-wilcox.test(cam_sin_in$cases,cam_sin_ws$cases)
wx_cam_sin_outxws <-wilcox.test(cam_sin_out$cases,cam_sin_ws$cases)
wx_cam_sol_inxout <-wilcox.test(cam_sol_in$cases,cam_sol_out$cases)

wx_mor_nei_inxout <-wilcox.test(mor_nei_in$cases,mor_nei_out$cases)
wx_mor_nei_inxws <-wilcox.test(mor_nei_in$cases,mor_nei_ws$cases)
wx_mor_nei_outxws <-wilcox.test(mor_nei_out$cases,mor_nei_ws$cases)
wx_mor_sin_inxout <-wilcox.test(mor_sin_in$cases,mor_sin_out$cases)
wx_mor_sin_inxws <-wilcox.test(mor_sin_in$cases,mor_sin_ws$cases)
wx_mor_sin_outxws <-wilcox.test(mor_sin_out$cases,mor_sin_ws$cases)
wx_mor_sol_inxout <-wilcox.test(mor_sol_in$cases,mor_sol_out$cases)

method <- rep(c("Caminade","Liu-Helmersson","Mordecai"), each = 7)
city <- rep(c(rep("Neiva",3),rep("Sincelejo",3),"Soledad"),3)
comparison1 <-rep(c(rep(c("in","in","out"),2),"in"),3)
comparison2 <- rep(c(rep(c("out","ws","ws"),2),"out"),3)
comparison <- paste(comparison1,comparison2,sep = "_vs_")
W <- c(wx_cam_nei_inxout$statistic,wx_cam_nei_inxws$statistic,wx_cam_nei_outxws$statistic,
       wx_cam_sin_inxout$statistic,wx_cam_sin_inxws$statistic,wx_cam_sin_outxws$statistic,
       wx_cam_sol_inxout$statistic,
       wx_liu_nei_inxout$statistic,wx_liu_nei_inxws$statistic,wx_liu_nei_outxws$statistic,
       wx_liu_sin_inxout$statistic,wx_liu_sin_inxws$statistic,wx_liu_sin_outxws$statistic,
       wx_liu_sol_inxout$statistic,
       wx_mor_nei_inxout$statistic,wx_mor_nei_inxws$statistic,wx_mor_nei_outxws$statistic,
       wx_mor_sin_inxout$statistic,wx_mor_sin_inxws$statistic,wx_mor_sin_outxws$statistic,
       wx_mor_sol_inxout$statistic)
p_value <- c(wx_cam_nei_inxout$p.value,wx_cam_nei_inxws$p.value,wx_cam_nei_outxws$p.value,
             wx_cam_sin_inxout$p.value,wx_cam_sin_inxws$p.value,wx_cam_sin_outxws$p.value,
             wx_cam_sol_inxout$p.value,
             wx_liu_nei_inxout$p.value,wx_liu_nei_inxws$p.value,wx_liu_nei_outxws$p.value,
             wx_liu_sin_inxout$p.value,wx_liu_sin_inxws$p.value,wx_liu_sin_outxws$p.value,
             wx_liu_sol_inxout$p.value,
             wx_mor_nei_inxout$p.value,wx_mor_nei_inxws$p.value,wx_mor_nei_outxws$p.value,
             wx_mor_sin_inxout$p.value,wx_mor_sin_inxws$p.value,wx_mor_sin_outxws$p.value,
             wx_mor_sol_inxout$p.value)
comparison_cases <- data.frame(Method = method, city = city, 
                               comparison = comparison, W = W, P_value = p_value)
comparison_cases


