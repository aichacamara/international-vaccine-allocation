library(tidyverse)
library(gridExtra)
library(ggpubr)

FILEPATH = "C:/Users/Abraham/miniconda3/envs/snowflakes/Scripts/Global_SEIR/Real Sims"
SHORTCUT = "global_0.75_max_3_iterations"
TN = 53
clean_and_plot_areas(180, 0.6, 180, 4)

sim_info_function <- function(i1, i2, i3, i4){
  #sim_info = paste("global_simulation_", i1, "_days_", i2,"_alpha_", i3, "_vax_", i4, "_areas", sep = "")
  sim_info = paste("donor_all_countries_global_simulation_", i2,"_alpha", sep = "")
  return(sim_info)
}

clean_csv <- function(T_days, alpha, vax_days, num_areas){
  sim_info = sim_info_function(T_days, alpha, vax_days, num_areas)
  sim_info = SHORTCUT
  print(paste(FILEPATH, "/simulation_data/", sim_info, ".csv", sep = ""))
  sim_data <- read.csv(paste(FILEPATH, "/simulation_data/", sim_info, ".csv", sep = ""), header = F)
  vax = sim_data[which(substr(sim_data[,1], 1, 3)=="vax"),]
  vax_data_clean <- vax %>%
    rename(cases = 2) %>%
    separate(col = 1, into = c("area", "day"), sep = ",") %>%
    mutate(day = as.integer(substr(day, start = 1, stop = nchar(day)-1))) %>%
    mutate(area = substr(area, start = nchar(area)-2, stop = nchar(area))) %>% #assumes fewer than 10 areas
    mutate(state = "V") %>%
    filter(area %in% c("USA", "PER", "IND", "ZMB")) %>%
    group_by(area) %>%
    mutate(cases = cumsum(cases))
  sim_data = sim_data[which(substr(sim_data[,1], 1, 5)=="state"),]
  sim_data_clean <- sim_data %>%
    rename(cases = 2) %>%
    separate(col = 1, into = c("area", "state", "day"), sep = ",") %>%
    mutate(day = as.integer(substr(day, start = 1, stop = nchar(day)-1))) %>%
    mutate(area = substr(area, start = nchar(area)-2, stop = nchar(area))) %>% #assumes fewer than 10 areas
    filter(area %in% c("USA", "PER", "IND", "ZMB"))
  
  
  df <- rbind(sim_data_clean, vax_data_clean)
  return(df)
  #return(sim_data_clean)
}
plot_simulation <- function(df, c1, c2, c3, c4){
  vline_df = data.frame(intercept = c(TN+20, TN+20, TN, TN+20), area = c("USA", "PER", "IND", "ZMB"))
  p1 <- df %>%
    #filter(day <= c1) %>%
    filter(area == "USA") %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[1,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("United States of America (donor)")+
    theme_light()+ theme(legend.position = "none")
  p2 <- df %>%
    filter(area == "PER") %>%
    
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[2,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("Peru")+
    theme_light()+ theme(legend.position = "none")
  p3 <- df %>%
    filter(area == "IND") %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[3,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("India")+
    theme_light()+ theme(legend.position = "none")
  p4 <- df %>%
    filter(area == "ZMB") %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[4,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("Zambia")+
    theme_light()+ theme(legend.position = "none")
  p <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
  return(p)
}
plot_simulation_low_alpha <- function(df, c1, c2, c3, c4){
  vline_df = data.frame(intercept = c(TN+20, TN+20, TN, TN+20), area = c("USA", "PER", "IND", "ZMB"))
  p1 <- df %>%
    filter(!(state %in% c("S", "SV", "R", "V"))) %>%
    filter(area == "USA") %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[1,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("United States of America (donor)")+
    theme_light()+ theme(legend.position = "none")
  p2 <- df %>%
    filter(!(state %in% c("S", "SV", "R", "V"))) %>%
    filter(area == "PER") %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[2,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("Peru")+
    theme_light()+ theme(legend.position = "none")
  p3 <- df %>%
    filter(!(state %in% c("S", "SV", "R", "V"))) %>%
    filter(area == "IND") %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[3,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("India")+
    theme_light()+ theme(legend.position = "none")
  p4 <- df %>%
    filter(!(state %in% c("S", "SV", "R", "V"))) %>%
    filter(area == "ZMB") %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df[4,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("Zambia")+
    theme_light()+ theme(legend.position = "none")
  p <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
  
  return(p)
}

plot_simulation_stacked <- function(df, c1, c2, c3, c4){
  #df$state[which(df$state %in% c("S", "SV"))] <- "S,SV"
  df$state[which(df$state %in% c("I", "IV",  "E", "EV"))] <- "I,E"
  df$state[which(df$state %in% c("H"))] <- "H"
  df <- df %>%
    group_by(area, state, day)%>%
    summarise(area, state, day, cases = sum(cases)) %>%
    ungroup() %>%
    unique()
  df$state <- factor(df$state, levels = c("V", "S", "SV", "I,E", "H", "R", "D"))
  vline_df = data.frame(intercept = c(TN+20, TN+20, TN, TN+20), area = c("USA", "PER", "IND", "ZMB"))
  usa_df <- df %>%
    filter(area == "USA")
  per_df <- df %>%
    filter(area == "PER")
  ind_df <- df %>%
    filter(area == "IND")
  zmb_df <- df %>%
    filter(area == "ZMB")
  p1 <-usa_df %>%
    ggplot()+
    geom_area(aes(x=day, y = cases, fill = state), data = subset(usa_df, state != "V"))+
    geom_line(aes(x=day, y = cases), data = subset(usa_df, state == "V"))+
    geom_vline(data = vline_df[1,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("United States of America (donor)")+
    theme_light()+ theme(legend.position = "none")
  p2 <-per_df %>%
    ggplot()+
    geom_area(aes(x=day, y = cases, fill = state), data = subset(per_df, state != "V"))+
    geom_line(aes(x=day, y = cases), data = subset(per_df, state == "V"))+
    geom_vline(data = vline_df[2,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("Peru")+
    theme_light() + theme(legend.position = "none")
  p3 <-ind_df %>%
    ggplot()+
    geom_area(aes(x=day, y = cases, fill = state), data = subset(ind_df, state != "V"))+
    geom_line(aes(x=day, y = cases), data = subset(ind_df, state == "V"))+
    geom_vline(data = vline_df[3,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("India")+
    theme_light()+ theme(legend.position = "none")
  p4 <-zmb_df %>%
    ggplot()+
    geom_area(aes(x=day, y = cases, fill = state), data = subset(zmb_df, state != "V"))+
    geom_line(aes(x=day, y = cases), data = subset(zmb_df, state == "V"))+
    geom_vline(data = vline_df[4,], mapping = aes(xintercept = intercept), lty = "dotted")+
    ggtitle("Zambia")+
    theme_light()+ theme(legend.position = "none")
  p <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
  
  return(p)
}
create_folder <- function(sim_inform){
  output_dir = paste(FILEPATH, "/plots/", sim_inform, "global_2", sep = "")
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print(paste("Dir /plots/", sim_inform, " already exists!", sep = ""))
  }
}



clean_and_plot_areas <- function(int1, int2, int3, int4){
  #int1 = T, int2 = alpha, int3 = vax_days
  #int1 = int1 - 10
  #int3 = int3 - 10 #Because we sim on 190 days but only show up to 180
  sim_information = sim_info_function(int1, int2, int3, int4)
  sim_information = SHORTCUT
  create_folder(sim_information)
  csv = clean_csv(int1, int2, int3, int4)
  p1 <- plot_simulation(csv, int1, int2, int3, int4)
  p2 <- plot_simulation_low_alpha(csv, int1, int2, int3, int4)
  p3 <- plot_simulation_only_vax(csv, int1, int2, int3, int4)
  p4 <- plot_simulation_stacked(csv, int1, int2, int3, int4)
  sub_dir = paste(FILEPATH, "/plots/", sim_information, "global_2", sep = "")
  ggsave(paste("global_default", int2, "_2.png", sep = ""), plot = p1, path = sub_dir)
  ggsave(paste("global_low_alpha", int2, "_2.png", sep = ""), plot = p2, path = sub_dir)
  ggsave(paste("global_stacked", int2, "_2.png", sep = ""), plot = p4, path = sub_dir)
}

#clean_and_plot(90, .75, 30)
#clean_and_plot(90, 1, 30)
#clean_and_plot(90, .4, 30)
#clean_and_plot(90, .55, 30)
#clean_and_plot(90, .65, 30)
#clean_and_plot_areas(90, 1, 90, 4)
#clean_and_plot_areas(300, .35, 1, 2)
#clean_and_plot_areas(90, 1, 0, 4)
#clean_and_plot_areas(90, 1, 20, 4)
#clean_and_plot_areas(TN0, .TN, 300, 4)
#clean_and_plot_areas(TN0, 3.5, TN0, 4)
#clean_and_plot_areas(120, 2, 30, 4)
#clean_and_plot_areas(120, 2, 120, 2)
#clean_and_plot_areas(55, 1.815, 0, 4)
#clean_and_plot_areas(120, 2, 100, 2)
#clean_and_plot_areas(120, 1.1, 100, 2)
#clean_and_plot_areas(120, .86643, 1, 2)
#clean_and_plot_areas(120, 1.1, 1, 2)
#clean_and_plot_areas(120, 1.03638, 1, 2)
#clean_and_plot_areas(120, 1.04, 1, 2)
#clean_and_plot_areas(120, 1.302, 120, 2)
#clean_and_plot_areas(120, 1.302, 120, 2)
clean_and_plot_areas(180, 0.6, 180, 4)


