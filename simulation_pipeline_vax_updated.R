library(tidyverse)
FILEPATH = "C:/Users/Abraham/miniconda3/envs/snowflakes/Scripts/Global_SEIR/Real Sims"
TN = 47
L = 20
SHORTCUT = "four_unrestrained_600"
clean_and_plot_areas(600, 0.7, 600, 4)

sim_info_function <- function(i1, i2, i3, i4){
  #sim_info = paste("global_simulation_", i1, "_days_", i2,"_alpha_", i3, "_vax_", i4, "_areas", sep = "")
  sim_info = paste("donor_four_countries_global_simulation_", i2,"_alpha", sep = "")
  return(sim_info)
}


clean_csv <- function(T_days, alpha, vax_days, num_areas){
  sim_info = sim_info_function(T_days, alpha, vax_days, num_areas)
  sim_info= SHORTCUT
  print(paste(FILEPATH, "/simulation_data/", sim_info, ".csv", sep = ""))
  sim_data <- read.csv(paste(FILEPATH, "/simulation_data/", sim_info, ".csv", sep = ""), header = F)
  vax = sim_data[which(substr(sim_data[,1], 1, 3)=="vax"),]
  vax_data_clean <- vax %>%
    rename(cases = 2) %>%
    separate(col = 1, into = c("area", "day"), sep = ",") %>%
    mutate(day = as.integer(substr(day, start = 1, stop = nchar(day)-1))) %>%
    mutate(area = substr(area, start = nchar(area)-4, stop = nchar(area))) %>% #assumes fewer than 10 areas
    mutate(state = "V") %>%
    group_by(area) %>%
    mutate(cases = cumsum(cases))
  sim_data = sim_data[which(substr(sim_data[,1], 1, 5)=="state"),]
  sim_data_clean <- sim_data %>%
    rename(cases = 2) %>%
    separate(col = 1, into = c("area", "state", "day"), sep = ",") %>%
    mutate(day = as.integer(substr(day, start = 1, stop = nchar(day)-1))) %>%
    mutate(area = substr(area, start = nchar(area)-4, stop = nchar(area))) #assumes fewer than 10 areas
  
  
  df <- rbind(sim_data_clean, vax_data_clean)
  return(df)
  #return(sim_data_clean)
}
plot_simulation <- function(df, c1, c2, c3, c4){
  vline_df = data.frame(intercept = c(TN+L, TN+L, TN, TN+L), area = c("area1", "area2", "area3", "area4"))
  p <-df %>%
    #filter(day <= c1) %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df, mapping = aes(xintercept = intercept), lty = "dotted")+
    #ggtitle(paste(c1, " day time horizon in ", c4, " areas with alpha = ",
    #              c2,"\nAnd ",
    #              c3, " vaccination days.", sep = ""))+
    facet_wrap(~area)+
    theme_light()
  return(p)
}
plot_simulation_low_alpha <- function(df, c1, c2, c3, c4){
  vline_df = data.frame(intercept = c(TN+L, TN+L, TN, TN+L), area = c("area1", "area2", "area3", "area4"))
  p <- df %>%
    filter(!(state %in% c("S", "SV", "R", "V"))) %>%
    #filter(day <= c1) %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df, mapping = aes(xintercept = intercept), lty = "dotted")+
    #ggtitle(paste(c1, " day time horizon in ", c4, " areas with alpha = ",
    #              c2,"\nAnd ",
    #              c3, " vaccination days, \nThere are no S or SV states.", sep = ""))+
    facet_wrap(~area)+
    theme_light()
  
  return(p)
}
plot_simulation_only_vax <- function(df, c1, c2, c3, c4){
  print(df[which(df$cases < 0),])
  vline_df = data.frame(intercept = c(TN+L, TN+L, TN, TN+L), area = c("area1", "area2", "area3", "area4"))
  p <- df %>%
    filter((state %in% c("EV", "SV", "IV"))) %>%
    #filter(day <= c1) %>%
    ggplot()+
    geom_line(aes(x=day, y = cases, color = state), lwd = 1.5)+
    geom_vline(data = vline_df, mapping = aes(xintercept = intercept), lty = "dotted")+
    #ggtitle(paste(c1, " day time horizon in ", c4, " areas with alpha = ",
    #              c2,"\nAnd ",
    #              c3, " vaccination days, \nOnly vax states.", sep = ""))+
    facet_wrap(~area)+
    theme_light()
  
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
  #df$state = factor(df$state, levels = c("V", "S", "SV", "I", "IV", "E", "EV", "H", "R", "D"))
  vline_df = data.frame(intercept = c(TN+L, TN+L, TN, TN+L), area = c("area1", "area2", "area3", "area4"))
  p <-df %>%
    #filter(day <= 110) %>%
    ggplot()+
    geom_area(aes(x=day, y = cases, fill = state), data = subset(df, state != "V"))+
    geom_line(aes(x=day, y = cases), data = subset(df, state == "V"))+
    geom_vline(data = vline_df, mapping = aes(xintercept = intercept), lty = "dotted")+
    #geom_vline(xintercept = TN, lty = "dotted")+
    #geom_vline(xintercept = TN+L, lty = "dotted")+
    #geom_hline(yintercept = 48960, lty = "dotted")+
    #geom_hline(yintercept = 86130, lty = "dotted")+
    #geom_hline(yintercept = 76410, lty = "dotted")+
    #geom_hline(yintercept = 93830, lty = "dotted")+
    #ggtitle(paste(c1, " day time horizon in ", c4, " areas with alpha = ",
    #              c2,"\nAnd ",
    #              c3, " vaccination days.", sep = ""))+
    facet_wrap(~area)+
    theme_light()
  return(p)
}
create_folder <- function(sim_inform){
  output_dir = paste(FILEPATH, "/plots/", sim_inform, "_editing", sep = "")
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
  sub_dir = paste(FILEPATH, "/plots/", sim_information, "_editing", sep = "")
  ggsave(paste("default", int2, ".png", sep = ""), plot = p1, path = sub_dir)
  ggsave(paste("low_alpha", int2, ".png", sep = ""), plot = p2, path = sub_dir)
  ggsave(paste("only_vax", int2, ".png", sep = ""), plot = p3, path = sub_dir)
  ggsave(paste("stacked", int2, ".png", sep = ""), plot = p4, path = sub_dir)
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
clean_and_plot_areas(180, 0.7, 180, 4)


