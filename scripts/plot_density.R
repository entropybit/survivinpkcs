library(tidyverse)
library(magrittr)
library(cowplot)
library(ggsci)
setwd("~/DNAPKcs")

data_a <- read_csv("density_data_hetero_a.csv")
data_d <- read_csv("density_data_hetero_d.csv")
data_h_a <- read_csv("density_data_frb_fat_a.csv")
data_h_d <- read_csv("density_data_frb_fat_d.csv")
data_k_a <- read_csv("density_data_kinase_a.csv")
data_k_d <- read_csv("density_data_kinase_d.csv")
data_a$X1 %<>% multiply_by(100 / 10000) 
data_d$X1 %<>% multiply_by(100 / 10000) 
data_h_a$X1 %<>% multiply_by(100 / 10000) 
data_h_d$X1 %<>% multiply_by(100 / 10000) 
data_k_a$X1 %<>% multiply_by(100 / 10000) 
data_k_d$X1 %<>% multiply_by(100 / 10000) 
data_a %<>% mutate("Kinase" = "A", "Type" = "Tetramer")
data_d %<>% mutate("Kinase" = "D", "Type" = "Tetramer")
data_h_a %<>% mutate("Kinase" = "A", "Type" = "Head dimer")
data_h_d %<>% mutate("Kinase" = "D", "Type" = "Head dimer")
data_k_a %<>% mutate("Kinase" = "A", "Type" = "Reference")
data_k_d %<>% mutate("Kinase" = "D", "Type" = "Reference")
data <-bind_rows(data_a,data_d,data_h_a,data_h_d,data_k_a,data_k_d)
rmsd_kin_left_hetero <- read_csv("heterodimer/rmsd_kinase_left.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "A") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Tetramer")
rmsd_kin_right_hetero <- read_csv("heterodimer/rmsd_kinase_right.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "D") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Tetramer")
rmsd_sur_left_hetero <- read_csv("heterodimer/rmsd_survivin_left.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "B") %>%
  mutate("Protein" = "Survivin") %>%
  mutate("Structure" = "Tetramer")
rmsd_sur_right_hetero <- read_csv("heterodimer/rmsd_survivin_right.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "C") %>%
  mutate("Protein" = "Survivin") %>%
  mutate("Structure" = "Tetramer")
radius_kin_left_hetero <- read_csv("heterodimer/radius_kinase_left.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>%
  mutate("Chain" = "A") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Tetramer")
radius_kin_right_hetero <- read_csv("heterodimer/radius_kinase_right.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>% 
  mutate("Chain" = "D") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Tetramer")
radius_sur_left_hetero <- read_csv("heterodimer/radius_survivin_left.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>% 
  mutate("Chain" = "B") %>%
  mutate("Protein" = "Survivin") %>%
  mutate("Structure" = "Tetramer")
radius_sur_right_hetero <- read_csv("heterodimer/radius_survivin_right.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>% 
  mutate("Chain" = "C") %>%
  mutate("Protein" = "Survivin") %>%
  mutate("Structure" = "Tetramer")
rmsd_kin_left_no_surv <- read_csv("no_surv/rmsd_kinase_left.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "A") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Head dimer")
rmsd_kin_right_no_surv <- read_csv("no_surv/rmsd_kinase_right.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "D") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Head dimer")
radius_kin_left_no_surv <- read_csv("no_surv/radius_kinase_left.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>%
  mutate("Chain" = "A") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Head dimer")
radius_kin_right_no_surv <- read_csv("no_surv/radius_kinase_right.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>% 
  mutate("Chain" = "D") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Head dimer")
rmsd_kin_left_nullmodel <- read_csv("nullmodel/nullmodel_rmsd_kinase_left.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "A") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Reference")
rmsd_kin_right_nullmodel <- read_csv("nullmodel/nullmodel_rmsd_kinase_right.csv") %>% 
  mutate("Type" = "RMSD") %>% 
  mutate("Chain" = "D") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Reference")
radius_kin_left_nullmodel <- read_csv("nullmodel/nullmodel_radius_kinase_left.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>%
  mutate("Chain" = "A") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Reference")
radius_kin_right_nullmodel <- read_csv("nullmodel/nullmodel_radius_kinase_right.csv") %>% 
  mutate("Type" = "Radius of Gyration") %>% 
  mutate("Chain" = "D") %>%
  mutate("Protein" = "Kinase") %>%
  mutate("Structure" = "Reference")

rmsds <- rbind.data.frame(rmsd_kin_left_hetero,
                          rmsd_kin_right_hetero,
                          rmsd_sur_left_hetero,
                          rmsd_sur_right_hetero,
                          rmsd_kin_left_no_surv,
                          rmsd_kin_right_no_surv,
                          rmsd_kin_left_nullmodel,
                          rmsd_kin_right_nullmodel)

radius <- rbind.data.frame(radius_kin_left_hetero,
                           radius_kin_right_hetero,
                           radius_sur_left_hetero,
                           radius_sur_right_hetero,
                           radius_kin_left_no_surv,
                           radius_kin_right_no_surv,
                           radius_kin_left_nullmodel,
                           radius_kin_right_nullmodel)

rmsds$X1 %<>% multiply_by(100 / 10000)
radius$X1 %<>% multiply_by(100 / 10000)


plt <- data %>%
  gather("Radius","Density",-X1,-Kinase,-Type)%>%
  filter(Radius %in% c("13")) %>%
  ggplot(aes(x=X1,y=Density,color=Kinase))+
  theme_classic() + 
  theme(legend.position = "bottom",legend.title = element_blank()) + 
  geom_line(alpha=0.3,size=0.2)+
  geom_smooth(size=1)+
  labs(x = expression("time "/" ns"),y = expression("Particle Density "/" N /"~ring(A)^3))+
  scale_color_manual(values=c("#E18727FF","#20854EFF"))+
  facet_grid(Type~.)

plot_rmsds <- rmsds %>%
  filter(Protein == "Kinase") %>%
  ggplot(aes(x=X1,y=`0`,color=Structure))+
  theme_classic() + 
  geom_line(alpha=0.3,size=0.3)+
  geom_smooth(size = 1)+
  scale_color_nejm()+ 
  theme(legend.position = "none") + 
  facet_grid(Chain~.)+
  labs(x=expression("time "/" ns"),y=expression("RMSD "/~ring(A)))

plot_rmsds_sur <- rmsds %>%
  filter(Protein == "Survivin") %>%
  ggplot(aes(x=X1,y=`0`,color=Chain))+
  theme_classic() + 
  geom_line(alpha=0.3,size=0.3)+
  geom_smooth(size = 1)+
  theme(legend.position = "none") + 
  scale_color_nejm()+
  scale_y_continuous(limits = c(0,10)) + 
  labs(x=expression("time "/" ns"),y=expression("RMSD "/~ring(A)))

plot_radius_kinase <- radius %>%
  filter(Protein == "Kinase" ) %>% 
  ggplot(aes(x=X1,y=`0`,color=Structure))+
  theme_classic() + 
  geom_line(alpha=0.3,size=0.3)+
  geom_smooth(size=1)+
  scale_color_nejm()+
  theme(legend.position = "bottom",legend.title = element_blank()) + 
  facet_grid(Chain~.)+
  labs(x=expression("time "/" ns"),y=expression("Radius of Gyration "/~ring(A)))

plot_radius_survivin <- radius %>%
  filter(Protein == "Survivin" ) %>% 
  ggplot(aes(x=X1,y=`0`,color=Chain))+
  theme_classic() + 
  geom_line(alpha=0.3,size=0.3)+
  geom_smooth(size=1)+
  theme(legend.position = "bottom",legend.title = element_blank()) + 
  scale_color_nejm()+
  labs(x=expression("time "/" ns"),y=expression("Radius of Gyration "/~ring(A)))

text_size <- theme(axis.title = element_text(size = 16),
      axis.text = element_text(size=12))

fig_5 <- ggdraw() +
  draw_plot(plot_rmsds + text_size,0,0.83,1,0.17) +
  draw_plot(plot_radius_kinase + text_size,0,0.62,1,0.21) +
  draw_plot(plot_rmsds_sur + text_size,0,0.485,1,0.146) + 
  draw_plot(plot_radius_survivin + text_size,0,0.28,1,0.211) + 
  draw_plot(plt + text_size,0,0,1,0.28)

save_plot(filename = "figure5.png",plot = fig_5,base_height = 12,base_width = 8)
