##################################################
# Read data
# Size and Power
##################################################
library(openxlsx)

for (ind in 1:2){
  if(ind==1){
    namei<-"_ind"
  } else{
    namei<-"_MAEE"
  }
  
  filename <- paste("binResults/final_bin", namei, ".xlsx", sep="")
  data <- read.xlsx(filename, colNames=TRUE)
  colnames(data)[16:20]<-c("MB.1","ROB.1","KC.1","MD.1","AVG.1")
  assign(paste("final_bin", namei, sep=""), data)
}



##################################################
# Figures of results for all variance estimators
# Size and Power
##################################################
library(reshape2)
library(ggplot2)
library(ggpubr)
library(scales)

df2 <- data.frame(sy1 = 0.036, sy2 = 0.064, py1 = -0.023, py2 = 0.023)


# Size

# MAEE
data_plot <- final_bin_MAEE[, -(9:15)]
data_plot <- data_plot[, -13]
data_plot <- data_plot[-(31:36), ]
colnames(data_plot)[9:12] <- c("MB", "ROB", "KC", "MD")
data_plot$model <- c(rep(rep(c("Model 1. NT-CTE", "Model 2. LT-CTE", "Model 3. CT-CTE", "Model 4. LT-TI", "Model 5. CT-TI"), each=6)))
data_plot$scenario <- seq(1, 30, 1)

data_plot <- melt(data_plot, id.vars = c("b", "alpha0", "alpha1", "alpha2", "N", "I", "K", "t", "model", "scenario"),
                  measure.vars = c("MB", "ROB", "KC", "MD"),
                  variable.name = "var", value.name = "size")
min(data_plot$size) # 0.033
max(data_plot$size) # 0.093

P_MAEE <- ggplot(data=data_plot, aes(scenario, size, colour = factor(var), shape = factor(var))) +
  geom_point(size = 3) + 
  geom_hline(data = df2, aes(yintercept = sy1, linetype = "Acceptable bound"), color = "darkgrey") +
  geom_hline(data = df2, aes(yintercept = sy2, linetype = "Acceptable bound"), color = "darkgrey") +
  scale_linetype_manual(values = 2, breaks = "Acceptable bound") +
  scale_shape_manual(values = c(8, 18, 16, 15)) +
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "darkgrey"),
                     breaks = c("MB", "ROB", "KC", "MD", "Acceptable bound")) +
  labs(title = "(a) GEE/MAEE",
       subtitle = "[                         No interaction effect                          ] [        With interaction effect(s)        ]",
       y = "Type I Error Rate", x = "Scenario",
       colour = expression(paste("Variance estimator")),
       linetype = " ") +
  facet_grid(. ~ model, scale="free_x", space="free_x") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) +
  scale_y_continuous(limits = c(0.025, 0.225), breaks = seq(0.025, 0.225, 0.025)) +
  scale_x_continuous(breaks = seq(1, 30, 2)) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = c(8, 18, 16, 15)))) +
  guides(shape = "none")

# ind
data_plot <- final_bin_ind[, -(9:15)]
data_plot <- data_plot[, -13]
data_plot <- data_plot[-(31:36), ]
colnames(data_plot)[9:12] <- c("MB", "ROB", "KC", "MD")
data_plot$model <- c(rep(rep(c("Model 1. NT-CTE", "Model 2. LT-CTE", "Model 3. CT-CTE", "Model 4. LT-TI", "Model 5. CT-TI"), each=6)))
data_plot$scenario <- seq(1, 30, 1)

data_plot <- melt(data_plot, id.vars = c("b", "alpha0", "alpha1", "alpha2", "N", "I", "K", "t", "model", "scenario"),
                  measure.vars = c("MB", "ROB", "KC", "MD"),
                  variable.name = "var", value.name = "size")
min(data_plot$size) # 0.034
max(data_plot$size) # 0.268

p_IND <- ggplot(data=data_plot, aes(scenario, size, colour = factor(var), shape = factor(var))) +
  geom_point(size = 3) + 
  geom_hline(data = df2, aes(yintercept = sy1, linetype = "Acceptable bound"), color = "darkgrey") +
  geom_hline(data = df2, aes(yintercept = sy2, linetype = "Acceptable bound"), color = "darkgrey") +
  scale_linetype_manual(values = 2, breaks = "Acceptable bound") +
  scale_shape_manual(values = c(8, 18, 16, 15)) +
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "darkgrey"),
                     breaks = c("MB", "ROB", "KC", "MD", "Acceptable bound")) +
  labs(title = "(b) GEE with working independence",
       subtitle = "[                         No interaction effect                          ] [        With interaction effect(s)        ]",
       y = "Type I Error Rate", x = "Scenario",
       colour = expression(paste("Variance estimator")),
       linetype = " ") +
  facet_grid(. ~ model, scale="free_x", space="free_x") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) +
  scale_y_continuous(limits = c(0.025, 0.225), breaks = seq(0.025, 0.225, 0.025)) +
  scale_x_continuous(breaks = seq(1, 30, 2)) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = c(8, 18, 16, 15)))) +
  guides(shape = "none")


p_size <- ggarrange(P_MAEE, p_IND, ncol = 2, nrow = 1)

ggsave("binResults/binSize.pdf", p_size, height = 6.5, width = 13.5)



# Power

# MAEE
data_plot <- final_bin_MAEE[, -(16:21)]
data_plot <- data_plot[, -14]
data_plot <- data_plot[-(31:36), ]
data_plot$model <- c(rep(rep(c("Model 1. NT-CTE", "Model 2. LT-CTE", "Model 3. CT-CTE", "Model 4. LT-TI", "Model 5. CT-TI"), each=6)))
data_plot$scenario <- seq(1, 30, 1)

data_plot <- melt(data_plot, id.vars = c("b", "alpha0", "alpha1", "alpha2", "N", "I", "K", "t", "model", "scenario", "a_power"),
                  measure.vars = c("MB", "ROB", "KC", "MD"),
                  variable.name = "var", value.name = "power")
data_plot$pd <- data_plot$power - data_plot$a_power
min(data_plot$pd) # -0.02652284
max(data_plot$pd) # 0.04364018

P_MAEE <- ggplot(data=data_plot, aes(scenario, pd, colour = factor(var), shape = factor(var))) +
  geom_point(size = 3) + 
  geom_hline(data = df2, aes(yintercept = py1, linetype = "Acceptable bound"), color = "darkgrey") +
  geom_hline(data = df2, aes(yintercept = py2, linetype = "Acceptable bound"), color = "darkgrey") +
  scale_linetype_manual(values = 2, breaks = "Acceptable bound") +
  scale_shape_manual(values = c(8, 18, 16, 15)) +
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "darkgrey"),
                     breaks = c("MB", "ROB", "KC", "MD", "Acceptable bound")) +
  labs(title = "(a) GEE/MAEE",
       subtitle = "[                         No interaction effect                         ] [        With interaction effect(s)       ]",
       y = "Power Difference (empirical - predicted)", x = "Scenario",
       colour = expression(paste("Variance estimator")),
       linetype = " ") +
  facet_grid(. ~ model, scale="free_x", space="free_x") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) +
  scale_y_continuous(limits = c(-0.029, 0.125), breaks = seq(-0.025, 0.125, 0.025), labels = comma) +
  scale_x_continuous(breaks = seq(1, 30, 2)) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = c(8, 18, 16, 15)))) +
  guides(shape = "none")

# ind
data_plot <- final_bin_ind[, -(16:21)]
data_plot <- data_plot[, -14]
data_plot <- data_plot[-(31:36), ]
data_plot$model <- c(rep(rep(c("Model 1. NT-CTE", "Model 2. LT-CTE", "Model 3. CT-CTE", "Model 4. LT-TI", "Model 5. CT-TI"), each=6)))
data_plot$scenario <- seq(1, 30, 1)

data_plot <- melt(data_plot, id.vars = c("b", "alpha0", "alpha1", "alpha2", "N", "I", "K", "t", "model", "scenario", "a_power"),
                  measure.vars = c("MB", "ROB", "KC", "MD"),
                  variable.name = "var", value.name = "power")
data_plot$pd <- data_plot$power - data_plot$a_power
min(data_plot$pd) # -0.02852284
max(data_plot$pd) # 0.1221302

p_IND <- ggplot(data=data_plot, aes(scenario, pd, colour = factor(var), shape = factor(var))) +
  geom_point(size = 3) + 
  geom_hline(data = df2, aes(yintercept = py1, linetype = "Acceptable bound"), color = "darkgrey") +
  geom_hline(data = df2, aes(yintercept = py2, linetype = "Acceptable bound"), color = "darkgrey") +
  scale_linetype_manual(values = 2, breaks = "Acceptable bound") +
  scale_shape_manual(values = c(8, 18, 16, 15)) +
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "darkgrey"),
                     breaks = c("MB", "ROB", "KC", "MD", "Acceptable bound")) +
  labs(title = "(b) GEE with working independence",
       subtitle = "[                         No interaction effect                         ] [        With interaction effect(s)       ]",
       y = "Power Difference (empirical - predicted)", x = "Scenario",
       colour = expression(paste("Variance estimator")),
       linetype = " ") +
  facet_grid(. ~ model, scale="free_x", space="free_x") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) +
  scale_y_continuous(limits = c(-0.029, 0.125), breaks = seq(-0.025, 0.125, 0.025), labels = comma) +
  scale_x_continuous(breaks = seq(1, 30, 2)) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = c(8, 18, 16, 15)))) +
  guides(shape = "none")


p_power <- ggarrange(P_MAEE, p_IND, ncol = 2, nrow = 1)

ggsave("binResults/binPower.pdf", p_power, height = 6.5, width = 13.5)


