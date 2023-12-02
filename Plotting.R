library(ggplot2)
library(gridExtra)


####### Replication Paper 3a, 3b

colnames(results_Table2[[2]][[1]]) <- c("qc1", "qh1", "qc2", "qh2")

Fig3a <- ggplot(data = results_Table2[[2]][[1]]) + 
  geom_point(aes(qh1, qc1), color = "red", shape = 1, size = 2.2) +
  geom_point(aes(qh2, qc2), color = "blue", shape = 1, size = 2.2) +
  labs(x = "qh", y = "qc") + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", linetype = "solid", colour = "black"),
    plot.title = element_text(face = "bold")
    ) +
  ggtitle("(a)")

colnames(results_Table3[[2]][[1]]) <- c("qc1", "qh1", "qc2", "qh2")

Fig3b <- ggplot(data = results_Table3[[2]][[1]]) + 
  geom_point(aes(qh1, qc1), color = "red", shape = 1, size = 2.2) +
  geom_point(aes(qh2, qc2), color = "blue", shape = 1, size = 2.2) +
  labs(x = "qh", y = "qc") + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", linetype = "solid", colour = "black"),
    plot.title = element_text(face = "bold")
  ) +
  ggtitle("(c)")

grid.arrange(Fig3a, Fig3b, nrow = 2, ncol = 1, heights = c(1,1))

