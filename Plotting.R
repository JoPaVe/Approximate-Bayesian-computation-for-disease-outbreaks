##Plots.R

```{r echo=FALSE}
library(ggplot2)

colnames(results_Table2[[2]][[1]]) <- c("qc1_post", "qh1_post", "qc2_post", "qh2_post")

ggplot(data = results_Table2[[2]][[1]]) +
  geom_point(aes(qh1_post, qc1_post, colour = "red")) +
  geom_point(aes(qh2_post, qc2_post, colour = "blue")) +
  labs(x = expression(q[c]), y = expression(q[h])) + 
  xlim(c(0, 1)) + ylim(c(0,1)) +
  ggtitle(label = "Figure 3(a)" , 
          subtitle = "ABC posterior distributions for parameters modeling 
  Supplementary Table 2")

colnames(results_Table3[[2]][[1]]) <- c("qc1_post", "qh1_post", "qc2_post", "qh2_post")

ggplot(data = results_Table3[[2]][[1]]) +
  geom_point(aes(qh1_post, qc1_post, colour = "red")) +
  geom_point(aes(qh2_post, qc2_post, colour = "blue")) +
  labs(x = expression(q[c]), y = expression(q[h])) + 
  xlim(c(0, 1)) + ylim(c(0,1)) +
  ggtitle(label = "Figure 3(b)" , 
          subtitle = "ABC posterior distributions for parameters modeling
  Supplementary Table 3")
```