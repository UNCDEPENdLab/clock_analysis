sum_data <- group_data %>% select(source_est, node, Time) %>% group_by(node, Time) %>% 
  summarise(mean_signal = mean(source_est),
            sd_signal = sd(source_est))
sum_data <- sum_data %>% filter(Time < 3.5) %>% mutate(
  side = str_extract(node, "L|R"), 
  lobe = str_remove(node, "_L|_R")
) 
ggplot(sum_data, aes(Time, mean_signal)) + geom_line(aes(color = side)) + 
  geom_ribbon(aes(ymin = mean_signal - sd_signal/100, ymax = mean_signal + sd_signal/100, fill = side)) + facet_wrap(~lobe)
ggplot(sum_data, aes(Time, mean_signal)) + geom_line(aes(color = side)) + facet_wrap(~lobe)
