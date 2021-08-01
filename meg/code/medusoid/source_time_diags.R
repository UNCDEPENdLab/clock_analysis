sum_data <- group_data %>% select(source_est, node, Time) %>% group_by(node, Time) %>% 
  summarise(mean_signal = mean(source_est),
            sd_signal = sd(source_est),
            n = n())
sum_data <- sum_data %>% filter(Time < 3.5) %>% mutate(
  side = str_extract(node, "L|R"), 
  lobe = str_remove(node, "_L|_R"),
  se = sd_signal/n^0.5
) 
ggplot(sum_data, aes(Time, mean_signal, color = side)) + geom_line() + 
  geom_ribbon(aes(ymin = mean_signal - se, ymax = mean_signal + se), fill = "grey70", alpha = .8) + facet_wrap(~lobe)

ggplot(sum_data, aes(Time, mean_signal)) + geom_line(aes(color = side)) + facet_wrap(~lobe)
