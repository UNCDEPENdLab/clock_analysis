
library(ggrepel)


power_grid <- tibble::tribble(
  ~scenario, ~std_beta, ~power,
  'Existing\ndata\n', 0.1, .63,
  'Existing\ndata\n', 0.125, .78,
  'Existing\ndata\n', 0.15, .91,
  'Existing\ndata\n', 0.2, .96,
  'Existing\ndata\n', 0.3, 1.0,
  'Existing\ndata\n', 0.4, 1.0,
  'Proposed\nfollow-up\n+ 95\nparticipants', 0.1, 0.69,
  'Proposed\nfollow-up\n+ 95\nparticipants', 0.125, 0.85,
  'Proposed\nfollow-up\n+ 95\nparticipants', 0.15, 0.97,
  'Proposed\nfollow-up\n+ 95\nparticipants', 0.2, .99,
  'Proposed\nfollow-up\n+ 95\nparticipants', 0.3, 1.0,
  'Proposed\nfollow-up\n+ 95\nparticipants', 0.4, 1.0
)

pdf("aim1_power_curve_n458.pdf", width=4.6, height=3)
ggplot(power_grid, aes(x=std_beta, y=power, lty=scenario)) +
  geom_point(size=3) + geom_line(size=1.2) + #stat_smooth(method="", formula=(y~exp(x)), se=F) +
ggtitle("Figure 7. Power to detect trait moderation \nof ideation-attempt pathway") +
  theme_minimal() + xlab("Standardized regression coefficient") + ylab("Statistical power") +
  annotate(x=0.386, y=0.805, label="Observed\ndis-\ninhibition\n effect in\nAllen et\n al.,\n2022\n->", geom="text", hjust=1, color="gray40") +
  geom_vline(xintercept=.39, linetype=5, color="gray40") +
  geom_hline(yintercept=0.8, color="gray50") +
  scale_color_brewer("", palette="Set2") +
  geom_vline(xintercept=.11, linetype=5, color="gray40") +
  geom_vline(xintercept=.29, linetype=5, color="gray40") + 
  annotate(x=0.16, y=0.69, label="25th\n%ile\n<--", geom="text", hjust=1, color="gray40") + 
  annotate(x=0.275, y=0.78, label="Gignac &\nSzodorai\n 2016\nmeta-analysis,\n50th\n%ile\n-->", geom="text", hjust=1, color="gray40") + 
  theme(legend.key.size = unit(.25, "in"))
  # geom_text_repel(data = power_grid %>% filter(std_beta == .125), aes(label = scenario), max.overlaps = 1) + theme(legend.position = "none")
  #theme(panel.grid = element_blank())
dev.off()

 