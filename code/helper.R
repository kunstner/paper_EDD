theme_alpha <- theme(axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.position =  "none",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     # axis.text.x = element_text(angle = 45, hjust = 1, size=12),
                     axis.text = element_text(size=12),
                     axis.text.y = element_text(size=12),
                     # axis.ticks.x = element_blank(),
                     strip.text = element_text(size=12, face = "bold"),
                     strip.background = element_rect(colour="white", fill="grey95",
                                                     linewidth=1.5, linetype="solid"))

colv3 <- c( 'ARCI' = "darkcyan",
            'Control' = "gold2",
            'IV' = 'darkorchid',
            'TRUE' = "gold2",
            'FALSE' = "steelblue",
            'Patient' = "steelblue",
            'Abdomen' = 'cornsilk3',
            'Arm left' = 'lightsteelblue', 'Arm_left' = 'lightsteelblue',
            'Arm right' = 'lightseagreen', 'Arm_right' = 'lightseagreen',
            'young' = 'seagreen1',
            'middle' = 'seashell1',
            'old' = 'skyblue1',
            'older' = 'skyblue1',
            'yes' = 'cornflowerblue',
            'no' = 'orange')

# Colors for taxonomic abundance plots
taxNumberPlot <- 19
cols <- c(pals::tableau20(taxNumberPlot+1), "magenta4", "grey35", "orange4", "grey75", "deepskyblue4")
