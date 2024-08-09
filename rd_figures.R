
#f3 heatmap

f3_realdata_plot=ggplot(data=f4sim_real_baf3_plot) +geom_tile(data=f4sim_real_baf3_plot,aes(x=A,y=reorder(B,-f3),fill=f3))+
theme_bw()+
  scale_fill_gradient2(low=lighten("#5FFBF1",0.5),mid =lighten("#86A8E7",0.3) ,high =lighten("#f47727",0.1) ,na.value = "white",midpoint = 0.1748)+
  scale_y_discrete(drop = TRUE, expand = c(0, 0)) + ylab("")+ xlab("\nEast Anatolia")+
  facet_grid(pop2~., scales = "free", space = "free_y", switch = "y") +
  theme(strip.placement = "outside",
        panel.spacing = unit(0, "in"),
        strip.text = element_text(size=14,face="bold"),
        strip.background.y = element_rect(fill = "white", color = "gray75",size=3),
        axis.title.x = element_text(face="bold",color="gray35",size = 14),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=10),
        plot.title = element_text(face = "bold",size = 15));f3_realdata_plot

#f4

f4_realdata_plot=ggplot(data=f4sim_real_baf4_red_red,) +
  stat_interval(aes(x=Y,y=-f4,col=ifelse(abs(Zscore)>3,"Z>3","Z<3")),width = 1.5)+
  geom_point(aes(x=Y,y=-f4),size=1,alpha=0.5)+
  theme_bw()+ stat_summary(aes(x=Y,y=-f4),
    geom = "point", fun = median, size = 3, shape = 8, stroke = 1.6, alpha=0.5
  ) + scale_color_manual(values = c("#e7c482","#f47727"))+
  geom_hline(yintercept = 0,linewidth=0.6,col="gray40")+
   ylab("f4(YRI, EAN; EAN, GRC)\n")+ xlab("\nEast Anatolia")+
  guides(col=guide_legend(title="Z Score"))+
  theme(
        axis.title.y = element_text(size = 13,face="bold"),
        axis.text.x = element_text(size=10),
        plot.title = element_text(face = "bold",size = 15),
        axis.title.x = element_text(size = 12,face = "bold",colour ="gray35")); f4_realdata_plot

#map

style_1 <- c(feature = 'road', element = 'all', visibility = 'off')
style_2 <- c('&style=', feature = 'administrative', element = 'all', visibility = 'off')
style_3 <- c('&style=', feature = 'administrative.province', element = 'all', visibility = 'off')
style_4 <- c('&style=', feature = 'administrative.administrative.neighborhood', element = 'all', visibility = 'on')
style_5 <- c('&style=', feature = 'administrative.administrative.land_parcel', element = 'all', visibility = 'off')
style_6 <- c('&style=', feature = 'administrative.locality', element = 'all', visibility = 'off')
style_7 <- c('&style=', feature = 'poi', element = 'all', visibility = 'off')
style_8 <- c('&style=', feature = 'transit', element = 'all', visibility = 'off')
style_9 <- c('&style=', feature = 'landscape.man_made', element = 'all', visibility = 'off')
style_10 <- c('&style=', feature = 'landscape.natural', element = 'labels', visibility = 'off')
style_11 <- c('&style=', feature = 'landscape.natural.landcover', element = 'all', visibility = 'on', color = '0xb9a3d7')
style_12 <- c('&style=', feature = 'water', element ="labels" , visibility = 'off')
style_13 <- c('&style=', feature = 'water', element ="geometry" , visibility = 'on', color="0xeee3fd")

install.packages("ggmap")
library(ggmap)


anatgrc_f4map_plot=ggmap(get_googlemap(center=c(32.5400,38.3524),
                    zoom = 5, scale = 2,
                    maptype ='terrain',
                    style = c(style_1,
                              style_2,style_3,style_4,style_5,style_6,
                              style_7,style_8,style_9,style_11,style_12,style_13,style_10)))+
  geom_label_repel(data=anatgrc_f4map,aes(x=Longitude, y=Latitude,label=Region),
                   size=4,hjust=0,box.padding = 0.6) + theme(
                     axis.text = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks = element_blank(),legend.position = "none") +
  geom_star(data= anatgrc_f4map,aes(x=Longitude, y=Latitude,
                                               fill=V5,color=V5),starshape=1,starstroke=0.8,size=6,show.legend = F)+
  scale_x_continuous(limits = c(20, 43), expand = c(0, 0))+
  scale_y_continuous(limits = c(33.5, 41.5), expand = c(0, 0)) +
  scale_fill_manual(values=rev(c("#A3B900","#C34F74")))+
  scale_color_manual(values=darken(rev(c("#A3B900","#C34F74")),0.4))+
  theme_bw(base_rect_size = 0.5)+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank());anatgrc_f4map_plot


#geographic distance versus genetic distance

geo_f3_dist_plot=ggplot()+
  geom_point(data =geo_f3_dists, aes(geodist,f3dist,col=pairtype,fill=pairtype),shape=21,size=3,stroke=1.5 )+
  theme_bw(base_line_size = 0.2,base_rect_size = 0.2) + scale_color_manual("Pair type",values = lighten(c(lighten("#A73088",0.3),"orange","#A3B900"),0.1))+
  scale_fill_manual("Pair type",values = darken(c("#A73088","orange","#A3B900"),0.3))+
  xlab("\nGeographic distance (m)") + ylab("Genetic distance (1-f3)\n")
  

#f3 mds

realdata_f3mds_plot=ggplot(data=realdata_f3mds,aes(X1,X2,fill=Population))+
  geom_encircle(alpha=0.2,color="white",expand = 0.06,s_shape=0.5)+
  geom_point(aes(X1,X2,fill=Population,col=Population),shape=21,size=3,stroke=1.1)+
  ylim(-0.4,0.4)+
  scale_fill_manual(values=rev(c("#A3B900","#A73088")))+
  scale_color_manual(values=darken(rev(c("#A3B900","#A73088")),0.4))+
  theme_bw(base_line_size = 0.2)+ xlab("Coordinate 1")+ylab("Coordinate 2")+
  theme(
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size=10),
    plot.title = element_text(face = "bold",size = 15),
    axis.title.x = element_text(size = 12),
    #plot.margin = unit(c(1,0.5,0.5,0.5), "cm")
    );realdata_f3mds_plot


#f3 diversity

realda_f3div_plot=ggplot(data=f4sim_real_baf3_div,aes(y=1-f3,x=pop1,col=pop1,fill=pop1))+
  geom_boxplot(alpha=0.9,width=0.5,show.legend = F)+
  geom_point(size=2,alpha=0.4,show.legend = F)+ theme_bw()+
  scale_fill_manual( values=c("#A3B900","#A73088"))+
  scale_color_manual(values=darken(c("#A3B900","#A73088"),0.4))+
  stat_summary(
    geom = "point",
    fun = median,
    size = 3,
    shape=8,
    stroke=1.6,show.legend = F
  ) + xlab("Population")+ ylab("Intra-regional Diversity (1-f3)")+
  theme(
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size=10),
    plot.title = element_text(face = "bold",size = 15),
    axis.title.x = element_text(size = 12),
    #plot.margin = unit(c(1,0.5,0.5,0.5), "cm")
  )

#f3 histogram

f3_dist_plot=ggplot()+
  geom_histogram(data=f3comp_mergedba,aes(f3,fill=comp,y=..density..),alpha=0.6)+
  geom_density(data=f3comp_mergedba,aes(f3,fill=comp),alpha=0.3,size=0.5)+
  scale_fill_manual("Pair type",values = c("#A73088","orange","#A3B900"))+
  theme_bw(base_rect_size = 0.2,base_line_size = 0.2)+ylab("Density\n")+xlab("\nf3")+
  theme(axis.title = element_text(size=12))
 
