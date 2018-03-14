X <- data.frame( V1 = LETTERS, V2 = runif( 26 ),
                 V3 = rep( c("F1", "F2"), each = 13 ) )

v <- ifelse( X$V1 %in% c( "D", "O", "T" ), "red", "black" )
g <- ggplot( X, aes( x = V1, y = V2 ) ) + geom_point() +
  theme( axis.text.x = element_text( color = v ) )

g + facet_wrap( ~V3 )
g + facet_wrap( ~V3, scales = "free" )

gt <- ggplotGrob( g + facet_wrap( ~V3, scales = "free" ) )
gt
gt$grobs[[6]]$children[[2]]$grobs[[2]]$children[[1]]$gp

pt$grobs[[6]]$children[[2]]$grobs[[2]]$children[[1]]$gp$col <- "orange"
pt$grobs[[7]]$children[[2]]$grobs[[2]]$children[[1]]$gp


pt <- ggplotGrob(p)

grid::grid.draw(pt)
