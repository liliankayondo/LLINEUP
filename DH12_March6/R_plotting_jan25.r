# A function to lighten the tone of the colour. Either "color" or "lightness" can be a vector, but not both.
lighten.col <- function(color, lightness = 1, alpha = alpha){
	if (length(color) > 1 & length(lightness) > 1)
		stop('At least one of "color" and "lightness" must be of length 1')
	col.rgb <- col2rgb(color)/255
	if (length(lightness) == 1)
		return(rgb(t(1-(1-col.rgb)*lightness), alpha = alpha))
	else
		return(rgb(t(1-(1-col.rgb) %*% lightness), alpha = alpha))
}
chrom.sizes <- c('2RL' = 102883511, '3RL' = 84636641, 'X' = 22264324)
gaps <- 5000000
# Write a function to draw the chromosomes on an existing plotting device
add.chromosomes <- function( gaps=gaps,chrom.sizes = chrom.sizes, gene.cex = 0.9, gene.col = 'grey20', point.cex = 1.2, point.col = 'grey30', chrom.col = NULL, chrom.cex = 1.4, chrom.offset = 0){
  #gaps <- 5000000
  #chrom.sizes <- c('2RL' = 102883511, '3RL' = 84636641, 'X' = 22264324)
  ce <- cumsum(chrom.sizes + c(0, gaps, gaps))
  cs <- ce - chrom.sizes 
  plot(c(cs[1], ce[3]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[3] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  # show the kdr region
  kdr.region.mean <- cs['3RL'] + mean(c(44105643, 44156624))
  points(kdr.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(kdr.region.mean, -1.7, 'Vgsc', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
  # show the Rdl region
  rdl.region.mean <- cs['3RL'] + mean(c(13537199, 13612792))
  points(rdl.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(rdl.region.mean, -1.7, 'Rdl', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)

  ## show the new region on 3RL 
  protein.O.region.mean <- cs ['3RL'] + mean(c(67925771, 67931656))
  points(protein.O.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(protein.O.region.mean, -1.7, 'Protein-O-TMTC1', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
  

  # show the Ace1 gene region
  ace1.region.mean <- cs['2RL'] + mean(c(22199782, 22201327))
  points(ace1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(ace1.region.mean, -1.7, 'Ace1', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
  # show the CYP6 region
  cyp6.region.mean <- cs['2RL'] + mean(c(8671966, 8676311))
  points(cyp6.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(cyp6.region.mean, -1.7, 'Cyp6p', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
  
  ## show new complexin gene/ region on 2R
  complexin.region.mean<- cs['2RL'] + mean(c(51600427, 51705880))
  points(complexin.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(complexin.region.mean, -1.7, 'Complexin', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
    
  
  # show the GST region
  gst.region.mean <- cs['2RL'] + mean(c(76403081, 76411803))
  points(gst.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(gst.region.mean, -1.7, 'Gste2', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
  
  
  
  
  
  
  # show the DGK region
  dgk.region.mean <- cs['X'] + mean(c(13590465, 13736044))
  points(dgk.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(dgk.region.mean, -1.7, 'Dgk', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
  # show the CYP9K1 region
  cyp9k1.region.mean <- cs['X'] + mean(c(8448477, 8450887))
  points(cyp9k1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
  text(cyp9k1.region.mean, -1.7, 'Cyp9k1', srt = 45, adj = 1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
  
  
  
  
  
  # Plot the outline of the chromosomes
  chrom.sizes_2 <- c('2R' = 57604655, '2L' = 45278856, '3R' = 45028052, '3L' = 39608589, 'X' = 22264324)
  ce2 <- cumsum(chrom.sizes_2 + c(0, 0, gaps, 0, gaps))
  cs2 <- ce2 - chrom.sizes_2
  chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
  chrom.y <- -8.5 - chrom.offset
  lines(c(ce2['2R'], ce2['2R'] - gaps/2, cs2['2R'], cs2['2R'], ce2['2R'] - gaps/2, ce2['2R']), 
        c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['2R'])
  text((cs2['2R'] + ce2['2R'])/2, chrom.y, '2R', adj = 0.5, xpd = NA, cex = chrom.cex)
  
  lines(c(cs2['2L'], cs2['2L'] + gaps/2, ce2['2L'], ce2['2L'], cs2['2L'] + gaps/2, cs2['2L']), 
        c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['2L'])
  text((cs2['2L'] + ce2['2L'])/2, chrom.y, '2L', adj = 0.5, xpd = NA, cex = chrom.cex)
  
  lines(c(ce2['3R'], ce2['3R'] - gaps/2, cs2['3R'], cs2['3R'], ce2['3R'] - gaps/2, ce2['3R']), 
        c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['3R'])
  text((cs2['3R'] + ce2['3R'])/2, chrom.y, '3R', adj = 0.5, xpd = NA, cex = chrom.cex)
  
  lines(c(cs2['3L'], cs2['3L'] + gaps/2, ce2['3L'], ce2['3L'], cs2['3L'] + gaps/2, cs2['3L']), 
        c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['3L'])
  text((cs2['3L'] + ce2['3L'])/2, chrom.y, '3L', adj = 0.5, xpd = NA, cex = chrom.cex)
  
  lines(c(cs2['X'], cs2['X'], ce2['X'] - gaps/2, ce2['X'], ce2['X'], ce2['X'] - gaps/2, cs2['X']), 
        c(-1, 1, 1, 0.2, -0.2, -1, -1), lwd = 2, col = chrom.col['X'])
  text((cs2['X'] + ce2['X'])/2, chrom.y, 'X', adj = 0.5, xpd = NA, cex = chrom.cex)
}


draw.gene.model <- function(start.pos, end.pos, gene.positions, exon.positions, include.gene.names = T, y = 0, text.cex = 0.5, gene.thickness.fraction = 15, lwd = 2){
	lines(c(start.pos, end.pos), c(y, y), lwd = lwd, col = 'grey20')
	if (nrow(gene.positions) > 0){
		gene.positions <- gene.positions[order(start)]
		gr <- par('usr')
		gene.thickness <- (gr[4] - gr[3])/gene.thickness.fraction
		# Draw the genes
		v.adj = ifelse(gene.positions$strand == '+', 0, -gene.thickness)
		rect(
			gene.positions[,start], 
			y + v.adj, 
			gene.positions[,end], 
			y + gene.thickness + v.adj, 
			col = 'grey95',
			border = 'black',
			lwd = lwd,
			xpd = NA
		)
		# Draw the exons
		exon.v.adj = ifelse(exon.positions$strand == '+', 0, -gene.thickness)
		rect(
			exon.positions[,start], 
			y + exon.v.adj, 
			exon.positions[,end], 
			y + gene.thickness + exon.v.adj, 
			col = 'black',
			border = NA,
			xpd = NA
		)
		# Label the genes
		if (include.gene.names){
			text(
				apply(gene.positions[,.(start, end)], 1, mean), 
				y - 1.6*gene.thickness, 
				gene.positions$gene.name, 
				adj = 1, 
				srt = 45,
				cex = text.cex
			)
		}
	}
}

