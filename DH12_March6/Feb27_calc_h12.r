library(stringi)
library(data.table)
library(magrittr)
library(future.apply)

plan(tweak(multisession, workers = 30))

# Load all of the H12 tables
setwd('/home/namulil/lstm_projects/funestus_llineup/notebooks/DH12_permutations')
h12.filenames <- list.files(pattern = '\\.tsv$', full.names = TRUE)# retrieve all tsv files in the DH12_permutation folder

print(h12.filenames)## checking if my tsv files have been loaded 


#study.pops <- setNames(nm = c('siaya'))
pops <- setNames(nm = c('NonPBO', 'PBO'))# create named vectors NonPBO and PBO (empty)
str(pops) # check nature of pop
# Get the randomisation ids. ## raed fisrt colum and get radomisation samples ids
randomisation.ids <- readLines(h12.filenames[1], n = 1) %>%## read first column of h12 files
                     {strsplit(., '\t')[[1]]} %>%## tab separate columns
                     grep('r\\d{4}', ., value = T)## read out only r' colums
str(randomisation.ids )## check if it has worked

# A function that looks for positive peaks by identifying windows more extreme than thrice the 
# 95% centile:
find.peaks <- function(values, centile = 0.95, multiplier = 3){
    center <- median(values)
	thresh1 <- center + (quantile(values, centile) - center) * multiplier
	thresh2 <- center + (quantile(values, 1-centile) - center) * multiplier
	(values > thresh1) | (values < thresh2)
}

# A function to load and combine all data for a given randomisation id
load.and.combine.data <- function(pop, find.peaks =T, calculate.P.values = T){
	
	## check and load/create empty files
	pre.filenames <- grep(paste('_', pop, '\\.pre', sep = ''), h12.filenames, value = T) %>%  ## in pop search for .pre files in h12 files and return matching file
 	                 setNames(., stri_extract_first_regex(., '(2RL|3RL|X)+(?=\\.tsv)'))## extract names out and label them according to chromosome
	post.filenames <- grep(paste('_', pop, '\\.post', sep = ''), h12.filenames, value = T) %>% 
	                  setNames(., stri_extract_first_regex(., '(2RL|3RL|X)+(?=\\.tsv)'))
	
	if (any(names(post.filenames) != names(post.filenames)))
		stop(paste('Different chromosomes have been found for pre and post samplesets for randomisation', rand.id))## check if pre and post chromosome name labelling differes
	

	print("Loaded population names:")
	#print(names(h12.tables))

	## load actual data oders it into pre and post data tables
	pre.data <- names(pre.filenames) %>%
	              lapply(function(chrom){ 
	                  fread(pre.filenames[chrom], sep = '\t') %>%## read pre.files
	                  .[, chromosome := ..chrom] %>%## add new column chromosome
					  setnames('h12', 'H12') %>% ## rename h12to H12
					  setcolorder(c('chromosome', 'startpoint', 'endpoint', 'midpoint')) ## set column names in the df to be saved
				  }) %>%
				  rbindlist()## combine all chromosome data into a single data table
	post.data <- names(post.filenames) %>%## do the same for post
	             lapply(function(chrom){ 
	                 fread(post.filenames[chrom], sep = '\t') %>%
	                 .[, chromosome := ..chrom] %>%##
					 setnames('h12', 'H12') %>%
					 setcolorder(c('chromosome', 'startpoint', 'endpoint', 'midpoint'))
				 }) %>%
				 rbindlist()
	## check if files were loaded correctly
	colnames(pre.data)
	colnames(post.data)
	# Check that the midpoints are identical
	if (!identical(pre.data$midpoint, post.data$midpoint))
		stop('Window positions do no match between pre and alive samples')
	
	## screen difference between pre and post data while maintaining chromosome position and cobine chromosome position of the calculated dif nto new data datable
	diff.data <- cbind(post.data[, .(chromosome, startpoint, endpoint, midpoint)],
	                   post.data[, c('H12', ..randomisation.ids)] - pre.data[, c('H12', ..randomisation.ids)]
	)
	
	# Look for peaks in the diff data
	if (find.peaks){
		diff.data[, is.peak := find.peaks(H12, centile = 0.98)]
		setcolorder(diff.data, c('chromosome', 'startpoint', 'endpoint', 'midpoint', 'H12', 'is.peak'))
	}
	
	## calculate pvalues for peak, by, subtracting observed h12 from randomised values, store output in pval
	if (calculate.P.values){
		if(find.peaks){
			diff.data[is.peak == T, pval := apply(.SD - H12, 
			                                      1, 
			                                      # get TWO-TAILED p-value
			                                      function(x) 1-(abs(sum(x >= 0)/length(x) -0.5)*2)
			                                ), 
			                                .SDcols = randomisation.ids
			] 
		}
		else {
			diff.data[, pval := apply(.SD - H12, 
			                          1, 
			                          function(x) 1-(abs(sum(x >= 0)/length(x) -0.5)*2)
			                    ), 
			                    .SDcols = randomisation.ids
			] 
		}
	}
	
	output <- list(pre = pre.data, post = post.data, diff = diff.data)
	output
}

cat('Loading data\n')
h12.tables <- lapply(pops, function(pop) {cat(pop, '\n'); load.and.combine.data(pop)})
#print(h12.tables) ## commenting out this as code below achives the same in a better way

## being chicky and trying to write out h12 files as csv files 


output_dir <- "."  

# Loop through populations and write CSVs
for (pop in names(h12.tables)) {
    cat("Saving:", pop, "\n")  # Print progress
    
    # Save each dataset: pre, post, and diff
    fwrite(h12.tables[[pop]]$pre, file = paste0(pop, "_pre.csv"), sep = ",")
    fwrite(h12.tables[[pop]]$post, file = paste0(pop, "_post.csv"), sep = ",")
    fwrite(h12.tables[[pop]]$diff, file = paste0(pop, "_diff.csv"), sep = ",")
}

cat("✅ All files have been saved in the current directory!\n")

## check h12 table names 
length(h12.tables)
str(h12.tables)
source("./Harun_R_plotting_jan25.r")

#chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Function to plot the H12 data and visualise h12 values
plot.h12.diff <- function(h12.table, 
                          filename = NULL, 
                          num.randomisations = NULL, 
                          p.thresh = 0.01, 
                          plot.title = '', 
                          gaps = 5000000, 
                          filter.name = 'is.peak'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 3.5
		if (grepl('\\.eps', filename))
			prescript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	# Get vectors of start and end points for each chromosome (ie: cumulative sizes + gaps)
	chrom.sizes <- c('2RL' = 102883511, '3RL' = 84636641, 'X' = 22264324)
	ce <- cumsum(chrom.sizes + c(0, gaps,gaps))# dtermine where each chromosme starts and ends
	cs <- ce - chrom.sizes
	layout(matrix(c(rep(1,4),rep(2,1)), nrow = 5, ncol = 1))## plot layout and colour details
	colours <- c(h12 = lighten.col('red', alpha = 0.8),
                 randomisations = lighten.col('#cccccc', alpha = 0.15))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), ## extract h12 values for radomisation ids
	                             length(randomisation.ids),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.ids[seq(1, num.randomisations, length.out = num.randomisations)]
	h12.columns <- c('H12', r.columns)

	max.y <- max(c(max(h12.table[, ..h12.columns, with = FALSE]), 0.05))
	min.y <- min(h12.table[, ..h12.columns,  with = FALSE])
	# Create the empty plot.
	plot(c(cs[1], ce[3]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[3] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	#plot(c(0, 219784476), c(min.y, max.y), xlim = c(0 + gaps/2, 219784476 - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h12.table$genome.pos <- h12.table$midpoint + cs[as.character(h12.table$chromosome)]## covert chromosome specific positions to genome wide positions
	# Add randomised data
	sapply(r.columns, function(x) h12.table[, lines(genome.pos, get(x), col = colours['randomisations'], lwd = 0.4), by = chromosome])
	# Add true data
	h12.table[, lines(genome.pos, H12, col = colours['h12'], lwd = 1.2), by = chromosome]
	# Add y axis
	h12.step.size <- ifelse(max.y - min.y > 0.2, 0.1, 0.05) ## aded -min.y to try and see if i will get the negative values on my scale
	#axis(2, at = seq(min.y, max.y, h12.step.size))
	y_ticks <- seq(min.y, max.y, by = h12.step.size)  # Generate tick marks## aded on origanal 
	axis(2, at = y_ticks, labels = formatC(y_ticks, format = "f", digits = 1))  # Format labels## added on original
	mtext('H12', 2, 2, cex = 0.8)#'H12 difference
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('orchid3', 'green', 'blue')
	filter.pass <- h12.table[[filter.name]]
	h12.table[filter.pass, points(genome.pos, H12, 
	                              pch = 21,
	                              bg = p.colours[(pval < p.thresh[1]) + (pval < p.thresh[2]) + 1], 
	                              col = colours['randomisations'], 
	                              cex = 1.1, lwd = .5)
	]
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(gaps, chrom.sizes, gene.cex = 0.7, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}

p.threshold = 0.01




for (pop in names(h12.tables))
	plot.h12.diff(h12.tables[[pop]]$post, #post
	              filename = paste(pop, 'peak_filter_plot_post.png', sep = '_'), ## removed './DH12_permutations/' since I was working inside the same folder as I want my results to be stored 
	              p.thresh = p.threshold,
	              plot.title = pop)

				  

for (pop in names(h12.tables))
	plot.h12.diff(h12.tables[[pop]]$pre, #pre, 
	              filename = paste(pop, 'peak_filter_plot_pre.png', sep = '_'), ## removed './DH12_permutations/' since I was working inside the same folder as I want my results to be stored 
	              p.thresh = p.threshold,
	              plot.title = pop)

#Diff data()
for (pop in names(h12.tables))
	plot.h12.diff(h12.tables[[pop]]$diff, # diff
	              filename = paste(pop, 'Dh12_plot5.png', sep = '_'), ## removed './DH12_permutations/' since I was working inside the same folder as I want my results to be stored 
	              p.thresh = p.threshold,
	              plot.title = pop)
