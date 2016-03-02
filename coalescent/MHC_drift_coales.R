cat('##############################\n')
cat('##        MHC mammoth        #\n')
cat('##  Coalescent simulations   #\n')
cat('##############################\n')
date()

options(scipen=999)

grid.size = 10
pop.size = round(10 * ((1000000/10) ^ (1/(grid.size-1))) ^ (0:(grid.size-1)))
bot.size = round(10 * ((1000/10) ^ (1/(grid.size-1))) ^ (0:(grid.size-1)))
anc.size = 100000
nsims = 1000

files = 'Locus1'
fsc = '~/software/fsc_linux64/fsc'
arlsumstats= 'arlsumstat3512_64bit'
def.path = '~/projects/mammoth/coalescent/'
system(paste0('mkdir ', files))

resmat = matrix(NA, ncol=length(bot.size), nrow=length(pop.size))

for (i in 1:length(pop.size)){
  for (j in 1:length(bot.size)){
    
    cat(paste0('>>> pop.size ', i, ' out of ', length(pop.size), ' | bot.size ', j,
               ' out of ', length(bot.size), '\n'))
    
BOTEND = bot.size[j]/pop.size[i]
BOTINIT = anc.size/bot.size[j]

parfile <- paste0('//Number of population samples (demes)
25
//Population effective sizes (number of genes)
', pop.size[i],
'0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
//Sample sizes
0
2 120
2 126
2 127
2 136
2 142
2 143
2 144
2 157
2 170
2 181
2 212
2 242
2 387
2 400
2 450
2 582
2 640
2 645
2 845
2 1226
2 1495
2 1529
2 1539
2 1581
//Growth rates	: negative growth implies population expansion
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
26  historical event
120 1 0 1 1 0 0
126 2 0 1 1 0 0
127 3 0 1 1 0 0
136 4 0 1 1 0 0
142 5 0 1 1 0 0
143 6 0 1 1 0 0
144 7 0 1 1 0 0
157 8 0 1 1 0 0
170 9 0 1 1 0 0
181 10 0 1 1 0 0
212 11 0 1 1 0 0
242 12 0 1 1 0 0
381 0 0 0 ', BOTEND, ' 0 0 //bottleneck ends
386 0 0 0 ', BOTINIT, ' 0 0 //bottleneck starts
387 13 0 1 1 0 0
400 14 0 1 1 0 0
450 15 0 1 1 0 0
582 16 0 1 1 0 0
640 17 0 1 1 0 0
645 18 0 1 1 0 0
845 19 0 1 1 0 0
1226 20 0 1 1 0 0
1495 21 0 1 1 0 0
1529 22 0 1 1 0 0
1539 23 0 1 1 0 0
1581 24 0 1 1 0 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA 99 0.00000 0.0000001 0.33')
write.table(parfile, file = paste0(files, '_', i, '_', j, '.par'), sep='\t',
            col.names=F, row.names=F, append=F, quote=F)
    
    cat(paste0('... simulating genealogies\n'))

    # launch fastsimcoal
    waste <- system(paste0(fsc, ' -i ', paste0(files, '_', i, '_', j, '.par'), ' -n', nsims, ' -c 8 -B 12'),
           intern = T)
    system(paste0('mv ', files, '_', i, '_', j, ' ', files))
    

    # Reshape .arp files with the correct groups
    ls.arp <- list.files(path=paste0(files, '/', files, '_', i, '_', j), pattern='.arp')
    for (d in 1:length(ls.arp)){
        x <- readLines(paste0(files, '/', files, '_', i, '_', j, '/', ls.arp[d]))
        y <- readLines('data/structure.txt')
        write.table(c(x[1:grep('StructureName', x)],y),
                    file = paste0(files, '/', files, '_', i, '_', j, '/', ls.arp[d]), 
                    sep='\t', col.names=F, row.names=F, append=F, quote=F)
    }
    
    cat(paste0('... sumarizing and analyzing output\n'))
    
    # execute ArlSumStat
    system(paste0('cp arlsumstat/* ', files, '/', files, '_', i, '_', j))
    setwd(paste0(files, '/', files, '_', i, '_', j))
    waste <- system(paste0('./', arlsumstats, ' ./', files, '_1_1.arp Stats_',
                  files, '_', i, '_', j, '.txt 0 2 run_silent'), intern = T)
    waste <- system(paste0('for file in *.arp; do ./', arlsumstats, ' ./$file Stats_',
                  files, '_', i, '_', j, '.txt 1 0 run_silent; echo "Processing file $file"; done'),
           intern = T)
    system('rm -r *.res')
    system(paste0('cp Stats_', files, '_', i, '_', j, '.txt ', def.path))
    setwd(def.path)
    
    # look at the stats and calculate prob. of observing at least 72% alleles in the
    # wrangler population v ancestral population
    stats <- read.table(paste0('Stats_', files, '_', i, '_', j, '.txt'), quote="\"", comment.char="")
    p <- sum(stats[,2]*0.72 <= stats[,1])/nsims
    resmat[i,j] = p
    
    #clean main folder
    system(paste0('mv *',files, '_', i, '_', j, '* ', files, '/', files, '_', i, '_', j))
    
    }
  }

write.table(resmat, file = paste0('resmat_', files,'.txt'), sep='\t', col.names=F, row.names=F, 
            append=F, quote=F)


cat('###########################\n')
cat('#     Plotting results    #\n')
cat('###########################\n')

# ## Black, blue, yellow, red
# resmat = as.matrix(resmat)
# matrix_0 = matrix(NA, nrow=grid.size, ncol=grid.size)
# matrix_0[resmat <= 0.001] = resmat[resmat <= 0.001]
# matrix_1 = matrix(NA, nrow=grid.size, ncol=grid.size)
# matrix_1[resmat <= 0.05 & resmat > 0.001] = resmat[resmat <= 0.05 & resmat > 0.001]
# matrix_2 = matrix(NA, nrow=grid.size, ncol=grid.size)
# matrix_2[resmat > 0.05] = resmat[resmat > 0.05]
# 
# pal.2 <- colorRampPalette(c('yellow', 'orange', 'red', 'darkred'))
# #pal.1 <- colorRampPalette(c('grey99', 'grey80'))
# pal.1 <- colorRampPalette(c('darkblue', 'blue'))
# #pal.1 <- colorRampPalette(c('black'))
# 
# pdf(paste0('resmat_', grid.size,'.pdf'), width = 10, height = 8)
# 
# par(fig=c(0,0.85,0,1), new=TRUE, mar=c(5.5,5.5,4,1)+0.1, mgp=c(4,1,0))
# 
# image(matrix_0, xaxt='n', yaxt='n', col='black', zlim=c(0, 0.001))
# image(matrix_1, xaxt='n', yaxt='n', col=pal.1(100), zlim=c(0.001, 0.05), add=T)
# image(matrix_2, xaxt='n', yaxt='n', col=pal.2(100), zlim=c(0.05, 1), add=T)

resmat = as.matrix(resmat)
matrix_1 = matrix(NA, nrow=grid.size, ncol=grid.size)
matrix_1[resmat <= 0.05] = resmat[resmat <= 0.05]
matrix_2 = matrix(NA, nrow=grid.size, ncol=grid.size)
matrix_2[resmat > 0.05] = resmat[resmat > 0.05]

pal.2 <- colorRampPalette(c('yellow', 'orange', 'red'))
pal.1 <- colorRampPalette(c('grey99', 'grey80'))

pdf(paste0('resmat_', files, '_', grid.size,'.pdf'), width = 10, height = 8)

par(fig=c(0,0.85,0,1), mar=c(5.5,5.5,4,1)+0.1, mgp=c(4,1,0))

image(matrix_1, xaxt='n', yaxt='n', col=pal.1(100), zlim=c(0, 0.05))
image(matrix_2, xaxt='n', yaxt='n', col=pal.2(100), zlim=c(0.05, 1), add=T)
axis(side=1, at=seq(0,1, length.out=grid.size), labels=round(pop.size), las=2, cex.axis=1)
axis(side=2, at=seq(0,1, length.out=grid.size), labels=round(bot.size), las=1, cex.axis=1)
title(main = paste0('Exon 2'), cex.main=1.5, ylab= "Population size at bottleneck",
      xlab= "Population size after bottleneck", cex.lab=0.8)
box()


# set place for legend
par(fig=c(0.80,0.96,0.05,0.95), new=TRUE)

image(t(as.matrix(100:1)), col=c(pal.1(19), pal.2(80)), xlab='', ylab='', xaxt='n', yaxt='n')
axis(2, cex.axis=0.9, las=2, at=seq(0, 1, length.out = 26),
     labels=c('1.0','','0.9', '','0.8', '','0.7', '','0.6','', '0.5','', '0.4', '','0.3','',
              '0.2', '','0.1','', '0.05', '', '0.03', '', '0.01', ''))
mtext(text='Probability of observing 72% of alleles in 12 Wrangler samples', side=4, line=1)
abline(a=0.802, b=0.0, lwd=1)
box()

dev.off()

cat('Analysis finished\n')
date()





