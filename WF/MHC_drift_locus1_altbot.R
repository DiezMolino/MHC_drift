##############################
##      MHC mammoth
## Genetic drif simulations
##############################

#################################
# Load packages
###############################

cat('Loading packages...\n')
# load required libraries
x <- c('foreach', 'doMC', 'data.table')
lapply(x, require, character.only=T)

# parallel backend
registerDoMC(7) # ************* WATCH OUT (12) ONLY FOR LEGION!!! (6-7 for MY laptop) **********

# change scientific notation in R
options(scipen=999)

#################################
# Parameters
###############################

gens = 260 #((12K-4K)/31y/gen)
freq.init = c(0.16667, 0.04167, 0.16667, 0.08333, 0.29167, 0.2, 0.16667, 0.2, 0.1, 0.25, 0.16667)
n.polymor = length(freq.init)
n.all.wrangler = 3
nsim = 100
grid.size = 30
path= 'locus1_altbot'

# pop.size <- seq(10,
#                 10000,
#                 length.out = grid.size)
pop.size = round(10 * ((1000000/10) ^ (1/(grid.size-1))) ^ (0:(grid.size-1)))

bot.size = round(seq(0.1,0.99, length.out = grid.size),2)

system(paste0('mkdir ', path, '_', grid.size))

#################################
# Simulations
#################################

resmat = matrix(NA, nrow=length(pop.size), ncol=length(bot.size))

for (i in 1:length(pop.size)){
    for (j in 1:length(bot.size)){
      
      # create a list of Ne values for each generation
      # EXPONENTIAL GROWTH
      #Ne.list = round(bot.size[y] * ((pop.size[i]/bot.size[y]) ^ (1/(gens-1))) ^ (0:(gens-1)))
      # Linear growth
      #Ne.list = round(seq(bot.size[y], pop.size[i], length.out = gens))
      
      freq.log = matrix(NA, nrow=nsim, ncol=length(freq.init))
      
      for (s in 1:nsim){
#     grid <- obsprop <- foreach(j = 1:nsim, .combine='c', .inorder=TRUE) %dopar%
#       {
        freq = rbinom(length(freq.init),
                      round(pop.size[i]*bot.size[j]),
                      freq.init)/round(pop.size[i]*bot.size[j])

        for (e in 1:gens){
          freq = rbinom(length(freq),
                        round(pop.size[i]),
                        freq)/round(pop.size[i])
          # TO-DO skip all generations for this simulation if all frequencies are equal to 0
          
         }
        freq.log[s, ] = freq
      }
      write.table(freq.log,
                  file=paste0(path, '_', grid.size,'/freq.log_', j, '_', i, '.txt'),
                  sep= '\t',
                  col.names = F,
                  row.names = F,
                  quote = F,
                  append = F)
      
      # Sample 24 chromosomes from each simulation and check polymorphism at each position
      alleles.log = matrix(NA, nrow=nsim, ncol=ncol(freq.log))
  
      for (e in 1:nsim){
        for (d in 1:ncol(freq.log)){
          alleles.log[e,d] = length(unique(rbinom(24,1,freq.log[e,d]))) == 2
        }
      }
      write.table(alleles.log,
                  file=paste0(path, '_', grid.size,'/alleles.log_', j, '_', i, '.txt'),
                  sep= '\t',
                  col.names = F,
                  row.names = F,
                  quote = F,
                  append = F)
      
      resmat[i,j] = sum(apply(alleles.log, 1, sum) >= 8)/nsim
      cat(paste0('... pop.size ', i, ' out of ', length(pop.size), ' | bot.size ', j, ' out of ', length(bot.size), '\n'))
    }
}
write.table(resmat,
            file=paste0('resmat_', grid.size,'.txt'),
            sep= '\t',
            col.names = F,
            row.names = F,
            quote = F,
            append = F)

#################################
# Plot
###############################

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

pdf(paste0('resmat_', grid.size,'.pdf'), width = 10, height = 8)

par(fig=c(0,0.85,0,1), mar=c(5.5,5.5,4,1)+0.1, mgp=c(4,1,0))

image(matrix_1, xaxt='n', yaxt='n', col=pal.1(100), zlim=c(0, 0.05))
image(matrix_2, xaxt='n', yaxt='n', col=pal.2(100), zlim=c(0.05, 1), add=T)
axis(side=1, at=seq(0,1, length.out=grid.size), labels=round(pop.size), las=2, cex.axis=1)
axis(side=2, at=seq(0,1, length.out=grid.size), labels=round(bot.size,2), las=1, cex.axis=1)
title(main ='Exon 2 | ', n.polymor,' polymorphic positions', cex.main=1.5, ylab= "Population size at bottleneck", xlab= "Population size after bottleneck", 
      cex.lab=0.8)
box()


# set place for legend
par(fig=c(0.80,0.96,0.05,0.95), new=TRUE)

image(t(as.matrix(100:1)), col=c(pal.1(19), pal.2(80)), xlab='', ylab='', xaxt='n', yaxt='n')
axis(2, cex.axis=0.9, las=2, at=seq(0, 1, length.out = 26),
     labels=c('1.0','','0.9', '','0.8', '','0.7', '','0.6','', '0.5','', '0.4', '','0.3','',
              '0.2', '','0.1','', '0.05', '', '0.03', '', '0.01', ''))
mtext(text='Probability of observing ', n.all.wrangler, ' or more polymorphisms', side=4, line=1)
abline(a=0.802, b=0.0, lwd=1)
box()

dev.off()













######### JUNK

# resmat <- matrix(0, ncol = grid.size, nrow = grid.size)
# for (i in 1:length(pop.size)){
#   for (j in 1:length(bot.size)){
#     alleles.log <- read.delim(paste0(path, '_',freq.init, '_', grid.size,'/alleles.log_', j, '_', i, '.txt'),
#                               header=FALSE)
#     prob.nofix <- sum(alleles.log==1)/nsim
#     resmat[i,j] <- 1-2*abs(0.5-prob.nofix)
#   }
# }


# 
# image(resmat, axes = F)
# axis(1, 
#      at=seq(0,1, length.out = grid.size),
#      labels = round(pop.size),
#      las=2,
#      cex.axis=0.3)
# axis(2,
#      at=seq(0,1, length.out = grid.size),
#      labels = round(bot.size),
#      las=2,
#      cex.axis=0.3)
# dev.off()



#           # save freq.trace
#           write.table(freq.trace,
#                       file=paste0('freq.trace/freq.trace_', pop.size[i], '.txt'),
#                       sep= '\t',
#                       col.names = F,
#                       row.names = F,
#                       quote = F,
#                       append = F)
#           # plot freq.trace
#           pdf(paste0('freq.trace/freq.trace_', pop.size[i], '.pdf'),
#               width = 7,
#               height = 5)
#           plot(freq.trace[,1],
#                col='white',
#                ylim=c(0,1), 
#                main=paste0('Frequency trajectories | Ne = ', pop.size[i]),
#                ylab = 'Allele frequency',
#                xlab = 'Generation')
#           for(u in 1:nsim){
#             lines(freq.trace[,u], col='grey')
#           }
#           dev.off()


# 
# plot(two.tail.p~pop.size, type='l', lwd=3, col='grey80', xlab=expression(N[E]),
#      ylab='Porb. observing allele', ylim=c(0,1), xaxt='n')
# axis(1,at=pop.size, labels = round(pop.size), las=2, cex.axis=0.8)
# points(two.tail.p~pop.size, pch=16)




