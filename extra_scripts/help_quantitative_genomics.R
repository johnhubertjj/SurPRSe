#set.seed(6147)
#set.seed(6154)
library(ggplot2)
library(gridExtra)
library(reshape2)

set.seed(6155)
m = 500                                  # number of SNPs
maf = runif(m, 0, .5)                    # random MAF for each SNP
# apologies if using "=" as the assignment operator in R makes your eyes hurt

x012 = rbinom(m, 2, maf)

n = 400                                    # number of individuals
x012 = t(replicate(n, rbinom(m, 2, maf)))  # n x m genotype matrix

x012 = t(replicate(n, rbinom(2*m, 2, c(maf, maf))))
polymorphic = apply(x012, 2, var) > 0
x012 = x012[,polymorphic][,1:m]
maf = c(maf, maf)[polymorphic][1:m]
round(maf[1:10], 2)
x012[1:5, 1:10]

x_o = scale(x012, scale=FALSE)               # mean 0
x = scale(x012, scale=TRUE)                  # mean 0, variance 1

maf_est = colMeans(x012)/2
qplot(maf, maf_est, col=maf) +
  geom_abline() + xlab('p') + ylab(expression(hat(p)))

ss = .02
sequ = seq(0, max(maf_est), ss)
bins = cut(maf_est, sequ)
binMeans = (sequ + ss/2)[-length(sequ)]
dat = data.frame(observed=tapply(colMeans(x012)/2 - maf, bins, var),
                 expected=binMeans*(1-binMeans)/(2*n),
                 binNum=1:length(binMeans))

ggplot(dat, aes(expected, observed)) +
  geom_text(aes(label=binNum, col=binMeans)) +
  geom_abline() +
  xlab(expression(paste(hat(SE)[hat(p)]^2, "  (",hat(var),"(", hat(p), " | p))"))) +
  ylab(expression(paste("var(",hat(p)," | p)")))

#Variance of a genotype
varx = apply(x012, 2, var)
p1 = qplot(2*maf*(1-maf), varx, col=maf) + geom_abline()
# genotype frequencies
p2 = ggplot(data.frame(x=c(0, 1)), aes(x)) +
  stat_function(fun=function(x) 2*x*(1-x), col='red') 
  stat_function(fun=function(x) x^2, col='green') +
  stat_function(fun=function(x) (1-x)^2, col='blue') +
  ylim(c(0,1)) +
  xlab('allele frequency') +
  ylab('genotype frequencies / variance of x (red)')
grid.arrange(p1, p2, ncol=2)

#LD matrix
ld1 = cor(x012)
ld2 = cor(x)
ld3 = cov(x)
ld4 = (t(x) %*% x) / n

# LD scores

ldscores_sample = colSums(ld1^2)
ldscores = (ldscores_sample*n - m) / (n + 1)

qplot(ldscores, ldscores_sample, col=ldscores) +
  xlim(c(min(ldscores), max(ldscores))) +
  ylim(min(ldscores), max(ldscores_sample)) +
  geom_abline() +
  scale_colour_continuous(low='black', high='green')

# Genetic relatedness matrix 

grm = (x %*% t(x))/m
grm2 = cov(t(x))

p1 = ggplot(melt(grm), aes(value)) + geom_histogram(bins=100) 
p2 = ggplot(melt(grm), aes(value)) + geom_histogram(bins=100) +
  coord_cartesian(ylim=c(0,n)) + theme(panel.background=element_blank())
grid.arrange(p1, p2, layout_matrix=matrix(c(1,1,2,1), 2))

xt = scale(t(x), scale=F) # transpose and mean center across individuals

principal_axes_eigen        = eigen(grm2)$vectors
principal_axes_svd          = svd(xt)$v
principal_axes_prcomp       = prcomp(xt)$rotation

principal_components_eigen  = xt %*% eigen(grm2)$vectors
principal_components_svd1   = xt %*% svd(xt)$v
principal_components_svd2   = svd(xt)$u %*% diag(svd(xt)$d)
principal_components_prcomp = prcomp(xt)$x

cor(principal_components_eigen[,1], principal_components_eigen[,2])

p1 = qplot(1:min(n,m), prcomp(xt)$sdev^2) + xlab('rank') + ylab('eigenvalues')
p2 = qplot(principal_components_eigen[,1], principal_components_eigen[,2]) +
  xlab('PC1') + ylab('PC2')
grid.arrange(p1, p2, ncol=2)

#Simulating genotypes with LD
x012ld = jitter(x012[,rep(1:m, 1:m)[1:m]], .03)
xld = scale(x012ld)
ldld = (t(xld) %*% xld)/n
grmld = (xld %*% t(xld))/m
ldscores_sampleld = colSums(ldld^2)
ldscoresld = (ldscores_sampleld*n - m) / (n + 1)
varxld = apply(x012ld, 2, var)

greens = colorRampPalette(c('black', 'green'))(12)
par(mfrow=c(1,2))
image(ld1, col=greens)
image(ldld, col=greens)
par(mfrow=c(1,1))

#Number of effective SNPs
(meld = m/mean(ldscoresld))
1/var(grmld[upper.tri(grmld)])

#Number of effective people
1/var(ld1[upper.tri(ld1)])
recomb = sqrt(2*m)/m
1/(var(ldld[upper.tri(ldld)]) * recomb)

