require('tseries')
source('measures2.R')
trn    <- read.matrix("TRNFILE")
trnout <- matrix(data=NA,nrow=nrow(trn),ncol=14)
for (i in 1:nrow(trn))
{
 i
 trnout[i,1] <- trn[i,1]
 trnout[i,-1] <- measures(trn[i,-1])
}
write(t(trnout), "TRNFILE.features",ncolumns=14)
