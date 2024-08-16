#Copyright (c) 2016 ZarLab

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

rm(list=ls())
library(gtools);
library(mvtnorm);
library(vegan);
############### define functions ###############
getp <- function(Y, x, p) {
    res = adonis(Y ~ x, perm=p)
    return(res$aov.tab$"Pr(>F)"[1])
}
getF <- function(Y, x, p) {
    res = adonis(Y ~ x, perm=p)
    return(res$aov.tab$F.Model[1])
}
gamma <- function(Y, x, max_itr=4) {
    for (i in 2:max_itr) {
        p = 10^i
        limit = 5/p
        pval = getp(Y, x, p)
        if (pval > limit) {
            break
        }
    }
    return(pval)
}

chol_solve <- function(K) {
    a = eigen(K)$vectors
    b = eigen(K)$values
    b[b<1e-13] = 1e-13
    b = 1/sqrt(b)
    return(a%*%diag(b)%*%t(a))
}

rotate <- function(Y, sigma) {
    U <- chol_solve(sigma)
    tU <-t(U)
    UY = tU%*%Y
    return(UY)
}

run_gamma <- function(Y, X) {
    Ng = dim(X)[2]
    pval = 1:Ng
    fval = 1:Ng
    newY = Y-min(Y)
    for (i in 1:Ng) {
        pval[i] = gamma(newY, X[,i])
        fval[i] = getF(newY,X[,i],1)
        cat(i,". f=",fval[i]," p=",pval[i],"\n")
    }
    return(list(pval, fval))
}
################## get input ###################
args=(commandArgs(TRUE))
if(length(args)!=5){
        print("Usage: R CMD BATCH --args -Xpath= -Kpath= -Ypath= -VCpath= -outputPath= GAMMA.R GAMMA.log")
        stop ()
}
for (i in 1:5){
        paramName =  strsplit(args[i],"=")[[1]][1]
        param =strsplit(args[i],"=")[[1]][2]
        if(paramName=="-Xpath")
                Xpath = param
        else if(paramName=="-outputPath")
                outputPath = param
        else if(paramName=="-Kpath")
        	Kpath = param
        else if(paramName=="-Ypath")
                Ypath = param
	else if(paramName=="-VCpath")
                VCpath = param
	else {
                cat("Error: Wrong parameter name ",paramName)
                stop()
        }
}
print(Xpath)
print(Kpath)
print(Ypath)
print(VCpath)
print(outputPath)
################## GAMMA  ##################
X = as.matrix(read.table(Xpath)) # read input
K = as.matrix(read.table(Kpath))
Y = as.matrix(read.table(Ypath))
VC = as.matrix(read.table(VCpath))
snpNum <- dim(X)[2]
indiNum <- dim(X)[1]
geneNum <- dim(Y)[2]
Vg = median(VC[,1])		# Variance components
Ve = median(VC[,2])
I = diag(indiNum)
sigma = Vg*K + Ve*I
UY = rotate(Y,sigma)		# Rotate genotypes and phenotypes
UX = rotate(X,sigma)
ps2 = run_gamma(UY, UX)		# GAMOVA
write.table(ps2[1], paste(outputPath, "/P.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(ps2[2], paste(outputPath, "/F.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
