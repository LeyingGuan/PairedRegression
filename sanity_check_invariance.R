library(PREGS)
set.seed(2024)
B = 2000
M = 2000
args = list()
args$D = "Paired"; args$E = "multinomial"; args$S=0; args$p=20;args$n=100

dat = data_gen(n=as.integer(args$n), p = as.integer(args$p), M = M, design =args$D, noise = args$E,seed = seeds[it])
if(as.numeric(args$S) <=0.0){
  dat$beta = 0.0
}else{
  dat$beta  = beta_gen(dat$x, dat$Z, dat$epsMat, power = as.numeric(args$S))
}
dat$y = y_gen(dat$x, dat$beta, dat$epsMat)


##Paired regression sequential
pval_PRegs_twosided_coef = rep(NA, ncol(dat$y))
pval_PRegs_onesided_coef = rep(NA, ncol(dat$y))
for(l in 1:ncol(dat$y)){
  print(l)
  tmp1coef <- PREG(y=dat$y[,l], x=dat$x, z=dat$Z, add_intercept = T, B = B, statType = "coef")$res
  pval_PRegs_twosided_coef[l] =  tmp1coef$unsigned
  pval_PRegs_onesided_coef[l] = 2*min(tmp1coef$pos,tmp1coef$neg)
}
par(mfrow  = c(2,1))
hist(pval_PRegs_twosided_coef, breaks = 100)
hist(pval_PRegs_onesided_coef, breaks = 100)
##Paired regression joint
pval_PRegs_twosided_coef2 = rep(NA, ncol(dat$y))
pval_PRegs_onesided_coef2 = rep(NA, ncol(dat$y))

for(l in 1:ncol(dat$y)){
  print(l)
  ww = cbind(dat$y[,l], dat$x)
  perm_idx_X =  sapply(1:B, function(i) sample(1:args$n, args$n, replace = F))
  bb = matrix(NA, ncol = 2, nrow = B)
  for(b in 1:B){
    wwb  = cbind(ww, dat$x[perm_idx_X[,b]])
    tmp1 =lm(wwb~dat$Z+dat$Z[perm_idx_X[,b],,drop = D])
    wwb = tmp1$residuals
    bb[b,1] = sum(wwb[,1]*wwb[,2])/sum(wwb[,2]^2)
    bb[b,2] = sum(wwb[,1]*wwb[,3])/sum(wwb[,3]^2)
  }
  pval_PRegs_twosided_coef2[l] =  (sum(abs(bb[,2])>=abs(bb[,1]))+1)/(B+1)
  pval_PRegs_onesided_coef2[l] = 2*min((sum(bb[,2]>= bb[,1])+1)/(B+1),(sum(bb[,2] <= bb[,1])+1)/(B+1))
}






##############################




eps = dat$epsMat[,1]

perm_idx = sapply(1:M, function(i) sample(1:args$n, args$n, replace = F))
perm_idx_inverse = apply(perm_idx, 2, order)
ii = 1:args$n

ii_check = sapply(1:B, function(i) ii[perm_idx[,i]][perm_idx_inverse[,i]])
mean(apply(ii_check,2, function(i) sum(i == ii) == args$n))


x = dat$x
z= dat$Z
Zsvd = svd(z)
d = Zsvd$d
U = Zsvd$u[,  d>=1E-10] # Hz = U U^T

pos_counts = matrix(NA, ncol = M, nrow = B+1)
neg_counts = matrix(NA, ncol = M, nrow = B+1)
abs_counts = matrix(NA, ncol = M, nrow = B+1)


pos_counts_inv = matrix(NA, ncol = M, nrow = B+1)
neg_counts_inv = matrix(NA, ncol = M, nrow = B+1)

check_equal = rep(NA, M)
for(m in 1:M){
  print(m)
  ##need to generate
  perm_idx_X = sapply(1:B, function(i) sample(1:args$n, args$n, replace = F))
  perm_idx_X_full = cbind(rep(1:args$n), perm_idx_X)
  #get (B+1) * (B+1) beta list
  eps_perm = eps[perm_idx[,m]]
  eps_perp = sapply(1:(B+1), function(i){eps_perm - U[perm_idx_X_full[,i],]%*%(t(U[perm_idx_X_full[,i],])%*%eps_perm)})
  beta_hat = matrix(0, B+1, B+1)
  for(b in 1:(B+1)){
    beta_hat[b,] =x[perm_idx_X_full[,b]]%*%(eps_perp-U[perm_idx_X_full[,b],]%*%(t(U[perm_idx_X_full[,b],])%*%eps_perp))
  }
  beta_hat1 = t(beta_hat)
  pos_counts[,m] = apply(beta_hat >= beta_hat1,1,sum)
  neg_counts[,m] = apply(beta_hat <= beta_hat1,1,sum)
  abs_counts[,m] =  apply(abs(beta_hat) <= abs(beta_hat1),1,sum)
  
  #get (B+1) * (B+1) beta list, inverted ordering
  perm_idx_X_full_inv =apply(perm_idx_X_full,2,function(i) i[perm_idx_inverse[,m]])
  eps_perp_inv = sapply(1:(B+1), function(i){eps - U[perm_idx_X_full_inv[,i],]%*%(t(U[perm_idx_X_full_inv[,i],])%*%eps)})
  beta_hat_inv = matrix(0, B+1, B+1)
  for(b in 1:(B+1)){
    beta_hat_inv[b,] =(x[perm_idx_X_full_inv[,b]]%*%eps_perp_inv - x[perm_idx_X_full_inv[,b]]%*%U[perm_idx_X_full_inv[,b],]%*%(t(U[perm_idx_X_full_inv[,b],])%*%eps_perp_inv))
  }
  beta_hat1_inv = t(beta_hat_inv)
  #correct after examining
  check_equal[m] = all.equal(beta_hat,beta_hat_inv)
  print(paste0("iter",m,":",check_equal[m]))

}
####examine if the first row is exchangeable with others
o_pos = apply(pos_counts,2, rank, ties.method = "random")
o_neg = apply(neg_counts,2, rank, ties.method = "random")
w = sample(1:(B+1), M, replace = T)
par(mfrow = c(2,2))
hist((o_pos[1,]/(B+1)), breaks = 20)
hist((o_neg[1,]/(B+1)), breaks = 20)
hist((w/(B+1)), breaks = 20)
hist(1-(w/(B+1)), breaks = 20)
alpha = 0.05
mean((o_pos[1,]/(B+1))<=alpha)
mean((o_neg[1,]/(B+1))<=alpha)
####check the number of counts smaller than the threshold
alpha = 0.05
small_pos = apply(pos_counts <= alpha * (B+1),2,sum)
small_neg = apply(neg_counts <= alpha * (B+1),2,sum)
small_abs = apply(abs_counts <= alpha * (B+1),2,sum)
print(paste0("pos:", max(small_pos),"neg:", max(small_neg),"abs:", max(small_abs)))
print(paste0("pos:", max(small_pos/(B+1)),"neg:", max(small_neg/(B+1)),"abs:", max(small_abs/(B+1))))

####constructed pvals
alpha = 0.05
pval_pos = pos_counts[1,]/(B+1)
pval_neg = neg_counts[1,]/(B+1)
pval_abs = abs_counts[1,]/(B+1)
mean(pval_pos<=alpha)
mean(pval_neg<=alpha)
mean(pval_abs<=alpha)


##Nothing is wrong? 

## PRegs coef: still wrong intercept = F.?!
pval_PRegs_twosided = rep(NA, M)
pval_PRegs_onesided = rep(NA, M)
pval_PRegs_pos = rep(NA, M)
pval_PRegs_neg = rep(NA, M)
for(l in 1:M){
  print(l)
  tmp1coef <- PREG(y= eps[perm_idx[,l]], x=dat$x, z=dat$Z, add_intercept = F, B = B, statType = "partial")$res
  pval_PRegs_twosided[l] =  tmp1coef$unsigned
  pval_PRegs_pos[l] = tmp1coef$pos
  pval_PRegs_neg[l] = tmp1coef$neg
  pval_PRegs_onesided[l] = 2*min(tmp1coef$pos,tmp1coef$neg)
}
alpha = 0.05
mean(pval_PRegs_twosided<=alpha)
mean(pval_PRegs_pos<=alpha)
mean(pval_PRegs_neg<=alpha)
mean(pval_PRegs_onesided<=alpha)

alpha = 0.01
mean(pval_PRegs_twosided<=alpha)
mean(pval_PRegs_pos<=alpha)
mean(pval_PRegs_neg<=alpha)
mean(pval_PRegs_onesided<=alpha)

alpha = 0.001
mean(pval_PRegs_twosided<=alpha)
mean(pval_PRegs_pos<=alpha)
mean(pval_PRegs_neg<=alpha)
mean(pval_PRegs_onesided<=alpha)
########
pval_PRegs_twosided = rep(NA, M)
pval_PRegs_onesided = rep(NA, M)
pval_PRegs_pos = rep(NA, M)
pval_PRegs_neg = rep(NA, M)
for(l in 1:M){
  print(l)
  eps_perm = eps[perm_idx[,l]]
  out0 = perm_prepare(eps_perm, x, z, add_intercept=F)
  U = out0$U; yfitted = out0$yfitted; yresid = out0$yresid; Xresid = out0$Xresid
  perm_idx_X = sapply(1:B, function(i) sample(1:args$n, args$n, replace = F))
  perm_idxC = perm_idx_X-1
  out = permutation_conformal_C(Xresid=Xresid, y=eps_perm, U=U, perm_idx=perm_idxC, type = "coef")
  
  pval_PRegs_twosided[l] =  out$unsigned
  pval_PRegs_pos[l] = out$pos
  pval_PRegs_neg[l] = out$neg
  pval_PRegs_onesided[l] = 2*min(pval_PRegs_pos[l], pval_PRegs_neg[l])
}
alpha = 0.1
mean(pval_PRegs_twosided<=alpha)
mean(pval_PRegs_onesided<=alpha)
mean(pval_PRegs_pos<=alpha)
mean(pval_PRegs_neg<=alpha)

hist(pval_PRegs_twosided, breaks = 100)
########
pval_PRegs_twosided = rep(NA, M)
pval_PRegs_onesided = rep(NA, M)
pval_PRegs_pos = rep(NA, M)
pval_PRegs_neg = rep(NA, M)
for(l in 1:M){
  print(l)
  eps_perm = eps[perm_idx[,l]]
  perm_idx_X = sapply(1:B, function(i) sample(1:args$n, args$n, replace = F))
  perm_idx_X_full = cbind(c(1:args$n), perm_idx_X)
  out0 = perm_prepare(eps_perm, x, z, add_intercept=F)
  U = out0$U; yfitted = out0$yfitted; yresid = out0$yresid; Xresid = out0$Xresid
  eps_perp = sapply(1:(B+1), function(i){eps_perm - U[perm_idx_X_full[,i],]%*%(t(U[perm_idx_X_full[,i],])%*%eps_perm)})
  beta_hat1 = matrix(0, B+1, B+1)
  beta_hat1[1,] =t(Xresid[,1])%*%eps_perp
  for(b in 1:(B+1)){
  beta_hat1[b,1] =Xresid[,1][perm_idx_X_full[,b]]%*%eps_perp[,1]
  }
  pval_PRegs_twosided[l] =  sum(abs(beta_hat1[,1])>=beta_hat1[1,])/(B+1)
  pval_PRegs_pos[l] = sum(beta_hat1[,1] <= beta_hat1[1,])/(B+1)
  pval_PRegs_neg[l] = sum(beta_hat1[,1] <= beta_hat1[1,])/(B+1)
  pval_PRegs_onesided[l] = 2*min(pval_PRegs_pos[l], pval_PRegs_neg[l])
}
alpha = 0.01
mean(pval_PRegs_twosided<=alpha)
mean(pval_PRegs_pos<=alpha)
mean(pval_PRegs_neg<=alpha)
mean(pval_PRegs_onesided<=alpha)
hist(pval_PRegs_twosided,breaks = 100)
hist(pval_PRegs_pos,breaks = 100)
l=1
eps_perm = eps[perm_idx[,l]]
perm_idx_X = sapply(1:B, function(i) sample(1:args$n, args$n, replace = F))
perm_idx_X_full = cbind(c(1:args$n), perm_idx_X)

eps_perp = sapply(1:(B+1), function(i){eps_perm - U[perm_idx_X_full[,i],]%*%(t(U[perm_idx_X_full[,i],])%*%eps_perm)})

beta_hat = matrix(0, B+1, B+1)
beta_hat[1,] = t(x)%*%(eps_perp-U%*%(t(U)%*%eps_perp))
for(b in 1:(B+1)){
  beta_hat[b,1] =t(x[perm_idx_X_full[,b]])%*%(eps_perp[,b]-U%*%(t(U)%*%eps_perp[,b]))
}
(sum(beta_hat[1,] <= beta_hat[,1]))/(B+1)
(sum(beta_hat[1,] >= beta_hat[,1]))/(B+1)
(sum(abs(beta_hat[1,]) >= abs(beta_hat[,1])))/(B+1)

out0 = perm_prepare(eps_perm, x, z, add_intercept=F)
U = out0$U; yfitted = out0$yfitted; yresid = out0$yresid; Xresid = out0$Xresid

beta_hat1 = matrix(0, B+1, B+1)
beta_hat1[1,] =t(Xresid[,1])%*%eps_perp
for(b in 1:(B+1)){
  beta_hat1[b,1] =Xresid[,1][perm_idx_X_full[,b]]%*%eps_perp[,1]
}

all.equal(beta_hat1[1,],beta_hat[1,])
all.equal(beta_hat1[,1],beta_hat[,1])


perm_idxC = perm_idx_X-1
out = permutation_conformal_C(Xresid=Xresid, yfitted=rep(0, length(eps_perm)),yresid=eps_perp[,1], U=U, perm_idx=perm_idxC, type = "coef")
out$bhat[1:10,1]
beta_hat1[1,2:11]
beta_hat[1,2:11]

out$bperm[1:10,1]
beta_hat[2:11,1]

(sum(beta_hat[1,] <= beta_hat[,1]))/(B+1)
(sum(beta_hat[1,] >= beta_hat[,1]))/(B+1)
(sum(abs(beta_hat[1,]) >= abs(beta_hat[,1])))/(B+1)


(sum(beta_hat1[1,] <= beta_hat1[,1]))/(B+1)
(sum(beta_hat1[1,] >= beta_hat1[,1]))/(B+1)
(sum(abs(beta_hat1[1,]) >= abs(beta_hat1[,1])))/(B+1)

(sum(out$bhat[,1] <= out$bperm[,1])+1)/(B+1)
(sum(out$bhat[,1] >= out$bperm[,1])+1)/(B+1)
(sum(abs(out$bhat[,1]) <= abs(out$bperm[,1]))+1)/(B+1)

print(out$neg)
print(out$pos)
print(out$unsigned)


#########
all.equal(Xresid, (diag(rep(1,n))-U%*%t(U))%*%x)

b=2

print(t(x)%*%(diag(rep(1,n))-U%*%t(U))%*%(diag(rep(1,n))-U[perm_idx_X_full[,b],]%*%t(U[perm_idx_X_full[,b],]))%*%eps_perm)
all.equal(as.numeric((diag(rep(1,n))-U[perm_idx_X_full[,b],]%*%t(U[perm_idx_X_full[,b],]))%*%eps_perm),eps_perp[,b])
all.equal(as.numeric(eps_perm-U[perm_idx_X_full[,b],]%*%(t(U[perm_idx_X_full[,b],])%*%eps_perm)),eps_perp[,b])
print(t(Xresid[,1])%*%eps_perp[,b])
print(beta_hat[1,b]);
print(beta_hat1[1,b]); 

print(t(x[perm_idx_X_full[,b]])%*%(diag(rep(1,n))-U[perm_idx_X_full[,b],]%*%t(U[perm_idx_X_full[,b],]))%*%(diag(rep(1,n))-U%*%t(U))%*%eps_perm)
print(beta_hat[b,1])
print(beta_hat1[b,1])
t(Xresid[,1][perm_idx_X_full[,b]])%*%eps_perp[,1]

all.equal(eps_perp[,1], eps_perm - (diag(rep(1,n))-U%*%t(U)))
all.equal(as.numeric((diag(rep(1,n))-U[perm_idx_X_full[,b],]%*%t(U[perm_idx_X_full[,b],]))%*%eps_perm),eps_perp[,b])


