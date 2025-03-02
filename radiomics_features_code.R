#### Pyradiomics ####

#### Radiomics Features ####
dt <- read.csv("radiomics_features-right_pallidum.csv", sep = ",")

str(dt)
table(dt$Group)

#### Lasso ####
library(glmnet)
cv_x <- as.matrix(dt[,c(2:1133)])
cv_y <- data.matrix(dt[,1])
set.seed(123)
# tow pictures for lasso
nocv_lasso <- glmnet(
  x = cv_x, y = cv_y,
  family = "binomial", alpha = 1,
)
par(font.lab = 4, mfrow = c(2, 1), mar = c(4.5, 5, 3, 2))
plot(nocv_lasso, xvar = "lambda", las = 1, lwd = 2, xlab = "log(lambda)") 
abline(v = log(nocv_lasso$lambda.min), lwd = 1, lty = 3, col = "black")

fitcv <- cv.glmnet(
  x = cv_x,
  y = cv_y, type.measure = "deviance", 
  family = "binomial", alpha = 1, nfolds = 1000
) # cross validation
fitcv
plot(x = fitcv, las = 1, xlab = "log(lambda)") 
coef(fitcv, s = "lambda.min") # min、1se

library(broom) # CRAN v1.0.5
tidy_df <- broom::tidy(nocv_lasso)
tidy_cvdf <- broom::tidy(fitcv)
head(tidy_cvdf)

#### Heatmap ####
library(corrplot)
library(ggcorrplot) # CRAN v0.1.4.1
colnames(dt) <- wrap_column_names(colnames(dt), width = 10)  

env.cor <- round(cor(as.matrix(dt), method = "pearson"), 3) 

env.p <- round(cor_pmat(as.matrix(dt), method = "pearson"), 3) 
env.p

env.uppCI <- round(cor.mtest(as.matrix(dt), conf.level = 0.95)$uppCI, 3)
env.uppCI 
env.lowCI <- round(cor.mtest(as.matrix(dt), conf.level = 0.95)$lowCI, 3)
env.lowCI 

m <- par(no.readonly = TRUE) 
pdf("Heatmap.pdf", width = 20, height = 20) 

par(mfrow = c(1, 1)) 

cor.plot <- corrplot(
  corr = env.cor, p.mat = env.p, type = "upper",
  tl.pos = "lt", tl.col = "black", tl.srt = 45, tl.offset = 1.0,
  insig = "label_sig", sig.level = c(.01, .05), number.cex = 1.8,
  pch.cex = 1.2, pch.col = "black", tl.cex = 1.6 # order = "AOE"
) 
cor.plot <- corrplot(
  corr = env.cor, type = "lower", add = TRUE, method = "number",
  tl.pos = "n", tl.col = "black", tl.cex = 0.82, lab_size = 2, 
  col = "black", diag = FALSE, cl.pos = "n", pch.col = "black",
  number.cex = 1.0 # order = "AOE"
)
dev.off()
par(m)

#### treeheatr ####
library(treeheatr)
treeheatr(dt, target_lab = 'Group')
