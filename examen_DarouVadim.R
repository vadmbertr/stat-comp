## --------------------------------------------
library(class)
library(dplyr)
library(ggpattern)
library(ggplot2)
library(kableExtra)
library(latex2exp)
library(patchwork)
library(pROC)
library(tibble)
library(tidyr)


## --------------------------------------------
set.seed(15)


## ---- Exercice1------------------------------
g <- function (x) { # définition de g(x)
  exp(-abs(x)) / 2
}


## ---- 1.1------------------------------------


## ----plot-1-g-density, fig.cap = "Aperçu de la fonction $g(x) = \\frac{1}{2} exp(-|x|)$, pour $x \\in [-5, 5]$"----
ggplot() +
  xlim(-5, 5) +
  geom_function(fun = g) +
  labs(x = "x", y = "y")


## ---- 1.2------------------------------------
rg <- function (n) { # retourne les réalisations de g
  Un <- runif(n, min = 0, max = 1) # tirage selon la loi U[0,1]
  Gn <- sapply(Un, function (u) { # tirage selon g via sa fonction quantile
    if (u > 0.5) -log(2*(1-u)) else log(2*u)
  })
}


## ---- 1.3------------------------------------


## ----plot-1-g-densities, fig.cap = "Comparaison de la densité théorique (en noir) et de la densité empirique (en bleu) de g"----
ggplot(data.frame(x = rg(1000)), aes(x)) +
  xlim(-5, 5) +
  geom_density(aes(color = "g_hat")) +
  geom_function(fun = g, aes(color = "g(x)")) +
  labs(x = "x", y = "Densité", color = "") +
  scale_color_manual(values = c("g(x)" = "black", "g_hat" = "blue"),
                     labels = unname(TeX(c("g(x)", "$\\hat{g}$"))))


## ---- 1.4------------------------------------
M <- sqrt(2 * exp(1) / pi)


## ----plot-1-maj, fig.cap = "Représentation de $\\frac{f(x)}{g(x)}$"----
ggplot() +
  geom_function(fun = function (x) {dnorm(x) / g(x)}, xlim = c(-5, 5)) +
  geom_hline(yintercept = M, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(x = "x", y = "y")


## ----plot-1-rej, fig.cap = "Illustration graphique de la procédure de rejet"----
ggplot() +
  xlim(-5, 5) +
  geom_function(fun = function (x) {g(x) * M}, aes(color = "M*g(x)")) +
  geom_function(fun = dnorm, aes(color = "f(x)")) +
  labs(x = "x", y = "y", color = "") +
  scale_color_manual(values = c("f(x)" = "black", "M*g(x)" = "blue"))


## ---- 1.5------------------------------------
rf <- function (n) { # retourne les réalisations de f et le taux de rejet
  Xn <- NULL
  rejected <- 0
  while (length(Xn) < n) { # on veut n réalisations
    X0 <- rg(n - length(Xn)) # génération du nombre de réalisations manquantes selon g
    U0 <- sapply(X0, function (x) runif(1, 0, M*g(x))) # tirages dans U[0, g(x0)]
    idx <- U0 <= dnorm(X0) # f=dnorm
    Xn <- c(Xn, X0[idx]) # on conserve les x0 inférieurs ou égaux à f(x0)
    rejected <- rejected + sum(1-idx) # maj du nombre de rejets
  }
  list(Xn = Xn, rate = rejected/(rejected+n))
}


## --------------------------------------------
rf.res <- rf(1000)


## ----plot-1-f-densities, fig.cap = "Comparaison de la densité théorique (en noir) et de la densité empirique (en bleu) de f"----
ggplot(data.frame(x = rf.res$Xn), aes(x)) +
  xlim(-5, 5) +
  geom_density(aes(color = "f_hat")) +
  geom_function(fun = dnorm, aes(color = "f(x)")) +
  labs(x = "x", y = "Densité", color = "") +
  scale_color_manual(values = c("f(x)" = "black", "f_hat" = "blue"),
                     labels = unname(TeX(c("f(x)", "$\\hat{f}$"))))


## ---- Exercice2------------------------------
M <- 10000
alpha <- .05
treshold <- qnorm(1-alpha) # one-sided test
P1 <- .95
P0 <- c(0.8, 0.85, 0.9, 0.92, 0.93)
N <- seq(50, 500, by = 50)


## ----plot-2-region, fig.cap = "Représentation de la région de rejet du test (en orange)"----
ggplot() +
  xlim(-5, 5) +
  geom_function(fun = dnorm) +
  stat_function(fun = dnorm, xlim = c(treshold, 5), fill = "orange", geom = "area") +
  geom_vline(xintercept = treshold, color = "orange", linetype = "dashed") +
  labs(x = "x", y = "Densité")


## --------------------------------------------
p.sd <- function (p, n) { # écart-type de p sous H0
  sqrt(p * (1 - p) / n)
}
p1.r <- function (m, n) { # génère m taux à partir de n-échantillons sous H1
  apply(matrix(data = rbinom(m*n, 1, P1), nrow = m), 1, mean)
}
p.T <- function (p, p0, n) { # statistique de test
  (p - p0) / p.sd(p0, n)
}

p.MC <- function (p0, n) { # procédure MC
  pi <- p1.r(M, n) # M taux sous H1
  ti <- p.T(pi, p0, n) # stat de test pour les M taux
  mean(ti > treshold) # puissance du test
}


## --------------------------------------------
mat <- do.call(cbind, lapply(N, function (n) { # pour chaque n
  p <- sapply(P0, function (p0) p.MC(p0, n)) # et pour chaque p0, on calcule la puissance du test
  names(p) <- P0
  p
}))


## ----plot-2-puissances, fig.cap = "Evolution de la puissance du test selon la taille de l'échantillon, pour différentes valeurs de $p_0$"----
df <- data.frame(data = mat, row.names = P0)
colnames(df) <- N
tbl <- as_tibble(df)
tbl$P0 <- as.factor(P0)
tbl <- tbl %>% pivot_longer(!P0, names_to = "N", values_to = "pow")
tbl$N <- as.numeric(tbl$N)

ggplot(data = tbl) +
  geom_line(aes(x = N, y = pow, color = P0)) +
  labs(x = "Taille de l'échantillon", y = "Puissance du test", color = "p0")


## ---- Exercice3------------------------------
load("exo3.Rdata")


## ---- 3.2------------------------------------
# courbes ROC
roc.1 <- roc(y, prob1, quiet = T)
roc.2 <- roc(y, prob2, quiet = T)


## --------------------------------------------
auc.comp <- function (cases, controls) {
  sapply(cases, function (c) {
    row <- as.numeric(controls < c)
    row[controls == c] <- 1/2 # les égalités sont réparties
    row
  })
}
auc <- function (cases, controls) { # calcul exact de l'AUC
  comp <- auc.comp(cases, controls)
  mean(comp)
}

cases.idx <- which(y == 1)
controls.idx <- which(y == -1)
cases.1 <- prob1[cases.idx]
cases.2 <- prob2[cases.idx]
controls.1 <- prob1[controls.idx]
controls.2 <- prob2[controls.idx]
auc.1 <- auc(cases.1, controls.1)
auc.2 <- auc(cases.2, controls.2)


## ----plot-3-roc, fig.cap = "Courbes ROC des deux classifieurs"----
ggroc(list(roc.1 = roc.1, roc.2 = roc.2)) +
  labs(x = "Spécificité", y = "Sensibilité", color = "Classifieur") +
  scale_color_discrete(labels = c(1, 2))


## ---- 3.3------------------------------------
B <- 1000

auc.b <- function (scores) {
  replicate(B, {
    cases <- sample(scores[cases.idx], replace = T)
    controls <- sample(scores[controls.idx], replace = T)
    auc(cases, controls)
  })
}

auc.df <- data.frame(c1 = auc.b(prob1), c2 = auc.b(prob2))
colnames(auc.df) <- 1:2


## --------------------------------------------
alpha <- .05

auc.conf.int.per <- do.call(
  cbind,
  apply(auc.df, 2, function (auc) {
    list(inf = quantile(auc, alpha/2)[[1]],
         sup = quantile(auc, 1-alpha/2)[[1]])
}))

auc.conf.int.piv <- do.call(
  cbind,
  apply(auc.df, 2, function (auc) {
    auc.hat <- mean(auc)
    list(inf = 2*auc.hat - quantile(auc, 1-alpha/2)[[1]],
         sup = 2*auc.hat - quantile(auc, alpha/2)[[1]])
}))


## ----plot-3-confint, fig.cap = "Distribution des AUC obtenus par bootstrap et les intervalles de confiance correspondants"----
auc.tbl <- as_tibble(t(auc.df))
auc.tbl$est <- as.factor(1:2)
auc.tbl <- auc.tbl %>% pivot_longer(!est, values_to = "auc")

auc.conf.int.per.df <- data.frame(data = t(auc.conf.int.per), row.names = 1:2)
colnames(auc.conf.int.per.df) <- c("inf", "sup")
auc.conf.int.per.tbl <- as_tibble(auc.conf.int.per.df)
auc.conf.int.per.tbl$est <- as.factor(1:2)
auc.conf.int.per.tbl <- auc.conf.int.per.tbl %>% pivot_longer(!est, values_to = "borne")
auc.conf.int.per.tbl$method <- "percentiles"
auc.conf.int.per.tbl$borne <- as.numeric(auc.conf.int.per.tbl$borne)

auc.conf.int.piv.df <- data.frame(data = t(auc.conf.int.piv), row.names = 1:2)
colnames(auc.conf.int.piv.df) <- c("inf", "sup")
auc.conf.int.piv.tbl <- as_tibble(auc.conf.int.piv.df)
auc.conf.int.piv.tbl$est <- as.factor(1:2)
auc.conf.int.piv.tbl <- auc.conf.int.piv.tbl %>% pivot_longer(!est, values_to = "borne")
auc.conf.int.piv.tbl$method <- "pivot"
auc.conf.int.piv.tbl$borne <- as.numeric(auc.conf.int.piv.tbl$borne)

auc.conf.int.tbl <- rbind(auc.conf.int.piv.tbl, auc.conf.int.per.tbl)

x.lim <- c(min(auc.tbl$auc), max(auc.tbl$auc))

gh <- ggplot(data = auc.tbl, aes(x = auc, fill = est)) +
  geom_histogram(alpha = 0.6, position = "identity") +
  geom_vline(data = auc.conf.int.tbl, aes(xintercept = borne, color = est, linetype = method)) +
  labs(x = "AUC", y = "Effectif", fill = "Classifieur", color = "Classifieur", linetype = "Intervalle") +
  xlim(x.lim) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 1, 0, "pt"))
gbp <- ggplot(data = auc.tbl, aes(x = auc, y = est, fill = est)) +
  geom_boxplot() +
  scale_fill_discrete(guide = "none") +
  labs(x = "AUC", y = "Classifieur") +
  xlim(x.lim) +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))

layout <- "
AAA
AAA
AAA
BBB
"
gh / gbp +
  plot_layout(design = layout, guides = "collect")


## ---- 3.4------------------------------------
P <- 1000

auc.s <- function (comp1, comp2) { # statistique de test
  return(abs(mean(comp1) - mean(comp2)))
}

auc.comp.1 <- sample(as.numeric(auc.comp(cases.1, controls.1))) # comparaisons "brutes" positifs / négatifs
auc.comp.2 <- sample(as.numeric(auc.comp(cases.2, controls.2)))
s0 <- auc.s(auc.comp.1, auc.comp.2)
auc.comp.Z <- c(auc.comp.1, auc.comp.2)
auc.diff.perm <- replicate(P, {
  auc.comp.i <- sample(length(auc.comp.Z), length(auc.comp.1)) # permutations des comparaisons "brutes"
  auc.comp.x <- auc.comp.Z[auc.comp.i]
  auc.comp.y <- auc.comp.Z[-auc.comp.i]
  auc.s(auc.comp.x, auc.comp.y)
})

auc.diff.perm <- c(auc.diff.perm, s0)
p.val.perm <- mean(auc.diff.perm >= s0)
p.val.boot <- roc.test(roc.1, roc.2, method = "bootstrap")$p.value


## ----auc-test--------------------------------
auc.test.df <- data.frame(permutation = p.val.perm, bootstrap = p.val.boot)
rownames(auc.test.df) <- "p-valeur"
kable(auc.test.df, caption = "p-valeurs des tests d'égalité des AUC",
      label="auc-test") %>%
  kable_styling(latex_options = "HOLD_position")


## ---- Exercice4------------------------------
load("exo4.Rdata")


## ----4.1-------------------------------------


## ----plot-4-data, fig.height = 6, fig.cap = "Aperçu des jeux d'entrainement et de test"----
train.df <- data.frame(X.train)
train.df$cat <- y.train
train.df$set <- "train"
test.df <- data.frame(X.test)
test.df$cat <- y.test
test.df$set <- "test"
df <- rbind(train.df, test.df)
df$cat <- as.factor(df$cat)

x1.lim <- c(min(df$V1), max(df$V1))
x2.lim <- c(min(df$V2), max(df$V2))

gp <- ggplot(data = df) +
  geom_point(aes(x = V1, y = V2, col = cat, shape = set)) +
  scale_shape_manual(values = c(3, 19)) +
  labs(color = "Catégorie", shape = "Jeu") +
  xlim(x1.lim) +
  ylim(x2.lim) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        plot.margin = margin(3, 3, 1, 1, "pt"))
gh11 <- ggplot(data = df[df$set == "test", ], aes(x = V1, fill = cat)) +
  geom_histogram_pattern(alpha = 0.6, position = "identity", pattern = "crosshatch", pattern_density = .0001,
                         pattern_alpha = .6) +
  scale_pattern_discrete(guide = "none") +
  scale_fill_discrete(guide = "none") +
  xlim(x1.lim) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        plot.margin = margin(0, 0, 0, 0, "pt"))
gh21 <- ggplot(data = df[df$set == "train", ], aes(y = V2, fill = cat)) +
  geom_histogram_pattern(alpha = 0.6, position = "identity", pattern = "circle", pattern_fill = "black") +
  scale_pattern_discrete(guide = "none") +
  scale_fill_discrete(guide = "none") +
  scale_x_reverse() +
  labs(y = "X2") +
  ylim(x2.lim) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"))
gh12 <- ggplot(data = df[df$set == "train", ], aes(x = V1, fill = cat)) +
  geom_histogram_pattern(alpha = 0.6, position = "identity", pattern = "circle", pattern_fill = "black") +
  scale_pattern_discrete(guide = "none") +
  scale_fill_discrete(guide = "none") +
  scale_y_reverse() +
  labs(x = "X1") +
  xlim(x1.lim) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"))
gh22 <- ggplot(data = df[df$set == "test", ], aes(y = V2, fill = cat)) +
  geom_histogram_pattern(alpha = 0.6, position = "identity", pattern = "crosshatch", pattern_density = .0001,
                         pattern_alpha = .6) +
  scale_pattern_discrete(guide = "none") +
  scale_fill_discrete(guide = "none") +
  ylim(x2.lim) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        plot.margin = margin(0, 0, 0, 0, "pt"))

layout <- "
#DDD#
ACCCB
ACCCB
ACCCB
#EEE#
"
gh21 + gh22 + gp + gh11 + gh12 +
  plot_layout(design = layout, guides = "collect")


## ----4.2-------------------------------------
K <- c(1, 3, 5, 7, 9, 11)
rate <- sapply(K, function (k) {
  pred.test <- knn(X.train, X.test, y.train, k)
  mean(pred.test == y.test)
})


## ----tab-4-ppv-------------------------------
knn.df <- t(data.frame(k = as.character(as.integer(K)), rate = rate))
rownames(knn.df) <- c("k", "taux")
kable(knn.df, caption = "Taux de bonne classification par k-ppv pour différents k",
      label="tab-4-ppv") %>%
  kable_styling(latex_options = "HOLD_position")


## ----4.3-------------------------------------
B <- 100

pred.b <- replicate(B, { # procédure de bagging
  idx <- sample(length(y.train), replace = T)
  sapply(K, function (k) {
    as.numeric(knn(X.train[idx,], X.test, y.train[idx], k)) - 1
  })
})
rate.b <- apply(pred.b, 2, function (p) { # résultat du bagging
  p.maj <- as.numeric(apply(p, 1, mean) >= .5 )
  mean(p.maj == y.test)
})
rate.ind <- apply(pred.b, 2, function (p) { # résultats individuels des classifieurs de bagging
  apply(p == y.test, 2, mean)
})


## ----tab-4-bagging---------------------------
bagging.df <- t(data.frame(k = as.character(as.integer(K)), rate = rate.b))
rownames(bagging.df) <- c("k", "taux")
kable(bagging.df, caption = "Taux de bonne classification par bagging pour différents k",
      label="tab-4-bagging") %>%
  kable_styling(latex_options = "HOLD_position")


## ----4.4-------------------------------------


## ----plot-4-all, fig.cap = "Taux de bonne classification des $k$-ppv, avec (en rouge) et sans (en bleu) bagging"----
rate.df <- data.frame(knn = rate, bagging = rate.b)
rate.tbl <- as_tibble(rate.df)
rate.tbl$K <- K
rate.tbl <- rate.tbl %>% pivot_longer(!K, names_to = "method", values_to = "rate")
rate.tbl$K <- factor(rate.tbl$K, levels = K)

rate.ind.df <- data.frame(t(rate.ind))
colnames(rate.ind.df) <- as.factor(1:B)
rate.ind.tbl <- as_tibble(rate.ind.df)
rate.ind.tbl$K <- K
rate.ind.tbl <- rate.ind.tbl %>% pivot_longer(!K, names_to = "B", values_to = "rate")
rate.ind.tbl$K <- factor(rate.ind.tbl$K, levels = K)

ggplot() +
  geom_boxplot(data = rate.ind.tbl, aes(x = K, y = rate), width = .25, alpha = 0,
               color = paste0("#F8766D", toupper(as.hexmode(round(.5*255))))) +
  geom_line(data = rate.ind.tbl, aes(x = K, y = rate, group = B, col = "bagging"), alpha = .1) +
  geom_line(data = rate.tbl, aes(x = K, y = rate, group = method, col = method)) +
  scale_x_discrete(breaks = K, labels = K) +
  scale_color_discrete(labels = c("bagging", "k-ppv")) +
  labs(x = "k", y = "Taux de bonne classification", col = "Méthode")

