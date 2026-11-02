############################################################
# TP2
############################################################

# Chargement des packages

library(actuar)
library(ggplot2)
library(dplyr)

## === EXERCICE 2.1 ========================================

# --- Paramètres ---
lambda <- 5
q <- 0.3
n_sim <- 1e4

# --- Simulation empirique ---
set.seed(123)
N <- rpois(n_sim, lambda)
S_sim <- sapply(N, function(n) sum(rbinom(n, 1, q)))

# --- Panjer manuel ---

# Paramètre
k_max <- 2  # calcul de P(S=0), P(S=1), P(S=2)

# Initialisation
p <- numeric(k_max + 1)
f <- c(1-q, q)  # Indices R : f[1] = P(X=0), f[2] = P(X=1)

# p0
p[1] <- exp(-lambda * q)  # P(S=0)

# Récursion de Panjer
for (k in 1:k_max) {
  p[k+1] <- (lambda * q / k) * p[k]  # formule de Panjer pour Poisson - Bernoulli
}

names(p) <- 0:k_max
p

# --- Panjer avec actuar ---
sev_probs <- c(1 - q, q)

Fs <- aggregateDist("recursive",
                    model.freq = "poisson", lambda = lambda,
                    model.sev = sev_probs,
                    x.scale = 1)

# --- Comparaisons graphiques ---
par(mfrow = c(1, 2))

# Densité (loi empirique vs Panjer)
hist(S_sim, freq = FALSE, breaks = seq(-0.5, max(S_sim)+0.5, 1),
     main = "Loi empirique vs théorique", xlab = "S")
points(0:8, diff(Fs(-1:8)), col = "red", pch = 19)
legend("topright", legend = c("Empirique", "Panjer"),
       col = c("grey", "red"), pch = 19)

# Fonction de répartition (empirique vs Panjer)
plot(ecdf(S_sim), do.points = FALSE, main = "Fonction de répartition",
     ylab = "F(x)", xlab = "x")
lines(0:8, Fs(0:8), type = "s", col = "red", lwd = 2)

# Revenir à un graphique par fenêtre pour les exercices suivants
par(mfrow = c(1, 1))

# --- Probabilités demandées ---
P0 <- Fs(0)
P1 <- Fs(1) - Fs(0)
P2 <- Fs(2) - Fs(1)
cat("P[S=0] =", P0, "\nP[S=1] =", P1, "\nP[S=2] =", P2, "\n")




## === EXERCICE 2.2 ========================================

# --- Paramètres ---
n <- 10
p <- 0.1
alpha <- 1/2
tau <- 1/3
set.seed(123)

# --- Simulation empirique ---
n_sim <- 1e5
N <- rbinom(n_sim, n, p)
S_sim <- sapply(N, function(m) sum(rgamma(m, shape = alpha, rate = tau)))

# --- Fonctions paramétrées ---
Fx_gamma  <- function(x) pgamma(x, shape = alpha, rate = tau)
lev_gamma <- function(x) levgamma(x, shape = alpha, rate = tau)

# --- Discrétisation des coûts gamma ---
step <- 0.05

# Borne alignée sur la grille : la probabilité au-delà est inférieure à 1e-10.
q_hi <- ceiling(
  qgamma(1 - 1e-10, shape = alpha, rate = tau) / step
) * step

fx_discr <- discretize(
  Fx_gamma,
  from = 0, to = q_hi, step = step,
  method = "unbiased", lev = lev_gamma
)

# Renormalisation après une troncature très faible.
# On conserve la masse en zéro produite par la discrétisation.
fx_discr <- fx_discr / sum(fx_discr)

# --- Loi de S_h : somme des coûts discrétisés ---
Fs <- aggregateDist(
  "recursive",
  model.freq = "binomial", size = n, prob = p,
  model.sev = fx_discr,
  x.scale = step,
  tol = 1e-10,
  maxit = 10000
)

# --- Comparaison des fonctions de répartition ---
plot(
  ecdf(S_sim), do.points = FALSE,
  main = "Fonctions de répartition : simulation et Panjer",
  xlab = "s", ylab = "Probabilité cumulée"
)

plot(
  Fs, add = TRUE,
  do.points = FALSE, verticals = TRUE,
  col = "red", lwd = 2
)

legend(
  "bottomright",
  legend = c("S : simulation", "S_h : Panjer"),
  col = c("black", "red"), lwd = 2
)

# --- Probabilités ponctuelles : distinguer S et S_h ---
# Pour S, les coûts gamma sont strictement positifs presque sûrement :
# P(S=0) = P(N=0), et P(S=x) = 0 pour tout x > 0.
# Pour S_h, les coûts discrétisés ont des masses ponctuelles, dont une en zéro.
# P(S_h=0) = (1 - p + p * fx_discr[1])^n peut donc différer de P(S=0).

# diff(Fs) fournit les masses aux points 0, step, 2*step, ...
mass_S_h <- diff(Fs)
valeurs <- c(0, 1, 2)
indices <- round(valeurs / step) + 1

res_point <- data.frame(
  Valeur = valeurs,
  "P(S = x), modèle initial" = c((1 - p)^n, 0, 0),
  "P(S_h = x), modèle discrétisé" = mass_S_h[indices],
  check.names = FALSE
)

cat("\nProbabilités ponctuelles : modèles initial et discrétisé\n")
print(res_point, row.names = FALSE)

# --- Complément : comparaison sur les mêmes intervalles ---
# Le premier intervalle inclut bien la masse en zéro.
cdf_aux_bornes <- cumsum(mass_S_h)[round((1:3) / step) + 1]

res_interv <- data.frame(
  Intervalle = c("[0,1]", "(1,2]", "(2,3]"),
  Panjer = c(cdf_aux_bornes[1], diff(cdf_aux_bornes)),
  Empirique = c(
    mean(S_sim <= 1),
    mean(S_sim > 1 & S_sim <= 2),
    mean(S_sim > 2 & S_sim <= 3)
  )
)

cat("\nProbabilités sur les mêmes intervalles :\n")
print(res_interv, row.names = FALSE)


## === EXERCICE 2.3 ========================================

# --- Loi exacte (par convolution discrète) ---
pS1 <- c(0.80, 0.15, 0.05)
pS2 <- c(0.90, 0.00, 0.10)
pS3 <- c(0.95, 0.05, 0.00)
conv <- function(p, q) stats::convolve(p, rev(q), type = "open")
pS <- conv(conv(pS1, pS2), pS3)
pS <- pS / sum(pS)
valeurs_S <- 0:(length(pS) - 1)
cdf_S <- cumsum(pS)

VaR_exact <- valeurs_S[min(which(cdf_S >= 0.95))]
E_stoploss_exact <- sum(pmax(valeurs_S - 2, 0) * pS)

cat("\n=== Valeurs exactes ===\n")
cat(sprintf("VaR[S; 0.95]    = %d\n", VaR_exact))
cat(sprintf("E[(S - 2)+]     = %.5f\n", E_stoploss_exact))

# --- Complément : approximation collective ---
set.seed(123)
Nsim <- 1e5
S_coll <- numeric(Nsim)

p_i <- c(0.80, 0.90, 0.95)

# --- Paramètres Poisson (approximation collective conservatrice) ---
lambda <- -log(p_i)
# lambda <- 1-p_i  # option plus précise pour la prime stop-loss (voir remarque slide 61 chapitre 5)
vals <- list(S1=c(1,2), S2=2, S3=1)
probs <- list(S1=c(0.15,0.05)/0.2, S2=1, S3=1)

# --- Simulation collective ---
for(i in 1:3){
  Ni <- rpois(Nsim, lambda[i])
  Xi <- numeric(Nsim)
  
  idx <- which(Ni >= 1)
  if(length(idx) > 0){
    if(length(vals[[i]]) == 1){
      Xi[idx] <- Ni[idx] * vals[[i]]
    } else {
      Xi[idx] <- sapply(Ni[idx], function(n)
        sum(sample(vals[[i]], n, replace=TRUE, prob=probs[[i]]))
      )
    }
  }
  S_coll <- S_coll + Xi
}

# --- Mesures de risque (approx. collective) ---
# type = 1 : inverse de la fonction de répartition empirique, comme au cours
VaR_coll <- quantile(S_coll, probs = 0.95, type = 1, names = FALSE)
E_stoploss_coll <- mean(pmax(S_coll - 2, 0))

cat("=== Approximation collective ===\n")
cat(sprintf("VaR[S; 0.95]    = %.3f\n", VaR_coll))
cat(sprintf("E[(S - 2)+]     = %.5f\n", E_stoploss_coll))

# --- Graphique comparatif ---
plot(ecdf(S_coll),
     main = "CDF : approximation collective (simulation) vs exacte",
     xlab = "S", ylab = "F(S)", col = "red", lwd = 2,
     do.points = FALSE)
lines(valeurs_S, cdf_S, type = "s", col = "black", lwd = 2)
legend("bottomright",
       legend = c("CDF exacte (convolution)", "CDF collective (simulation)"),
       col = c("black", "red"),
       lty = c(1, 1), lwd = 2)


## === EXERCICE 2.4 ========================================

# Paramètres
alpha <- 0.99       # niveau de risque
theta <- 1         # paramètre de l’exponentielle
n <- 1e6            # taille de l'échantillon simulé

# ---------------------------
# 1. Valeurs théoriques
# ---------------------------
VaR_th <- -log(1 - alpha) / theta
TVaR_th <- VaR_th + 1/theta   # pour l’exponentielle, TVaR = VaR + 1/theta
CTE_th <- TVaR_th            # égalité pour ce modèle exponentiel continu

# ---------------------------
# 2. Simulation Monte-Carlo
# ---------------------------
set.seed(123)
X <- rexp(n, rate = theta)

# VaR empirique (quantile empirique)
VaR_emp <- quantile(X, probs = alpha, type = 1, names = FALSE)

# Estimateur Monte-Carlo commun de la TVaR et de la CTE,
# qui sont égales pour le modèle exponentiel continu étudié.
# Cette égalité ne se généralise pas automatiquement aux lois discrètes.
TVaR_emp <- mean(X[X > VaR_emp])
CTE_emp <- TVaR_emp

# ---------------------------
# 3. Comparaison
# ---------------------------
results <- data.frame(
  Theorique = c(VaR_th, TVaR_th, CTE_th),
  Empirique = c(VaR_emp, TVaR_emp, CTE_emp),
  row.names = c("VaR", "TVaR", "CTE")
)

print("Comparaison entre valeurs théoriques et empiriques :")
print(results)

# ---------------------------
# 4. Visualisation
# ---------------------------
hist(X, breaks = 100, probability = TRUE, col = "lightblue",
     main = "Loi exponentielle : VaR, TVaR et CTE",
     xlab = "Valeur")
curve(dexp(x, rate = theta), col = "red", lwd = 2, add = TRUE)

# Ajout de lignes pour la VaR et la valeur commune TVaR = CTE
abline(v = VaR_th, col = "blue", lwd = 2, lty = 2)
abline(v = TVaR_th, col = "darkgreen", lwd = 2, lty = 2)

legend("topright",
       legend = c("Densité théorique", "VaR théorique", "TVaR = CTE théoriques"),
       col = c("red", "blue", "darkgreen"), lty = c(1, 2, 2), lwd = 2, bty = "n")


###################################

### Convergence VaR, TVaR et CTE pour une exponentielle ###

library(ggplot2)
library(dplyr)

# Paramètres
alpha <- 0.99
theta <- 1
N <- 1e5  # taille max de l'échantillon

# ---------------------------
# 1. Valeurs théoriques
# ---------------------------
VaR_th <- -log(1 - alpha) / theta
TVaR_th <- VaR_th + 1/theta
CTE_th <- TVaR_th

# ---------------------------
# 2. Simulation
# ---------------------------
set.seed(123)
X <- rexp(N, rate = theta)

# Tailles d’échantillons intermédiaires
n_seq <- unique(floor(exp(seq(log(100), log(N), length.out = 50))))

# Fonction pour calculer les estimations Monte-Carlo
# TVaR et CTE : même estimateur du risque exponentiel continu.
compute_risks <- function(n) {
  x <- X[1:n]
  VaR_emp <- quantile(x, probs = alpha, type = 1, names = FALSE)
  TVaR_emp <- mean(x[x > VaR_emp])
  data.frame(n = n, VaR = as.numeric(VaR_emp), TVaR = TVaR_emp)
}

# Calcul pour toutes les tailles
res <- do.call(rbind, lapply(n_seq, compute_risks))

# Passage au format long pour ggplot
res_long <- res %>%
  tidyr::pivot_longer(cols = c("VaR", "TVaR"), names_to = "Mesure", values_to = "Valeur")
res_long$Mesure[res_long$Mesure == "TVaR"] <- "TVaR = CTE"

# ---------------------------
# 3. Visualisation ggplot
# ---------------------------
ggplot(res_long, aes(x = n, y = Valeur, color = Mesure)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = VaR_th, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = TVaR_th, color = "darkgreen", linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c("VaR" = "blue", "TVaR = CTE" = "darkgreen")) +
  labs(title = "Convergence de la VaR et de la TVaR/CTE (exponentielle)",
       x = "Taille de l'échantillon (log)", y = "Valeur",
       caption = paste("Paramètre theta =", theta, " - niveau alpha =", alpha)) +
  theme_minimal(base_size = 14)


## === EXERCICE 2.5 ========================================
library(ggplot2)
set.seed(123)

# ---- Simulation ----
n <- 1e6
X1 <- rnorm(n)
Y1 <- rnorm(n)   # indépendant
V  <- sample(c(-1,1), n, replace = TRUE)
X2 <- X1
Y2 <- V * X1

# ---- Sommes ----
S1 <- X1 + Y1
S2 <- X2 + Y2

# ---- Corrélations ----
# En théorie, les quatre marginales sont N(0,1) et les deux corrélations sont nulles.
# cor() fournit ici leurs estimations sur les échantillons simulés.
rho1 <- cor(X1, Y1)
rho2 <- cor(X2, Y2)

cat("Corrélations :\n")
cat(sprintf("rho(X1, Y1) = %.3f\n", rho1))
cat(sprintf("rho(X2, Y2) = %.3f\n", rho2))

# ---- Niveaux de probabilité ----
p <- seq(0.8, 0.999, by = 0.001)

# ---- VaR théoriques ----
# S2 = (1 + V) * X1 : masse 1/2 en zéro et composante N(0,4) de poids 1/2.
# La CDF saute de 0.25 à 0.75 en zéro, d'où le plateau de la fonction quantile.
qS2 <- function(p) {
  resultat <- numeric(length(p))
  bas <- p < 0.25
  haut <- p > 0.75

  resultat[bas] <- 2 * qnorm(2 * p[bas])
  resultat[haut] <- 2 * qnorm(2 * p[haut] - 1)
  resultat
}

VaR1_theo <- sqrt(2) * qnorm(p)
VaR2_theo <- qS2(p)

# ---- Données pour ggplot ----
df_var <- data.frame(
  p = rep(p, 2),
  VaR = c(VaR1_theo, VaR2_theo),
  Portefeuille = factor(rep(c("X1+Y1", "X2+Y2"), each = length(p)))
)
# ---- Graphique ggplot ----
ggplot(df_var, aes(x = p, y = VaR, color = Portefeuille, linetype = Portefeuille)) +
  geom_line(linewidth = 1) +
  labs(title = "VaR des portefeuilles X1+Y1 et X2+Y2",
       x = "Niveau de probabilité p",
       y = "VaR(p)") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme(legend.position = "top")
