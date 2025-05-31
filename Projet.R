### Projet de séries temporelles linéaires ###
require(zoo)
library(zoo)
require(tseries)
library(tseries)
library(fUnitRoots)


## Récupération des données et cleaning ##

#Indice CVS-CJO de la production industrielle (base 100 en 2021) - Extraction de pétrole brut (NAF rév. 2, niveau groupe, poste 06.1) 
datafile <- "C:/Users/Antony/Desktop/Time series/valeurs_mensuelles.csv"
df <- read.csv(datafile, sep=';')

#On garde seulement les lignes où libellé est au format YYYY-MM
df <- df[grepl("^\\d{4}-\\d{2}$", df$Libellé), ]
df <- df[, 1:2]
colnames(df)[1:2] <- c("date", "indice") 


df$date <- as.yearmon(df$date, format = "%Y-%m") 
df <- df[df$date >= as.yearmon("2000-01"), ]
df$indice <- as.numeric(df$indice)

#reset les indices du dataframe
rownames(df) <- NULL 

## PARTIE 1 : Les données

# Utilisez directement les dates du dataframe pour tracer le graphique
plot(df$date, df$indice, type = "l", xlab = "Date", ylab = "Indice", main = "Indice de la production industrielle")

indice <- zoo(df$indice, order.by = df$date)
dindice <- diff(indice,1)

plot(indice)
plot(dindice)
#sans différenciation, on observe une trend décroissante
#avec différenciation de premier ordre, la série semble stationnaire, 
#on observe cependant quelques pics de volatilité, notamment, de manière prévisible, après 2020

lm <- lm(indice~dates)
summary(lm)
#Avant de procéder aux tests de racine unitaire, il convient de vérifier s’il y a une constante et/ou une tendance
#linéaire non nulle. La représentation graphique de indice a montré que la tendance est plutôt linéaire et décroissante
#Vérifions cela par des tests statistiques 

#Test de racine unitaire (pour tester la stationnarité d'une série)
adf <- adfTest(log_indice, type="ct", lags=0)  # tendance et constante

#Avant d’interpréter le test, vérifions que les résidus du modèle de régression sont bien non autocorrélés, 
#sans quoi le test ne serait pas valide.

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
#Comme la série est mensuelle, testons l’autocorrélation des résidus jusqu’à l’ordre 24 (deux ans), sans oublier de corriger les degrés de libertés du nombre de régresseurs.
Qtests(adf@test$lm$residuals,24,length(adf@test$lm$coefficients))

#L’absence d’autocorrélation des résidus est rejetée au moins une fois, le test ADF avec aucun retard n’est donc pas valide 
#Ajoutons des retards de ∆Xt jusqu’à ce que les résidus ne soient plus autocorrélés.

adfTest_valid <- function(series,kmax,type){ 
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("nope \n")
    k <- k + 1
  }
  return(adf)
}
adf <- adfTest_valid(indice,24,"ct")

#Il a fallu considérer 9 retards au test ADF pour supprimer l’autocorrélation des résidus.

adf #affichage des résultats du test valide maintenu
#l'hypothèse de présence de racine unitaire (non-stationnarité) est acceptée, la série n'est pas stationnaire


#Testons maintenant la racine unitaire pour la série différenciée dindince
#La représentation graphique précédente semble montrer l’absence de constante et de tendance non nulle
summary(lm(dindice ~ dates[-1])) #sans la première date car on a différencié la série
#on a bien ni trend ni constante significative
adf <- adfTest_valid(dindice,24, type="nc")
# on doit inclure 5 lags pour rendre le test valide 
adf
#on a une p-value de 0.01, la série est donc stationnaire
#la série indice est donc I(1)

pp.test(dindice)
#le pp test confirme notre résultat, p-value inférieure à 0.01


## PARTIE 2 : Modèle ARMA

#choisir un modèle ARIMA(p,d,q) sur indice revient à choisir un ARMA(p,q) pour dindice

par(mfrow=c(1,1))
acf(dindice, 24)
pacf(dindice, 24)

x <- dindice
pmax <- 5
qmax <- 2

mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1)
rownames(mat) <- paste0("p=",0:pmax)
colnames(mat) <- paste0("q=",0:qmax)
AICs <- mat
BICs <- mat
pqs <- expand.grid(0:pmax,0:qmax)

for (row in 1:nrow(pqs)){
  p <- pqs[row,1]
  q <- pqs[row,2]
  estim <- try(arima(x, c(p,0,q), include.mean = FALSE))
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim)
}

AICs
BICs
AICs==min(AICs, na.rm = TRUE) #ARMA(2,2)
BICs==min(BICs, na.rm = TRUE) #MA(2)


signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1 - pnorm(abs(t))) * 2
  return(rbind(coef, se, pval))
}

arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals, 24, length(estim$coef))
  pvals <- matrix(apply(matrix(1:24, nrow=6), 2, function(c) round(pvals[c,],3)), nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"), 4)
  cat("tests de nullité des coefficients :\n")
  print(adjust)
  cat("\ntests d’absence d’autocorrélation des résidus : \n")
  print(pvals)
}

# MA(2)
ma2 <- arima(dindice, c(0,0,2), include.mean = FALSE)
cat("MA(2)\n")
arimafit(ma2)
#modèle non valide

# ARMA(2,2)
ar2ma2 <- arima(train, c(2,0,2), include.mean = FALSE)
cat("ARMA(2,2)\n")
arimafit(ar2ma2)
#Le lag 5 a une p-value de 0.093, le modèle n'est pas valide

#ARMA(1,2)
ar1ma2 <- arima(train, c(1,0,2), include.mean = FALSE)
cat("ARMA(1,2)\n")
arimafit(ar1ma2)
#modèle valide, coef significatifs

#ARMA(2,1)
ar2ma1 <- arima(train, c(2,0,1), include.mean = FALSE)
cat("ARMA(2,1)\n")
arimafit(ar2ma1)

#Le AIC et le BIC du modèle ARMA(2,1) sont légèrement inférieurs à ceux du modèle ARMA(1,2)
#On retient donc finalement le modèle ARMA(2,1)

### Partie 3 : Prédiction ###
library(forecast)
library(ggplot2)

forecast_values <- forecast(ar2ma1, h = 2)

forecast_df <- data.frame(
  Date = seq(max(time(dindice)) + 1/12, by = 1/12, length.out = 2),
  Forecast = as.numeric(forecast_values$mean),
  Lower = as.numeric(forecast_values$lower[, 1]),
  Upper = as.numeric(forecast_values$upper[, 1])
)

historical_df <- data.frame(
  Date = tail(time(dindice), 24), # Considérer les 2 dernières années pour la visibilité
  Value = as.numeric(tail(dindice, 24))
)

# Tracer les données historiques et les prévisions
ggplot() +
  geom_line(data = historical_df, aes(x = Date, y = Value), color = "black") +
  geom_point(data = forecast_df, aes(x = Date, y = Forecast), color = "blue", size = 3) +
  geom_ribbon(data = forecast_df, aes(x = Date, ymin = Lower, ymax = Upper), fill = "grey", alpha = 0.5) +
  labs(title = "Prévision à l'horizon T+1 et T+2",
       x = "Temps",
       y = "Valeur") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

