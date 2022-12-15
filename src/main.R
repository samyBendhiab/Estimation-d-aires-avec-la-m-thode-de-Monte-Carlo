#-------------------------------------------------------------#
# Description: Estimation de pi par la méthode de Monte Carlo
#-------------------------------------------------------------#
# Author: Samy Ben Dhiab
#-------------------------------------------------------------#

#Comment tirer uniformément au hasard n points dans un rectangle ? Indication: utiliser la primitive runif, qui permet de tirer uniformément au hasard un ou plusieurs points dans un intervalle.

#runif(n, min = 0, max = 1)

#Comment faire pour savoir si un point de coordonnées (x,y) appartient ou non au cercle de centre (0,0) et de rayon R ?
#(x,y)€ cercle => x²+y²<=R²


#Définir une fonction mc.pi qui prend en argument un entier n, et renvoie une valeur approchée de π, obtenue à l’aide de la méthode de Monte Carlo, et avec n points tirés uniformément au hasard.

mc.pi<-function(n){
  #point alea
  x<-runif(n, min = 0, max = 1)
  y<-runif(n, min = 0, max = 1)
  R<-1 #rayon du cercle
  N<-0 #nombre de points dans le cercle
  for (i in 1:n){
    if (x[i]^2+y[i]^2<=R^2){N<-N+1} #si le point est dans le cercle
  }
  return(4*N/n)
}

#function test grandeur de mc.pi
testmc.pi<-function(n){
  s<-1:n
  s<-10**s
  res<-sapply(s,mc.pi)
  return(res)
}

testmc.pi(8)

#Estimations
#Dans cette section, vous allez définir une matrice PIE de taille t*p (avec t=50 et p=7), contenant des estimations de π. Plus précisément, la coordonnée (i,j) de PIE (avec 1<=i<=t et 1<=j<=p) doit contenir une estimation de π effectuée avec n=10**j points. Ainsi, la j-ième colonne de PIE contient t estimations de π, toutes effectuées avec n=10**j points. Toutes les estimations seront faites de façon indépendante les unes des autres.
t <- 50
p <- 7
PIE <- matrix(0, t, p)
for (i in 1:t) {
  for (j in 1:p) {
    PIE[i, j] <- mc.pi(10**j)
  }
}
#Dans cette section, vous allez définir un vecteur tE de taille p contenant le temps moyen mis pour obtenir de telles estimations. Plus précisément, la j-ième coordonnée du vecteur tE (avec 1<=j<=p) doit contenir le temps moyen mis pour effectuer une estimation de π, chacune de ces estimations étant effectuée avec n=10**j points.
#Indication: utiliser la primitive system.time, qui renvoie le temps mis pour évaluer un expression donnée et la primitive replicate qui répète l’évaluation d’une expression.

tE <- vector("numeric", p)
for (j in 1:p) {
  tE[j] <- system.time(replicate(t, mc.pi(10**j)))[3]
}
tE

#Erreur relative
#Quelle est la formule pour l’erreur relative entre une valeur et son estimation ? Définir une matrice ERR de taille t*p dont la coordonnée (i,j) est l’erreur relative entre π et son estimation PIE[i,j].

ERR <- matrix(0, t, p)
for (i in 1:t) {
  for (j in 1:p) {
    ERR[i, j] <- abs(PIE[i, j] - pi) / pi
  }
}

ERR <- abs(PIE/pi - 1)
par(mfrow=c(1,2),mar=c(4,4,2,2)+0.1)
boxplot(ERR, main='Erreur relative sur PI', log='y', xlab='#points', ylab='Rel. Error')
plot(10 ** seq(1:p), tE, type='b', main='Temps moyen d\'une simulation', log='x', xlab='#points', ylab='Time')

#Conclusion
#Quelle est la relation entre l’erreur relative et le temps mis pour obtenir une estimation ?
#La relation est linéaire, plus le nombre de points est grand, plus le temps de calcul est long et plus l'erreur est petite.
#Quelle est la relation entre l’erreur relative et le nombre de points utilisés pour obtenir une estimation ?
#Plus le nombre de points est grand, plus l'erreur est petite.
#Quelle est la relation entre le temps mis pour obtenir une estimation et le nombre de points utilisés pour obtenir une estimation ?
#Plus le nombre de points est grand, plus le temps de calcul est long.


# Definition d'une fonction très utile
creer_polygone <- function (x,y) {
  matrix(c(x, x[1], y, y[1]), ncol=2,dimnames=list(c(), c("x","y")))
}

carre <- creer_polygone(c(10,10,90,90), c(30, 70, 70, 30))
# Une permutation cyclique des points donne le même polygone
carre <- creer_polygone(c(10,90,90,10), c(70, 70, 30, 30))
# En revanche, le code suivant ne définit pas un rectangle,
# mais un polygone dont les arêtes se croisent.
papillon <- creer_polygone(c(10,90,10,90), c(30,70,70,30))
# pour finir, voici un losange.
losange <- creer_polygone(c(50,10,50,90),c(30,50,70,50))

#affiche un plot du carre
#+BEGIN_SRC R :exports both :file act07/dessin_poly.jpg :width 300 :height 300 :session poly plot(carre, type=’l’) lines(papillon -1, type=’b’, col=’firebrick’) lines(losange, type=’l’, col=’darkblue’) #+END_SRC R
dessin_polynome<-function(polynome){

  plot(polynome, type='l')
  lines(polynome, type="l", col="darkblue")
}
dessin_polynome(carre)





#Définir une fonction reg_poly <- function(n, r=1) { ... } qui prend en argument un entier n, un réel strictement positif r (de valeur 1 par défaut), et qui renvoie un polygone p vérifiant:
#le polygone p a n côtés,
#il est inscrit dans un cercle de centre (0,0) et de rayon r.

reg_poly<-function(n, r=1){
  #Tous les sommets d’un polygone régulier peuvent être engendrés à partir d’un seul sommet et les images successives d’une rotation: quel est l’angle de cette rotation ? Donnez-le en radians et en fonction du nombre n de côtés du polygone.
  angle<-2*pi/n

  poly<-matrix(0,n+1,2)

  for (i in 1:n) {
    poly[i,1]<-r*cos(angle*i)
    poly[i,2]<-r*sin(angle*i)
  }
  poly[n+1,1]<-r*cos(angle*1)
  poly[n+1,2]<-r*sin(angle*1)

  #Quelles sont les coordonnées cartésiennes (x,y) d’un point dont les coordonnées polaires sont (r,theta)?
  #x=r*cos(theta)
  #y=r*sin(theta)

  return(poly)
}
#dessin des polynomes de 5 à 10 côtés
par(mfrow=c(2,3),mar=c(4,4,2,2)+0.1)
for (i in 5:10) {
  dessin_polynome(reg_poly(i))
}

#4.2 Polygone surprise
x <- c(0,0,9,11,11,9,8,11,9,6,3,3,8,9,9,8,2,2)
y <- c(0,12,12,10,7,5,5,0,0,5,5,7,7,8,9,10,10,0)
surprise <- creer_polygone(x,y)
plot(surprise,col="black", type="l")
# donne le logo de R

#5 Approximation de l’aire d’un polygone simple
#5.1 Tirer un ou plusieurs points uniformément au hasard dans un rectangle
#Définir une fonction boite qui prend en argument un polygone donné sous la forme d’une matrice (comme précédemment) et renvoie une matrice contenant l’abscisse minimale du plus petit rectangle contenant ce polygone, son abscisse maximale, son ordonnée minimale et son ordonnée maximale. Par exemple:

boite<-function(poly){

  boite<-matrix(0,2,2)

  boite[1,1]<-min(poly[,1])
  boite[1,2]<-min(poly[,2])
  boite[2,1]<-max(poly[,1])
  boite[2,2]<-max(poly[,2])
  return(boite)
}


bo <- boite(losange)
#losange
#x  y
#[1,] 50 30
#[2,] 10 50
#[3,] 50 70
#[4,] 90 50
#[5,] 50 30

#bo
#x  y
#min 10 30
#max 90 70


#Tirage de points uniformément aléatoirement dans un rectangle
#Définir une fonction points_aleatoires qui prend en argument un couple (n, bo), où n est un entier et bo une boîte rectangulaire, et renvoie une matrice M contenant n points tirés uniformément au hasard dans le rectangle r=[xmin;xmax]*[ymin;ymax]. Plus précisément, la matrice M est de taille n*2 et chaque colonne contient un point tiré uniformément au hasard dans le rectangle r. Par exemple:

points_aleatoires<-function(n, bo){

  points<-matrix(0,n,2)

  for (i in 1:n) {
    points[i,1]<-runif(1,bo[1,1],bo[2,1])
    points[i,2]<-runif(1,bo[1,2],bo[2,2])
  }
  return(points)
}
bo <- matrix(c(3, 5, 6, 8),nrow=2, dimnames=list(c("min","max"), c("x","y")))
pts <- points_aleatoires(5, bo)

cat("bo=",bo,"\t","pts=",pts,"\n")


#5.2 Un point donné est-il à l’intérieur ou à l’extérieur du polygone ?

#Définir une fonction appartient(points, polygone) qui prend en arguments des points et un polygone, et renvoie pour chaque point TRUE si le point est à l’intérieur du polygone, FALSE sinon. On pourra s’appuyer sur une fonction auxiliaire appartient_poly(point, polygone) qui prend en arguments un point et un polygone, et renvoie TRUE si le point est à l’intérieur du polygone, FALSE sinon. On ne se souciera pas du résultat renvoyé par la fonction dans le cas où le point appartient à un des côtés du polygone. En théorie, les points seront tirés uniformément aléatoirement et la probabilité qu’un tel cas se produise sera nulle.

#Indication: Lorsqu’un point est à l’intérieur d’un polygone, toute demi-droite partant de ce point possède un nombre impair d’intersections avec les côtés du polygone. Lorsqu’il est à l’extérieur, elle en possède nécessairement un nombre pair.


appartient_poly<-function (point,polygone){
  # Verifier si le point est à l'intérieur du polygone
  # On va utiliser la methodes du raytracing
  # Le point est à l'intérieur du polygone si :
  # le nombre d'intersections de la demi-droite partant du point vers la droite est impair
  # le x du point est plus grand que le x min du polygone
  # le x du point est plus petit que le x max du polygone
  # le y du point est plus grand que le y min du polygone
  # le y du point est plus petit que le y max du polygone
  resultat<-FALSE
  n<-nrow(polygone)
  for (i in 1:n) {
    if (i==n) {
      j<-1
    } else {
      j<-i+1
    }
    if (point[2]>min(polygone[i,2],polygone[j,2]) && point[2]<=max(polygone[i,2],polygone[j,2])) {
      if (point[1]<=max(polygone[i,1],polygone[j,1])) {
        if (polygone[i,2]!=polygone[j,2]) {

          # y = ax+b
          # x = (y-b)/a
          # x = (y - y1) * (x2 - x1) / (y2 - y1) + x1
          # on evalue x pour y=point[2]
          # x = (point[2] - y1) * (x2 - x1) / (y2 - y1) + x1

          x1<-polygone[i,1]
          x2<-polygone[j,1]
          y1<-polygone[i,2]
          y2<-polygone[j,2]
          xinters<-((point[2]-y1)*(x2-x1)/(y2-y1)+x1)
          if (polygone[i,1]==polygone[j,1] || point[1]<=xinters) {
            resultat<-!resultat
          }
        }
      }
    }
  }
  return  (resultat)
}


appartient<-function(points,polynome){
  #qui prend en arguments des points et un polygone, et renvoie pour chaque point TRUE si le point est à l’intérieur du polygone, FALSE sinon.
  resultat<-rep(FALSE,nrow(points))
  for (i in 1:nrow(points)) {
    resultat[i]<-appartient_poly(points[i,],polynome)
  }
  return(resultat)
}
## Réaliser un test de la fonction
carre <- creer_polygone(c(0, 0, 1, 1), c(0, 1, 1, 0))
cc <- seq(from=-0.25,to=1.25,by=0.25)
points <- do.call(rbind,lapply(cc, FUN=cbind, cc,deparse.level = 0))
pin <- appartient(points,carre);

## Dessiner le résultat du test
par(mar=c(2,2,3,2)+0.1)
plot(carre, type='l', main="Test de la fonction appartient", xlim=range(carre[,1],points[,1]), ylim=range(carre[,2],points[,2]))
points(points[pin,1], points[pin,2], col='firebrick', pch=20)
points(points[!pin,1], points[!pin,2], col='darkblue', pch=20)


#fonction qui calcul l'aire d'un polygone
aire<-function(bo){
  diffx<-bo[2,1]-bo[1,1]
  diffy<-bo[2,2]-bo[1,2]
  aire<-diffx*diffy
  return(abs(aire))
}


dessin_points<-function(points,pin){
  #dessiner les points
  #si le point est in le point est rouge
  if (pin) { points(points[,1], points[,2], col='firebrick', pch=20)
  }
    #sinon le point est bleu
  else { points(points[,1], points[,2], col='darkblue', pch=20)}

}
#Définir une fonction mc.poly qui prend en argument un entier n correspondant au nombre de points à tirer au hasard et un polygone, et qui renvoie une valeur approchée de l'aire du ~polygone par la méthode de Monte Carlo.
mc.poly<-function(n,polygone,DRAW=FALSE){
  #faire un carre de taille max des dimensions du polygone
  bo<-boite(polygone)
  #tirer n points aleatoires
  pts <- points_aleatoires(n, bo)
  #calculer le nombre de points dans le polygone
  nbpts<-sum(appartient(pts,polygone))
  #portion de points dans le polygone
  p<-nbpts/n
  if(DRAW){
    #dessine le polygone et les points
    dessin_polynome(polygone)
    #points dans le polygone
    dessin_points(pts[appartient(pts,polygone),],TRUE)
    #points hors du polygone
    dessin_points(pts[!appartient(pts,polygone),],FALSE)
  }

  return(p*aire(bo))


}
#affiche les polynomes et les points sur un graphique coupe en 4
par(mfrow=c(2,2))
cat("approximation de l'aire du polygone par la methode de Monte Carlo\n")
cat("Aire du losange avec 10 points",mc.poly(10,losange,TRUE),"\n")
cat("Aire du losange avec 100 points",mc.poly(100,losange,TRUE),"\n")
cat("Aire du losange avec 1000 points",mc.poly(1000,losange,TRUE),"\n")
cat("Aire du losange avec 10000 points",mc.poly(10000,losange,TRUE),"\n")

## l'aire du losange et du carré est de 1 c'est ce que j'ai obtenu a la main , le resultat est satisfaisant

# Test avec affichage de surprises car le polygone n'est regulier
print(mc.poly(10000, surprise,TRUE))

#Définir une fonction aire.poly qui prend en argument un polygone et calcule son aire exacte.
aire.poly<-function(polygone){
  #calculer l'aire du polygone avec la somme des triangles
  aire<-0
  for (i in 1:(nrow(polygone)-1)) {
    aire<-aire+(polygone[i+1,1]+polygone[i,1])*(polygone[i+1,2]-polygone[i,2])
  }
  aire<-abs(aire/2)
  return(aire)

}

print(aire.poly(surprise))#-71
print(mc.poly(5000,surprise))

poly3<-reg_poly(3,1)
poly10<-reg_poly(10,1)

#calcul des aires exactes
aire3<-aire.poly(poly3)
aire10<-aire.poly(poly10)
airelosange<-aire.poly(losange)
airecarre<-aire.poly(carre)
airesurprise<-aire.poly(surprise)

#calcul des aires approchees en fonction du nombre de points
nf<-5#10^n = plus grand nombre de points teste
nd<-3 #10^nd = plus petit nombre de points teste
echelle<-seq(nd,nf)
n<-length(echelle) #10 puissance cmb de points
nrep<-7 #nombre de simulation
aire3mc<-numeric(n)
aire10mc<-numeric(n)
airelosangemc<-numeric(n)
airecarremc<-numeric(n)
airesurprisemc<-numeric(n)



for (i in 1:n){
  aire3mc[i]<-mean(replicate(nrep,mc.poly(10**echelle[i],poly3)))
    aire10mc[i]<-mean(replicate(nrep,mc.poly(10**echelle[i],poly10)))
    airelosangemc[i]<-mean(replicate(nrep,mc.poly(10**echelle[i],losange)))
    airecarremc[i]<-mean(replicate(nrep,mc.poly(10**echelle[i],carre)))
    airesurprisemc[i]<-mean(replicate(nrep,mc.poly(10**echelle[i],surprise)))
}

#calcul des erreurs
erreur3<-abs(aire3-aire3mc)
erreur10<-abs(aire10-aire10mc)
erreurlosange<-abs(airelosange-airelosangemc)
erreurcarre<-abs(airecarre-airecarremc)
erreursurprise<-abs(airesurprise-airesurprisemc)

#on va tracer les courbes

par(new=TRUE)
plot(10^echelle,erreur3,type="l",col="darkblue",xlab="Nombre de points",ylab="Erreur",main="Erreur en fonction du nombre de points",log="xy")
legend("top",legend=c("polygone a 3 cote","polygone a 10 cote","losange","carre","surprise"),col=c("darkblue","firebrick","darkgreen","darkorange","darkviolet"),lty=1)
par(new=TRUE)
lines(10^echelle,erreur10,type="l",col="firebrick")
par(new=TRUE)
lines(10^echelle,erreurlosange,type="l",col="darkgreen")
par(new=TRUE)
lines(10^echelle,erreurcarre,type="l",col="darkorange")
par(new=TRUE)
lines(10^echelle,erreursurprise,type="l",col="darkviolet")

#on voit que plus le nombre de points est grand plus l'erreur est petite
#on voit que pour le polygone a 3 cote l'erreur est plus grande que pour le polygone a 10 cote
#on voit que pour le polygone a 3 cote l'erreur est plus grande que pour le losange
#on voit que pour le polygone a 3 cote l'erreur est plus grande que pour le carre
#on voit que pour le polygone a 3 cote l'erreur est plus grande que pour le surprise
#on voit que pour le carre et le lonsange l'erreur reste constante

#Conclusion:
# on voit que plus le nombre de points est grand plus l'erreur est petite