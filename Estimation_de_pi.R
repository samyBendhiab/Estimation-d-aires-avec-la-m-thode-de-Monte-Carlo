#-------------------------------------------------------------#
# Description: Estimation de pi par la méthode de Monte Carlo
#-------------------------------------------------------------#

#-------------------------------------------------------------#
#Parti booléen pour les partie

section_pi<-FALSE


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

#testmc.pi(8)

#Estimations
#Dans cette section, vous allez définir une matrice PIE de taille t*p (avec t=50 et p=7), contenant des estimations de π. Plus précisément, la coordonnée (i,j) de PIE (avec 1<=i<=t et 1<=j<=p) doit contenir une estimation de π effectuée avec n=10**j points. Ainsi, la j-ième colonne de PIE contient t estimations de π, toutes effectuées avec n=10**j points. Toutes les estimations seront faites de façon indépendante les unes des autres.
if(section_pi){ t <- 50
  p <- 7
  PIE <- matrix(0, t, p)
  for (i in 1:t) {
    for (j in 1:p) {
      PIE[i, j] <- mc.pi(10**j)
    }
  }
  #Dans cette section, vous allez définir un vecteur tE de taille p contenant le temps moyen mis pour obtenir de telles estimations. Plus précisément, la j-ième coordonnée du vecteur tE (avec 1<=j<=p) doit contenir le temps moyen mis pour effectuer une estimation de π, chacune de ces estimations étant effectuée avec n=10**j points.
  #Indication: utiliser la primitive system.time, qui renvoie le temps mis pour évaluer un expression donnée et la primitive replicate qui répète l’évaluation d’une expression.

  tE <- system.time(replicate(p, mc.pi(10**j)))
  tE


  #Erreur relative
  #Quelle est la formule pour l’erreur relative entre une valeur et son estimation ? Définir une matrice ERR de taille t*p dont la coordonnée (i,j) est l’erreur relative entre π et son estimation PIE[i,j].

  #la valeur absolue
  ERR <- abs(PIE / pi - 1)

  par(mfrow = c(1, 2), mar = c(4, 4, 2, 2) + 0.1)
  boxplot(ERR, main = 'Erreur relative sur PI', log = 'y', xlab = '#points', ylab = 'Rel. Error')
  plot(10^(1:p), tE, type = 'b', main = 'Temps moyen d\'une simulation', log = 'x', xlab = '#points', ylab = 'Time') }





## Definition d'une fonction très utile
creer_polygone <- function (x,y) {
  matrix(c(x, x[1], y, y[1]), ncol=2,dimnames=list(c(), c("x","y")))
}

carre <- creer_polygone(c(10,10,90,90), c(30, 70, 70, 30))
## Une permutation cyclique des points donne le même polygone
carre <- creer_polygone(c(10,90,90,10), c(70, 70, 30, 30))
## En revanche, le code suivant ne définit pas un rectangle,
## mais un polygone dont les arêtes se croisent.
papillon <- creer_polygone(c(10,90,10,90), c(30,70,70,30))
## pour finir, voici un losange.
losange <- creer_polygone(c(50,10,50,90),c(30,50,70,50))

print(carre)
#affiche un plot du carre
#+BEGIN_SRC R :exports both :file act07/dessin_poly.jpg :width 300 :height 300 :session poly plot(carre, type=’l’) lines(papillon -1, type=’b’, col=’firebrick’) lines(losange, type=’l’, col=’darkblue’) #+END_SRC R

dessin_polynome<-function(polynome){
    plot(polynome, type="l")
    lines(polynome -1, type="b", col="firebrick")
    lines(polynome, type="l", col="darkblue")
}

dessin_polynome(carre)

#Définir une fonction reg_poly <- function(n, r=1) { ... } qui prend en argument un entier n, un réel strictement positif r (de valeur 1 par défaut), et qui renvoie un polygone p vérifiant:
#le polygone p a n côtés,
#il est inscrit dans un cercle de centre (0,0) et de rayon r.

reg_poly<-function(n, r=1){
  #Tous les sommets d’un polygone régulier peuvent être engendrés à partir d’un seul sommet et les images successives d’une rotation: quel est l’angle de cette rotation ? Donnez-le en radians et en fonction du nombre n de côtés du polygone.
  angle=2*pi/n
  #Définir une matrice de taille n*2, dont la j-ième ligne contient les coordonnées du j-ième sommet du polygone p.
  poly<-matrix(0,n+1,2)
  #Définir une matrice de taille n*2, dont la j-ième ligne contient les coordonnées du j-ième sommet du polygone p.
  for (i in 1:n) {
    poly[i,1]<-r*cos(angle*i)
    poly[i,2]<-r*sin(angle*i)
  }
  poly[n+1,1]<-r*cos(angle*1)
  poly[n+1,2]<-r*sin(angle*1)

  #Quelles sont les coordonnées cartésiennes (x,y) d’un point dont les coordonnées polaires sont (r,theta)?
  #x=r*cos(theta)
  #y=r*sin(theta)
  #Tester votre fonction reg_poly en dessinant un exemple de polygone régulier (vous pouvez utiliser la fonction dessin_polygone définie précédemment).
  return(poly)
}

dessin_polynome(reg_poly(5, 1))

#4.2 Polygone surprise
x <- c(0,0,9,11,11,9,8,11,9,6,3,3,8,9,9,8,2,2)
y <- c(0,12,12,10,7,5,5,0,0,5,5,7,7,8,9,10,10,0)
surprise <- creer_polygone(x,y)
plot(surprise,col="black", type="l")

# donne le logo de R

#5 Approximation de l’aire d’un polygone simple
#5.1 Tirer un ou plusieurs points uniformément au hasard dans un rectangle
#Définir une fonction boite qui prend en argument un polygone donné sous la forme d’une matrice (comme précédemment) et renvoie une matrice contenant l’abscisse minimale du plus petit rectangle contenant ce polygone, son abscisse maximale, son ordonnée minimale et son ordonnée maximale. Par exemple:

#print(losange)
#bo <- boite(losange)
#print(bo)
#x  y
#[1,] 50 30
#[2,] 10 50
#[3,] 50 70
#[4,] 90 50
#[5,] 50 30
#x  y
#min 10 30
#max 90 70


boite<-function(poly){
  #Définir une matrice de taille 2*2, dont la j-ième ligne contient les coordonnées du j-ième sommet du polygone p.
    boite<-matrix(0,2,2)

    boite[1,1]<-min(poly[,1])
    boite[1,2]<-min(poly[,2])
    boite[2,1]<-max(poly[,1])
    boite[2,2]<-max(poly[,2])
  return(boite)
}

print(losange)
bo <- boite(losange)
print(bo)

#Tirage de points uniformément aléatoirement dans un rectangle
#Définir une fonction points_aleatoires qui prend en argument un couple (n, bo), où n est un entier et bo une boîte rectangulaire, et renvoie une matrice M contenant n points tirés uniformément au hasard dans le rectangle r=[xmin;xmax]*[ymin;ymax]. Plus précisément, la matrice M est de taille n*2 et chaque colonne contient un point tiré uniformément au hasard dans le rectangle r. Par exemple:

points_aleatoires<-function(n, bo){
  #Définir une matrice de taille n*2, dont la j-ième ligne contient les coordonnées du j-ième sommet du polygone p.
  points<-matrix(0,n,2)

  for (i in 1:n) {
    points[i,1]<-runif(1,bo[1,1],bo[2,1])
    points[i,2]<-runif(1,bo[1,2],bo[2,2])
  }
  return(points)
}
bo <- matrix(c(3, 5, 6, 8),nrow=2, dimnames=list(c("min","max"), c("x","y")))
pts <- points_aleatoires(5, bo)
print(pts)

#5.2 Un point donné est-il à l’intérieur ou à l’extérieur du polygone ?

#Définir une fonction appartient(points, polygone) qui prend en arguments des points et un polygone, et renvoie pour chaque point TRUE si le point est à l’intérieur du polygone, FALSE sinon. On pourra s’appuyer sur une fonction auxiliaire appartient_poly(point, polygone) qui prend en arguments un point et un polygone, et renvoie TRUE si le point est à l’intérieur du polygone, FALSE sinon. On ne se souciera pas du résultat renvoyé par la fonction dans le cas où le point appartient à un des côtés du polygone. En théorie, les points seront tirés uniformément aléatoirement et la probabilité qu’un tel cas se produise sera nulle.

#Indication: Lorsqu’un point est à l’intérieur d’un polygone, toute demi-droite partant de ce point possède un nombre impair d’intersections avec les côtés du polygone. Lorsqu’il est à l’extérieur, elle en possède nécessairement un nombre pair.


appartient_poly<-function (point,polygone){
  lineaire<-function (x1,y1,x2,y2){
    a<-x1-x2
    b<-y2-y1
    a<-a/b
    return(c(a,b))
  }
  #tout les fonction lineaire de chaque segment
  l_segment_x<-numeric(nrow(polygone)-1)
  l_segment_y<-numeric(nrow(polygone)-1)
  l_segment_z<-numeric(nrow(polygone)-1)
    for (i in 1:(nrow(polygone)-1)){
        l_segment_x[i]<-lineaire(polygone[i,1],polygone[i,2],polygone[i+1,1],polygone[i+1,2])[1]
        l_segment_y[i]<-lineaire(polygone[i,1],polygone[i,2],polygone[i+1,1],polygone[i+1,2])[2]



    }
    #on calcule le produit scalaire de chaque segment avec le vecteur point
    cat('\npoint x\n')
    cat(l_segment_x)
    cat('\npoint y\n')
    cat(l_segment_y)


}
x<-c(0,-2,-1.25,1.25,2)
y<-c(6,4,2,2,4)
poly<-creer_polygone(x,y)
poly
appartient_poly(c(0,0),poly)

  resultat<-0
  #prendre tous les points du polygone
  for(i in 1:(nrow(polygone)-1)){}
  #pour chaque point on trace une demie droite imaginaires horizontale

  #on compte le nombre de fois ou cette demie droite coupe le polygone
  #si le nombre de fois est pair alors le point est à l'extérieur
  #si le nombre de fois est impair alors le point est à l'intérieur
  #on renvoie le resultat

    return(vecteur)

}
appartient<-function(points,polynome){
  #Définir une matrice de taille n*2, dont la j-ième ligne contient les coordonnées du j-ième sommet du polygone p.
  n<-nrow(points)
  appartient<-matrix(0,n,1)

  for (i in 1:n) {
    appartient[i,1]<-appartient_poly(points[i,],polynome)
  }
  return(appartient)
}