library(rgdal)
w <- readOGR('~/shp/TM_WORLD_BORDERS-0.3') # https://thematicmapping.org/downloads/world_borders.php
b <- readOGR('~/shp/Biomas_250mil/lm_bioma_250') # https://geoftp.ibge.gov.br/informacoes_ambientais/estudos_ambientais/biomas/vetores/
p <- read.csv('refs.csv',sep='\t')
p = p[!is.na(p$lat),] # remove pontos sem coordenadas
head(p)
pV <- p[is.na(p$soApis),] # outras espécies (também)
pA <- p[!is.na(p$soApis),] # só Apis

Rjitter <- function(x,y,amount=0.5) { # função para embaralhar coordenadas
  ang <- runif(length(x),0,2*pi)
  raio <- runif(length(x),0,amount)
  cbind(x+raio*cos(ang),y+raio*sin(ang))
}

155+125 # 280
72+56 # 128

png('map_metabarcodingA.png',2800,1280)
par(mar=c(0,0,0,0))
plot(w,xlim=c(-125,155),ylim=c(-56,72),asp=1,border=NA,xaxs='i',yaxs='i')
abline(h=0,lwd=1.5)
abline(h=c(-23.43657,23.43657),lty='FF',lwd=1.5)
plot(w,lwd=1,col='#dddddd',border='#808080',add=T)
plot(b[which(b$Bioma == 'Cerrado'),],col='orange',border=NA,add=T)
points(Rjitter(pV$lon,pV$lat),bg='red',pch=21,cex=2.5)
points(Rjitter(pA$lon,pA$lat),bg='yellow',pch=21,cex=2.5)
dev.off()

png('map_metabarcodingV.png',2800,1280)
par(mar=c(0,0,0,0))
plot(w,xlim=c(-125,155),ylim=c(-56,72),asp=1,border=NA,xaxs='i',yaxs='i')
abline(h=0,lwd=1.5)
abline(h=c(-23.43657,23.43657),lty='FF',lwd=1.5)
plot(w,lwd=1,col='#dddddd',border='#808080',add=T)
plot(b[which(b$Bioma == 'Cerrado'),],col='orange',border=NA,add=T)
points(Rjitter(pA$lon,pA$lat),bg='yellow',pch=21,cex=2.5)
points(Rjitter(pV$lon,pV$lat),bg='red',pch=21,cex=2.5)
dev.off()
