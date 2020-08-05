# Análisis del canto de un grillo con R
# www.overfitting.net
# https://www.overfitting.net/2018/07/analisis-del-canto-de-un-grillo-con-r.html

library(tuneR)
library(corrplot)


# FUNCIONES AUXILIARES

# Conversiones tiempo <-> samples. t=0 <-> n=1
time2sample=function(t, fs=44100) return(round(t*fs+1))
sample2time=function(n, fs=44100) return((n-1)/fs)


# ANÁLISIS GRILLO.WAV
                                     
grillo=readWave("grillo.wav")
play(grillo)
grillo

fs=as.numeric(grillo@samp.rate)
waveform=grillo@left

dft=abs(fft(waveform))
N=round(length(dft)/2)  # Primera mitad de la FFT
maxfreq=grillo@samp.rate/2/1000  # Máx. frecuencia FFT en kHz
plot(seq(from=0, to=maxfreq, len=N),
    dft[1:N]/max(dft), main='FFT "grillo.wav"',
    xlab='Frecuencia (kHz)', ylab='Amplitud (Lin.)', col='red', type='l')
axis(side=1, at=c(0:maxfreq))

fgrillo=which( round(dft)==max(round(dft)) )[1] * fs/length(dft)
Tgrillo=1/fgrillo


# DETECCIÓN DE PULSOS INDIVIDUALES Y TRENES DE PULSOS

# Detección de envolvente y normalización
waveformabs=abs(waveform)  # Rectificamos
envelope=waveformabs*0
WINDOW=time2sample(Tgrillo)  # Abarcamos dos semiciclos
delta=round((WINDOW-1)/2)

LEN=length(waveform)
for (i in ((1+delta):(LEN-delta))) {
    envelope[i]=max(waveformabs[(i-delta):(i+delta)])
}

MAX=max(envelope)
envelope=envelope/MAX
waveform=waveform/MAX


# Detección de pulsos individuales
ThUP=0.2  # Umbrales
ThDOWN=0.2  # Recomendable ThUP>=ThDOWN
ThON=time2sample(0.05/4*0.8)  # s
ThOFF=time2sample(0.05/6*0.8)  # s

STATUS=0
pulse=as.data.frame(rbind(c(1, STATUS)))
colnames(pulse)=c('n','STATUS')  # n entero [1..LEN], STATUS=0/1

i=1
while (i<=LEN) {
    # (posible) subida
    if (STATUS==0 & envelope[i]>=ThUP) {
        iCandidato=i  # Instante del posible cambio
        iCount=1
        while(envelope[i+1]>ThDOWN & i+1<=LEN) {
            iCount=iCount+1
            i=i+1
        }
        if (iCount>=ThON) {
            STATUS=1
            pulse=rbind(pulse, c(iCandidato, STATUS)) 
        }
    }
    
    # (posible) bajada
    if (STATUS==1 & envelope[i]<=ThDOWN) {
        iCandidato=i  # Instante del posible cambio
        iCount=1
        while(envelope[i+1]<ThUP & i+1<=LEN) {
            iCount=iCount+1
            i=i+1
        }
        if (iCount>=ThOFF) {
            STATUS=0
            pulse=rbind(pulse, c(iCandidato, STATUS)) 
        }      
    }
  
    i=i+1
}
pulse=rbind(pulse, c(LEN, 0))  # Cerramos la secuencia con


# Detección de trenes de pulsos
cri=as.data.frame(rbind(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)))
colnames(cri)=c('n', 'N', 'D', 'Dpre', 'Dpost', 'T', 'dp', 'f', 'W',
    'duty', 'E', 'dutypre')
                  
pulsos=pulse[2:(nrow(pulse)-1),]  # Ignoramos los '0' extremos
LEN=nrow(pulsos)  # Siempre par
NUM=LEN/2  # Número de pulsos
ThINTER=3500  # Separación máxima de pulsos para considerarlos del mismo tren

iTren=1
i=1
while (i<=NUM-1) {
    N=1
    n=pulsos$n[i*2-1]  # Inicio del tren de pulsos
    noduty=0
    energy=sum( waveform[n:(pulsos$n[i*2])] ^ 2 )
    while(pulsos$n[(i+1)*2-1]-pulsos$n[i*2]<ThINTER & i<=NUM-1) {
        noduty=noduty+pulsos$n[(i+1)*2-1]-pulsos$n[i*2]
        energy=energy+
            sum( waveform[(pulsos$n[(i+1)*2-1]):(pulsos$n[(i+1)*2])] ^ 2 )
        N=N+1
        i=i+1
    }
    D=pulsos$n[i*2]-n
    dutypre=(D-noduty)/D
    W=energy/(D*dutypre)  # Promedio de energía emitida durante los pulsos
    dft=abs(fft(waveform[n:(n+D)]))
    f=which( round(dft)==max(round(dft)) )[1] * fs/length(dft)
    
    cri[iTren,]=c(n, N, D, NA, NA, NA, NA, f, W, NA, NA, dutypre)
    
    iTren=iTren+1
    i=i+1
}

# Dpre y Dpost
for (i in 1:(nrow(cri)-1)) {
    delay=cri$n[i+1]-(cri$n[i]+cri$D[i])
    cri$Dpost[i]=delay
    cri$Dpre[i+1]=delay
}

# Otras variables derivadas
cri$dp=cri$D*cri$dutypre/cri$N
cri$dsilence=cri$D*(1-cri$dutypre)/(cri$N-1)  # Var. auxiliar
cri$duty=cri$dp/(cri$dp+cri$dsilence)  # Versión de dc indepte. de N
cri$T=cri$dp/cri$duty
cri$E=cri$W*cri$dp
cri=subset(cri, select = -c(dutypre, dsilence))  # Eliminamos vars. auxiliares

# Descartamos trenes con algún NA (primero y último): 99 cri's -> 97 cri's
cri=cri[2:(nrow(cri)-1),]


# CORRELACIONES

# f en el t
plot(sample2time(cri$n), cri$f, col='red', xlab='Time (s)', ylab='f')
abline(lm(f ~ sample2time(cri$n), data=cri), lty='dotted')

# E y f
scatter.smooth(x=cri$E, y=cri$f, main="f ~ E", col='red', xlab='E', ylab='f (Hz)')

# Correlaciones entre todas las variables
col1=colorRampPalette(c("red", "darkred", "white", "darkgreen", "green"))
M=cor(cri)
M=cor(subset(cri, select=-c(N, D, Dpost, T, dp, duty, W)))
corrplot(M, order = "AOE", col = col1(200), addCoef.col = "black",
    tl.cex=1.5, diag=F)


# PREMIO DE CONSOLACIÓN

# Temperatura ambiente
# https://es.wikipedia.org/wiki/Gryllidae (grillo campestre)
Temperature=((60/sample2time(mean(cri$D+cri$Dpost))-40)/4+18)/1.8  # 21.5 ºC


# Visualización ad hoc
MEDIADpost=mean(cri$Dpost)
MEDIANAf=median(cri$f)
MINW=min(cri$W)

LADO=5000
img=NewBitmap(4000, 10000)
imgr=img
imgb=img
imgl=img

x=500
y=500
alpha=pi/2
i=1
R=(cri$W[i]-MINW)*1500+5
img =DrawCircle(img , x, y, r=R, fill=T)
imgr=DrawCircle(imgr, x, y, r=R, fill=T, val= max(0,(cri$f[i]-MEDIANAf)))
imgb=DrawCircle(imgb, x, y, r=R, fill=T, val=-min(0,(cri$f[i]-MEDIANAf)))
for (i in 2:nrow(cri)) {
    xpre=x
    ypre=y
    
    L=(cri$Dpost[i-1]-MEDIADpost)/100+100
    x=x+L*cos(alpha)
    y=y+L*sin(alpha)
    imgl=DrawLine(imgl, xpre, ypre, x, y)
    
    R=(cri$W[i]-MINW)*1500+5
    img= DrawCircle(img , x, y, r=R, fill=T)
    imgr=DrawCircle(imgr, x, y, r=R, fill=T, val= max(0,(cri$f[i]-MEDIANAf)))
    imgb=DrawCircle(imgb, x, y, r=R, fill=T, val=-min(0,(cri$f[i]-MEDIANAf)))
    
    alpha=alpha-(cri$dp[i]-cri$dp[i-1])/250*(2*pi)
}

img= SaveBitmap(img , gamma=2.2, trunc=F, "grillocamino")
imgr=SaveBitmap(imgr, gamma=2.2, trunc=F, "grillocaminor")
imgb=SaveBitmap(imgb, gamma=2.2, trunc=F, "grillocaminob")
imgb=SaveBitmap(imgl, gamma=2.2, trunc=T, "grillocaminol")
