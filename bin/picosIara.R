rm(list = ls())

#calcula a primeira derivada do gráfico
fderivada = function(df){
  derivada <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(derivada)<-c("x","y")
  for(i in seq(1,length(df$x)-1,1)){
    derivada[i,2] <- (df[i,2]-df[i+1,2])/
      (df[i,1]-df[i+1,1])
    derivada[i,1] <- df[i,1]
  }
  return(derivada)
}


#localiza picos do gráfico
achaPico = function (derivada){
  #compara os sinais da derivada. Havendo diferença houve inflexão no gráfico
  sinal<- sign(derivada$y)
  #variavel de resultado
  resultado<-data.frame(matrix(ncol = 2, nrow = 0))
  tmp<-data.frame(matrix(ncol = 2, nrow = 1))
  for(i in seq(1,length(sinal)-1,1)){
    if(sinal[i] != sinal[i+1]){
      #em tipo picos recebem 1, vales recebem -1
      if(sinal[i] > sinal[i+1]){
        tmp[1,2]<-1
      }else{
        tmp[1,2]<- -1
      }
      tmp[1,1]<-derivada$x[i]
      
      resultado<- rbind(resultado,tmp)
    }
  }
  colnames(resultado)<-c("x","tipo")
  return(resultado)
}

#corta picos e vales que não ultrapassem o threshold
aplicaTH= function (df, picos, limiarUp){
  #Cancelar picos e vales menores q limiarDw
  #Deve haver N picos e N+1 vales. Acrescentar no inicio e fim se preciso
  # Como valores são normalizados entre 0 e 1 colocar primeiro vale em -0.1 e ultimo em 1
  #Deve haver alternancia entre picos e vales tb. caso não haja, acusar erro de processamento
  #testar se pico é excede os vales adjacntes em limiarUp
  #caso contrário, cortar pico e vale posterior
  
  #cria variavel de retorno
  resultado<-data.frame(matrix(ncol = 3, nrow = 0))
  #df temporario juntando tabela original e picos
  tmp <- merge(df,picos, by= 1)
  colnames(tmp)<- c("x","y","tipo")
  #erro se tmp ficar vazio
  if(nrow(tmp)<=1){
    cat(paste("Nenhum pico e vale a ser processado... \n"))
    #print("Erro: ")
    colnames(resultado)<-c("x","y","tipo")
    return(resultado)
  }
  
  #acrescenta o primeiro vale, caso este não exista
  if(tmp$tipo[1] == 1){
    tmp<-rbind(c(-0.1,-1,-1),tmp)
  }
  #acrescenta o ultimo vale, caso este não exista
  if(tmp$tipo[nrow(tmp)] == 1){
    tmp<-rbind(tmp,c(max((tmp$x)+1000),-1,-1))
  }
  #verifica se existe alternancia entre picos e vales
  for(i in seq(1,nrow(tmp)-1,2)){
    #primeiro vale depois pico
    if(tmp$tipo[i]== 1 | tmp$tipo[i+1] == -1 ){
      cat(paste("Inconsistência no número de picos e vales... \n"))
      #      print("Erro: ")
      colnames(resultado)<-c("x","y","tipo")
      return(resultado)
    }
  }

  while(nrow(tmp) >= 3){
    #realiza verificação onde vale/pico/vales > limiarUp
    valeAnt<-tmp$y[1]+1
    picoAnt<-tmp$y[2]+1
    valePos<-tmp$y[3]+1
    if(nrow(tmp) >= 4){
      picoPos<-tmp$y[4]+1
    }else{
      picoPos<- 1
    }
    if(limiarUp < picoAnt - valeAnt & limiarUp < picoAnt - valePos){
      resultado <- rbind(resultado, tmp[2,],tmp[3,])
    }
    if(limiarUp < picoAnt - valePos & limiarUp < picoPos - valePos){
      resultado <- rbind(resultado, tmp[3,],tmp[4,])
    }
    tmp <- tmp[-c(1,2),]
  }
  colnames(resultado)<-c("x","y","tipo")
  resultado<-unique(resultado)
  return(resultado)
}


# Principal ----

  load(file="exemplo.RData")
  limiarUp = 0.025
  derivada<-fderivada(df)
  picos<-achaPico(derivada)
  picos<-aplicaTH(df = df,picos = picos,limiarUp = limiarUp )
  
  library(ggplot2)
  ggplot()+theme_bw()+
    geom_line(data = df, aes(x,y), col="blue")+
    geom_line(data = derivada, aes(x,y*100), col="red",lty=3)
