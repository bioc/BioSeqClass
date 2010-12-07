elements <- function(ele.type){
## Define the basic elements of biological sequence.
##
## ele.type - a string for the type of biological sequence.
##        "rnaBase" - basic elements of RNA.
##        "dnaBase" - basic elements of DNA.
##        "aminoacid" - 20 amino acides.
##        "aminoacid2" - 20 amino acides and 1 pseudo amino acid "O". Unknown or
##                       uncomplete amino acides will be substituted by pseudo 
##                       amino acid.
##
    e = switch(ele.type,
           "rnaBase" = list("A","T","C","G"),
           "dnaBase" = list("A","U","C","G"),
           "aminoacid" = list("R","K","E","D","Q","N","W","G","A","S","T","P",
                         "H","Y","C","V","L","I","M","F"),
           "aminoacid2" = list("R","K","E","D","Q","N","W","G","A","S","T","P",
                         "H","Y","C","V","L","I","M","F","O"),
           stop("parameter 'ele.type' has to be: rnaBase, dnaBase, aminoacid or aminoacid2")
        )
   names(e) = unlist(e)    
   return(e)                       
}
                       
aaClass <- function(aa.type){
## Group amino acids depend on their physical-chemical properties. 
##
## aa.type="aaH"
## Polar                        RKEDQN
## Neutral                      GASTPHY
## Hydrophobic                  CVLIMFW
##                              
## aa.type="aaV"                             
## Small                        GASCTPD
## Medium                       NVEQIL
## Large                        MHKFRYW
## 
## aa.type="aaZ"
## Low polarizability	        GASDT
## Medium polarizability	CPNVEQIL
## High polarizability	        KMHFRYW
## 
## aa.type="aaP"
## Low polarity	                LIFWCMVY
## Neutral polarity	        PATGS
## High polarity	        HQRKNED
## 
## aa.type="aaF"
## Acidic	                DE
## Basic	                HKR
## Polar	                CGNQSTY
## Nonpolar	                AFILMPVW
## 
## aa.type="aaS"
## Acidic	                DE
## Basic	                HKR
## Aromatic	                FWY
## Amide	                NQ
## Small hydroxyl	        ST
## Sulfur-containing	        CM
## Aliphatic	                AGPILV
## 
## aa.type="aaE"
## Acidic	                DE
## Basic	                HKR
## Aromatic	                FWY
## Amide	                NQ
## Small hydroxyl	        ST
## Sulfur-containing	        CM
## Aliphatic 1	                AGP
## Aliphatic 2	                ILV
##
## Ref: "Yu CS, Chen YC, Lu CH, Hwang JK. Prediction of protein subcellular localization. Proteins. 2006 Aug 15;64(3):643-51."
  class = switch(aa.type,
  	             "aaH" = list("polar"=c("R","K","E","D","Q","N"), 
  	                 "neutral"=c("G","A","S","T","P","H","Y"), 
  	                 "hydrophobic"=c("C","V","L","I","M","F","W")),
  	             "aaV" = list("small"=c("G","A","S","C","T","P","D"),
  	                 "medium"=c("N","V","E","Q","I","L"),
  	                 "large"=c("M","H","K","F","R","Y","W")),
  	             "aaZ" = list("lowPolarizability"=c("G","A","S","D","T"),
  	                 "mediumPolarizability"=c("C","P","N","V","E","Q","I","L"),
  	                 "highPolarizability"=c("K","M","H","F","R","Y","W")),
  	             "aaP" = list("lowPolarity"=c("L","I","F","W","C","M","V","Y"),
  	                 "neutralPolarity"=c("P","A","T","G","S"),
  	                 "highPolarity"=c("H","Q","R","K","N","E","D")),
  	             "aaF" = list("acidic"=c("D","E"),
  	                 "basic"=c("H","K","R"),
  	                 "polar"=c("C","G","N","Q","S","T","Y"),
  	                 "nonpolar"=c("A","F","I","L","M","P","V","W")),
  	             "aaS" = list("acidic"=c("D","E"),
  	                 "basic"=c("H","K","R"),
  	                 "aromatic"=c("F","W","Y"),
  	                 "amide"=c("N","Q"),
  	                 "smallHydroxyl"=c("S","T"),
  	                 "sulfurContaining"=c("C","M"),
  	                 "aliphatic"=c("A","G","P","I","L","V")),
  	             "aaE" = list("acidic"=c("D","E"),
  	                 "basic"=c("H","K","R"),
  	                 "aromatic"=c("F","W","Y"),
  	                 "amide"=c("N","Q"),
  	                 "smallHydroxyl"=c("S","T"),
  	                 "sulfurContaining"=c("C","M"),
  	                 "aliphatic1"=c("A","G","P"),
  	                 "aliphatic2"=c("I","L","V")),
  	)
  class	       
}   

featureBinary <- function(seq,class=elements("aminoacid")){
## Use 0-1 vector to code sequence.
##
## "aminoacid", each amino acid is represented by a 20-dimensional binary vector.
## "rnaBase", each base is represented by a 4-dimensional binary vector.
## "dnaBase", each base is represented by a 4-dimensional binary vector.
    L = unique(sapply(seq,function(x){length(unlist(strsplit(x,split="")))}))
    if( length(L)>1 ){
      stop("Sequences in seq must have equal length")
    }
  
    binary <- sapply(0:(length(class)-1), function(x){paste(c(rep(0,x), 1, 
                    rep(0,length(class)-x-1)),collapse="")} )
    names(binary) <- names(class)
    binary2 <- rep(binary,sapply(class,length))
    names(binary2) <- unlist(class)
    code <- t(sapply(seq,function(x){ as.numeric(
     unlist(sapply(binary2[unlist(strsplit(x,split=""))],strsplit,split=""))) }))
    colnames(code) <- paste("BIN:",paste(rep(1:L,each=length(class)),names(class),sep="_"),sep="")
    code
}

featureFragmentComposition <- function(seq,k,class=elements("aminoacid")){
## Cacluate the the frequency of a k-length sequence fragment.
##
## Ref: "Yu CS, Lin CJ, Hwang JK. Predicting subcellular localization of proteins for Gram-negative bacteria by support vector machines based on n-peptide compositions. Protein Sci. 2004 May;13(5):1402-6."
##      "Luo, R. Y., Feng, Z. P., and Liu, J. K. Prediction of protein structural class by amino acid and polypeptide composition. Eur J Biochem, 2002, 269 (17): 4219-4225."
##    
## k>=1
    e <- names(class)
    if(k>=2){
      for(i in (2:k) ){
        e <- paste(rep(e,each=length(class)),rep(names(class),length(class)))
      }
    }
    binary2 <- rep(names(class),sapply(class,length))
    names(binary2) <- unlist(class)
    composition <- sapply(seq, function(x){
                         nTotal = length(unlist(strsplit(x,split="")))-k+1
                         fragment = sapply(1:nTotal,function(i){substr(x,i,i+k-1)})                        
                         compo = rep(0,length(e))
                         names(compo) = e
                         tmp = table(sapply(fragment,function(y){
                               paste(binary2[unlist(strsplit(y,split=""))],collapse=" ")}))
                         compo[names(tmp)] = tmp
                         compo/nTotal
                         })
    composition <- t(composition)
    colnames(composition) <- sapply(colnames(composition), function(x){
                               paste("FC:",gsub(" ","_",x),sep="") 
                             })
    return(composition)
}   	                          

featureGapPairComposition <- function(seq,g,class=elements("aminoacid")){
## Cacluate the the frequency of a g-gap pair aminoacids/base.
##
## Ref: "Yu CS, Chen YC, Lu CH, Hwang JK. Prediction of protein subcellular localization. Proteins. 2006 Aug 15;64(3):643-51."
##
## g>=0
##
    e = paste(rep(names(class),each=length(class)),rep(names(class),length(class)))
    binary2 = rep(names(class),sapply(class,length))
    names(binary2) = unlist(class)
    composition = sapply(seq, function(x){
                         nTotal = length(unlist(strsplit(x,split="")))-g-1                         
                         compo = rep(0,length(e))
                         names(compo) = e
                         tmp = table(sapply(1:nTotal,function(i){
                         paste(binary2[c(substr(x,i,i),substr(x,i+g+1,i+g+1))],collapse=" ")}))
                         compo[names(tmp)] = tmp
                         compo/nTotal
                         })
    composition <- t(composition)
    tmp <- rep("X",g)
    colnames(composition) <- sapply(colnames(composition), function(x){
                               y = unlist(strsplit(x,split=" "))
                               paste("GPC:", paste(c(y[1],tmp,y[2]),collapse="_"), sep="") 
                             })
    return(composition)
}

featureHydro <- function(seq, hydro.method="SARAH1"){
## "kpm": The "hydrophobic effect" is believed to play a fundamental role in the
##        spontaneous folding of proteins. It can be expressed as the
##        free energy (kilocalories per mole) of transfer of amino
##        acid side chains from cyclohexane to water. The amino
##        acids with positive values of free energy in transferring
##        cyclohexane to water are hydrophobic and the ones with
##        negative values are hydrophilic.
##        I 4.92 
##        L 4.92 
##        V 4.04 
##        P 4.04 
##        F 2.98 
##        M 2.35 
##        W 2.33 
##        A 1.81 
##        C 1.28 
##        G 0.94
##        Y -0.14
##        T -2.57
##        S -3.40
##        H -4.66
##        Q -5.54
##        K -5.55
##        N -6.64
##        E -6.81
##        D -8.72
##        R -14.92
##        Ref: Radzicka A, Wolfenden R: Comparing the polarities of the
##             amino acids: Side-chain distribution coefficients between the
##             vapor phase, cyclohexane, 1-Octanol, and neutral aqueous
##             solution. Biochemistry 1988, 27:1664-1670.

## "SARAH1": SARAH1 assigns each amino acid a unique five-bit signed
##           code where exactly two bits are non-zero. SARAH1 ranks
##           twenty possible amino acids according to the Rose hydrophobicity
##           scale (Table 8). Each amino acid is assigned a
##           five bit code in descending order of the binary value of the
##           corresponding code. The ten most hydrophobic residues are positive,
##           and the ten least hydrophobic residues are negative.
##           C 1,1,0,0,0         
##           F 1,0,1,0,0
##           I 1,0,0,1,0
##           V 1,0,0,0,1
##           L 0,1,1,0,0
##           W 0,1,0,1,0
##           M 0,1,0,0,1
##           H 0,0,1,1,0
##           Y 0,0,1,0,1
##           A 0,0,0,1,1
##           G 0,0,0,-1,-1
##           T 0,0,-1,0,-1
##           S 0,0,-1,-1,0
##           R 0,-1,0,0,-1
##           P 0,-1,0,-1,0
##           N 0,-1,-1,0,0
##           D -1,0,0,0,-1
##           Q -1,0,0,-1,0
##           E -1,0,-1,0,0
##           K -1,-1,0,0,0
##           Ref: "Korenberg MJ, David R, Hunter IW, Solomon JE: Automatic 
##           classification of protein sequences into structure/function groups 
##           via parallel cascade identification: a feasibility study. Ann Biomed 
##           Eng 2000, 28(7):803-811."
    L = unique(sapply(seq,function(x){length(unlist(strsplit(x,split="")))}))
    if( length(L)>1 ){
      stop("Sequences in seq must have equal length")
    }
    
    if(hydro.method=="kpm"){
        H = c(4.92,4.92,4.04,4.04,2.98,2.35,2.33,1.81,1.28,0.94,-0.14,-2.57,-3.40,-4.66,-5.54,-5.55,-6.64,-6.81,-8.72,-14.92)
        names(H) = c("I","L","V","P","F","M","W","A","C","G","Y","T","S","H","Q","K","N","E","D","R")
        hydro = sapply(seq,function(x){
                       H[unlist(strsplit(x,split=""))]
                })
        rownames(hydro) = paste("H:",1:L,sep="")
    }
    if(hydro.method=="SARAH1"){
        H = list(c(1,1,0,0,0),c(1,0,1,0,0),c(1,0,0,1,0),c(1,0,0,0,1),c(0,1,1,0,0),c(0,1,0,1,0),c(0,1,0,0,1),c(0,0,1,1,0),c(0,0,1,0,1),c(0,0,0,1,1),c(0,0,0,-1,-1),c(0,0,-1,0,-1),c(0,0,-1,-1,0),c(0,-1,0,0,-1),c(0,-1,0,-1,0),c(0,-1,-1,0,0),c(-1,0,0,0,-1),c(-1,0,0,-1,0),c(-1,0,-1,0,0),c(-1,-1,0,0,0))
        names(H) = c("C","F","I","V","L","W","M","H","Y","A","G","T","S","R","P","N","D","Q","E","K")
        hydro = sapply(seq,function(x){
                       unlist(H[unlist(strsplit(x,split=""))])
                })    
        rownames(hydro) = paste("H:",paste(rep(1:L,each=5),1:5,sep="_"),sep="")
    }
    t(hydro)
}

featureACH <- function(seq,hydro.index="hydroE" ){
## Average cumulative hydrophobicity
## Ref: 
## Accurate sequence-based prediction of catalytic residues. Bioinformatics. 2008. 
## Kurgan,L. et al. (2007) Novel scales based on hydrophobicity indices for secondary
## protein structure. J. Theor. Biol., 248, 354每366.

  L = unique(sapply(seq,function(x){length(unlist(strsplit(x,split="")))}))
  if( length(L)>1 ){
    stop("Sequences in seq must have equal length")
  }
  if( !(L-1)%%2==0 ){
    stop("Sequences in seq must have be odd length")
  }
    
  hydro = switch(hydro.index,
                 # Eisenberg's hydrophobicity indices
                 "hydroE"=c(0.62,0.29,-0.9,-0.74,1.19,0.48,-0.4,1.38,-1.5,1.06,0.64,-0.78,0.12,-0.85,-2.53,-0.18,-0.05,1.08,0.81,0.26),                 
                 # Fauchere每Pliska's hydrophobicity indices
                 "hydroF"=c(0.42,1.34,-1.05,-0.87,2.44,0.00,0.18,2.46,-1.35,2.32,1.68,0.82,0.98,-0.30,-1.37,-0.05,0.35,1.66,3.07,1.31), 
                 # Cid's hydrophobicity indices
                 "hydroC"=c(0.17,1.24,-1.07,-1.19,1.29,-0.57,-0.25,2.06,-0.62,0.96,0.6,-0.9,-0.21,-1.2,-0.7,-0.83,-0.62,1.21,1.51,0.66), 
  )
  names(hydro) = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  ach <- t(sapply(seq,function(x){
    y=unlist(strsplit(x,split=""))
    sapply(1:((L-1)/2),function(w){mean(hydro[y[((L+1)/2-w):((L+1)/2+w)]])})
  }))
  colnames(ach) <- paste("ACH:",1:((L-1)/2),sep="")
  ach
}


featureCKSAAP <- function(seq,g,class=elements("aminoacid")){
## CKSAAP: the composition of k-spaced element pairs (0<=k<=g).
##
## g>=0
## 
## Ref: "Yong-Zi Chen, Yu-Rong Tang, Zhi-Ya Sheng and Ziding Zhang. Prediction 
##        of mucin-type O-glycosylation sites in mammalian proteins using the 
##        composition of k-spaced amino acid pairs."
  class2 = rep(names(class),sapply(class,length) )
  names(class2) = unlist(class)
  couple = paste(rep(names(class),each=length(class)*(g+1)), 
           rep(0:g,each=length(class),),
           names(class), sep="")

  coupleCount = t(sapply(seq,function(x){
                       y=unlist(strsplit(x,split=""))
                       pair=unlist(sapply( 1:(length(y)-1), function(i){
                              j=0:min(length(y)-i-1,g) ; 
                              paste(class2[y[i]],j,class2[y[i+j+1]],sep="") 
                            }))
                       tmp = table(pair)     
                       c = rep(0,length(couple))
                       names(c) = couple
                       c[names(tmp)] = tmp
                       c
                }))
  colnames(coupleCount) = paste("CKSAAP:", colnames(coupleCount), sep="")
  coupleCount
}

featureAAindex <- function(seq,aaindex.name="all"){
## AAindex (http://www.genome.jp/aaindex) is a database of numerical indices 
## representing various physicochemical and biochemical properties of amino acids 
## and pairs of amino acids. AAindex consists of three sections now: 
## AAindex1 for the amino acid index of 20 numerical values, 
## AAindex2 for the amino acid mutation matrix and 
## AAindex3 for the statistical protein contact potentials.
##
## Ref: "Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., and Kanehisa, M.; AAindex: amino acid index database, progress report 2008. Nucleic Acids Res. 36, D202-D205 (2008)."
##    
    L = unique(sapply(seq,function(x){length(unlist(strsplit(x,split="")))}))
    if( length(L)>1 ){
      stop("Sequences in seq must have equal length")
    }
    
    data(aa.index)
    if(aaindex.name=="all"){
      index = sapply(aa.index, function(x){x$I})
      index = index[,apply(index,2,function(x){sum(is.na(x))==0})]
      name = colnames(index)
      indexFeature = sapply(seq, function(x){
              	            x2 = unlist(strsplit(x,split=""))
                            x3 = as.vector(apply(index,2,function(y){ y[x2] }) )          
                            names(x3) = paste(rep(name,each=length(x2)), rep(1:length(x2),ncol(index)),sep="_")
                            x3
                     })
    }else{
      index = aa.index[[aaindex.name]]$I
      if( sum(is.na(index))>0 ){
        stop(paste("aaindex.name",aaindex.name,"is not the property name in 
        AAindex or has NA value"))
      }
      indexFeature = sapply(seq, function(x){
          	            x2 = unlist(strsplit(x,split=""))
                            x3 = index[x2]
                            names(x3) = paste(aaindex.name, 1:length(x2),sep="_")
                            x3
                     })
    }
    indexFeature = t(indexFeature)
    colnames(indexFeature) = paste("AAindex:",colnames(indexFeature),sep="")
    indexFeature
}

featureACI <- function(seq,aaindex.name="all"){
  L = unique(sapply(seq,function(x){length(unlist(strsplit(x,split="")))}))
  if( length(L)>1 ){
    stop("Sequences in seq must have equal length")
  }
  if( !(L-1)%%2==0 ){
    stop("Sequences in seq must have be odd length")
  }
  data(aa.index)
  
  if(aaindex.name=="all"){
    index = sapply(aa.index, function(x){x$I})
    index = index[,apply(index,2,function(x){sum(is.na(x))==0})]
    name = colnames(index)
    aci = t(sapply(seq, function(x){
             y = unlist(strsplit(x,split=""))
              z= as.vector(apply(index,2,function(i){sapply(1:((L-1)/2),function(w){mean(i[y[((L+1)/2-w):((L+1)/2+w)]])}) }))
              names(z) = paste(rep(name,each=(L-1)/2),1:((L-1)/2),sep="_")
              z
    }))
  }else{
    index = aa.index[[aaindex.name]]$I
    if( sum(is.na(index))>0 ){
      stop(paste("aaindex.name",aaindex.name,"is not the property name in 
      AAindex or has NA value"))
    }
    aci <- t(sapply(seq,function(x){
      y=unlist(strsplit(x,split=""))        
      z=sapply(1:((L-1)/2),function(w){mean(index[y[((L+1)/2-w):((L+1)/2+w)]])})
      names(z)=paste(aaindex.name, 1:((L-1)/2),sep="_")
      z
    }))
  }  
  colnames(aci) <- paste("ACI:",colnames(aci),sep="")
  aci
}

featureACF <- function(seq,n,aaindex.name="all"){
## L-2>=n>=1, L is the the length of sequence X.
## Auto-Correlation between fragment X(1)...X(L-m) and X(m+1)...X(L) (1<=m<=n).
##
## Ref: "Zhang, C. T., Lin, Z. S., Zhang, Z., and Yan, M. Prediction of the 
##      helix/strand content of globular proteins based on their primary 
##      sequences. Protein Eng, 1998, 11 (11): 971-979."
##      "Liu, Wei-min and Chou, Kuo-Chen. Prediction of Protein Structural 
##      Classes by Modified Mahalanobis Discriminant Algorithm. Journal of Protein 
##      Chemistry, 1998, 17 (3): 209-217."
##      
  data(aa.index)
  if(aaindex.name=="all"){
    index = sapply(aa.index, function(x){x$I})
    index = index[,apply(index,2,function(x){sum(is.na(x))==0})]
    name = colnames(index)
    indexFeature = sapply(seq, function(x){         
       		          x2 = unlist(strsplit(x,split=""))
       		          L = length(x2)
                          x3 = as.vector(apply(index,2,function(y){ 
                                         sapply(1:n, function(m){                                             
                                                mean(y[x2[1:(L-m)]]*y[x2[(m+1):L]])
                                         })
                               }))         
                          x3
                   })
    indexFeature = t(indexFeature)
    colnames(indexFeature) = paste("ACF:", paste(rep(name,each=n),1:n,sep="_" ), sep="")
  }else{
    index = aa.index[[aaindex.name]]$I       
    indexFeature = sapply(seq, function(x){         
       		          x2 = unlist(strsplit(x,split=""))
       		          L = length(x2)
                          x3 = sapply(1:n, function(m){                                             
                                      mean(index[x2[1:(L-m)]]*index[x2[(m+1):L]])
                               })                                      
                          x3
                   })    
    if( is.null(nrow(indexFeature)) ){
      indexFeature = matrix(indexFeature,nrow=1)
    }
    indexFeature = t(indexFeature)
    colnames(indexFeature) = paste("ACF:", paste(rep(aaindex.name,each=n),1:n,sep="_" ), sep="")
  } 
  indexFeature           
}

##featureCenterACF <- function(seq,name="all"){
## L is the the length of sequence X, and it is even.
## Auto-Correlation between fragment X(1)...X(L/2) and X(L)...X(L/2+1).
##  
##  data(aa.index)
##  if(name=="all"){
##    index = sapply(aa.index, function(x){x$I})
##    index = index[,apply(index,2,function(x){sum(is.na(x))==0})]
##    name = colnames(index)
##    indexFeature = sapply(seq, function(x){         
##       		          x2 = unlist(strsplit(x,split=""))
##                          L = length(x2)
##                          x3 = apply(index,2,function(y){                                
##                                      mean(y[x2[1:(L/2)]]*y[x2[L:(L/2+1)]])  
##                                      #sum((y[x2[1:(L/2)]]-y[x2[L:(L/2+1)]])^2)
##                                      #y[x2[1:(L/2)]]-y[x2[L:(L/2+1)]]                             
##                               })    
##                          #names(x3) = name
##                          x3
##                   })
##  }else{
##    index = aa.index[[name]]$I
##    indexFeature = sapply(seq, function(x){         
##    	          x2 = unlist(strsplit(x,split=""))
##                      L = length(x2)
##                      x3 = mean(index[x2[1:(L/2)]]*index[x2[L:(L/2+1)]])
##                      names(x3) = name
##                      x3
##               })  
##  }
##  t(indexFeature)
##}

featurePseudoAAComp <- function(seq,d,w=0.05){
## pseudo-amino acid composition: 
## Ref: "Chou, K. C. Prediction of protein cellular attributes using pseudo-amino 
##       acid composition. Proteins, 2001, 43 (3): 246-255."
##
## d>=1
##
    # H1 is the hydrophobicity value.
    # Ref: "Tanford C. Contribution of hydrphobic interactions to the stability of the globular conformation of proteins. J Am Chem Soc 1962; 84: 4240-4274."
    # http://www.expasy.ch/tools/pscale/Hphob.Tanford.html
    H1 = c(0.620,-2.530,-0.780,-0.090,0.290,-0.850,-0.740,0.480,-0.400,1.380,1.530,-1.500,0.640,1.190,0.120,-0.180,-0.050,0.810,0.260,1.800)   
    names(H1) = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")    
    # H2 is the hydrophilicity value.
    # Ref: "Hopp TP, Woods KR. Prediction of protein antigenic determinants from amino acid sequences. Proc Natl Acad Sci USA 1981; 78: 3824-3828."
    # AAindex: HOPT810101 Hydrophilicity value (Hopp-Woods, 1981)
    H2 = c(-0.5,3.0,0.2,3.0,-1.0,0.2,3.0,0.0,-0.5,-1.8,-1.8,3.0,-1.3,-2.5,0.0,0.3,-0.4,-3.4,-2.3,-1.5)
    names(H2) = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V") 
    # M is the the mass of amino acid side chain. (Here, M is the mass of residue).
    M = c(71.03711,156.10111,114.04293,115.02694,103.00919,128.05858,129.04259,57.02146,137.05891,113.08406,113.08406,128.09496,131.04049,147.06841,97.05276,87.03203,101.04768,186.07931,163.06333,99.06841)
    names(M) = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V") 
    H1 = (H1-mean(H1))/sd(H1)
    H2 = (H2-mean(H2))/sd(H2)
    M = (M-mean(M))/sd(M)
    R <- function(aa1,aa2){
    	((H1[aa1]-H1[aa2])^2+(H2[aa1]-H2[aa2])^2+(M[aa1]-M[aa2])^2)/3
    }

    compo = featureFragmentComposition(seq,1,elements("aminoacid"))
    pseudo = sapply(seq, function(x){           
                    y = unlist(strsplit(x,split=""))
                    nTotal = length(y) 
                    sapply(1:d, function(i){
                           mean( sapply(1:(nTotal-i), function(j){
                                  R(y[j],y[j+i])
                                 }) )
                    })
             })
    pseudo = cbind(compo,matrix(unlist(pseudo),ncol=d))
    pseudo = t(apply(pseudo,1,function(x){ y=sum(x[1:20])+w*sum(x[21:length(x)]); 
                                           x[1:20] = x[1:20]/y;
                                           x[21:length(x)] = w*x[21:length(x)]/y;
                                           x
              }))
    colnames(pseudo) = c(sub("FC:", "PAC:", colnames(pseudo)[1:20]), paste("PAC:",1:d,sep=""))
    pseudo
}

getDSSP <- function(pdb){
## DSSP is a database of secondary structure assignments (and much more) for all
## protein entries in the Protein Data Bank (PDB)
## Download: rsync -avz rsync://rsync.cmbi.ru.nl/dssp/ /path_to_dssp/
## Ftp: ftp://ftp.cmbi.kun.nl/pub/molbio/data/
##
## The very short explanation of the DSSP result: 
##    H = alpha helix
##    B = residue in isolated beta-bridge
##    E = extended strand, participates in beta ladder
##    G = 3-helix (3/10 helix)
##    I = 5 helix (pi helix)
##    T = hydrogen bonded turn
##    S = bend 
##    A blank in the DSSP secondary structure determination stands for loop or 
##    irregular. Loops and irregular elements are often, very incorrectly, called 
##    "random coil" or "coil". Here, we replace this blank by a C.
##
## Ref: "Dictionary of protein secondary structure: pattern recognition of 
##      hydrogen-bonded and geometrical features. Biopolymers. 1983 Dec;22(12):
##      2577-637. PMID: 6667333; UI: 84128824. "
##

  data(dssp.ss)  
  tmp = is.na(dssp.ss[pdb])
  if( sum(tmp)>0 ){
    warning(paste(pdb[tmp],"do not have secondary structure in DSSP"))
  }
  pdb = pdb[!tmp]
  ss = dssp.ss[pdb]
  ss
}

##predict.ss <- function(seq, ss.method,
##  psipred.program="psipred", psipred.filter="complex", 
##  proteus2.organism="euk"
##){
## Predict secondary structure. 
##  
## seq - a chracter for protein sequence.
## 
## ss.method - a character for secondary structure prediction method:
##           "PSIPRED": http://bioinf.cs.ucl.ac.uk/psipred/psiform.html
##                      Ref: "Bryson K, McGuffin LJ, Marsden RL, Ward JJ, Sodhi JS, 
##                      Jones DT: Protein structure prediction servers at University 
##                      College London. Nucl. Acids Res. 2005, 33(Web Server issue):W36-38."
##           "JNET": http://www.compbio.dundee.ac.uk/~www-jpred/index.html 
##                   Ref: "Cuff JA, Barton GJ: Application of multiple sequence 
##                   alignment profiles to improve protein secondary structure
##                   prediction. Proteins 2000, 15:502-11."
##           "PROTEUS2": http://wks16338.biology.ualberta.ca/proteus2/index.jsp
##                       Ref: "Montgomerie S, Sundararaj S, Gallin WJ, Wishart DS: 
##                       Improving the accuracy of protein secondary structure 
##                       prediction using structural alignment. BMC Bioinformatics 2006, 14:301."
##
## read.ss returns a list for secondary structure prediction results. It contains
## three elements: "sequence" for the protein sequence; "secondaryStructure"
## for the protein secondary structure prediction  (H = Helix, E = Beta Strand, 
## C = Coil); "score" for the confidence score (0-9, 0 = low, 9 = high).
##             
##
## if(is.null(names(seq))){
##    names(seq) <- paste("seq",1:length(seq),sep="_")
##  }
##  ssPredict <- switch(ss.method,
##        "PSIPRED" = return(predict.psipred(seq,psipred.program,psipred.filter)),
##        "JNET" = return(predict.jnet(seq)), 
##        "PROTEUS2" = return(predictPROTEUS(seq,proteus2.organism)), 
##        stop(paste("ss.method", ss.method," is not supported. Has to be: 
##        PSIPRED, JNET, or PROTEUS2"))
##  )
##}
##
##predict.psipred <- function(seq,psipred.program="psipred",psipred.filter="complex"){
##  OS <- .Platform$OS.type
##  ssPredict <- lapply(names(seq), function(x){
##    tmp <- tempfile()  
##    perlName <- paste(tmp, "pl", sep=".")
##    outName <- paste(tmp, "psipred", sep=".")
##    .pathPerl(perlName,.Platform$OS.type) 
##    write(paste("$Program='",psipred.program,"';",sep=""), file = perlName, append = TRUE)
##    write(paste("$Filter='",psipred.filter,"';",sep=""), file = perlName, append = TRUE)
##    write(paste("$Sequence='",seq[x],"';",sep=""), file = perlName, append = TRUE)    
##    write(paste("$name='",x,"';",sep=""), file = perlName, append = TRUE)  
##    write(paste("open(OUT,'>",outName,"');",sep=""), file = perlName, append = TRUE)
##    write(readLines( file.path(.path.package("BioSeqClass"),"scripts","psipred") ), 
##      file = perlName, append = TRUE) 
##    .callPerl(perlName, OS)
##    as.list(read.csv(outName,header=T,sep="\t"))    
##  })
##}
##
##predict.jnet <- function(seq){
##  OS <- .Platform$OS.type
##  ssPredict <- lapply(seq, function(x){
##    tmp <- tempfile()  
##    perlName <- paste(tmp, "pl", sep=".")
##    outName <- paste(tmp, "jnet", sep=".")
##    .pathPerl(perlName,.Platform$OS.type)    
##    write(paste("$seq='",x,"';",sep=""), file = perlName, append = TRUE)    
##    write(paste("open(OUT,'>",outName,"');",sep=""), file = perlName, append = TRUE)
##    write(readLines( file.path(.path.package("BioSeqClass"),"scripts","jnet") ), 
##      file = perlName, append = TRUE) 
##    .callPerl(perlName, OS)
##    as.list(read.csv(outName,header=T,sep="\t"))    
##  })
##}

predictPROTEUS <- function(seq,proteus2.organism="euk"){
## proteus2.organism="gram-" for "Gram negative prokaryote"
## proteus2.organism="gram+" for "Gram positive prokaryote"
## proteus2.organism="euk" for "Eukaryote"
  OS <- .Platform$OS.type
  ssPredict <- lapply(seq, function(x){
    tmp <- tempfile()  
    perlName <- paste(tmp, "pl", sep=".")
    outName <- paste(tmp, "proteus2", sep=".")
    .pathPerl(perlName,.Platform$OS.type)  
    write(paste("$organism='",proteus2.organism,"';",sep=""), file = perlName, append = TRUE)
    write(paste("$sequence='",x,"';",sep=""), file = perlName, append = TRUE)    
    write(paste("open(OUT,'>",outName,"');",sep=""), file = perlName, append = TRUE)
    write(readLines( file.path(.path.package("BioSeqClass"),"scripts","proteus2") ), 
      file = perlName, append = TRUE) 
    .callPerl(perlName, OS)    
    data <-as.matrix(read.csv(outName,header=F,sep="\t",row.names=1))
    result <- list()
    tmpR <- unlist(strsplit(data["Residue",1],split="")) 
    for( i in rownames(data)[-1] ){
      tmpN <- unlist(strsplit(i,split="_")) 
      result[[tmpN[1]]][[tmpN[2]]] <- unlist(strsplit(data[i,1],split=""))       
      result[[tmpN[1]]][["Residue"]] <- tmpR
    }
    result
  })
  ssPredict
}

featureSSC <- function(secondaryStructure, confidenceScore){
## Coding for for the secondary structure of peptides. It is 
## suitable for peptides with odd residues and the central residue has important 
## role.
##
## Coding scheme for secondary structure. Four scheme are supported: 
## "Binary": binary value denoting secondary structure for the central residue.
## "SCORE": normalized confidence score. confidence value (division by 10 results 
##            in normalization to a unit interval) obtained for the central residue 
##            by the SSP prediction method
## "PATTERN": binary value denoting a specific configuration of the secondary 
##              structure for the central and the two adjacent residues.   
## Ref: "Ce Zheng, Lukasz Kurgan. Prediction of beta-turns at over 80% accuracy 
##       based on an ensemble of predicted secondary structures and multiple 
##       alignments. BMC Bioinformatics 2008, 9:430"
##
## seq - a named string vecotr for the peptide sequences.
## ss - a list for secondary structure prediction results. It contains
##             three elements: "Residue" for the protein sequence; "SecondaryStructure"
##             for the protein secondary structure prediction  (H = Helix, 
##             E = Beta Strand, C = Coil); "confidenceScore" for the confidence score 
##             (0-9, 0 = low, 9 = high).
##             
##  
##  'secondaryStructure' is a string vector for the protein secondary structure.
##  It is consisted of three kinds of secondary structures: H = Helix, 
##  E = Beta Strand, C = Coil.        
## 'confidenceScore' is a string vector for the confidence score of secondary 
## structure prediction (0-9, 0 = low, 9 = high).
 
   name <- c( paste("SS: Binary_central_", c("C","H","E"), sep=""),
              paste("SS: PATTERN_central_",unlist(sapply( c("C","H","E"), function(y){ 
                    c(paste(y,y,y,sep=""),paste(y,y,"X",sep=""),
                    paste("X",y,y,sep=""),paste("X",y,"X",sep="")) } ) ), sep="") 
           ) 
   tmp1 <- sapply(1:length(secondaryStructure),function(x){
               tmpS <- unlist(strsplit(secondaryStructure[x],split=""))
               z = rep(0,length(name))
               names(z) = name       
               middle = (length(tmpS)+1)/2
               z[paste("SS: Binary_central_", tmpS[middle], sep="")] <- 1         
               tmpP <- tmpS[(middle-1):(middle+1)]
               tmpP[tmpP!=tmpS[middle+1]] <- "X"
               z[paste("SS: PATTERN_central_", paste(tmpP,collapse=""), sep="")] <- 1
               z
           })
   if( is.null(confidenceScore) ){
     ssFeature <- t(tmp1)
   }else{
      if( max(as.numeric(as.vector((sapply(confidenceScore,function(x){unlist(strsplit(x,split=""))})))))>9 ){
        stop(paste("Parameter confidenceScore should ranged from 0 to 9"))
      }
      if( length(secondaryStructure)!=length(confidenceScore) ){
        stop(paste("Parameter secondaryStructure and confidenceScore should have 
        equal length"))
      }
     tmp2 <- sapply(1:length(secondaryStructure),function(x){
                 tmpC <- unlist(strsplit(confidenceScore[x],split=""))
                 middle = (length(tmpC)+1)/2
                 as.numeric(tmpC[middle])/10 
             })
     ssFeature <- cbind(t(tmp1),"SS: SCORE_central"=tmp2)        
   }  
  ssFeature
}

##predict.sasa <- function(seq,method,asap.model="1"){
## predict solvent accessible surface area (SASA). 
## 
## seq - a character for protein sequence.
## 
## method - a character for SASA prediction program:
##          "Naccess": a stand alone program that calculates the accessible area 
##                     of a molecule from a PDB (Protein Data Bank) format file.  
##                     Link: http://www.bioinf.manchester.ac.uk/naccess/
##                     Ref: "Lee, B. and Richards, F.M. 1971. The interpretation 
##                          of protein structures: Estimation of static 
##                          accessibility. J. Mol. Biol. 55: 379每400."
##          "ASAP": Protein Solvent Accessible Surface Area Predictor.
##                     Link: http://ccb.imb.uq.edu.au/ASAP/
##                     Ref: 
## asap.model - a character string for model used by "ASAP" method.
## asap.model="1" for "Soluble Protein"
## asap.model="2" for "Transmembrane Helix Protein"
## asap.model="4" for "Transmembrane	Beta-barrel Protein"
##
##
## if(is.null(names(seq))){
##   names(seq) <- paste("seq",1:length(seq),sep="_")
## }
## sasaPredict <- switch(method,
##                  "ASAP" = return(predictASAP(seq ,asap.model)),
##                  stop(paste("method", method," is not supported. Has to be: ASAP"))
##                )
##}


predictASAP <- function(seq, asap.model="1"){
  OS <- .Platform$OS.type
  sasaPredict <- lapply(seq, function(x){
    tmp <- tempfile()  
    perlName <- paste(tmp, "pl", sep=".")
    outName <- paste(tmp, "asap", sep=".")
    .pathPerl(perlName,.Platform$OS.type)  
    write(paste("$model='",asap.model,"';",sep=""), file = perlName, append = TRUE)
    write(paste("$seq='",x,"';",sep=""), file = perlName, append = TRUE)    
    write(paste("open(out,'>",outName,"');",sep=""), file = perlName, append = TRUE)
    write(readLines( file.path(.path.package("BioSeqClass"),"scripts","asap") ), 
      file = perlName, append = TRUE) 
    .callPerl(perlName, OS)
    as.list(read.csv(outName,header=T,sep="\t"))
  }) 
  sasaPredict 
}


##predict.ep <- function(seq,method){
## predict electrostatic potential of proteins.
##
## method - a character for electrostatic potential prediction method.
## "PCE": Protein Continuum Electrostatics performs the calculation of the electrostatic 
##         potentials for a protein by solving numerically the Poisson-Boltzmann 
##         equation (the Finite Difference Poisson- Boltzmann method, FDPB). It is 
##         a server adaptation of the MEAD potential program (MEAD: Macroscopic 
##         Electrostatics with Atomic Detail, D. Bashford). 
##         Link: http://bioserv.rpbs.jussieu.fr/cgi-bin/PCE-Pot
##         Ref: "Miteva, M.A., Tuffe∩ry, P., and Villoutreix, B.O. 2005. PCE: Web 
##              tools to compute protein continuum electrostatics. Nucleic Acids Res. 
##              33: W372每W375. doi: 10.1093/nar/gki365.".
##}

##predict.pKa <- function(seq,method){
## predict protein pKa.
##
## method - a character for protein pKa prediction software: 
## "PROPKA": Link: http://propka.ki.ku.dk/~drogers/.
##           Ref: "Li, H., Robertson, A.D., and Jensen, J.H. 2005. Very fast 
##                 empirical prediction and rationalization of protein pKa values. 
##                 Proteins 61: 704每721."

##}


featurePSSM <- function(seq, start.pos, stop.pos, psiblast.path, database.path){
##
## PSSM: a position-specific scoring matrix generated by 
##       PSI-BLAST(Position-Specific Iterated BLAST), representing 
##       the likelihood of that particular residue substitution at that position. 
##       It is a matrix of 20℅N elements, where N is the size of the window.  
##       Highly conserved positions receive high scores and weakly conserved 
##       positions receive scores near zero.
##
## WOP: a frequency distribution of 20 amino acids at that sequence position.
##
## PSI-BLAST:
##           Ref:"Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, 
##           Lipman DJ: Gapped BLAST and PSI-BLAST: a new generation of protein 
##           database search programs. Nucleic Acids Res 1997, 25:3389-3402."
##
## Ref: "Qidong Zhang, Sukjoon Yoon and William J. Welsh. Improved method for 
##      predicting beta-turn using support vector machine. Bioinformatics 2005 
##      21(10):2370-2374."
##      "Chung-Tsai Su, Chien-Yu Chen and Yu-Yen Ou. Protein disorder prediction 
##      by condensed PSSM considering propensity for order or disorder. BMC 
##      Bioinformatics 2006, 7:319"
##      "Ahmad S, Sarai A. PSSM-based prediction of DNA binding sites in proteins. 
##      BMC Bioinformatics. 2005 Feb 19;6:33."
##
  if( any(length(seq)!=length(start.pos),length(stop.pos)!=length(start.pos)) ){
    stop("Parameters seq, start.pos, and stop.pos must have equal length")
  }
  N = unique(stop.pos-start.pos+1)
  if( length(N)>1 ){
    stop("Sequence fragments determined by 'start.pos' and 'stop.pos' must have equal length")
  }
  if( !file.exists(paste(database.path,".pin",sep="")) ){
    stop(paste("Do not have formated database file:",paste(database.path,".pin",sep="")))
  }
  result = t(sapply(1:length(seq),function(x){
                  tmp = tempfile("tmp")
                  inFasta = paste(tmp,"fasta",sep=".") ;                
                  out1 = paste(tmp,"psiout",sep=".") ;     
                  out2 = paste(tmp,"pssm",sep=".") ;     
                  outputSeq = list()
                  outputSeq[[1]]=list("desc"=names(seq)[x],"seq"=seq[x])
                  writeFASTA(outputSeq,inFasta)  
                  system(paste(psiblast.path, "-d", database.path, "-i",inFasta,"-j 3 -o ", out1, "-Q", out2))
                  data = readLines(out2)                
                  pssm = sapply(start.pos[x]:stop.pos[x],function(y){                   
                           z = unlist(strsplit(data[y+3],split="\\s+"))                  
                           # PSSM will be normalized to <0, 1> interval.
                           1/(1+exp(-as.numeric(z[4:23])))
                         })
                  
                  ent = sapply(start.pos[x]:stop.pos[x],function(y){       
                          z = unlist(strsplit(data[y+3],split="\\s+"))   
                          p = as.numeric(z[24:43])
                          p = as.vector(table(p))/length(p)
                          p = p[p>0]
                          sum(-p*log(p))
                        })            
                  c(as.vector(pssm),as.vector(ent))
  }))
  aa = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")    
  colnames(result) = c(paste("PSSM:",paste(rep(1:N,each=20),aa,sep="_"),sep=""), 
    paste("EntWOP:",1:N,sep="") )
  result  
}

predictPFAM <- function(seq, hmmpfam.path, pfam.path, Evalue=10^-5){
  if( !file.exists(hmmpfam.path) ){
    stop(paste("hmmpfam","is not existed in",hmmpfam.path))
  }
  if( !file.exists(pfam.path) ){
    stop(paste("Pfam database","is not existed in",pfam.path))
  }
  perlName <- paste(tempfile("tempPerl"), "pl", sep=".")
  .pathPerl(perlName,.Platform$OS.type)    
  inFasta = tempfile("tempFasta") ;
  outFile = paste(inFasta,"pfam",sep=".") ;  
  outputSeq = list()                  
  for (x in 1:length(seq) ){outputSeq[[x]]=list("desc"=names(seq)[x],"seq"=seq[x]) }
  writeFASTA(outputSeq,inFasta)    
  write(paste("$inFasta='",inFasta,"';",sep=""), file = perlName, append = TRUE)     
  write(paste("$cutoff=",Evalue,";",sep=""), file = perlName, append = TRUE)       
  write(paste("$path='",hmmpfam.path,"';",sep=""), file = perlName, append = TRUE)   
  write(paste("$database='",pfam.path,"';",sep=""), file = perlName, append = TRUE)   
  write(readLines( file.path(.path.package("BioSeqClass"),"scripts",
    "hmmpfam") ), file = perlName, append = TRUE) 
  .callPerl(perlName, .Platform$OS.type) 
  tmp = as.matrix(read.csv(outFile,header=T,sep="\t"))
  tmp = matrix(tmp,ncol=5)
  domain = as.list(tapply(1:nrow(tmp),tmp[,1],function(x){tmp[x,4]}))
  domain = sapply(domain,function(x){y=table(x); z=as.vector(y); names(z)=names(y); z}) 
  
  unlink(perlName)
  unlink(file.path(dirname(inFasta),dir(path=dirname(inFasta),pattern=basename(inFasta))) ) 
  domain
}

featureDOMAIN <- function(domain){
  cols <- unique(unlist(sapply(domain,names)))
  m <- matrix(0,nrow=length(domain),ncol=length(cols))  
  rownames(m) <- names(domain)
  colnames(m) <- cols
  for( i in rownames(m) ){ m[i, names(domain[[i]]) ]=domain[[i]] }
  colnames(m) <- paste("DOMAIN",colnames(m),sep=":")
  m
}

featureCTD <- function(seq, class=elements("aminoacid") ){
# Three descriptors, composition (C), transition (T) and
# distribution (D), are used to describe global composition of
# each of these properties (33). C is the number of amino acids
# of a particular property (such as hydrophobicity) divided by
# the total number of amino acids. T characterizes the percent
# frequency with which amino acids of a particular property is
# followed by amino acids of a different property. D measures
# the chain length within which the first, 25, 50, 75 and 100% of
# the amino acids of a particular property is located respectively.

# Ref: "C.Z. Cai1,, L.Y. Han1, Z.L. Ji, X. Chen and Y.Z. Chen. SVM-Prot: web-based support vector machine
# software for functional classification of a protein
# from its primary sequence. Nucleic Acids Research, 2003, Vol. 31, No. 13"

  k = 0 ;
  pair = vector();
  for(i in 1:(length(class)-1) ){
    for(j in (i+1):length(class) ){
      k = k+1
      pair[k] <- paste(names(class)[i],names(class)[j],sep="_")
    }
  }            
  name <- c( paste("CTD:C",names(class),sep="_"), 
           paste("CTD:T",pair,sep="_"),
           paste("CTD:D",rep(names(class),each=5),c("1st","25%","50%","75%","100%"),sep="_") )

  binary2 <- rep(names(class),sapply(class,length))
  names(binary2) <- unlist(class)
  ctd <- sapply(seq,function(x){
           y = binary2[unlist(strsplit(x,split=""))]
           z = rep(0,length=length(name))
           names(z) = name
           tmp = table(y)/length(y)
           z[paste("CTD:C",names(tmp),sep="_")] = tmp
           tmp1 = table(sapply(1:(length(y)-1),function(i){paste(y[i],y[i+1],sep="_")}))
           tmp2 = table(sapply(length(y):2,function(i){paste(y[i],y[i-1],sep="_")}))
           if( length(intersect(pair,names(tmp1)))>0 ){
             z[paste("CTD:T",intersect(pair,names(tmp1)),sep="_")] = tmp1[intersect(pair,names(tmp1))]/sum(tmp1)
           }
           if( length(intersect(pair,names(tmp2)))>0 ){
             z[paste("CTD:T",intersect(pair,names(tmp2)),sep="_")] = z[paste("CTD:T",intersect(pair,names(tmp2)),sep="_")] + tmp2[intersect(pair,names(tmp2))]/sum(tmp2)
           }
           for( c in unique(y) ) { 
             tmp = (1:length(y))[y==c]
             z[paste("CTD:D",rep(c,each=5),c("1st","25%","50%","75%","100%"),sep="_")] = c(tmp[1],tmp[ceiling(length(tmp)*c(0.25,0.5,0.75,1))])/length(y) 
           }
           z
         })
  t(ctd)
}

featureBDNAVIDEO <- function(seq){
## B-DNA-VIDEO: Ref: http://wwwmgs.bionet.nsc.ru/mgs/systems/bdnavideo/
## "PROPERTY" in B-DNA-VIDEO : A distributed and intelligent database for the activities of the functional sites in DNA and RNA.  The current release has 38 entries and was indexed 22-Feb-2007. 

  data(PROPERTY) 
  result <- t(sapply(toupper(seq),function(x){
                z <- sapply(1:(length(unlist(strsplit(x,split="")))-1),function(i){substr(x,i,i+1)})
                sapply(names(PROPERTY),function(y){
                  tmp <- PROPERTY[[y]][["DINUCLEOTIDE"]]
                  mean(tmp[z])
                })
              }))
  colnames(result) <- paste("BDNAVIDEO",colnames(result),sep=":")
  result
}

featureDIPRODB <- function(seq, na.type="all", na.strand="all", 
  diprodb.method="all", diprodb.type="all"){
## DiProDB: Ref: http://diprodb.fli-leibniz.de
## DiProDB is a database of conformational and thermodynamic dinucleotide properties. It includes datasets both for DNA and RNA, as well as for single and double strands.
## They were collected from the literature and are classified according to nucleic acid type (DNA and RNA), strand information (double or single), how the data were obtained (experimental, theoretical/calculated) and also according to the general type of the dinucleotide property: thermodynamical (e.g. free energy), conformational (e.g. twist) or letter-based (e.g. GC content). 
  
## na.type - a string for nucleic acid type. It can be "DNA", "DNA/RNA", "RNA", 
##           or "all".  
## na.strand - a string for strand information. It can be "double", "single", 
##             or "all".
## diprodb.method - a string for mode of property determination. It can be 
##                  "experimental", "calculated", or "all".
## diprodb.type - a string for property type. It can be "physicochemical", 
##                "conformational", "letter based", or "all".

  if( !any(na.type==c("DNA", "DNA/RNA", "RNA", "all")) ){
    stop("Parameter na.type is not correct, must be DNA, DNA/RNA, RNA or all")
  }
  if( !any(na.strand==c("double", "single", "all")) ){
    stop("Parameter na.strand is not correct, must be double, single or all")
  }
  if( !any(diprodb.method==c("experimental", "calculated", "all")) ){
    stop("Parameter diprodb.method is not correct, must be experimental, calculated or all")
  }
  if( !any(diprodb.type==c("physicochemical", "conformational", "letter based", "all")) ){
    stop("Parameter diprodb.type is not correct, must be physicochemical, conformational, letter based or all")
  }
  
  data(DiProDB)   
  tmp <- names(DiProDB)
  if(na.type!="all"){
    tmp <- names(unlist(sapply(tmp, function(x){ grep(na.type,DiProDB[[x]]$NucleicAcid)})))
  }
  if(na.strand!="all"){
    tmp <- names(unlist(sapply(tmp, function(x){ grep(na.strand,DiProDB[[x]]$Strand)})))
  }
  if(diprodb.method!="all"){
    tmp <- names(unlist(sapply(tmp, function(x){ grep(diprodb.method,DiProDB[[x]]$HowCreated)})))
  }
  if(diprodb.type!="all"){
    tmp <- names(unlist(sapply(tmp, function(x){ grep(diprodb.type,DiProDB[[x]]$Type)})))
  }
  result <- t(sapply(toupper(seq),function(x){
                z <- sapply(1:(length(unlist(strsplit(x,split="")))-1),function(i){substr(x,i,i+1)})
                sapply(tmp,function(y){                  
                  mean(unlist(DiProDB[[y]][z]))
                })
              }))
  colnames(result) <- paste("DIPRODB",colnames(result),sep=":")
  result
}


# phastCons
# PhastCons is a program for identifying evolutionarily conserved elements in a multiple alignment, given a phylogenetic tree.
# Link: http://compgen.bscb.cornell.edu/phast/phastCons-HOWTO.html
# Ref: J. Felsenstein and G. Churchill, 1996. A hidden Markov model approach to variation among sites in rate of evolution. Mol Biol Evol, 13:93-104.

# UNAFold
# Partition functions can be computed to derive base pair probabilities and stochastic samples of foldings or hybridizations. Energy minimization methods compute minimum free energy foldings or hybridizations, and can also compute suboptimal foldings that mimic the performance of the famous mfold software.
# Link: http://dinamelt.bioinfo.rpi.edu/unafold/
# Ref: Markham, N. R. & Zuker, M. (2008) UNAFold: software for nucleic acid folding and hybriziation. In Keith, J. M., editor, Bioinformatics, Volume II. Structure, Functions and Applications, number 453 in Methods in Molecular Biology, chapter 1, pages 3每31. Humana Press, Totowa, NJ. ISBN 978-1-60327-428-9.

# Vienna RNA Package
# RNAfold -- predict minimum energy secondary structures and pair probabilities 
# RNAdistance -- compare secondary structures 
# Link: http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-1.7.2.tar.gz
# Ref: Ivo L. Hofacker. Ref: Vienna RNA secondary structure server. Nucleic Acids Research, 2003, Vol. 31, No. 13 3429-3431. 

# mfold
