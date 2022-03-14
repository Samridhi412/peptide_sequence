library(Peptides)
library(peptider)
library(log4r)
# print(getwd())
log_file <- file("101916086.log")
args = commandArgs(trailingOnly=TRUE)
if (length(args)>1) {
  cat("pass only one argument(file name)", file = log_con)
  stop("only file name argument must be supplied (input file).n", call.=FALSE)
}
tryCatch(
  expr = {
    filename <- args[1]
    df<-read.csv(args[1])
    df<-read.csv("input.csv")
    # head(df)
    len<-c()
    s<-df$Peptide.Sequence
    t<-df$Target
    lengthpep<-c()
    aIndex<-c()
    boman<-c()
    instaIndex<-c()
    ppeptide<-c()
    hmoment_alpha<-c()
    hmoment_beta<-c()
    mw_false<-c()
    mw_true<-c()
    massshift_1<-c()
    massshift_2<-c()
    mz_1<-c()
    mz_2<-c()
    mz_3<-c()
    kiderafactors1<-c()
    kiderafactors2<-c()
    kiderafactors3<-c()
    kiderafactors4<-c()
    kiderafactors5<-c()
    kiderafactors6<-c()
    kiderafactors7<-c()
    kiderafactors8<-c()
    kiderafactors9<-c()
    kiderafactors10<-c()
    hydro1<-c()
    hydro2<-c()
    hydro3<-c()
    hydro4<-c()
    hydro5<-c()
    hydro6<-c()
    hydro7<-c()
    hydro8<-c()
    hydro9<-c()
    hydro10<-c()
    hydro11<-c()
    hydro12<-c()
    hydro13<-c()
    hydro14<-c()
    hydro15<-c()
    hydro16<-c()
    hydro17<-c()
    hydro18<-c()
    hydro19<-c()
    hydro20<-c()
    hydro21<-c()
    hydro22<-c()
    hydro23<-c()
    hydro24<-c()
    charge1<-c()
    charge2<-c()
    charge3<-c()
    charge4<-c()
    charge5<-c()
    charge6<-c()
    charge7<-c()
    charge8<-c()
    charge9<-c()
    pp1<-c()
    pp2<-c()
    pp3<-c()
    pp4<-c()
    pp5<-c()
    pp6<-c()
    pp7<-c()
    pp8<-c()
    pp9<-c()
    crucianiProperties1=c()
    crucianiProperties2=c()
    crucianiProperties3=c()
    fasgaiVectors1=c()
    fasgaiVectors2=c()
    fasgaiVectors3=c()
    fasgaiVectors4=c()
    fasgaiVectors5=c()
    fasgaiVectors6=c()
    MSWHIM1<-c()
    MSWHIM2<-c()
    MSWHIM3<-c()
    ProtFP1<-c()
    ProtFP2<-c()
    ProtFP3<-c()
    ProtFP4<-c()
    ProtFP5<-c()
    ProtFP6<-c()
    ProtFP7<-c()
    ProtFP8<-c()
    st1<-c()
    st2<-c()
    st3<-c()
    st4<-c()
    st5<-c()
    st6<-c()
    st7<-c()
    st8<-c()
    t1<-c()
    t2<-c()
    t3<-c()
    t4<-c()
    t5<-c()
    vhse1<-c()
    vhse2<-c()
    vhse3<-c()
    vhse4<-c()
    vhse5<-c()
    vhse6<-c()
    vhse7<-c()
    vhse8<-c()
    z1<-c()
    z2<-c()
    z3<-c()
    z4<-c()
    z5<-c()
    b1<-c()
    b2<-c()
    b3<-c()
    b4<-c()
    b5<-c()
    b6<-c()
    b7<-c()
    b8<-c()
    b9<-c()
    b10<-c()
    # print("here")
    j=0
    for(i in df$Peptide.Sequence){
      # print(j)
      j=j+1
      if(i=="")
      {
        d=i
        e=print(paste(d,"warning empty sequence found "))
        cat(e, file = log_file)
        next
      }
        len<-append(len,nchar(i))
      aIndex<-append(aIndex,aIndex(seq=i))
      boman<-append(boman,boman(seq=i))
      a=instaIndex(seq = i)
      instaIndex<-append(instaIndex,a)
      ppeptide<-append(ppeptide,ppeptide(i, libscheme="NNK", N=10^8))
      hmoment_alpha<-append(hmoment_alpha,hmoment(seq=i, angle = 100, window = 11))
      hmoment_beta<-append(hmoment_beta,hmoment(seq=i, angle = 160, window = 11))
      mw_false=append(mw_false,mw(seq=i,monoisotopic = FALSE))
      mw_true=append(mw_true,mw(seq=i,monoisotopic = TRUE))
      massshift_1=append(massshift_1,massShift(seq=i, label = "silac_13c"))
      # massshift_2=append(massshift_2,massShift(seq=i, aaShift = c(K = 6.020129, R = 6.020129)))
      mz_1=append(mz_1,mz(seq=i))
      # mz_2=append(mz_2,mz(seq=i, aaShift = c(K = 6.020129, R = 6.020129)))
      # mz_3=append(mz_3,mz(seq=i, label = "silac_13c", cysteins = 58.005479))
      lengthpep=append(lengthpep,lengthpep(seq=i))
      y <- as.numeric(unlist(kideraFactors(seq =i)))
      kiderafactors1<-append( kiderafactors1,y[1])
      kiderafactors2<-append( kiderafactors2,y[2])
      kiderafactors3<-append( kiderafactors2,y[3])
      kiderafactors4<-append( kiderafactors4,y[4])
      kiderafactors5<-append( kiderafactors5,y[5])
      kiderafactors6<-append( kiderafactors6,y[6])
      kiderafactors7<-append( kiderafactors7,y[7])
      kiderafactors8<-append( kiderafactors8,y[8])
      kiderafactors9<-append( kiderafactors9,y[9])
      kiderafactors10<-append( kiderafactors10,y[10])
      x_num <- as.numeric(unlist(stScales(seq =i)))
      st1<-append(st1,x_num[1])
      st2<-append(st2,x_num[2])
      st3<-append(st3,x_num[3])
      st4<-append(st4,x_num[4])
      st5<-append(st5,x_num[5])
      st6<-append(st6,x_num[6])
      st7<-append(st7,x_num[7])
      st8<-append(st8,x_num[8])
      x_num <- as.numeric(unlist(mswhimScores(seq =i)))
      MSWHIM1<-append(MSWHIM1,x_num[1])
      MSWHIM2<-append(MSWHIM2,x_num[2])
      MSWHIM3<-append(MSWHIM3,x_num[3])
      x_num <- as.numeric(unlist(crucianiProperties(seq=i)))
      crucianiProperties1<-append(crucianiProperties1,x_num[1])
      crucianiProperties2<-append(crucianiProperties2,x_num[2])
      crucianiProperties3<-append(crucianiProperties3,x_num[3])
      x_num <- as.numeric(unlist(blosumIndices(seq=i)))
      b1<-append(b1,x_num[1])
      b2<-append(b2,x_num[2])
      b3<-append(b3,x_num[3])
      b4<-append(b4,x_num[4])
      b5<-append(b5,x_num[5])
      b6<-append(b6,x_num[6])
      b7<-append(b7,x_num[7])
      b8<-append(b8,x_num[8])
      b9<-append(b9,x_num[9])
      b10<-append(b10,x_num[10])
      x_num <- as.numeric(unlist(zScales(seq=i)))
      z1<-append(z1,x_num[1])
      z2<-append(z2,x_num[2])
      z3<-append(z3,x_num[3])
      z4<-append(z4,x_num[4])
      z5<-append(z5,x_num[5])
      x_num <- as.numeric(unlist(vhseScales(seq =i)))
      vhse1<-append(vhse1,x_num[1])
      vhse2<-append(vhse2,x_num[2])
      vhse3<-append(vhse3,x_num[3])
      vhse4<-append(vhse4,x_num[4])
      vhse5<-append(vhse5,x_num[5])
      vhse6<-append(vhse6,x_num[6])
      vhse7<-append(vhse7,x_num[7])
      vhse8<-append(vhse8,x_num[8])
      x_num <- as.numeric(unlist(tScales(seq = i)))
      t1<-append(t1,x_num[1])
      t2<-append(t2,x_num[2])
      t3<-append(t3,x_num[3])
      t4<-append(t4,x_num[4])
      t5<-append(t5,x_num[5])
      #####
      #######
      # a=hydrophobicity(seq = i,scale = "Aboderin")
      # hydro1<-append(hydro1,a)
      # a=hydrophobicity(seq = i,scale = "AbrahamLeo")
      # hydro2<-append(hydro2,a)
      # a=hydrophobicity(seq = i,scale = "Argos")
      # hydro3<-append(hydro3,a)
      # a=hydrophobicity(seq = i,scale = "BlackMould")
      # hydro4<-append(hydro4,a)
      # a=hydrophobicity(seq = i,scale = "BullBreese")
      # hydro5<-append(hydro5,a)
      # a=hydrophobicity(seq = i,scale = "Casari")
      # hydro6<-append(hydro6,a)
      # a=hydrophobicity(seq = i,scale = "Chothia")
      # hydro7<-append(hydro7,a)
      # a=hydrophobicity(seq = i,scale = "Cid")
      # hydro8<-append(hydro8,a)
      # a=hydrophobicity(seq = i,scale = "Cowan3.4")
      # hydro9<-append(hydro9,a)
      # a=hydrophobicity(seq = i,scale = "Cowan7.5")
      # hydro10<-append(hydro10,a)
      # a=hydrophobicity(seq = i,scale = "Eisenberg")
      # hydro11<-append(hydro11,a)
      # a=hydrophobicity(seq =i,scale = "Engelman")
      # hydro12<-append(hydro12,a)
      # a=hydrophobicity(seq = i,scale = "Fasman")
      # hydro13<-append(hydro13,a)
      # a=hydrophobicity(seq = i,scale = "Fauchere")
      # hydro14<-append(hydro14,a)
      # a=hydrophobicity(seq = i,scale = "Goldsack")
      # hydro15<-append(hydro15,a)
      # a=hydrophobicity(seq = i,scale = "Guy")
      # hydro16<-append(hydro16,a)
      # a=hydrophobicity(seq = i,scale = "HoppWoods")
      # hydro17<-append(hydro17,a)
      # a=hydrophobicity(seq = i,scale = "Janin")
      # hydro18<-append(hydro18,a)
      # a=hydrophobicity(seq = i,scale = "Jones")
      # hydro19<-append(hydro19,a)
      # a=hydrophobicity(seq = i,scale = "Juretic")
      # hydro20<-append(hydro20,a)
      # a=hydrophobicity(seq =i,scale = "Kidera")
      # hydro21<-append(hydro21,a)
      # a=hydrophobicity(seq = i,scale = "kuhn")
      # hydro22<-append(hydro22,a)
      # a=hydrophobicity(seq = i,scale = "KyteDoolittle")
      # hydro23<-append(hydro23,a)
      # a=hydrophobicity(seq = i,scale = "Levitt")
      # hydro24<-append(hydro24,a)
      # 
      a=hydrophobicity(seq = i,scale = "Aboderin")
      hydro1<-append(hydro1,a)
      
      a=hydrophobicity(seq = i,scale = "AbrahamLeo")
      hydro2<-append(hydro2,a)
      
      a=hydrophobicity(seq = i,scale = "Argos")
      hydro3<-append(hydro3,a)
      
      a=hydrophobicity(seq = i,scale = "BlackMould")
      hydro4<-append(hydro4,a)
      
      a=hydrophobicity(seq = i,scale = "BullBreese")
      hydro5<-append(hydro5,a)
      
      a=hydrophobicity(seq = i,scale = "Casari")
      hydro6<-append(hydro6,a)
      
      a=hydrophobicity(seq = i,scale = "Chothia")
      hydro7<-append(hydro7,a)
      
      a=hydrophobicity(seq = i,scale = "Cid")
      hydro8<-append(hydro8,a)
      ##########################
      a=hydrophobicity(seq = i,scale = "Cowan3.4")
      hydro9<-append(hydro9,a)
      
      a=hydrophobicity(seq = i,scale = "Cowan7.5")
      hydro10<-append(hydro10,a)
      
      a=hydrophobicity(seq = i,scale = "Eisenberg")
      hydro11<-append(hydro11,a)
      
      a=hydrophobicity(seq = i,scale = "Engelman")
      hydro12<-append(hydro12,a)
      
      a=hydrophobicity(seq = i,scale = "Fasman")
      hydro13<-append(hydro13,a)
      
      a=hydrophobicity(seq = i,scale = "Fauchere")
      hydro14<-append(hydro14,a)
      
      a=hydrophobicity(seq = i,scale = "Goldsack")
      hydro15<-append(hydro15,a)
      
      a=hydrophobicity(seq = i,scale = "Guy")
      hydro16<-append(hydro16,a)
      #################
      a=hydrophobicity(seq = i,scale = "HoppWoods")
      hydro17<-append(hydro17,a)
      
      a=hydrophobicity(seq = i,scale = "Janin")
      hydro18<-append(hydro18,a)
      
      a=hydrophobicity(seq = i,scale = "Jones")
      hydro19<-append(hydro19,a)
      
      a=hydrophobicity(seq = i,scale = "Juretic")
      hydro20<-append(hydro20,a)
      
      a=hydrophobicity(seq = i,scale = "Kidera")
      hydro21<-append(hydro21,a)
      
      a=hydrophobicity(seq = i,scale = "Kuhn")
      hydro22<-append(hydro22,a)
      
      a=hydrophobicity(seq = i,scale = "KyteDoolittle")
      hydro23<-append(hydro23,a)
      
      a=hydrophobicity(seq = i,scale = "Levitt")
      hydro24<-append(hydro24,a)
      a=charge(seq= i,pH= 5, pKscale= "Bjellqvist")
      charge1<-c(charge1,a)
      
      a=charge(seq= i,pH= 5, pKscale= "EMBOSS")
      charge2<-c(charge2,a)
      
      a=charge(seq= i,pH= 5, pKscale= "Murray")
      charge3<-c(charge3,a)
      
      a=charge(seq= i,pH= 5, pKscale= "Sillero")
      charge4<-c(charge4,a)
      
      a=charge(seq= i,pH= 5, pKscale= "Solomon")
      charge5<-c(charge5,a)
      
      a=charge(seq= i,pH= 5, pKscale= "Stryer")
      charge6<-c(charge6,a)
      
      a=charge(seq= i,pH= 5, pKscale= "Lehninger")
      charge7<-c(charge7,a)
      
      a=charge(seq= i,pH= 5, pKscale= "Dawson")
      charge8<-c(charge8,a)
      
      a=charge(seq= i,pH= 5, pKscale= "Rodwell")
      charge9<-c(charge9,a)

      # ####################
      # a=charge(seq= i,pH= 5, pKscale= "Bjellqvist")
      # charge1_p5<-c(charge1_p5,a)
      # 
      # a=charge(seq= i,pH= 5, pKscale= "EMBOSS")
      # charge2_p5<-c(charge2_p5,a)
      # 
      # a=charge(seq= i,pH= 5, pKscale= "Murray")
      # charge3_p5<-c(charge3_p5,a)
      # 
      # a=charge(seq= i,pH= 5, pKscale= "Sillero")
      # charge4_p5<-c(charge4_p5,a)
      # 
      # a=charge(seq= seq,pH= 5, pKscale= "Solomon")
      # charge5_p5<-c(charge5_p5,a)
      # 
      # a=charge(seq= i,pH= 5, pKscale= "Stryer")
      # charge6_p5<-c(charge6_p5,a)
      # 
      # a=charge(seq= i,pH= 5, pKscale= "Lehninger")
      # charge7_p5<-c(charge7_p5,a)
      # 
      # a=charge(seq= i,pH= 5, pKscale= "Dawson")
      # charge8_p5<-c(charge8_p5,a)
      # 
      # a=charge(seq= i,pH= 5, pKscale= "Rodwell")
      # charge9_p5<-c(charge9_p5,a)
      # a=charge(seq= i,pH= 9, pKscale= "Bjellqvist")
      # charge1_p9<-c(charge1_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "EMBOSS")
      # charge2_p9<-c(charge2_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "Murray")
      # charge3_p9<-c(charge3_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "Sillero")
      # charge4_p9<-c(charge4_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "Solomon")
      # charge5_p9<-c(charge5_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "Stryer")
      # charge6_p9<-c(charge6_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "Lehninger")
      # charge7_p9<-c(charge7_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "Dawson")
      # charge8_p9<-c(charge8_p9,a)
      # 
      # a=charge(seq= i,pH= 9, pKscale= "Rodwell")
      # charge9_p9<-c(charge9_p9,a)
      a=pI(seq= i,pKscale= "Bjellqvist")
      pp1<-append(pp1,a)
      a=pI(seq=i,pKscale= "EMBOSS")
      pp2<-append(pp2,a)
      
      a=pI(seq= i,pKscale= "Solomon")
      pp3<-append(pp3,a)
      
      a=pI(seq= i,pKscale= "Stryer")
      pp4<-append(pp4,a)
      
      a=pI(seq= i,pKscale= "Lehninger")
      pp5<-append(pp5,a)
      
      a=pI(seq= i,pKscale= "Dawson")
      pp6<-append(pp6,a)
      
      a=pI(seq= i,pKscale= "Rodwell")
      pp7<-append(pp7,a)
      a=pI(seq= i,pKscale= "Murray")
      pp8<-append(pp8,a)
      a=pI(seq= i,pKscale= "Sillero")
      pp9<-append(pp9,a)
      x_num <- as.numeric(unlist(fasgaiVectors(seq = i)))
      fasgaiVectors1<-append(fasgaiVectors1,x_num[1])
      fasgaiVectors2<-append(fasgaiVectors2,x_num[2])
      fasgaiVectors3<-append(fasgaiVectors3,x_num[3])
      fasgaiVectors4<-append(fasgaiVectors4,x_num[4])
      fasgaiVectors5<-append(fasgaiVectors5,x_num[5])
      fasgaiVectors6<-append(fasgaiVectors6,x_num[6])
      x_num <- as.numeric(unlist(protFP(seq = i)))
      ProtFP1<-append(ProtFP1,x_num[1])
      ProtFP2<-append(ProtFP2,x_num[2])
      ProtFP3<-append(ProtFP3,x_num[3])
      ProtFP4<-append(ProtFP4,x_num[4])
      ProtFP5<-append(ProtFP5,x_num[5])
      ProtFP6<-append(ProtFP6,x_num[6])
      ProtFP7<-append(ProtFP7,x_num[7])
      ProtFP8<-append(ProtFP8,x_num[8])
    }
    print(length(s))
    print(length(aIndex))
    print(length(boman))
    print(length(len))
    print(j)
      final<-data.frame(
                       "peptideSequence"=s,
                       "length_sequence"=len,
                       "aIndex"=aIndex,
                       "boman"=boman,
                       "instaIndex"=instaIndex,
                       "lengthpep"=lengthpep,

                       "ppeptide"=ppeptide,

                       "hmoment_alpha"=hmoment_alpha,
                       "hmoment_beta"= hmoment_beta,
                       "mw_false"= mw_false,
                       "mw_true"=mw_true,
                       "mz_1"= mz_1,
                       "mz_2"= mz_2,
                       "mz_3"=  mz_3,
                       "massshift_1"= massshift_1,
                       "massshift_2"= massshift_2,
                       "kiderafactors(KF1)"=kiderafactors1,
                       "kiderafactors(KF2)"=kiderafactors2,
                       "kiderafactors(KF3)"= kiderafactors3,
                       "kiderafactors(KF4)"=kiderafactors4,
                       "kiderafactors(KF5)"= kiderafactors5,
                       "kiderafactors(KF6)"= kiderafactors6,
                       "kiderafactors(KF7)"=kiderafactors7,
                       "kiderafactors(KF8)"=kiderafactors8,
                       "kiderafactors(KF9)"=kiderafactors9,
                       "kiderafactors(KF10)"=kiderafactors10,
                       "hydrophobicity1"=hydro1,
                       "hydrophobicity2"=hydro2,
                       "hydrophobicity3"=hydro3,"hydrophobicity4"=hydro4,"hydrophobicity5"=hydro5,"hydrophobicity6"=hydro6,"hydrophobicity7"=hydro7,"hydrophobicity8"=hydro8,"hydrophobicity9"=hydro9,"hydrophobicity10"=hydro10,"hydrophobicity11"=hydro11,"hydrophobicity12"=hydro12,"hydrophobicity13"=hydro13,"hydrophobicity14"=hydro14,"hydrophobicity15"=hydro15,"hydrophobicity16"=hydro16,"hydrophobicity17"=hydro17,"hydrophobicity18"=hydro18,"hydrophobicity19"=hydro19,"hydrophobicity20"=hydro20,"hydrophobicity21"=hydro21,"hydrophobicity22"=hydro22,"hydrophobicity23"=hydro23,"hydrophobicity24"=hydro24,
                      "charge1"=charge1,
                      "charge2"=charge2,
                      "charge3"=charge3,
                      "charge4"=charge4,
                      "charge5"=charge5,
                      "charge6"=charge6,
                      "charge7"=charge7,
                      "charge8"=charge8,
                      "charge9"=charge9,

                       ###########
                      "MSWHIM1"=MSWHIM1,
                      "MSWHIM2"= MSWHIM2,
                      "MSWHIM3"=MSWHIM3,
                      "pI1"=pp1,
                      "pI2"=pp2,
                      "pI3"=pp3,
                      "pI4"=pp4,
                      "pI5"=pp5,
                      "pI6"=pp6,
                      "pI7"=pp7,
                      "pI8"=pp8,
                      "pI9"=pp9,
                      ###########
                      #############
                      "crucianiProperties1"= crucianiProperties1,
                      "crucianiProperties2"=crucianiProperties2,
                      "crucianiProperties3"=crucianiProperties3,
                      "fasgaiVectors1"=fasgaiVectors1,
                      "fasgaiVectors2"=fasgaiVectors2,
                      "fasgaiVectors3"=fasgaiVectors3,
                      "fasgaiVectors4"=fasgaiVectors4,
                      "fasgaiVectors5"=fasgaiVectors5,
                      "fasgaiVectors6"=fasgaiVectors6,
                      "ProtFP1"= ProtFP1,
                      "ProtFP2"= ProtFP2,
                      "ProtFP3"= ProtFP3,
                      "ProtFP4"= ProtFP4,
                      "ProtFP5"= ProtFP5,
                      "ProtFP6"= ProtFP6,
                      "ProtFP7"= ProtFP7,
                      "ProtFP8"= ProtFP8,
                      ##########
                      #########
                      "st_scales1"=st1,
                      "st_scales2"=st2,
                      "st_scales3"=st3,
                      "st_scales4"=st4,
                      "st_scales5"=st5,
                      "st_scales6"=st6,
                      "st_scales7"=st7,
                      "st_scales8"=st8,

                      "t_scales1"=t1,
                      "t_scales2"=t2,
                      "t_scales3"=t3,
                      "t_scales4"=t4,
                      "t_scales5"=t5,
                      #######

                      #############
                      "z_scales1"=z1,
                      "z_scales2"=z2,
                      "z_scales3"=z3,
                      "z_scales4"=z4,
                      "z_scales5"=z5,

                      ###########
                      "vhse1"= vhse1,
                      "vhse2"= vhse2,
                      "vhse3"= vhse3,
                      "vhse4"= vhse4,
                      "vhse5"= vhse5,
                      "vhse6"= vhse6,
                      "vhse7"= vhse7,
                      "vhse8"= vhse8,
                      "blosumIndices1"=b1,
                      "blosumIndices2"=b2,
                      "blosumIndices3"=b3,
                      "blosumIndices4"=b4,
                      "blosumIndices5"=b5,
                      "blosumIndices6"=b6,
                      "blosumIndices7"=b7,
                      "blosumIndices8"=b8,
                      "blosumIndices9"=b9,
                      "blosumIndices10"=b10,
                      "target"=t
                 )
      write.csv(final,"101916086_ouput.csv",row.names = FALSE)
  },
  error = function(e){
    cat("File not found", file = log_file)
    print("File not found error")
  }
)


