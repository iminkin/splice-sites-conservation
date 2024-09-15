#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = splice-sites-conservation
PYTHON_INTERPRETER = python
CONFIGURATION = release

ifeq ($(CONFIGURATION),debug)
complete_chromosomes := chr21 chr22
rand_limit := 20000
chromosomes := chr21 chr22
ucsc_genomes := panTro6 panPan3
zoo_genomes := HLmacFus1 HLallNig1
ncbi_genomes := HLhylMol2 HLsemEnt1
else
rand_limit := 180000
chromosomes = chr19_KI270868v1_alt chr3_KI270784v1_alt chr11_KI270829v1_alt chr14_KZ208919v1_alt chr8_KI270819v1_alt chr14_KI270722v1_random chr7_KI270803v1_alt chr1_KI270760v1_alt chr7_KZ208912v1_fix chr19 chr10_KQ090020v1_alt chrUn_KI270330v1 chr19_KV575253v1_alt chr19_GL383574v1_alt chr2_KI270716v1_random chr19_KI270920v1_alt chr2_KI270894v1_alt chr19_KI270885v1_alt chr21_GL383578v2_alt chr19_GL949749v2_alt chr11_KZ559109v1_fix chr5_KI270791v1_alt chr17_KI270862v1_alt chrUn_KI270753v1 chr8_KI270818v1_alt chr7_KI270899v1_alt chr16_KQ031390v1_alt chr10 chr4_KI270925v1_alt chr15_KQ031389v1_alt chr1_KI270761v1_alt chr8_KI270926v1_alt chr22_KN196485v1_alt chr19_KI270884v1_alt chr19_KI270921v1_alt chr19_KV575252v1_alt chr1_KZ208906v1_fix chr17 chr19_GL949748v2_alt chr11_KZ559108v1_fix chrX chr14_GL000225v1_random chr21_GL383579v2_alt chr14_KI270723v1_random chr9_KI270823v1_alt chr1_KI270762v1_alt chr19_GL949752v1_alt chr3_KZ208909v1_alt chrUn_KI270749v1 chr22 chr5_KI270793v1_alt chr17_KI270860v1_alt chr10_KN196480v1_fix chr17_GL000205v2_random chr9_KI270718v1_random chr19_GL383576v1_alt chr19_KV575251v1_alt chr19_KI270922v1_alt chr19_KI270887v1_alt chr22_KN196486v1_alt chr17_KV766197v1_alt chr12_KI270833v1_alt chr1_KI270763v1_alt chr9_KI270719v1_random chrY_KZ208924v1_fix chr18_GL383567v1_alt chrUn_GL000220v1 chr6_KV766194v1_fix chr17_GL383564v2_alt chr19_KQ458386v1_fix chrUn_KI270745v1 chr17_KI270861v1_alt chrM chr12_KN538369v1_fix chr5_KI270792v1_alt chr19_KI270886v1_alt chr19_KI270923v1_alt chr19_KV575250v1_alt chr5_KV575243v1_alt chr1_KN538360v1_fix chr9_GL383539v1_alt chr5_KI270795v1_alt chr14_KI270725v1_random chr19_KV575257v1_alt chr7_KI270807v1_alt chr1_KI270764v1_alt chr12_KI270834v1_alt chr11_GL383547v1_alt chr4_KI270785v1_alt chr20_KI270869v1_alt chr17_KI270730v1_random chrY_KZ208923v1_fix chr22_KQ458388v1_alt chr6_KQ031387v1_fix chr7_KV880765v1_fix chr10_KN538367v1_fix chr3_KI270780v1_alt chrUn_KI270579v1 chr7_KZ208913v1_alt chr10_KQ090021v1_fix chr22_KI270928v1_alt chr5_KI270794v1_alt chr1_KN538361v1_fix chr19_KV575256v1_alt chr17_KI270907v1_alt chr5 chrUn_KI270752v1 chr1_KQ031383v1_fix chr2 chr7_KI270806v1_alt chr12_KI270835v1_alt chr1_KI270765v1_alt chr3_KI270924v1_alt chrUn_KI270438v1 chr14_KI270724v1_random chr7_KV880764v1_fix chr10_KN538366v1_fix chr3_KI270781v1_alt chr19_KI270883v1_alt chrUn_KI270748v1 chr19_KV575255v1_alt chr5_KV575244v1_fix chr5_KN196477v1_alt chr17_KV766196v1_fix chr1_KZ208904v1_alt chr5_KI270897v1_alt chr5_GL949742v1_alt chr14_KI270726v1_random chr3_KI270782v1_alt chr10_KN538365v1_fix chr4_GL383528v1_alt chr4_KI270787v1_alt chr12_KI270836v1_alt chr1_KI270766v1_alt chr1_KI270707v1_random chrUn_KI270521v1 chr7_KI270805v1_alt chr19_GL383573v1_alt chr19_KV575254v1_alt chr2_KI270893v1_alt chr19_KI270882v1_alt chr1_KZ208905v1_alt chr1_KI270706v1_random chr5_KI270796v1_alt chr12_KZ559112v1_alt chr3_KI270783v1_alt chrUn_KI270744v1 chr12_KI270837v1_alt chr7_KI270804v1_alt chr4_KI270786v1_alt chr16_KI270853v1_alt chr3_KV766192v1_fix chr19_KI270929v1_alt chrUn_GL000219v1 chr18_KI270911v1_alt chr5_KI270898v1_alt chr21_KI270874v1_alt chr12_KN196482v1_fix chr17_GL383563v3_alt chr3_KN538364v1_fix chr9_GL383542v1_alt chr8_KI270810v1_alt chr2_KI270770v1_alt chr4_GL383527v1_alt chr13_KN538373v1_fix chr4_KI270788v1_alt chr19_KN196484v1_fix chr9_KQ090019v1_alt chr4_KV766193v1_alt chr2_KI270771v1_alt chr12_GL383552v1_alt chr8_KI270811v1_alt chr7_KZ559106v1_alt chrUn_KI270751v1 chr9_KQ090018v1_alt chr13_KN538372v1_fix chr4_KI270789v1_alt chr1_KQ458383v1_alt chr15_KI270848v1_alt chrX_KI270913v1_alt chr17_KI270909v1_alt chr19_KV575258v1_alt chr16 chrY chr13_KN538371v1_fix chr7_KI270808v1_alt chr11 chr15_GL383554v1_alt chr8_KI270812v1_alt chr12_GL383551v1_alt chr9_GL383540v1_alt chr2_KI270772v1_alt chr22_KQ458387v1_alt chr18_KI270912v1_alt chrUn_KI270747v1 chr15_KI270849v1_alt chr1_KQ458382v1_alt chr3_KQ031385v1_fix chr12_KQ090023v1_alt chr9_KI270717v1_random chr19_KV575259v1_alt chr17_KI270908v1_alt chr7_KI270809v1_alt chrUn_KI270442v1 chrY_KN196487v1_fix chr18 chr2_KI270773v1_alt chr12_KI270904v1_alt chr8_KI270813v1_alt chr9_GL383541v1_alt chr15_KI270906v1_alt chr17_KZ559114v1_alt chr2_KQ983256v1_alt chr19_KI270865v1_alt chr4_GL000008v2_random chr3_GL383526v1_alt chrUn_GL000218v1 chr2_KI270774v1_alt chr18_GL383569v1_alt chr8_KI270814v1_alt chr17_KV766198v1_alt chr19_KI270888v1_alt chr14_GL000194v1_random chr12_GL877876v1_alt chr18_GL383568v1_alt chr8_KI270815v1_alt chr2_KI270775v1_alt chr1_KQ458384v1_alt chr5_GL383532v1_alt chrUn_GL000214v1 chr19_KI270889v1_alt chrUn_KI270750v1 chr16_KI270856v1_alt chr11_KI270826v1_alt chr2_KI270776v1_alt chr8_KI270816v1_alt chr3 chr12_GL877875v1_alt chr22_GL383583v2_alt chr19_KI270867v1_alt chr8_KI270900v1_alt chr11_JH159137v1_alt chr1_KI270709v1_random chr16_KI270855v1_alt chr5_GL383531v1_alt chr21_KI270872v1_alt chr4 chr22_GL383582v2_alt chr19_KI270866v1_alt chr1_KI270708v1_random chr8_KI270817v1_alt chr11_KI270827v1_alt chr15_KI270905v1_alt chrUn_KI270746v1 chr16_KI270854v1_alt chr11_JH159136v1_alt chr8_KI270901v1_alt chr15_KI270727v1_random chr9_KN196479v1_fix chrUn_GL000195v1 chr5_KZ208910v1_alt chr19_GL949747v2_alt chr16_KV880768v1_fix chr21_KI270873v1_alt chr5_GL383530v1_alt chr4_GL000257v2_alt chr4_KQ983257v1_fix chrUn_GL000213v1 chr10_KI270825v1_alt chr1_KN196474v1_fix chr17_JH159147v1_alt chr2_KZ208907v1_alt chr22_KI270878v1_alt chr18_KI270864v1_alt chr11_KZ559110v1_alt chr20_KI270871v1_alt chr22_KI270736v1_random chr8_KV880767v1_fix chr19_KI270914v1_alt chr13_KQ090025v1_alt chr8_KZ208914v1_fix chr17_GL000258v2_alt chr22_KI270737v1_random chr11_KZ559111v1_alt chr22_KI270879v1_alt chr17_JH159146v1_alt chr10_KI270824v1_alt chr19_KI270915v1_alt chr8_KV880766v1_fix chr20_KI270870v1_alt chr6_KQ090017v1_alt chr8_KZ208915v1_fix chrUn_KI270581v1 chr13_KQ090024v1_alt chr17_KI270857v1_alt chrUn_KI270528v1 chr12_KZ208917v1_fix chrUn_GL000224v1 chr2_KN538363v1_fix chr1_KZ559100v1_fix chr3_KI270777v1_alt chr22_KI270735v1_random chr11_KN538368v1_alt chr7 chr13_KN538373v1_fix chr12_KQ759760v1_fix chr2_KI270767v1_alt chr19_KI270916v1_alt chr3_JH636055v2_alt chrUn_KI270741v1 chr16_KZ559113v1_fix chr2_KN538362v1_fix chr12_KZ208916v1_fix chr6_KI270797v1_alt chr1_KI270892v1_alt chrUn_KI270448v1 chr9 chr19_KI270917v1_alt chr20_GL383577v2_alt chr22_KI270734v1_random chr11_KI270831v1_alt chr19_GL383575v2_alt chr19_GL000209v2_alt chrUn_KI270756v1 chr17_GL383566v1_alt chr11_KI270927v1_alt chr19_KI270938v1_alt chr13_KN196483v1_fix chr22_KI270731v1_random chr13_KI270838v1_alt chr19_GL949751v2_alt chr17_KV575245v1_fix chr16_KZ208921v1_alt chr11_KI270830v1_alt chr6_KQ090016v1_fix chr20 chr2_GL383522v1_alt chr22_KB663609v1_alt chr19_GL949750v2_alt chr13_KI270839v1_alt chr12 chr22_KI270732v1_random chr19_GL949753v2_alt chr12_KN538370v1_fix chr17_KI270729v1_random chr7_GL383534v2_alt chr2_GL383521v1_alt chr19_KV575249v1_alt chr15 chr1_KN196472v1_fix chr1_GL383520v2_alt chr19_KV575260v1_alt chr17_GL383565v1_alt chr11_KI270832v1_alt chr2_KQ031384v1_fix chrY_KI270740v1_random chr2_GL582966v2_alt chr1_KN196473v1_fix chr19_KV575248v1_alt chr22_KI270733v1_random chr15_KN538374v1_fix chr18_KI270863v1_alt chr7_KQ031388v1_fix chr21_GL383580v2_alt chr14_KI270846v1_alt chr11_KI270902v1_alt chr3_KI270934v1_alt chr18_GL383571v1_alt chr8 chrX_KI270880v1_alt chr18_KZ559116v1_alt chr13_KI270843v1_alt chr18_KQ458385v1_alt chr6_KN196478v1_fix chr10_GL383546v1_alt chr6_KI270758v1_alt chr19_KV575246v1_alt chr22_KI270738v1_random chrUn_KI270755v1 chr19_KI270890v1_alt chr18_KZ208922v1_fix chr3_KN196475v1_fix chr11_KI270903v1_alt chr1 chr14_KI270847v1_alt chr21_GL383581v2_alt chr13_KI270842v1_alt chrX_KI270881v1_alt chr6_GL000256v2_alt chr18_GL383570v1_alt chrUn_KI270583v1 chr3_KI270935v1_alt chr12_KZ208918v1_alt chrUn_GL000216v2 chr6_GL383533v1_alt chr14_KZ208920v1_fix chr6 chr19_KI270891v1_alt chr19_KV575247v1_alt chr3_GL000221v1_random chr13_KI270841v1_alt chr3_KZ559105v1_alt chr1_KQ983255v1_alt chr14_GL000009v2_random chr1_KI270714v1_random chr3_KI270936v1_alt chr6_GL000255v2_alt chrUn_KI270743v1 chr3_KN196476v1_fix chr4_KI270896v1_alt chr5_GL000208v1_random chr14_KI270844v1_alt chr5_GL339449v2_alt chr11_KQ090022v1_fix chr11_KI270721v1_random chr4_KQ090014v1_alt chr6_GL000254v2_alt chr18_GL383572v1_alt chr3_KI270937v1_alt chr13_KI270840v1_alt chr14_KI270845v1_alt chr19_GL949746v1_alt chr11_KQ759759v1_fix chr4_KQ983258v1_alt chrX_KV766199v1_alt chr4_KQ090015v1_alt chr22_KQ759762v1_fix chr10_GL383545v1_alt chr3_KI270778v1_alt chr19_KI270931v1_alt chr1_KI270759v1_alt chr6_KI270799v1_alt chr22_KI270875v1_alt chr15_KI270852v1_alt chr6_KI270802v1_alt chr8_KI270820v1_alt chr16_GL383557v1_alt chr6_GL000253v2_alt chrUn_KI270754v1 chr2_KI270768v1_alt chr3_KI270895v1_alt chr3_KZ559103v1_alt chr4_KI270790v1_alt chr1_KI270710v1_random chr19_KI270919v1_alt chr6_KI270798v1_alt chrUn_KI270582v1 chr12_GL383553v2_alt chr22_KQ759761v1_alt chr19_KI270930v1_alt chr3_KI270779v1_alt chr8_KI270821v1_alt chr14 chr4_KQ090013v1_alt chr1_KI270711v1_random chr3_KZ559102v1_alt chr2_KI270769v1_alt chr13 chr16_GL383556v1_alt chr11_KV766195v1_fix chr6_GL000252v2_alt chr19_KI270918v1_alt chr1_KI270713v1_random chrUn_KI270742v1 chr15_GL383555v2_alt chr8_KI270822v1_alt chr6_KI270800v1_alt chr6_KZ208911v1_fix chr15_KI270850v1_alt chr16_KQ090026v1_alt chr22_KI270877v1_alt chr6_KB021644v2_alt chr19_KI270933v1_alt chr12_GL383550v2_alt chr17_JH159148v1_alt chr2_KZ208908v1_alt chr18_KZ559115v1_fix chr3_KZ559104v1_fix chr8_KZ559107v1_alt chr11_KN196481v1_fix chr1_GL383519v1_alt chr3_KZ559101v1_alt chr17_KI270859v1_alt chr6_GL000251v2_alt chr12_GL383549v1_alt chr21 chr15_KI270851v1_alt chr6_KI270801v1_alt chr9_KI270720v1_random chr1_KV880763v1_alt chr19_KI270932v1_alt chr17_KI270910v1_alt chr16_KQ090027v1_alt chr22_KI270876v1_alt chr1_KI270712v1_random chr16_KI270728v1_random chr18_KQ090028v1_fix chr6_GL000250v2_alt chr17_KI270858v1_alt chr1_GL383518v1_alt
ucsc_genomes := canFam4 tarSyr2 panPan3 cerSim1 equCab3 panTro6 canFam5 hetGla2 macEug2 myoLuc2 felCat9 tupBel1 speTri2 turTru2 otoGar3 cavPor3 mm39 nasLar1 chlSab2 susScr11 saiBol1 triMan1 dasNov3 vicPac2 neoSch1 gorGor6 ochPri3 rn6 manPen1 oryCun2 bosTau9 bisBis1 balAcu1 eriEur2 rheMac10 ponAbe3 enhLutNer1 sorAra2 mm10 echTel2 monDom5
zoo_genomes := HLbalEde1 HLconTau2 HLmirAng2 HLcasCan3 HLursAme2 HLeulMon1 HLgirCam2 HLpseCup1 HLmyoSep1 HLbasSum1 HLmyoLuc1 HLlamGla1 HLturAdu2 HLchoHof3 HLmacFus1 HLpseOcc1 HLeidDup1 HLodoRos1 HLeubGla1 HLneoNeb1 HLoryDam1 HLsylBac1 HLokaJoh2 HLlycPic3 HLeulFla1 HLpepEle1 HLeriBar1 HLeidHel2 HLnasNar1 HLphoPho2 HLequQuaBoe1 HLpanOnc2 HLallNig1 HLtapInd2 HLpseCor1 HLpteRuf1 HLcryFer2 HLeleMax1 HLcalPym1 HLrouMad1 HLailFul2 HLodoVir2 HLpotFla1 HLperNas1 HLcoePre1 HLsaiBol1 HLrhiUni1 HLperCri1 HLpteBra2
ncbi_genomes := HLnomLeu4 HLhylMol2 HLtheGel1 HLmacFas6 HLcerMon1 HLpilTep2 HLpapAnu5 HLmanSph1 HLsemEnt1 HLtraFra1 HLrhiRox2 HLpygNem1 HLcerNeg1 HLeryPat1 HLpitPit1 HLateGeo1 HLpleDon1 HLaloPal1 HLsagImp1 HLcalJac4 HLsapApe1 HLcebAlb1 HLdauMad1 HLindInd1 HLeulFul1 HLlemCat1 HLproSim1 HLcheMed1 HLmirZaz1 HLmirCoq1 HLmicSpe31 HLmicTav1 HLnycCou1 HLgalVar2 HLpumCon1 HLpanPar1 HLpanOnc1 HLphoVit1 HLaciJub2 HLhalGry1 HLpriBen1 HLursThi1 HLlynPar1 HLmirLeo1 HLlynCan1 HLpanLeo1 HLcalUrs1 HLailMel2 HLzalCal1 HLeumJub1 HLpumYag1 HLursArc1 HLarcGaz2 HLursAme1 HLfelNig1 HLvulVul1 HLvulLag1 HLcanLupDin1 HLlycPic2 HLcroCro1 HLhyaHya1 HLparHer1 HLmarZib1 HLgulGul1 HLspiGra1 HLtaxTax1 HLmelCap1 HLsurSur2 HLsurSur1 HLmunMug1 HLhelPar1 HLlonCan1 HLpteBra1 HLlutLut1 HLmusErm1 HLneoVis1 HLmusPut1 HLmusFur2 HLtapTer1 HLtapInd1 HLdicBic1 HLdicSum1 HLcerSimCot1 HLequAsi1 HLequAsiAsi2 HLmanTri1 HLmanPen2 HLmanJav1 HLmanJav2 HLsolPar1 HLscaAqu1 HLtalOcc1 HLuroGra1 HLbalBon1 HLescRob1 HLphyCat2 HLdelLeu2 HLhipAmp3 HLhipAmp1 HLphoPho1 HLcamFer3 HLcamBac1 HLsouChi1 HLcamDro2 HLturTru4 HLturAdu1 HLturTru3 HLgirTip1 HLbubBub2 HLbosInd2 HLbosMut2 HLtraStr1 HLammLer1 HLcapHir2 HLcapAeg1 HLodoVir3 HLoviCan2 HLoviAri5 HLelaDav1 HLoviAmm1 HLeubJap1 HLmegNov1 HLbalPhy1 HLbalMus1 HLmesBid1 HLplaMin1 HLzipCav1 HLlniGeo1 HLmonMon1 HLneoAsi1 HLkogBre1 HLphoSin1 HLponBla1 HLvicPacHua3 HLgloMel1 HLlagObl1 HLvicVicMen1 HLlamGuaCac1 HLlamGlaCha1 HLcatWag1 HLgirCam1 HLsynCaf1 HLmosBer1 HLmosMos1 HLmosChr1 HLcerHanYar1 HLbosGau1 HLbosFro1 HLprzAlb1 HLhipEqu1 HLcepHar1 HLhipNig1 HLbosGru1 HLsylGri1 HLphiMax1 HLoryGaz1 HLantAme1 HLmunRee1 HLmunCri1 HLcerEla1 HLtraImb1 HLtraScr1 HLkobEll1 HLmunMun1 HLdamLun1 HLoviCan1 HLkobLecLec1 HLcapPyg1 HLalcAlc1 HLbeaHun1 HLaepMel1 HLodoHem1 HLredRed1 HLcapSib1 HLranTarGra2 HLranTar1 HLoreOre1 HLhydIne1 HLoviNivLyd1 HLneoMos1 HLodoVir1 HLhemHyl1 HLoviOri1 HLneoPyg1 HLnanGra1 HLproPrz1 HLrapCam1 HLeudTho1 HLantMar1 HLlitWal1 HLmadKir1 HLaxiPor1 HLtraJav1 HLtraKan1 HLoreAme1 HLsaiTat1 HLcapIbe1 HLchoDid2 HLchoDid1 HLmyrTri1 HLtamTet1 HLtolMat1 HLrhiSin1 HLhipArm1 HLptePse1 HLtadBra1 HLpteVam2 HLpteGig1 HLhipGal1 HLcynBra1 HLeonSpe1 HLrouLes1 HLmegLyr2 HLmacSob1 HLminSch1 HLminNat1 HLcraTho1 HLmorBla1 HLmacCal1 HLmicHir1 HLanoCau1 HLmurAurFea1 HLnocLep1 HLcarPer3 HLtonSau1 HLartJam1 HLartJam2 HLstuHon1 HLaeoCin1 HLantPal1 HLnycHum2 HLlasBor1 HLpipPip1 HLpipPip2 HLsciCar1 HLsciVul1 HLxerIna1 HLaplRuf1 HLmarFla1 HLmarMar1 HLmarVan1 HLmarHim1 HLspeDau1 HLmarMon1 HLmarMon2 HLuroPar1 HLcynGun1 HLgliGli1 HLpedCap1 HLgraMur1 HLlepTim1 HLlepAme1 HLoryCunCun4 HLoryCun3 HLhysCri1 HLereDor1 HLmusAve1 HLfukDam2 HLdasPun1 HLcteGun1 HLallBul1 HLdipSte1 HLrhiPru1 HLdinBra1 HLzapHud1 HLdolPat1 HLperLonPac1 HLhydHyd1 HLpetTyp1 HLcavTsc1 HLthrSwi1 HLcriGam1 HLneoLep1 HLcteSoc1 HLmyoCoy1 HLcriGri3 HLondZib1 HLperCal2 HLperEre1 HLonyTor1 HLperLeu1 HLellTal1 HLperPol1 HLperManBai2 HLsigHis1 HLellLut1 HLmyoGla2 HLarvAmp1 HLpsaObe1 HLacoRus1 HLgraSur1 HLarvNil1 HLmicOec1 HLmicAgr2 HLmicFor1 HLacoCah1 HLmicArv1 HLrhoOpi1 HLmasCou1 HLmerUng1 HLratRat7 HLratNor7 HLmusPah1 HLmusCar1 HLmusSpi1 HLmusSpr1 HLapoSyl1 HLdugDug1 HLhydGig1 HLhetBru1 HLproCap3 HLmicTal1 HLvomUrs1 HLphaCin1 HLgraAgi1 HLtriVul1 HLgymLea1 HLthyCyn1 HLantFla1 HLsarHar2 HLornAna3 HLtacAcu1
endif

gencode_version := 45
refseq_version := 110
mane_version := 1.3
chess_version := 3.1.0
random_version := 1.0

#################################################################################
# RAW FILES                                                                     #
#################################################################################
human_version := GCF_000001405.40_GRCh38.p14
human_genome := data/raw/human/$(human_version)_genomic.fna.gz
human_stats := data/raw/human/$(human_version)_assembly_report.txt
phast_wig := $(foreach chr,$(chromosomes),data/raw/phast/$(chr).phastCons470way.wigFix.gz)
maf_files := $(foreach chr,$(chromosomes),data/raw/maf/$(chr).maf)

gnomad_files := $(foreach chr,$(complete_chromosomes),data/raw/gnomad/gnomad.genomes.v4.0.sites.$(chr).vcf.bgz)

genomes_dir := data/raw/genomes
ucsc_genomes_files := $(foreach genome,$(ucsc_genomes),$(genomes_dir)/$(genome).fa.gz)
zoo_genomes_files := $(foreach genome,$(zoo_genomes),$(genomes_dir)/$(genome).fa.gz)
ncbi_genomes_files := $(foreach genome,$(ncbi_genomes),$(genomes_dir)/$(genome).fa.gz)
all_genomes := $(ucsc_genomes) $(zoo_genomes) $(ncbi_genomes)

gencode_gtf := data/raw/annotation/gencode.v$(gencode_version).annotation.gtf.gz
refseq_filename := $(human_version)_genomic.gtf.gz
refseq_gtf := data/raw/annotation/$(refseq_filename)
mane_gtf := data/raw/annotation/MANE.GRCh38.v$(mane_version).ensembl_genomic.gtf.gz
chess_gtf := data/raw/annotation/chess$(chess_version).GRCh38.gtf.gz
random_gtf := data/raw/annotation/random.gtf.gz

#################################################################################
# INTERIM FILES                                                                 #
#################################################################################

clinvar_tracks := data/interim/clinvar/clinvar.csv
gnomad_tracks := $(foreach chr,$(complete_chromosomes),data/interim/gnomad/$(chr).csv)

gencode_dir := data/interim/gencode/$(gencode_version)
gencode_introns := $(foreach chr,$(chromosomes),$(gencode_dir)/introns_all/$(chr))
gencode_query := $(foreach chr,$(chromosomes),$(gencode_dir)/query_all/$(chr))

refseq_dir := data/interim/refseq/$(refseq_version)
refseq_introns := $(foreach chr,$(chromosomes),$(refseq_dir)/introns_all/$(chr))
refseq_query := $(foreach chr,$(chromosomes),$(refseq_dir)/query_all/$(chr))

mane_dir := data/interim/mane/$(mane_version)
mane_introns := $(foreach chr,$(chromosomes),$(mane_dir)/introns_all/$(chr))
mane_query := $(foreach chr,$(chromosomes),$(mane_dir)/query_all/$(chr))

chess_dir := data/interim/chess/$(chess_version)
chess_introns := $(foreach chr,$(chromosomes),$(chess_dir)/introns_all/$(chr))
chess_query := $(foreach chr,$(chromosomes),$(chess_dir)/query_all/$(chr))

random_dir := data/interim/random/1.0
random_introns := $(foreach chr,$(chromosomes),$(random_dir)/introns_all/$(chr))
random_query := $(foreach chr,$(chromosomes),$(random_dir)/query_all/$(chr))

realignment_dir := data/interim/realignment
unique_exons := $(foreach chr,$(chromosomes),$(realignment_dir)/unique_exons/$(chr))
mapped_exons := $(foreach chr,$(chromosomes),$(realignment_dir)/mapped_exons/$(chr))
alignment_jobs := $(foreach genome,$(all_genomes),$(realignment_dir)/jobs/$(genome))
alignment_results := $(foreach genome,$(all_genomes),$(realignment_dir)/results/$(genome))
alignment_reports := $(foreach genome,$(all_genomes),$(realignment_dir)/extra_cons/$(genome))
identity_results := $(foreach genome,$(all_genomes),$(realignment_dir)/identity/$(genome))

#################################################################################
# FEATURES                                                                      #
#################################################################################

all_csv := data/processed/splice_sites.csv
model_csv := data/processed/model_out.csv
mane_csv := data/processed/mane/$(mane_version).csv
gencode_csv := data/processed/gencode/$(gencode_version).csv
refseq_csv := data/processed/refseq/$(refseq_version).csv
chess_csv := data/processed/chess/$(chess_version).csv
random_csv := data/processed/random/$(random_version).csv

#################################################################################
# COMMANDS                                                                      #
#################################################################################

#all: requirements $(maf_files) $(gnomad_tracks) $(ucsc_genomes_files) $(zoo_genomes_files) $(ncbi_genomes_files) $(gencode_query) $(refseq_query) $(mane_query) $(chess_query) $(random_query) data/interim/clinvar/clinvar.csv


#all: requirements $(maf_files) $(gnomad_tracks) $(ucsc_genomes_files) $(zoo_genomes_files) $(ncbi_genomes_files) $(gencode_query) $(refseq_query) $(mane_query) $(chess_query) $(random_query) data/interim/clinvar/clinvar.csv $(alignment_reports)


all: $(model_csv)

clean:
	rm -rf data/raw/annotation/*
	rm -rf $(realignment_dir)/unique_exons/*
	rm -rf $(realignment_dir)/mapped_exons/*
	rm -rf $(realignment_dir)/jobs/*
	rm -rf $(realignment_dir)/alignment_reports/*
	rm -rf $(realignment_dir)/results/*
	rm -rf $(realignment_dir)/extra_cons/*
	rm -rf $(gencode_dir)
	rm -rf $(chess_dir)
	rm -rf $(mane_dir)
	rm -rf $(refseq_dir)
	rm -rf $(random_dir)
	rm -f $(mane_csv) $(gencode_csv) $(refseq_csv) $(chess_csv) $(random_csv)
	rm -f $(all_csv) $(model_csv)
	rm -f data/processed/roc/*

## Install requirements

requirements:

## Run the model

$(model_csv): $(all_csv)
	$(PYTHON_INTERPRETER) src/models/model.py $(all_csv) $(model_csv) data/processed/roc

## Generate feature tables

$(all_csv): $(mane_csv) $(gencode_csv) $(refseq_csv) $(chess_csv) $(random_csv)
	$(PYTHON_INTERPRETER) src/features/merge.py $(mane_csv) $(gencode_csv) $(refseq_csv) $(chess_csv) $(random_csv) $(all_csv)

$(mane_csv): $(mane_gtf) $(mane_query) $(alignment_reports) $(gnomad_tracks) $(clinvar_tracks) $(phast_wig)
	$(PYTHON_INTERPRETER) src/features/generate_table.py $(mane_gtf) MANE $(mane_dir)/introns $(mane_dir)/query $(mane_dir)/introns data/interim/gnomad/ $(clinvar_tracks) $(realignment_dir)/extra_cons data/raw/phast  "$(all_genomes)" > $(mane_csv)

$(gencode_csv): $(gencode_gtf) $(gencode_query) $(alignment_reports) $(gnomad_tracks) $(clinvar_tracks) $(phast_wig)
	$(PYTHON_INTERPRETER) src/features/generate_table.py $(gencode_gtf) GENCODE $(gencode_dir)/introns $(gencode_dir)/query $(mane_dir)/introns data/interim/gnomad/ $(clinvar_tracks) $(realignment_dir)/extra_cons data/raw/phast  "$(all_genomes)" > $(gencode_csv)

$(refseq_csv): $(refseq_gtf) $(refseq_query) $(alignment_reports) $(gnomad_tracks) $(clinvar_tracks) $(phast_wig)
	$(PYTHON_INTERPRETER) src/features/generate_table.py $(refseq_gtf) RefSeq $(refseq_dir)/introns $(refseq_dir)/query $(mane_dir)/introns data/interim/gnomad/ $(clinvar_tracks) $(realignment_dir)/extra_cons data/raw/phast  "$(all_genomes)" > $(refseq_csv)

$(chess_csv): $(chess_gtf) $(chess_query) $(alignment_reports) $(gnomad_tracks) $(clinvar_tracks) $(phast_wig)
	$(PYTHON_INTERPRETER) src/features/generate_table.py $(chess_gtf) CHESS $(chess_dir)/introns $(chess_dir)/query $(mane_dir)/introns data/interim/gnomad/ $(clinvar_tracks) $(realignment_dir)/extra_cons data/raw/phast  "$(all_genomes)" > $(chess_csv)

$(random_csv): $(random_gtf) $(random_query) $(alignment_reports) $(gnomad_tracks) $(clinvar_tracks) $(phast_wig)
	$(PYTHON_INTERPRETER) src/features/generate_table.py $(random_gtf) Random $(random_dir)/introns $(random_dir)/query $(mane_dir)/introns data/interim/gnomad/ $(clinvar_tracks) $(realignment_dir)/extra_cons data/raw/phast  "$(all_genomes)" > $(random_csv)


## Realign exons to other genomes

$(alignment_reports): $(alignment_results) $(identity_results)
	$(PYTHON_INTERPRETER) src/data/parse_alignment_results.py $(human_genome) $(human_stats) $(realignment_dir)/results/$(@F) $(realignment_dir)/identity/$(@F) $(@)

$(alignment_results): $(alignment_jobs)
	$(PYTHON_INTERPRETER) src/data/run_alignment_jobs.py $(realignment_dir)/jobs/$(@F) $(genomes_dir)/$(@F).fa.gz $(@)

$(identity_results): $(mapped_exons)
	$(PYTHON_INTERPRETER) src/data/parse_identity_results.py $(realignment_dir)/mapped_exons $(@)

$(alignment_jobs): $(mapped_exons)
	$(PYTHON_INTERPRETER) src/data/generate_alignment_jobs.py $(human_genome) $(human_stats) $(mane_dir)/exons $(realignment_dir)/mapped_exons $(@)

$(mapped_exons): $(unique_exons) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/map_exons.py data/raw/maf/$(@F).maf $(realignment_dir)/unique_exons/$(@F) $(realignment_dir)/mapped_exons/$(@F)

## Generate unique exons from all datasets:
$(unique_exons): $(random_gtf) $(gencode_gtf) $(refseq_gtf) $(chess_gtf)
	$(PYTHON_INTERPRETER) src/data/unique_exons.py $(gencode_dir)/exons/$(@F) $(refseq_dir)/exons/$(@F) $(chess_dir)/exons/$(@F) $(mane_dir)/exons/$(@F) $(random_dir)/exons/$(@F) > $(realignment_dir)/unique_exons/$(@F)

## Get GRCh38:

$(human_stats):
	wget --directory-prefix=data/raw/human https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/$(human_version)/$(human_version)_assembly_report.txt

$(human_genome):
	wget --directory-prefix=data/raw/human https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/$(human_version)/$(human_version)_genomic.fna.gz

## Query random dataset

$(random_query): $(random_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/random/1.0 data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/random/1.0/ $(@F)

$(random_introns): $(random_gtf)
	src/data/make_introns.sh data/interim/random/1.0 $(@F)

## Generate random dataset

$(random_gtf): $(human_genome) $(human_stats) $(mane_introns)
	$(PYTHON_INTERPRETER) src/data/generate_random.py data/interim/mane/$(mane_version)/introns $(human_stats) $(human_genome) $(rand_limit) > $(basename $(random_gtf))
	gzip $(basename $(random_gtf))
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(random_gtf) data/interim/random/1.0/ "$(chromosomes)"

## Query GENCODE

$(gencode_query): $(gencode_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/gencode/$(gencode_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/gencode/$(gencode_version)/ $(@F)

$(gencode_introns): $(gencode_gtf)
	src/data/make_introns.sh data/interim/gencode/$(gencode_version) $(@F)

## Get GENCODE

$(gencode_gtf):
	wget --directory-prefix=data/raw/annotation https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$(gencode_version)/gencode.v$(gencode_version).annotation.gtf.gz
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(gencode_gtf) data/interim/gencode/$(gencode_version)/ "$(chromosomes)"

## Query RefSeq
$(refseq_query): $(refseq_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/refseq/$(refseq_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/refseq/$(refseq_version)/ $(@F)

$(refseq_introns): $(refseq_gtf)
	src/data/make_introns.sh data/interim/refseq/$(refseq_version) $(@F)

## Get RefSeq:

$(refseq_gtf): $(human_stats)
	wget --directory-prefix=data/raw/annotation https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/110/$(human_version)/$(refseq_filename)
	$(PYTHON_INTERPRETER) src/data/gtf_replace.py $(human_stats) $(refseq_gtf) data/raw/annotation/temp.gtf
	rm $(refseq_gtf)
	gzip data/raw/annotation/temp.gtf
	mv data/raw/annotation/temp.gtf.gz data/raw/annotation/$(refseq_filename)
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(refseq_gtf) data/interim/refseq/$(refseq_version)/ "$(chromosomes)"


## Query MANE
$(mane_query): $(mane_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/mane/$(mane_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/mane/$(mane_version)/ $(@F)

$(mane_introns): $(mane_gtf)
	src/data/make_introns.sh data/interim/mane/$(mane_version) $(@F)

## Get MANE:

$(mane_gtf):
	wget --directory-prefix=data/raw/annotation https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_$(mane_version)/MANE.GRCh38.v$(mane_version).ensembl_genomic.gtf.gz
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(mane_gtf) data/interim/mane/$(mane_version)/ "$(chromosomes)"


## Query CHESS
$(chess_query): $(chess_introns) $(maf_files)
	$(PYTHON_INTERPRETER) src/data/query_chr.py data/interim/chess/$(chess_version) data/raw/maf/$(@F).maf $(@F)
	$(PYTHON_INTERPRETER) src/data/collect_queries.py data/interim/chess/$(chess_version)/ $(@F)

$(chess_introns): $(chess_gtf)
	src/data/make_introns.sh data/interim/chess/$(chess_version) $(@F)

## Get CHESS:

$(chess_gtf):
	wget --directory-prefix=data/raw/annotation https://github.com/chess-genome/chess/releases/download/v.$(chess_version)/chess$(chess_version).GRCh38.gtf.gz
	$(PYTHON_INTERPRETER) src/data/split_gtf.py $(chess_gtf) data/interim/chess/$(chess_version)/ "$(chromosomes)"


## Get and parse the clinvar data

$(clinvar_tracks): data/raw/clinvar/clinvar.vcf.gz $(human_stats)
	$(PYTHON_INTERPRETER) src/data/parse_clinvar.py src/data/clinvar_pathogenic.txt $(human_stats) data/raw/clinvar/clinvar.vcf.gz $(clinvar_trracks)

data/raw/clinvar/clinvar.vcf.gz:
	wget --directory-prefix=data/raw/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

## Get the alignmnet

$(maf_files):
	wget --directory-prefix=data/raw/maf https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf/$(@F)

## Get and parse the gnomAD data

$(gnomad_tracks): $(gnomad_files)
	$(PYTHON_INTERPRETER) src/data/parse_gnomad.py data/raw/gnomad $@

$(gnomad_files):
	wget --content-on-error --directory-prefix=data/raw/gnomad https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/$(@F)

## Get the phastCons data

$(phast_wig):
	wget --directory-prefix=data/raw/phast https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons470way/hg38.470way.phastCons/$(@F)

## Get the genomes from the alignment hosted by UCSC

$(ucsc_genomes_files):
	wget --directory-prefix=data/raw/genomes https://hgdownload.soe.ucsc.edu/goldenPath/$(basename $(basename $(@F)))/bigZips/$(@F)

## Get the genomes from the alignment hosted by DNA Zoo

$(zoo_genomes_files):
	src/data/download_zoo.sh src/data/dna_zoo.txt $(@F) data/raw/genomes

## Get the genomes from the alignment hosted by NCBI

$(ncbi_genomes_files):
	src/data/download_ncbi.sh src/data/ncbi.txt $(@F) data/raw/genomes data/tmp
