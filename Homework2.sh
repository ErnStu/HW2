#!/bin/sh

#failų atsisiuntimas
prefetch ERR204044 SRR15131330 SRR18214264
fastq-dump --split-files ERR204044 SRR15131330 SRR18214264

threads=6

for i in ../inputs/*_1.fastq
do
  R1=${i}
  R2="../inputs/"$(basename ${i} _1.fastq)"_2.fastq"
  #fastqc raw failams
  fastqc -t ${threads} ${R1} ${R2} -o ../outputs/fastqc
  trim_galore -j ${threads} --paired --quality 20 --length 20 --illumina -o ../outputs/ ${R1} ${R2}
done

for i in ../outputs/*_1_val_1.fq
do
  R1=${i}
  ID=$(basename $i _1_val_1.fq)
  R2="../outputs/"$(basename ${i} _1_val_1.fq)"_2_val_2.fq"
  #fastqc trimmed failams
  fastqc -t ${threads} ${R1} ${R2} -o ../outputs/fastqc

  #Genome assambly naudojant spades
  spades.py -t 6 -o ../outputs/Spades/${ID}/ --pe1-1 ${R1} --pe1-2 ${R2}

  #Mapping
  bwa index ../outputs/Spades/${ID}/contigs.fasta
  bwa mem -t 6 ../outputs/Spades/${ID}/contigs.fasta ${R1} ${R2} -o ../outputs/Mapping/${ID}/${ID}.sam
  samtools view -bS ../outputs/Mapping/${ID}/${ID}.sam -@ 6 | samtools sort -@ 6 - -o ../outputs/Mapping/${ID}/${ID}.bam

  #Mapping fraction
  samtools flagstat ../outputs/Mapping/${ID}/${ID}.bam -@ 6
  #Genome coverage
  samtools depth -a ../outputs/Mapping/${ID}/${ID}.bam | awk '{sum += $3} END {print sum / NR}'

  #Blast duombazių sukūrimas
  makeblastdb -in ../outputs/Spades/${ID}/contigs.fasta -dbtype nucl -out ../outputs/BlastDB/${ID}/${ID}_db

  #Nucleotide blast
  blastn -query ../CP015498_nucl.fasta -db ../outputs/BlastDB/${ID}/${ID}_db -out ../outputs/BlastResults/${ID}/nuclResults.txt -outfmt "6 qseqid sseqid pident length mismatch 
        gapopen qstart qend sstart send evalue bitscore"
  blastn -query ../CP015498_nucl.fasta -db ../outputs/BlastDB/${ID}/${ID}_db -out ../outputs/BlastResults/${ID}/nuclResults_nonformatted.txt

  #Protein blast
  tblastn -query ../CP015498_prot.fasta -db ../outputs/BlastDB/${ID}/${ID}_db -out ../outputs/BlastResults/${ID}/protResults.txt -outfmt "6 qseqid sseqid pident length mismatch 
        gapopen qstart qend sstart send evalue bitscore"
  tblastn -query ../CP015498_prot.fasta -db ../outputs/BlastDB/${ID}/${ID}_db -out ../outputs/BlastResults/${ID}/protResults_nonformatted.txt
done

#multiqc raw ir trimmintiems failams
multiqc -o ../outputs/fastqc ../outputs/fastqc

#alternatyvus assembly naudojant abyss
abyss-pe name=ERR204044 k=64 in='../ERR204044_1_val_1.fq ../ERR204044_2_val_2.fq'
abyss-pe name=SRR15131330 k=64 in='../SRR15131330_1_val_1.fq ../SRR15131330_2_val_2.fq'
abyss-pe name=SRR18214264 k=64 in='../SRR18214264_1_val_1.fq ../SRR18214264_2_val_2.fq'

#Viskas kita buvo atlika naudojantis UseGalaxy, RAST, GeneMark, SeaView.

#ATSAKYMAI

#Data QA/QC

#Perform fastqc to evaluate your data quality. In your code write a short comment about your data quality.
#Matome, kad yra daug pasikartojančių sekų SRR15131330, tačiau jų quality yra labai aukštas. 
#Su ERR204044 viskas atrodo pakankamai gerai, pasikartojančių sekų procentas atrodo normalus, taip pat ir quality yra geras.
#Su SRR18214264 pasikartojančių sekų procentas normalus, tačiau quality galėtų būti kiek geresnis, matome, jog kai kur jis yra apie 34, 
#kai kur apie 27, ir net kartais mažesnis nei 20, kas nėra gerai.

#Repeat fastqc and evaluate if there were any changes. In your code write a short comment about your results.
#Po trimmingo, SRR15131330 vistiek išliko didelis pasikartojančių sekų procentas, quality toliau išliko labai aukštas.
#Su ERR204044 viskas išliko daug maž taip pat, o su SRR18214264 šiek tiek pagerėjo quality, tačiau galėtų būti ir geriau.
#taip pat matome, su visais mėginiais po trimminimo "per base sequence content" grafe linijos yra šiek tiek iškraipytos, nors 
#turėtų būti daug maž paralelios. Taip pat nukirpti visi adapteriai ir po to, nuskaitimai, kurie buvo trumpesni už 20bp buvo pašalinti.

#Genome assembly

#In your code describe your results: do all assemblies look the same? What are the main results?
#Atlikę trijų genomų surinkimus dviem skirtingais būdais matome, kad yra šiokių tokių panašumų, tačiau ir skirtumų.
#Visi genomų rinkiniai turi panašų fraction procentą, kad reiškia, jog daug maž tokia pati genomo dalis buvo apimta.
#Taip pat visų rinkinių duplikacijų koeficientas buvo artimas 1, tai reiškia, jog buvo minimalus kiekis pasikartojančių sekų. 
#Largest alignment ir total aligned lenght tarp visų buvo labai panašus, kas reiškia, kad panašūs genetiniai regionai buvo apimti. 
#Contigų kiekis taip pat yra panašus, matome, kas skiriasi yra tai, jog SRR15131330 mėginyje largest contig yra dvigubai mažesnis nei kitų mėginių(~50000bp ir ~100000bp), 
#tačiau jis yra artimesnis largest alignment, manau tai galėtų reikšti, kad su šiuo mėginiu atliktas genomo surinkimas yra tikslesnis mūsų reference genomui. Taip pat
#matome, kad šiame mėginyje yra mažiausias misassemblies kiekis.

#Using quast results select the best/better assembly for each sample. In your code explain why you chose a specific assembly.
#Palyginimai:
#ERR204044 - atliktame surinkime su Spades fraction yra geresnis nei su Abyss (77% ir 75%), duplikacijų koeficientas taip pat (1.003 ir 1.015), 
#o tiek largest alignment, tiek total lenght aligned yra beveik identiški. Pasirinkau Spades.
#SRR15131330 - atliktame surinkime su Spades fraction yra geresnis nei su Abyss (82% ir 76%), tačiau duplikacijų koeficientas yra geresnis Abyss(1.032 ir 1.016)
#o tiek largest alignment, tiek total lenght aligned yra beveik identiški. Pasirinkau Spades, dėl geresnio Fraction procento.
#SRR18214264 - atliktame surinkime su Spades fraction yra geresnis nei su Abyss (77% ir 73%), tačiau duplikacijų koeficientas yra geresnis Abyss(1.026 ir 1.017),
#largest alignment ir aligment lenght buvo didesnis atlient Spades programa. Pasirinkau Spades.

#Using appropriate mapper, map original reads to you assemblies. Evaluate mapping fraction as well 
#as genome coverage from mapped reads (and as in other questions, provide a comment on your results).

#ERR204044 99.72% 293
#SRR15131330 99.84% 2194.12
#SRR18214264 99.73% 276.485

#Visi mėginiai ERR204044, SRR15131330, SRR18214264 turėjo labai aukštą mapping procentą, atitinkamai: 99.72%, 99.84% ir 99.73%, kas parodo, 
#jog genomų rinkiniai buvo surinkti sėkmingai atsižvelgiant į reference genomą. Coverage atitinkamai buvo toks: 293, 2194 ir 276, 
#įtarimų sukėlė SRR15131330 mėginio covarage skaičius, nes jis atrodo labai didelis ir gerokai išsiskiria iš kitų.

#Genome analysis and annotation

#Using Gepard tool create dotplots to show similarities/dissimilarities between your samples. 
#Describe, your results (in your code). In the last question you will have to upload dotplots, so save them.

#Gepard programa susikūrus ir palyginus visus gautus dotplots, visose nuotraukose gavau gražią istrižą liniją, kad parodo, jog visi mėginiai yra panašūs.
#Vieniniteliai skirtumai, kuriuos pastebėjau yra tai, jog dotplots tarp ERR204044 ir SRR15131330 linijoje yra mažiau trūkių, 
#nei tarp SRR15131330 ir SRR18214264 bei SRR18214264 ir ERR204044, tačiau nemanau, jog tai galėtų kažką labai reikšti. 
#Taip pat linija tarp ERR204044 ir SRR18214264 mėginių yra šiek tiek tiesesnė ir ryškesnė, nei kituose dotplot, kas manau taip pat neturi labai didelės įtakos.

#Using BUSCO analysis tool, evaluate your assemblies. Provide a short comment on BUSCO results.

#Visų trijų surinkimų matomi BUSCO rezultatai parodo, didelį completeness procentą (99%), taip pat fragmented ir missing procentai yra labai maži(po 0,5%), 
#todėl manau galima teigti, kad genomų surinkimai buvo atlikti gerai.

#At the moment you should have 3 gene predictions. Compare and describe them (you should compare number of predicted genes and genes overlap. 
#You don't have to include functional annotations).

#ERR204044 - Rast: 2564, GeneMark: 2310, Blast: 1813
#SRR15131330 - Rast: 2774, GeneMark: 2552, Blast: 1909
#SRR18214264 - Rast: 2518, GeneMark: 2330, Blast: 1811

#Matome pasikartojančią tendenciją, tai, jog daugiausiai genų surandama su RAST, o mažiausiai su Blast. Mano nuomone, kodėl taip galėtų būti yra tai, 
#jog su Blast mes randame tik tuos genus, kurie yra mūsų reference genome, o su RAST visus įmanomus.
#Taip pat matome, kad genų skaičius ERR204044 ir SRR18214264 yra gana panašus, o SRR15131330 mėginyje šiek tiek skiriasi nuo jų.

#Compare your phylogenetic trees. do they look the same? Do they show same/identical clusters?

#Abu medžiai atrodo panašiai, matome pasikartojančią tendenciją, kad SRR15131330 mėginys šiek tiek skiriasi nuo ERR204044 ir SRR18214264 mėginių, 
#tai labai aiškiai matome 16S sekų medyje. ERR204044 ir SRR18214264 mėginiai yra suportuoti, o SRR15131330 yra kartu su reference genomu. Multi-gene 
#medyje mano nuomone yra šiek tiek netikslumų. Nors ir mūsų surinktų genomų baltymai yra tuose pačiuose cluteriuose, tačiau kai kur referentinio
#genomo baltymai yra šiek tiek atsiskyrę: tai matome su DNA topoisomerase I, DNA topoisomerase IV subunit A ir su phosphonate ABC transporter ATP-binding protein. 
#Outgroup abiejuose medžiuose yra atsiskyręs nuo mūsų mėginių cluterio, tačiau yra toje pačioje atšakoje. 

#Using all data you got, can you identify if any of you genomes are more similar to each other than to the third one (or reference genome)? Explain your ideas.
#Remdamasis visais turimais duomenimis, manau, kad ERR204044 ir SRR18214264 mėginiai yra panašesni vienas į kitą ir skiriasi nuo SRR15131330 mėginio. Taip sprendžiu, nes 
#pradedant nuo genomo surinkimo, SRR15131330 mėginio largest contig (~50000bp) skyrėsi nuo kitų dviejų mėginių (abu po ~100000bp), taip pat susikuriant 
#Gepard programa dotplots galima buvo įžvelgti, kad linija tarp ERR204044 ir SRR18214264 buvo kiek tiesesnė ir ryškesnė. Suskaičiavus genus RAST, GeneMark ir Blast įrankiais, 
#taip pat matome, kad genų kiekis tarp ERR204044 ir SRR18214264 buvo vos ne identiškas, o nuo SRR15131330 mėginio skyrėsi per ~100.
#Iš RAST sukurtų diagramų matome, kad ERR204044 ir SRR18214264 panašumas yra tikrai didelis, vos ne visa diagrama padengta tamsesne mėlyna spalva, 
#kas reiškia arti 100% esantį atitikimą. Kitose dviejose diagramose: tarp ERR204044 ir SRR15131330 bei SRR18214264 ir SRR15131330, matome, 
#jog žiedinės diagramos yra labiau žalsvos su protarpiais pasitaikančia mėlyna spalva, kad reiškia, jog atitikimas yra arčiau 95%. Ir galiausiai pažiūrėję į
#susikurtus medžius taip pat galime pastebėti, kad daug kur pasitaikė, kad ERR204044 ir SRR18214264 mėginiai yra suportuoti, o SRR15131330 yra kartu su reference genomu.