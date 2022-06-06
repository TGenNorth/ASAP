#!/bin/sh

#$1 name of .xml from asap
#$2 path to first output transform
#$3 path to second output transform

export PYTHONHASHSEED=1

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

pathToTranforms=${SCRIPTPATH}/../output_transforms

formatOutputName=$(echo "$1" | sed 's/_analysis.xml/.html/g')

reformattedName=$(echo "$1" | sed 's/.xml/_Reformatted.xml/g')

reReformattedName=$(echo $reformattedName | sed 's/.xml/_2.xml/g')

formatOutput -s ${pathToTranforms}/new_Tick_AmpSeq_Output_1.xsl  -x $1 -o ./$formatOutputName

echo "first output transform done"

reformatXML -x $1

echo "reformatting done"

python ${SCRIPTPATH}/tick_formatter.py -x $reformattedName -o $reReformattedName

echo "re-reformatting done"

formatOutput -s ${pathToTranforms}/new_Tick_AmpSeq_Output_2.xsl -x $reReformattedName

mkdir ./fasta

formatOutput -s ${pathToTranforms}/new_Tick_AmpSeq_Output_3.xsl -x $reReformattedName

echo "generating fastas"

assays=(BbovSBP2I_SBP2I01 BbovSBP2I_SBP2I02 BbovSBP2I_SBP2I03 BbovSBP2I_SBP2I04 BbovSBP2I_SBP2I05 BbovSBP2I_SBP2I06 BbovSBP2I_SBP2I07 BbovSBP2II_SBP2II01 BbovSBP2II_SBP2II02 BbovSBP2II_SBP2II03 BbovSBP2II_SBP2II04 BbovSBP2II_SBP2II05 BbovSBP2II_SBP2II06 BbovSBP2II_SBP2II07 AnaplasmamarginaleMSP1beta_Amarginale01 AnaplasmamarginaleMSP1beta_Amarginale02 BabesiabovisVESA1a_Bbov01  BabesiabovisVESA1a_Bbov02 BabesiabovisVESA1a_Bbov03 BabesiabovisVESA1a_Bbov04 BabesiabovisVESA1a_Bbov05 BabesiabovisVESA1a_Bbov06 BabesiabovisVESA1a_Bbov07 BabesiabovisVESA1a_Bbov08 BabesiabovisVESA1a_Bbov09 BabesiabovisVESA1a_Bbov10 BabesiabovisVESA1a_Bbov11 BabesiabovisVESA1a_Bbov12 BabesiabovisVESA1a_Bbov13 BabesiabovisVESA1a_Bbov14 BabesiabovisVESA1a_Bbov15 BabesiabovisVESA1a_Bbov16 BabesiabovisVESA1a_Bbov17 BabesiabovisVESA1a_Bbov18 BabesiabovisVESA1a_Bbov19 BabesiabovisVESA1a_Bbov20 BabesiabovisVESA1a_Bbov21 BabesiabovisVESA1a_Bbov22 BabesiabovisVESA1a_Bbov23 BabesiabovisVESA1a_Bbov24)

cd ./bowtie2

samples=(*bam)

out_dir=../fasta/temp

mkdir ${out_dir}

#in outdir create a bam file for each sample and region of interest that contains all the reads that were in that sample and region of interest
for assay in "${assays[@]}"; do
  for sample in "${samples[@]}"; do
    samtools view -h $sample ${assay} > ${out_dir}/${sample}_${assay}
  done
done


#for each of the files created above turn it into a fasta
for assay in "${assays[@]}"; do
  for sample in "${samples[@]}"; do
    {
    samtools fasta ${out_dir}/${sample}_${assay} > ${out_dir}/${sample}_${assay}.fasta
    } 1>/dev/null 2>&1
  done
    #call python script that parses these newly created fasta files
    python ${SCRIPTPATH}/tick_formatter.py --switch -q ${out_dir} -a ${assay} -s "${samples[@]}"
    #/${sample}_${assay}.fasta
done

cd ../fasta

rm -r ./temp

echo "done"
