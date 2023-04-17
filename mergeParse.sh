#!/bin/bash
usage(){
        cat <<EOF

Usage: bash joinBam.sh [Options]  or ./joinBam.sh [Options]  or ./place/with/the/joinBam.sh [Options] or ./path/where/is/the/joinBam.sh [Options]

Options:

  -h                 Show this help message and exit.

  -b BAM_FILES       Input bam files separated by comma WITHOUT any space.
                     example_1: </folder_path/file1.bam>, [</folder_path/file2.bam>, ... ,</folder_path/fileN.bam>]
                     example_2: bam1.bam,bam2.bam,bam3.bam  ...

  -g genome_name     Name that you put in the flag --genome_name <inserted_name>
                     example:Parse Bioscience
                          split-pipe \\
                          --mode mkref \\
                          --genome_name hg38 \\
                          ...

  -t threads         Number of threads to use in samtools
                    
  -s suffix          Suffixes will be added (__1 to __7) in the same order that was specified for each sublibrary.
                     The {first sublibrary} will be <s1>, The {second sublibrary} will be <s2>, ...
                     example: s1,s2,s3,...

  -m cell_metadata   Input <cell_metadata.csv> will be transformed in format barcode.tsv            
                     This file is located in <ref_name>/subs_all_comb/all-well/DGE_filtered/cell_metadata.csv


EOF
}


while getopts "hb:g:t:s:m:" opt; do
        case $opt in
                h)      # Call the usage function to print the help and how to use
                        usage
                        exit
                                ;;
 
                b)       # $OPTARG is the value after the -b flag, separated by commas without space between them
                        BAMFILES=$OPTARG
                                ;;                              

                g)       # $OPTARG is the value after the -g flag
                        GENOME_NAME=$OPTARG
                                ;;                               
                     
                t)       # $OPTARG is the value after the -t flag
                        THREADS=$OPTARG
                                ;;

                s)       # $OPTARG is the value after the -s flag
                        SUFFIX=$OPTARG
                                ;;

                m)       # $OPTARG is the value after the -m flag
                        METADATA=$OPTARG
                                ;;                

                *)       # Putting other parameters
                        usage
                        exit
                                ;;                              
        esac
done

#Non-option arguments
if [ $OPTIND -eq 1 ];then
    echo "NO OPTIONS WERE PASSED! See below for Usage:"
    usage
    exit
fi

#Required parameter GENOME_NAME
if [[ -z $GENOME_NAME ]];then
    echo "I NEED -g FLAG! See below for Usage:"
    usage
    exit
fi

#Required parameter SUFFIX
if [[ -z $SUFFIX ]];then
    echo "I NEED -s FLAG! See below for Usage:"
    usage
    exit
fi

#Required parameter METADATA
if [[ -z $METADATA ]];then
    echo "I NEED -m FLAG! See below for Usage:"
    usage
    exit
fi

#Default for parameter threads
if [[ -z $THREADS ]];then
    THREADS=8
fi

echo "Getting all the inputs..."

processing() {
        #Using my variables "b:g:t:s:m" from while getopts
        #Because the argument in the flag -b is with comma, use IFS=,
        #The parameter, read -a mylist, sets up a list with the values separated by comma

        IFS=, read -a mylist <<<"${BAMFILES}"
        
        IFS=, read -a mysuffix <<<"${SUFFIX}"

        #Find the lenght
        x=${#mylist[@]}
        echo "Your inputs.bam are: $x files"

        #Convert to a sam file and get up the header
        samtools view -H ${mylist[0]} | sed -e "s/${GENOME_NAME}_//g" > barcode_headAligned.merged.sam

        #Loop for every bam files
        # ${mylist[@]} has all the item of my list

        index=0

        for i in ${mylist[@]};do
                echo "Modifiying data of sublibrary $i ..."
                subindex=${mysuffix[$index]}
                samtools view $i | sed -r -e "s@(CB:.:.._.._..)@\1__$subindex@g" -e "s@${GENOME_NAME}_@@g" -e 's@pN:@UB:@g' >> barcode_headAligned.merged.sam       
                index=$((index+1))
        done

        echo "Doing the convertion to bam file"
        #Convert to a bam file
        samtools view -@ ${THREADS} -b -o barcode_headAligned.master.bam barcode_headAligned.merged.sam

        #Remove the file.sam
        rm barcode_headAligned.merged.sam

        #Sort my file.bam
        samtools sort -@ ${THREADS} -o headAligned.master.bam barcode_headAligned.master.bam

	#Remove the file.bam
        rm barcode_headAligned.master.bam

        echo "Doing the barcode-Parse-test.tsv"
        #Doing my specified variable
        index="__${mysuffix[0]}"
        for j in ${mysuffix[@]:1};do
                index="$index|__$j"
        done

        #Create the tsv file
        awk -F ',' '{print $1}' ${METADATA} | grep -E ${index} > barcode-Parse-test.tsv

        echo "Successfully finished!"
}

processing

