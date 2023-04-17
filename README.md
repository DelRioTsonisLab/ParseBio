# mergeParse.sh
This shell script is **_mergeParse.sh_**, which can join multiple BAM files from ParseBioscience outputs and transform them into a new file.bam to enter as imput in SouporCell tool.

## Installation can be accomplished through clone

**git clone** https://github.com/DelRioTsonisLab/ParseBio.git

## The options for using mergeParse.sh are:
```
bash  mergeParse.sh   -b  </folder_path/file1.bam>, [</folder_path/file2.bam>, ... ,</folder_path/fileN.bam>]    -g  <inserted_name>   -t   N    -m  <inserted_name>    -s <s1,s2,s3...>
```
For details of all arguments type: bash `mergeParse.sh  -h`
