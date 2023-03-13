for BC in bc2001 bc2002 bc2004 bc2007 bc2011 bc2019 bc2022 bc2015
do
    wget -O "$BC".bam https://downloads.pacbcloud.com/public/dataset/2021-11-Microbial-96plex/demultiplexed-reads/m64004_210929_143746."${BC}".bam
    samtools fastq -0 "$BC".fastq "$BC".bam
    gzip "$BC".fastq
done

