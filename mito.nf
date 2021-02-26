#!/usr/bin/env nextflow

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.bam), file(row.bai)) }
	.into { mutserve_bam; eklipse_bam }


// run mutationserver on MT BAM files
process run_mutserve {
    cpus 4
    memory '16 GB'
    time '15m'
    
    input:
        set group, id, file(mito_bam), file(bai) from mutserve_bam

    output:
        set group, id, file("${mito_bam.baseName}.ms.vcf") into ms_vcfs_1, ms_vcfs_2
    
    """
    java -Xmx16G -Xms16G -jar /opt/bin/mutserve.jar call \
        --output ${mito_bam.baseName}.ms.vcf \
        --reference $params.rCRS_fasta \
        --deletions \
        --insertions \
        --threads ${task.cpus} \
        ${mito_bam}
    """

}

// split and left-align variants
process run_bcftools {


    input:
        set group, id, file(ms_vcf) from ms_vcfs_1

    output:
        set group, id, file("${ms_vcf.baseName}.adjusted.vcf") into adj_vcfs

    shell:
    '''
    sed 's/^chrM/M/' !{ms_vcf} | bgzip > !{ms_vcf}.tmp
    tabix -p vcf !{ms_vcf}.tmp
    bcftools norm -m-both -o !{ms_vcf}.split !{ms_vcf}.tmp
    perl -na -e 'chomp; @c=split(/\t/); if($c[4] eq "*"){ print "$c[0]\t".($c[1]-2)."\t".($c[1]-1)."\t$c[1]:$c[3]:$c[4]\n"; }' !{ms_vcf}.split > !{ms_vcf}.bed
    bedtools getfasta -fi /fs1/resources/ref/hg38/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna -name -bed !{ms_vcf}.bed -tab > !{ms_vcf}.prebase_del
    perl -na -e '$i++; @col=split(/\t/); if( $i == 1 ){ @del = `cat !{ms_vcf}.prebase_del`; \
    for $row(@del){ chomp($row); ($id,$base)=split/\t/,$row; $del{$id}=$base; } } \
    $id="$col[1]:$col[3]:$col[4]"; if($del{$id}){ $col[1]--; $col[3]=$del{$id}.$col[3]; $col[4]=$del{$id}; \
    print"$col[0]"; for $i(1..$#col){ print"\t$col[$i]"; } }else{ print; } ' !{ms_vcf}.split | bcftools sort | bgzip > !{ms_vcf}.split.prebase
    tabix -p vcf !{ms_vcf}.split.prebase
    bcftools norm -f /fs1/resources/ref/hg38/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna -o !{ms_vcf.baseName}.adjusted.vcf !{ms_vcf}.split.prebase    
    '''

}


// use python tool HmtNote for annotating vcf
process run_hmtnote {
    cpus 1
    memory '5GB'
    time '15m'


    input:
        set group, id, file(adj_vcf) from adj_vcfs

    output:
        set group, file("${adj_vcf.baseName}.hmtnote") into hmtnote_vcfs

    """
    source activate tools
    hmtnote annotate ${adj_vcf} ${adj_vcf.baseName}.hmtnote --offline
    """

}

//    hmtnote annotate /data/bnf/dev/paul/Hans/MITO/simulate_vars.vcf ${adj_vcf.baseName}.hmtnote --offline

process annotate_vep {
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	cpus 2
	tag "$group"
	memory '20 GB'
	time '5h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, file(vcf) from hmtnote_vcfs

	output:
		set group, file("${group}.vep.vcf") into vep

	"""
	vep \\
		-i ${vcf} \\
		-o ${group}.vep.vcf \\
		--offline \\
		--everything \\
		--merged \\
		--vcf \\
		--no_stats \\
		--fork ${task.cpus} \\
		--force_overwrite \\
		--plugin CADD,$params.CADD \\
		--plugin LoFtool \\
		--plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS \\
		--fasta $params.VEP_FASTA \\
		--dir_cache $params.VEP_CACHE \\
		--dir_plugins $params.VEP_CACHE/Plugins \\
		--distance 200 \\
		-cache \\
		-custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_popmax,AF,popmax \\
		-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_popmax,AF,popmax \\
		-custom $params.PHYLOP \\
		-custom $params.PHASTCONS
	"""
}

// run haplogrep 2 on resulting vcf
process run_haplogrep {

    memory '16 GB'

    input:
        set group, id, file(ms_vcf) from ms_vcfs_2

    output:
        set group, id, file("${ms_vcf.baseName}.hg2.vcf.dot") into hg2_dots

    """
    java  -Xmx16G -Xms16G -jar /opt/bin/haplogrep.jar classify \
        --in ${ms_vcf} \
        --out ${ms_vcf.baseName}.hg2.vcf \
        --format vcf \
        --lineage 1
    """

}

// create a pdf from dot file (created with haplogrep --lineage flag)
process run_dot2pdf {

    input:
        set group, id, file(hg2_dot) from hg2_dots

    output:
        set group, id, file("${hg2_dot}.png") into hg_pngs

    """
    dot ${hg2_dot} -Tps2 > ${hg2_dot}.ps
    gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -r1200 -dDownScaleFactor=3 -sOutputFile=${hg2_dot}.png ${hg2_dot}.ps
    """

}

//    gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -r1000 -sDEVICE=pngalpha -sOutputFile=${hg2_dot}.png ${hg2_dot}.ps


// // use eKLIPse for detecting mitochondrial deletions
// process run_eklipse {

//     container = '/data/bnf/dev/paul/NextFlow/mito_sv/container/mitoSV_2020-10-01.sif'

//     input:
//     file mito_bam from mito_bams2
//     set sampleId, file(bam) from bam_ch2

//     """
//     source activate htslib10
//     echo "${mito_bam}\tsample" > infile.txt
//     python /eKLIPse/eKLIPse.py \
//     -in infile.txt \
//     -ref /eKLIPse/data/NC_012920.1.gb
//     mv eKLIPse_*/eKLIPse_deletions.csv ./${sampleId}_deletions.csv
//     mv eKLIPse_*/eKLIPse_genes.csv ./${sampleId}_genes.csv
//     mv eKLIPse_*/eKLIPse_sample.png ./${sampleId}_plot.png
//     perl /data/bnf/scripts/hetplasmid_frequency_eKLIPse.pl --bam ${mito_bam} --in ${sampleId}_deletions.csv
//     """

// }
