#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils qw(first_index);
use File::Basename;
use lib dirname (__FILE__);
use vcf2 qw( parse_vcf );
use ggModel qw(gm);

my %opt = ();
GetOptions( \%opt, 'sv=s', 'snv=s', 'ped=s', 'annotsv=s', 'osv=s' );
my @files = checkoptions(\%opt);
my $svfile = $opt{sv};
my $snvfile = $opt{snv};
my $outsv = $opt{osv};

if (!defined $outsv ) {
	$outsv = "outsv.vcf";
}

open (OSV, '>', $outsv) or die $!;
##### READ PED #####################################
my $pedfile = $opt{ped};
my ($PED, $proband, $pedsize) = read_ped($pedfile);
my $father;
my $mother;
if ($pedsize > 2) {
	$mother = $PED->{$proband}->{MOTHER};
	$father = $PED->{$proband}->{FATHER};
}
####################################################
##### READ SV VCF ##################################
my $vcf = CMD::vcf2->new('file'=>$svfile );
my $ref = readSV($vcf);
my %SV = %$ref;

####################################################


my $desc;
my $annot;	
if ($opt{annotsv}) {
	($desc, $annot) = annotsv($opt{annotsv});
}

## Print header to STDOUT
## Add info fields
my @header = split/\n/,$vcf->{header_str};
my $c = 0;
foreach my $line (@header) {
	$c++;
	print OSV $line,"\n";
	if ($c == scalar(@header)-1) {
		print OSV "##INFO=<ID=INHER,Number=.,Type=Integer,Description=\"Inheritance, de-novo or not\">\n";
		print OSV "##INFO=<ID=Omim,Number=.,Type=Integer,Description=\"Reported Omim gene\">\n";
		print OSV "##INFO=<ID=GQC,Number=.,Type=Integer,Description=\"gatk genotype quality score\">\n";
		print OSV "##INFO=<ID=MANTAPENALTY,Number=.,Type=Integer,Description=\"big manta events nto supported by gatk/cnvnator\">\n";
		print OSV "##INFO=<ID=GeneticModel,Number=.,Type=String,Description=\"Genetic model for variant\">\n";
		print OSV "##INFO=<ID=QZERO,Number=.,Type=Float,Description=\"Fraction of reads mapped with MAPQ=0 in variants called only by CNVnator\">\n";
		print OSV "##INFO=<ID=RD,Number=.,Type=Float,Description=\"Estimated read depth from variants only called by CNVnator\">\n";
		print OSV "##INFO=<ID=LOWP,Number=.,Type=Float,Description=\"lowest p-value from variants only called by CNVnator\">\n";
		print OSV "##INFO=<ID=HOMHEM,Number=.,Type=String,Description=\"If variant follows geneticmodel and is a homo- or hemi-zygous deletion\">\n";
		if ($opt{annotsv}) {
			print OSV join"", @$desc;
		}
	}
}
## Print each original VCF entry with new annotations
foreach my $chrom (keys %SV) {
	foreach my $var (keys %{ $SV{$chrom} }) {
		my @vcf_split = split/\t/,$SV{$chrom}->{$var}->{vcf_str};
		my @info_field = split/;/,$vcf_split[7];
		my $compound = $SV{$chrom}->{$var}->{COMPOUND};
		my $omim = $SV{$chrom}->{$var}->{OMIM_GENES};
		my $inher = $SV{$chrom}->{$var}->{GENETIC_MODEL};
		my $qzero = $SV{$chrom}->{$var}->{QZERO};
		my $rd = $SV{$chrom}->{$var}->{RD};
		my $lowp = $SV{$chrom}->{$var}->{LOWP};
		my $homhem = $SV{$chrom}->{$var}->{HOMHEM};
		my $GQC = $SV{$chrom}->{$var}->{GQC};
		my $MANTAPENALTY = $SV{$chrom}->{$var}->{MANTAPENALTY};
		
		my @add_info;

		## If proband is female dont print Y-chromosome variants, CNVnator false positives
		if ($chrom eq 'Y') {
			if ($PED->{$proband}->{SEX} == 2) {
				next;
			}
		}
		if (defined $GQC) {
			push @add_info, "GQC=".$GQC;
		}
		if (defined $MANTAPENALTY) {
			push @add_info, "MANTAPENALTY=".$MANTAPENALTY;
		}
		if (defined $compound) {
			$inher="AR_comp";
		}
		if (defined $omim) {
			push @add_info, "Omim=Found";
		}
		if (defined $inher) {
			#unless ($inher eq "NA") {
				push @add_info, "GeneticModel=".$inher;
		#	} 
			
		}
		if (defined $qzero) {
			push @add_info, "QZERO=".$qzero;
			push @add_info, "RD=".$rd;
			push @add_info, "LOWP=".$lowp;
		}
		if (defined $homhem) {
			if (defined $inher) {
				if ($inher ne 'NA') {
					push @add_info, "HOMHEM=".$homhem;
				}
			} 
			else {
				push @add_info, "HOMHEM=".$homhem;
			}
		}
		## If AnnotSV is used
		if ($opt{annotsv}) {
			## morbidGenes
			if ( $annot->{$chrom}->{$var}->{morbidGenes} ) {
				push @add_info, "morbidGenes=".$annot->{$chrom}->{$var}->{morbidGenes};
			}
			## HI_DDDpercent
			if ( $annot->{$chrom}->{$var}->{HI_DDDpercent} ) {
				push @add_info, "HI_DDDpercent=".$annot->{$chrom}->{$var}->{HI_DDDpercent};
			}
			## DGV_CNV_Frequency
			if ( $annot->{$chrom}->{$var}->{DGV_CNV_Frequency} ) {
				push @add_info, "DGV_CNV_Frequency=".$annot->{$chrom}->{$var}->{DGV_CNV_Frequency};
			}
			## DDD_CNV_Frequency
			if ( $annot->{$chrom}->{$var}->{DDD_CNV_Frequency} ) {
				push @add_info, "DDD_CNV_Frequency=".$annot->{$chrom}->{$var}->{DDD_CNV_Frequency};
			}
			## IMH_AF
			if ( $annot->{$chrom}->{$var}->{IMH_AF} ) {
				push @add_info, "IMH_AF=".$annot->{$chrom}->{$var}->{IMH_AF};
			}
			## dbVar_status
			if ( $annot->{$chrom}->{$var}->{dbVar_status} ) {
				push @add_info, "dbVar_status=".$annot->{$chrom}->{$var}->{dbVar_status};
			}
			## AnnotSVrank
			if ( $annot->{$chrom}->{$var}->{AnnotSVrank} ) {
				push @add_info, "AnnotSVrank=".$annot->{$chrom}->{$var}->{AnnotSVrank};
			}
			## GnomadMax
			my $gnomadf = $annot->{$chrom}->{$var}->{gnomadmax};
			if ( $gnomadf) {
				unless ($gnomadf == -1) {
					push @add_info, "gnomad_svAF=".$annot->{$chrom}->{$var}->{gnomadmax};
				}
			}
		}
		print OSV join "\t",@vcf_split[0..6];
		print OSV "\t";
		push @info_field,@add_info;
		print OSV join ";", @info_field;
        print OSV "\t";
        #print everything after info field
        print OSV join "\t", @vcf_split[8..$#vcf_split];
		print OSV "\n";
	}
	
}

sub checkoptions {
	my %opt = %{ $_[0] };

	help() unless ($opt{sv} && $opt{annotsv});

}

sub help {
	my $in = shift;

	print "./compound_finder.pl --sv --snv --ped > sv_vcf_with_snvcompounds.vcf\n\n";
	print "--sv\t\tVCF containing structural variants from MANTA/CNVnator/TIDDIT/LUMPY REQUIRED\n";
	print "--snv\t\tVCF containg single nucleotide variants and indels REQUIRED\n";
	print "--ped\t\tPED-file containing all individuals in SNV/SV VCF REQUIRED\n";
	print "--annotsv\t\t.tsv file from AnnotSV, will be annotated into SV-VCF\n";
	exit;
}

sub read_ped {
	my $pedfile = shift;

	open (PED, $pedfile) or die $!;
	
	my $proband;
	my $mother;
	my $father;
	my %PED;
	
	while ( <PED> ) {
		my @line = split/\t/,$_;
		my %ind;
		
		$ind{FATHER} = $line[2];
		$ind{MOTHER} = $line[3];
		$ind{SEX} = $line[4];
		$ind{PHENO} = $line[5];
		$PED{$line[1]} = \%ind;

		unless ($line[2] eq "0" && $line[3] eq "0") {
			$proband = $line[1];
			$father = $line[2];
			$mother = $line[3];
		}
	}
	my $count = keys %PED;
	## if single sample, proband is obvious.
	if ($count == 1 ) {
		foreach my $ind (keys %PED) {
			$proband = $ind;
		}
	}
	elsif ($count == 2) {
		foreach my $ind (keys %PED) {
			unless ($PED{$ind}->{FATHER} eq 0 && $PED{$ind}->{MOTHER} eq 0) {
				$proband = $ind;
			}
		}
	}
	return \%PED, $proband, $count;
}

sub annotsv {
	my $tsv = shift;

	open (TSV, $tsv) or die $!;

	my @header;
	my $c;
	my @descriptions;
	my %av_annot;

	while ( <TSV> ) {
		$c++;
		chomp;
		my %add_anno;
		if ( $c == 1) {
			@header = split/\t/;
			push @descriptions, "##INFO=<ID=morbidGenes,Number=1,Type=String,Description=\"In morbidGenes?\">\n";
			push @descriptions, "##INFO=<ID=HI_DDDpercent,Number=1,Type=Float,Description=\"0-10% indicate a gene is more likely to exhibit haploinsufficiency\">\n";
			push @descriptions, "##INFO=<ID=DGV_CNV_Frequency,Number=1,Type=Float,Description=\"DGV_CNV_Frequency of matched CNV type and IDs with CNV\">\n";
			push @descriptions, "##INFO=<ID=DDD_CNV_Frequency,Number=1,Type=Float,Description=\"DDD_CNV_Frequency of matched CNV type\">\n";
			push @descriptions, "##INFO=<ID=IMH_AF,Number=1,Type=Float,Description=\"Ira M. Hall’s Allele frequency of matched CNV type\">\n";
			push @descriptions, "##INFO=<ID=dbVar_status,Number=1,Type=String,Description=\"dbVar_status pathogenesis\">\n";
			push @descriptions, "##INFO=<ID=AnnotSVrank,Number=1,Type=Integer,Description=\"AnnotSV ranking 1-5\">\n";
			push @descriptions, "##INFO=<ID=gnomad_svAF,Number=1,Type=Float,Description=\"gnomad popmax AF\">\n";
			
		}
		
		else {
			my @line = split/\t/;
			my @isdupdel = split/_/,$line[0];
			#### CHROM POS ####
			my $chrom = $line[first_index { $_ eq "SV chrom" } @header];
			my $start = $line[first_index { $_ eq "SV start" } @header];
			my $end = $line[first_index { $_ eq "SV end" } @header];	
			my $ref = $line[first_index { $_ eq "REF" } @header];
			my $alt = $line[first_index { $_ eq "ALT" } @header];
			my $ID = $line[first_index { $_ eq "ID" } @header];
			if (!defined $end) { $end = $alt; }
			#### MORBID GENE ####
			$add_anno{morbidGenes} = $line[first_index { $_ eq "morbidGenes" } @header];
			#### HAPLOTYPE INSUFFICIENCY ####
			$add_anno{HI_DDDpercent} = $line[first_index { $_ eq "HI_DDDpercent" } @header];
			#### DGV ####
			if ($isdupdel[3] eq 'DEL') {
				#my $s = $line[first_index { $_ eq "DGV_LOSS_IDs" } @header];
				$add_anno{DGV_CNV_Frequency} = $line[first_index { $_ eq "DGV_LOSS_Frequency" } @header];	
			}
			elsif ($isdupdel[3] eq 'DUP') {
				#my $s = $line[first_index { $_ eq "DGV_GAIN_IDs" } @header];
				$add_anno{DGV_CNV_Frequency} = $line[first_index { $_ eq "DGV_GAIN_Frequency" } @header];
			}
			#### DDD ####
			if ($isdupdel[3] eq 'DEL') {
				$add_anno{DDD_CNV_Frequency} = $line[first_index { $_ eq "DDD_DEL_Frequency" } @header];	
			}
			elsif ($isdupdel[3] eq 'DUP') {
				$add_anno{DDD_CNV_Frequency} = $line[first_index { $_ eq "DDD_DUP_Frequency" } @header];
			}
			#### IMH ####
			my @tmp = split/_/,$line[first_index { $_ eq "IMH_ID" } @header];
			if (defined $tmp[3]) {
				if ($isdupdel[3] eq $tmp[3] ) {
					$add_anno{IMH_AF} = $line[first_index { $_ eq "IMH_AF" } @header];	
				}
			}
			#### DBVAR STATUS ####
			my @dbvar =  split/;/,$line[first_index { $_ eq "dbVar_status" } @header];
			my $dbvar = join(',',@dbvar);
			$add_anno{dbVar_status} = $dbvar;
			#### ANNOTSV RANKING ####
			$add_anno{AnnotSVrank} = $line[first_index { $_ eq "AnnotSV ranking" } @header];
			$add_anno{gnomadmax} = $line[first_index { $_ eq "GD_POPMAX_AF"} @header];

			$av_annot{$chrom}{"$start\_$end"} = \%add_anno;
		}
		
	}
	return \@descriptions, \%av_annot;
}

sub readSV {
	my $vcf = shift;
	
	######
	###### SAVE all structural variants in hash. Name pos chrom gene which individuals: TODO MOVE TO SUBROUTINE
	######
	my %SV;
	while ( my $A = $vcf->next_var() ) {
		my %INFO;
		#print Dumper($A);
		my $ref = $A->{REF};
		my $alt = $A->{ALT};
		my $pos = $A->{POS};
		my $chrom = $A->{CHROM};
		my $id = $A->{ID};
		my $end = $A->{INFO}->{END};
		if (!defined $end) { $end = $alt; }

		## VCF STRING ##
		$INFO{ vcf_str } = $A->{vcf_str};
		## POS ##
		$INFO{ POS } = $pos;
		## END ##
		$INFO{ END } = $A->{INFO}->{END};
		## SVLEN ##
		$INFO{ SVLEN } = $A->{INFO}->{SVLEN};
		## GENE ##
		$INFO{ GENE } = $A->{INFO}->{CSQ}->[0]->{SYMBOL};
		## OMIM GENES
		$INFO{ OMIM_GENES } = $A->{INFO}->{OMIM_GENES};
		## TYPE ##
		$INFO{ TYPE } = $A->{INFO}->{SVTYPE};
		## GT ##
		my $gt = $A->{GT};
		for my $ind ( 0..scalar(@$gt)-1 ) {
			my $GT = $A->{GT}->[$ind]->{GT};
			$GT =~ s/\./0/g;
			my @sum_GT = split/\//,$GT;
			my $sum = 0;
			foreach (@sum_GT) { $sum = $_ + $sum;  }
			$INFO{ $A->{GT}->[$ind]->{_sample_id} } = $sum;
		}
		## ignore proband 0 variants
		if ($INFO{$proband} == 0 ) {next;}
		## GENETIC MODEL ##
		if ($pedsize > 2) {
			my $gm = model(\%INFO, $chrom);
			$INFO{ GENETIC_MODEL } = $gm;
		}
		## Find homozygous and hemizygous deletions
		if ($A->{INFO}->{SVTYPE} eq 'DEL') {
			if ($chrom eq 'X') {
				if ($PED->{$proband}->{SEX} == 1 && $INFO{$proband} > 0) {
					$INFO{HOMHEM} = "LOSS";
				}
				elsif ($PED->{$proband}->{SEX} == 2 && $INFO{$proband} > 1) {
					$INFO{HOMHEM} = "LOSS";
				}
			}
			elsif ($INFO{$proband} == 2) {
				$INFO{HOMHEM} = "LOSS";
			}
		}
		## Callers ##
		my $callers = callers($A->{vcf_str});

		## if only in manta, large (>100,000) mark as susp
		if ($A->{INFO}->{SVLEN} and $A->{INFO}->{SVTYPE} ne 'INS') { ##could be bnd, need to check svlen defined
			if (abs($A->{INFO}->{SVLEN}) > 100000) {
				if ($$callers[3] == 0 && $$callers[1] == 0 && $$callers[2] == 1 && $$callers[3] == 0 ) {
					$INFO{ MANTAPENALTY } = "1";
				}
			}
		}
		
		## if found only in gatk
		
		if ($$callers[1] == 0 && $$callers[2] == 0 && $$callers[3] == 1) {
			my $gatkqc = 0; my $gatk_count = 0;
			foreach my $gcsc (@{ $A->{GT} }) {
				next if ($gcsc->{QS} eq ".");
				$gatkqc = $gatkqc + $gcsc->{QS};
				$gatk_count++;
			}
			$gatkqc = ($gatkqc/$gatk_count);
			$INFO{ GQC } = int($gatkqc);
		}
		## if only found in cnvnator
		if ($$callers[0] == 1 && $$callers[1] == 0 && $$callers[2] == 0) {
			my @QZERO = split/,/,$A->{INFO}->{natorQ0};
			my @RD = split/,/,$A->{INFO}->{natorRD};
			## mean of merged CNVnator outputs. See mergeCNVnator.pl ~ author Björnieboy
			my $rd = cnvnatorfix(@RD);
			my $qzero = cnvnatorfix(@QZERO);
			$INFO{ QZERO } = $qzero;
			$INFO{ RD } = $rd;
			my @lowp;
			my @P1 = split/,/,$A->{INFO}->{natorP1};
			my @P2 = split/,/,$A->{INFO}->{natorP2};
			my @P3 = split/,/,$A->{INFO}->{natorP3};
			my @P4 = split/,/,$A->{INFO}->{natorP4};
			my $p1 = cnvnatorfix(@P1);
			my $p2 = cnvnatorfix(@P2);
			my $p3 = cnvnatorfix(@P3);
			my $p4 = cnvnatorfix(@P4);
			push @lowp,$p1,$p2,$p3,$p4;
			my $lowestP = findlow(@lowp);
			$INFO{ LOWP } = $lowestP;
		}
		## CHROMOSOME -> VARIANT -> ABOVE HASH
		$SV{$chrom}{"$pos\_$end"} = \%INFO;
		#print Dumper($A);
	}
	return \%SV;
}

sub model {
	my ($INFO,$chrom) = @_;

	my $m_GT = $INFO->{$mother}; 
	my $f_GT = $INFO->{$father}; 
	my $p_GT = $INFO->{$proband}; 
	my $type = $INFO->{TYPE};
	my $m_p = $PED->{$mother}->{PHENO};
	my $f_p = $PED->{$father}->{PHENO};
	chomp $f_p;
	chomp $m_p;
	my $x = 0;
	$x = 1 if $chrom eq 'X';
	my $gm = gm($PED->{$proband}->{SEX},$p_GT,$x,$f_GT,$f_p,$m_GT,$m_p);
	#print $PED->{$proband}->{SEX},$p_GT,$x,$f_GT,$f_p,$m_GT,$m_p,"  $gm \n";
	## at the moment rediculous genetic patterns including homozygous dominant traits are marked *, rank model version 1.0 cannot handle this.
	$gm =~ s/\*//;
	return $gm;
}

sub callers {
	my $vcfstr = shift;
	my @a = split/;/,$vcfstr;

	my @b = grep /^set.+$/, @a;
	my $cnvnator =  grep /cnvnator/, @b;
	my $tiddit =  grep /tiddit/, @b;
	my $manta =  grep /manta/, @b;
	my $gatk = grep /gatk/, @b;
	my @callers;
	push @callers, $cnvnator, $tiddit, $manta, $gatk;
	if ($cnvnator == 1 && $manta == 0 && $tiddit == 0) {

	}
	return \@callers;

}

sub findlow {
	my @in = @_;

	my $low = 10000000;
	
	foreach my $val (@in) {
		if ($val < $low) {
			$low = $val;
		}
	}
	return $low;
}


sub cnvnatorfix {
	my @in = @_;

	my $sum = 0;

	foreach my $val (@in) {

		$sum = $sum + $val;
	}
	my $mean = $sum / scalar(@in);
	return $mean;
}
