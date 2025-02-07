#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils qw(first_index);
use File::Basename;
use lib dirname (__FILE__);
use vcf2 qw( parse_vcf );
use ggModel qw(gm);

my %opt = ();
GetOptions( \%opt, 'sv=s', 'snv=s', 'ped=s', 'skipsv', 'osnv=s', 'osv=s' );
my @files = checkoptions(\%opt);
my $svfile = $opt{sv};
my $snvfile = $opt{snv};
my $outsnv = $opt{osnv};
my $outsv = $opt{osv};
open (LOG, '>', "compound.log") or die $!;
##
if (!defined $outsnv ) {
	$outsnv = "outsnv.vcf";
}

if (!defined $outsv ) {
	$outsv = "outsv.vcf";
}
open (OSNV, '>', $outsnv) or die $!;
##### READ PED #####################################
my $pedfile = $opt{ped};
my ($PED, $proband, $pedsize) = read_ped($pedfile);
my $father;
my $mother;
if ($pedsize > 2) {
	$mother = $PED->{$proband}->{MOTHER};
	$father = $PED->{$proband}->{FATHER};
}
else {
	system( "zcat ".$svfile." >".$outsv);
	system( "zcat ".$snvfile." >".$outsnv);
	print LOG "duo or single will not calculate compounds\n";
	exit;
}
####################################################
##### READ SV VCF ##################################
my $vcf = CMD::vcf2->new('file'=>$svfile );
my $ref = readSV($vcf);
my %SV = %$ref;
####################################################




#####
##### Loop over all SNVs variants, save compounding SNVs into SV-hash
#####

my $vcf2 = CMD::vcf2->new('file'=>$snvfile );
## Print header of altered SNV-vcf
print OSNV $vcf2->{header_str};
my $rankresult_meta = $vcf2->{meta}->{INFO}->{RankResult}->{Description};
$rankresult_meta =~ s/\|/:/g;
print LOG "SNV recalc initiated\n";
while ( my $b = $vcf2->next_var() ) {
	my $ref = $b->{REF};
	my $alt = $b->{ALT};
	my $pos = $b->{POS};
	my $chrom = $b->{CHROM};
	my $svref = \%SV;
	my $gene = $b->{INFO}->{Annotation};
	my @rs = split/:/,$b->{INFO}->{RankScore};
	# Updated genmod outputs float RankScore values
	my $rankscore = int($rs[1]);
	my $rankresult = $b->{INFO}->{RankResult};
	my $svc = $svref->{$chrom};
	my $snvGT = $b->{GT};


	## check every SV-variant on the same chromosome
	## if gene-match check compounding
	## else check overlap
	## else not confounding
	foreach my $SVvar (keys %{$svc}) {
		my $svPOS = $svc->{$SVvar}->{POS};
		my $svEND = $svc->{$SVvar}->{END};
		my $svlen = $svc->{$SVvar}->{SVLEN};
		my $svGENE = $svc->{$SVvar}->{GENE};
		my $SV_GT = $svc->{$SVvar};
		my @rankscore_sv = split/:/,$svc->{$SVvar}->{RankScore};
		my $rank_sv = $rankscore_sv[1];
		
		## ATM BNDs are not checked for compounds, have neither END nor SVLENsvGENE
		if (!defined $svEND) {next;}
		## If trio
		if ($pedsize > 2) {
			## Find GENE match, if array of SNV-genes match array of SV-genes
			if (defined $svGENE && defined $gene ) {
				my @SVgene = split/,/,$svGENE;
				my @SNVgene = split/,/,$gene;
				foreach my $n (@SNVgene) {
					if (grep /^$n$/, @SVgene) {
						print LOG "Gene match found SNV:$chrom:$pos SV:$chrom:$SVvar\n";
						my $check = compoundsolver($SV_GT, $snvGT);
						## compatible as compounds SV + SNV
						if ($check == 1 ) {
							## Check whether SNV needs to be rescored due to lowscoring/lack of SNV compound. Return new values
							my ($sub,$total) = snv_score_analyzer($rankscore,$rankresult,$rankresult_meta,$rank_sv);
							my $new_rankresult = join'|',@$sub;
							#print LOG "$new_rankresult\t$total\n";
							$b->{INFO}->{RankResult} = $new_rankresult;
							$b->{INFO}->{RankScore} = "$rs[0]:".$total;
							#$b->{INFO}->{Compounds} = $b->{INFO}->{Compounds}."$chrom:$svPOS\_SV";
							## SNV compounds to SVs
							if ($rankscore >= 12) {
								$svref->{$chrom}->{$SVvar}->{COMPOUND} = $check;
							}
						}
						# else {
						# 	print LOG "no compound eligable SNV-SV partner\n";
						# }
						last;
					}
				}
			}
			## Else if no gene match, match overlaps
			elsif ($pos >= $svPOS && $pos <= $svEND) {
				print LOG "Overlap match found SNV:$chrom:$pos SV:$chrom:$SVvar\n";
				my $check = compoundsolver($SV_GT, $snvGT);
				if ($check == 1 ) {
					## Check whether SNV needs to be rescored due to lowscoring/lack of SNV compound. Return new values
					my ($sub,$total) = snv_score_analyzer($rankscore,$rankresult,$rankresult_meta,$rank_sv);
					my $new_rankresult = join'|',@$sub;
					#print LOG "$new_rankresult\t$total\n";
					$b->{INFO}->{RankResult} = $new_rankresult;
					$b->{INFO}->{RankScore} = "$rs[0]:".$total;
					if ($rankscore >= 12 ) {
						$svref->{$chrom}->{$SVvar}->{COMPOUND} = $check;
					}
				}
			}
		}
		## Else single sample, never used
		else {
			print LOG "Overlap match found SNV:$chrom:$pos SV:$chrom:$SVvar\n";
			## Only compounds if overlapping SV
			if ($pos >= $svPOS && $pos <= $svEND) {
				if ($rankscore >= 12 ) {
					$svref->{$chrom}->{$SVvar}->{COMPOUND} = 1;
				}
			}
		}
	}
	my $tot_str = vcfstr($b, []);
	print OSNV $tot_str;
}
print LOG "SNV recalc done\n\n\n\n";
if ($opt{skipsv}) {
	print LOG "skipping SV recalc\n";
}
else {
	open (OSV, '>', $outsv) or die $!;
	## Print SV Header to OSV
	print OSV $vcf->{header_str};
	## Rescore SV based on SNV compounds ##
	#######################################
	my $sv_vcf = CMD::vcf2->new('file'=>$svfile );
	## Get meta-info RankResults to be able to modify inheritance models
	my $rankresult_meta_sv = $sv_vcf->{meta}->{INFO}->{RankResult}->{Description};
	$rankresult_meta_sv =~ s/\|/:/g;
	print LOG "SV recalc initiated\n";
	while ( my $c = $sv_vcf->next_var() ) {
			#print Dumper($c);
			my $ref = $c->{REF};
			my $alt = $c->{ALT};
			my $pos = $c->{POS};
			my $chrom = $c->{CHROM};
			my @meta = split':',$rankresult_meta_sv;
			my $index = first_index { $_ eq "inheritance_models" } @meta;
			my $sub = $c->{INFO}->{RankResult};
			my @rs = split/:/,$c->{INFO}->{RankScore};
			my $rankscore = $rs[1];
			$sub =~ s/\|/:/g;
			my @sub = split':',$sub;
			## If SV has a high scoring SNV-compound
			if ($SV{$chrom}{"$pos\_$ref\_$alt"}->{COMPOUND}) {
				$c->{INFO}->{GeneticModel} = "AR_comp";
				print LOG "SNV compound was found for SV $chrom:$pos\_$ref\_$alt, increase score by 12 for SV\n";
				## recalc Inheritance_Models add 12 again
				my $inher_score = $sub[$index] + 12;
				splice @sub, $index, 1, $inher_score;
				$rankscore = $rankscore+13;
				$c->{INFO}->{RankScore} = "$rs[0]:".$rankscore;
				$c->{INFO}->{RankResult} = join'|',@sub;
			}

			my $tot_str = vcfstr($c, []);
			print OSV $tot_str;

	}
	print LOG "SV recalc done\n";
}

close LOG;
sub compoundsolver {
	my ($SV, $SNV) = @_;
	my @snv_GT = @$SNV;

	my $check = 0;

	my %ID_VAR;
	foreach my $id (@snv_GT) {
		$ID_VAR{ $id->{_sample_id} } = $id->{GT};
	}
	## SV GENOTYPES
	my $sv_proband;
	if ($SV->{$proband} == 1 ) { $sv_proband = $SV->{$proband}; } 
	my $sv_mother;
	if ($SV->{$mother} == 1 ) { $sv_mother = $SV->{$mother}; }
	my $sv_father;
	if ($SV->{$father} == 1 ) { $sv_father = $SV->{$father}; }
	## Is SNV in proband as Heterozygous?
	if ($ID_VAR{$proband} eq "0/1" ) {
		## if above is true, is SV variant in proband?
		if (defined $sv_proband) {
			#print "$ID_VAR{$proband}\t$SV->{$proband}\n";
			## if father heterozygous and mother homozygous for ref
			if ($ID_VAR{$father} eq "0/1" && $ID_VAR{$mother} eq "0/0") {
				## does mother carry SV-variant in heterozygous form and father not?
				if (defined $sv_mother && !defined $sv_father) {
					#print "$proband:$ID_VAR{$proband}|0/$sv_proband\t$father:$ID_VAR{$father}|0/0\t$mother:$ID_VAR{$mother}|0/1\t";
					$check = 1;
				}
			}
			## if mother heterozygous and father homozygous for ref
			elsif ($ID_VAR{$father} eq "0/0" && $ID_VAR{$mother} eq "0/1") {
				## does mother carry SV-variant in heterozygous form and father not?
				if (!defined $sv_mother && defined $sv_father) {
					#print "$proband:$ID_VAR{$proband}|0/$sv_proband\t$father:$ID_VAR{$father}|0/1\t$mother:$ID_VAR{$mother}|0/0\t";
					$check = 1;
				}
			}
			## else neither parents are heterozygous or both are
			else {
				$check = 0;
				
			}
		}
		## else, variant not in proband, not of interest.
		else {
			$check = 0;
		}
	}
	return $check;
}

sub checkoptions {
	my %opt = %{ $_[0] };

	help() unless ($opt{sv} && $opt{snv});

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
	if ($count <= 2 ) {
		foreach my $ind (keys %PED) {
			$proband = $ind;
		}
	}
	return \%PED, $proband, $count;
}

sub readSV {
	my $vcf = shift;
	
	######
	###### SAVE all structural variants in hash. Name pos chrom gene which individuals: TODO MOVE TO SUBROUTINE
	######
	my %SV;
	while ( my $A = $vcf->next_var() ) {
		my %INFO;
		
		my $ref = $A->{REF};
		my $alt = $A->{ALT};
		my $pos = $A->{POS};
		my $chrom = $A->{CHROM};

		## VCF STRING ##
		$INFO{ vcf_str } = $A->{vcf_str};
		## POS ##
		$INFO{ POS } = $pos;
		## END ##
		$INFO{ END } = $A->{INFO}->{END};
		## GENE ##
		# my @gene;
		# foreach my $transcript ( @{ $A->{INFO}->{CSQ} } ) {
		# 	push @gene,$transcript->{SYMBOL};
		# }
		# $INFO{ GENE } = join(',',@gene);
		$INFO{ GENE } = $A->{INFO}->{CSQ}->[0]->{SYMBOL};
		## TYPE ##
		$INFO{ TYPE } = $A->{INFO}->{SVTYPE};
		## GeneticModel ##
		$INFO{ GeneticModel } = $A->{INFO}->{GeneticModel};
		## GeneticModel ##
		$INFO{ RankScore } = $A->{INFO}->{RankScore};
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

		## CHROMOSOME -> VARIANT -> ABOVE HASH
		$SV{$chrom}{"$pos\_$ref\_$alt"} = \%INFO;
		#print Dumper($A);
	}
	return \%SV;
}

sub snv_score_analyzer {
	my ($total, $sub, $meta, $rank_sv) = @_;
	#print $sub,"\n";
	$sub =~ s/\|/:/g;
	my @sub = split':',$sub;
	my @meta = split':',$meta;
	my $sum = 0;
	my $index = first_index { $_ eq "Inheritance_Models" } @meta;
	#print "$index\n";
	foreach my $res (@sub) {
		$sum = $res + $sum;
	}
	my $check = 0;
	## If no SNV got no penalty for low scoring Compound
	if ($sum == $total) {
		## If SNV had no compounds and Inheritance model was penalized
		if ($sub[$index] == -12) {
			$check = 1;
			print LOG "SNV has no SNV partner and inheritance was penalized \n";
		}
	}
	## else, SNV got penalized for having no high scoring SNV buddy
	else {
		$check = 2;
		print LOG "SNV has no high scoring partner\n";
		
	}
	#print LOG "$total\t$sum\t$sub[$index]\t$check\n";

	## If diff_res == 2 SNV variant has no highscoring SNV-compound, if svrank is above 12 increase SNV by 8 again
	if ($check == 2 && $rank_sv >=12 ) {
		$total = $total+8;
		print LOG "SNV had no high-scoring compound & SV was high scoring, increasing SNV score by 8 \n";
	}
	## If diff_res == 1 SNV variant has no SNV-compound at all, if svrank is above 12 increase SNV by 13 (-12 vs 1)
	elsif ($check == 1 && $rank_sv >=12 ) {
		splice @sub, $index, 1, "1";
		print LOG "SNV had no compound & SV was high scoring, increasing SNV score by 13 \n";
		$total = $total+13;
	}
	else {
		print LOG "SNV had no high scoring compound, SV score was below 12 = no score change for SNV\n";
	}

	return \@sub, $total;
}

sub vcfstr {
	my( $v, $sample_order ) = @_;
	
	my @all_info;
	my $tot_str = $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

	# Generate and print INFO field
	for my $info_key (@{$v->{INFO_order}}) {
		if($info_key eq "CSQ") {
			push @all_info, $info_key."=".$v->{_CSQstr};
		}
		else {
			push @all_info, $info_key."=".$v->{INFO}->{$info_key};
		}
	}
	$tot_str = $tot_str.join(";", @all_info)."\t";

	# Print FORMAT field
	$tot_str = $tot_str.join(":", @{$v->{FORMAT}})."\t";


	my %order;
	my $i=0;
	if( @$sample_order > 0 ) {
		$order{$_} = $i++ foreach @{$sample_order};
	}
	else {
		$order{$_->{_sample_id}} = $i++ foreach @{$v->{GT}};
	}
	# Print GT fields for all samples
	for my $gt ( sort {$order{$a->{_sample_id}} <=> $order{$b->{_sample_id}}} @{$v->{GT}}) {
		my @all_gt;
		for my $key ( @{$v->{FORMAT}} ) {
			push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
		}
		$tot_str = $tot_str.join(":", @all_gt)."\t";
	}
	$tot_str = $tot_str."\n";
	return $tot_str;
}
