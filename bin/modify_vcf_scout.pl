#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use List::Util qw( min max );
# use List::MoreUtils qw(first_index);
my %vcf_meta;
my @vcf_data;
my @head;
open(VEP, $ARGV[0]);


my %maxentscan = (
    "MES-NCSS_downstream_acceptor" => 1,"MES-NCSS_downstream_donor"=> 1,"MES-NCSS_upstream_acceptor"=> 1,
    "MES-NCSS_upstream_donor"=> 1,"MES-SWA_acceptor_alt"=> 1,"MES-SWA_acceptor_diff"=> 1,"MES-SWA_acceptor_ref"=> 1,
    "MES-SWA_acceptor_ref_comp"=> 1,"MES-SWA_donor_alt"=> 1,"MES-SWA_donor_diff"=> 1,"MES-SWA_donor_ref"=> 1,
    "MES-SWA_donor_ref_comp"=> 1,"MaxEntScan_alt"=> 1,"MaxEntScan_diff"=> 1,"MaxEntScan_ref"=> 1
);

my %clinmod = (
    "Pathogenic" => "_5_",
    "Likely_pathogenic" => "_4_",
    "Likely_benign" => "_3_",
    "Benign" => "_2_",
    "Uncertain_significance" => "_0_",
    "not_provided" => "_1_", 
    "drug_response" => "_6_"
);


my %rank = (
    'transcript_ablation' => 1,
    'initiator_codon_variant' => 2,
    'frameshift_variant' => 3,
    'stop_gained' => 4,
    'start_lost' => 5,
    'stop_lost' => 6,
    'splice_acceptor_variant' => 7,
    'splice_donor_variant' => 8,
    'inframe_deletion' => 9,
    'transcript_amplification' => 10,
    'splice_donor_5th_base_variant' => 11,
    'splice_region_variant' => 12,
    'splice_donor_region_variant' => 13,
    'splice_polypyrimidine_tract_variant' => 14,
    'missense_variant' => 15,
    'protein_altering_variant' => 16,
    'inframe_insertion' => 17,
    'incomplete_terminal_codon_variant' => 18,
    'non_coding_transcript_exon_variant' => 19,
    'synonymous_variant' => 20,
    'mature_mirna_variant' => 21,
    'non_coding_transcript_variant' => 22,
    'regulatory_region_variant' => 23,
    'upstream_gene_variant' => 24,
    'regulatory_region_amplification' => 25,
    'tfbs_amplification' => 26,
    '5_prime_utr_variant' => 27,
    'intron_variant' => 28,
    '3_prime_utr_variant' => 29,
    'feature_truncation' => 30,
    'coding_transcript_variant' => 31,
    'tf_binding_site_variant' => 32,
    'start_retained_variant' => 33,
    'stop_retained_variant' => 34,
    'feature_elongation' => 35,
    'regulatory_region_ablation' => 36,
    'tfbs_ablation' => 37,
    'coding_sequence_variant' => 38,
    'downstream_gene_variant' => 39,
    'nmd_transcript_variant' => 40,
    'intergenic_variant' => 41,
    'sequence_variant' => 42
    );

my $vep_csq;

while( <VEP>) {
    ## Print and store Meta-info
    if( /^##/ ) {
        print;
        my( $type, $meta ) = parse_metainfo( $_ );
	    $vcf_meta{$type}->{$meta->{ID}} = $meta if defined $type;
        if ( /^##INFO=<ID=CSQ,Number=/) {$vep_csq = $_;}
    }
    # Print and store header
    elsif( /^#/ ) {
        print "##INFO\=<ID=GNOMADAF\,Number=1\,Type=Float,Description=\"Average AF GnomAD\">\n";
        print "##INFO=<ID=GNOMADAF_MAX,Number=1,Type=Float,Description=\"Highest reported AF in gnomAD\">\n";
        print "##INFO=<ID=GNOMADPOP_MAX,Number=1,Type=Float,Description=\"Population of highest AF\">\n";
        print "##INFO=<ID=dbNSFP_GERP___RS,Number=1,Type=Float,Description=\"GERP score\">\n";
        print "##INFO=<ID=dbNSFP_phyloP100way_vertebrate,Number=1,Type=Float,Description=\"phyloP100 score\">\n";
        print "##INFO=<ID=dbNSFP_phastCons100way_vertebrate,Number=1,Type=Float,Description=\"phastcons score\">\n";
        print "##INFO=<ID=CLNSIG_MOD,Number=.,Type=String,Description=\"Modified Variant Clinical Significance, for genmod score _0_ - Uncertain significance, _1_ - not provided, _2_ - Benign, _3_ - Likely benign, _4_ - Likely pathogenic, _5_ - Pathogenic, _6_ - drug response, _7_ - histocompatibility, _255_ - other\">\n";
        print "##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description=\"Most severe genomic consequence.\">\n";
        print "##INFO=<ID=CADD,Number=.,Type=String,Description=\"CADD phred score\">\n";
        print "##INFO=<ID=nhomalt,Number=.,Type=Integer,Description=\"number of alt allele homozygous individuals in gnomad\">\n";
	    print;
        $_ =~ s/^#//;
	    @head = split /\t/;
    }
    # Print and store variant information, add gnomadg and conservation scores
    # to info-field.
    else {
        my $doobi = parse_variant( $_, \@head, \%vcf_meta );
        my @add_info_field;
        #print Dumper($doobi);
        my @VARIANTS = split /\t/;
        my @info_field = split/;/,$VARIANTS[7];

        if ($doobi->{CHROM} =~ /^M/) {
            $vep_csq =~ /Consequence annotations from Ensembl VEP\. Format: (.*?)$/;
            my @field_names = split(/\|/, $1);
            my $info_field_mt = "";
            #print join("|",@field_names)."\n";
            my $trans_c = 0;
            #print Dumper($doobi->{INFO}->{CSQ}->[0]);
            foreach my $trans ( @{ $doobi->{INFO}->{CSQ} }) {
                my @csq_mt;
                foreach my $key ( @field_names) {
                   # print $key."  =>  ".$doobi->{INFO}->{CSQ}->[$trans_c]->{$key}."\n";
                    if ($maxentscan{$key}) {
                        push @csq_mt,"";
                    }
                    elsif ($key eq 'Consequence') {
                        push @csq_mt, join('&', @{$doobi->{INFO}->{CSQ}->[$trans_c]->{$key}});
                    }
                    else {
                        push @csq_mt,$doobi->{INFO}->{CSQ}->[$trans_c]->{$key};
                    }
                }
                my $csq_trans =  join("|",@csq_mt);
                $csq_trans =~ s/^\|//;
                #print $csq_trans."\n";
                if ($trans_c == 0) {
                    $info_field_mt = $info_field_mt.$csq_trans;
                }
                else {
                    $info_field_mt = $info_field_mt.",".$csq_trans;
                }
                $trans_c++; ##next transcript
            }
            my @tmpinfo;
            foreach my $info (@info_field) {
                if ($info =~ /CSQ/) {
                    push @tmpinfo, "CSQ=".$info_field_mt;
                }
                else {
                    push @tmpinfo,$info;
                }
            }
            @info_field = @tmpinfo;
            push @info_field,"GeneticModels=mt";
        }
        
        print join "\t", @VARIANTS[0..6];
        
        print "\t";
        ## GNOMAD 
        ### OVERALL
        my $gAF = $doobi->{INFO}->{CSQ}->[0]->{gnomADg_AF};
        if ($gAF) {
             push @add_info_field,"GNOMADAF=$gAF";
        }
	
        ### AF MAX POPULATION
        my $max = $doobi->{INFO}->{CSQ}->[0]->{gnomADg_AF_popmax};
        # my @max = split '&', $max;
        # $max = findmax(@max);
        # my $index = first_index {$_ eq $max } @max;
        if ($max) {
            push @add_info_field,"GNOMADAF_MAX=$max";
        }
        ### POPULATION WITH MAX
        my $max_pop = $doobi->{INFO}->{CSQ}->[0]->{gnomADg_popmax};
        if ($max_pop) {
            push @add_info_field,"GNOMADPOP_MAX=$max_pop";
        }
        ## GERP
        my $GERP = $doobi->{INFO}->{CSQ}->[0]->{GERP};
        if ($GERP) {
            push @add_info_field,"dbNSFP_GERP___RS=$GERP";
        }
        ## PHASTCONS
        my $pC = $doobi->{INFO}->{CSQ}->[0]->{phastCons};
        if ($pC) {
            push @add_info_field,"dbNSFP_phastCons100way_vertebrate=$pC";

        }
        ## PHYLOP
        my $pP = $doobi->{INFO}->{CSQ}->[0]->{phyloP100way};
        if ($pP) {
            push @add_info_field,"dbNSFP_phyloP100way_vertebrate=$pP";
        }
        ## CADD
        my $CADD = $doobi->{INFO}->{CSQ}->[0]->{CADD_PHRED};
        if ($CADD) {
            push @add_info_field,"CADD=$CADD";
        }
        # ## HomAltCount
        # my $hac = $doobi->{INFO}->{CSQ}->[0]->{gnomADg_nhomalt};
        # my @hac = split '&', $hac;
        # if ($hac) {
        #     if ($max && $max <= 0.02) {
        #         push @add_info_field,"nhomalt=$hac[$index]";
        #         #print STDERR "$max => $hac[$index]  $hac\n";
        #     }
            
        # }
        ## CLINSIG MODIFY
        my $csM = $doobi->{INFO}->{CLNSIG};
        my @mods;
        if (defined $csM) {
            my @CSm = split/,/,$csM;
            foreach my $entry (@CSm) {
                my @split_slash = (split '/', $entry);
                foreach my $entry_slash (@split_slash) {
                    if ($clinmod{$entry_slash}) {
                        push @mods,$clinmod{$entry_slash};
                    }
                    else {
                        push @mods,"_255_";
                    }
                }

            }
            push @add_info_field, "CLNSIG_MOD=".join('|',@mods);
        }
        ## MOST SEVERE CONSEQUENCE
        my $csq_ref = $doobi->{INFO}->{CSQ};
        my $m_s_c = CSQ($csq_ref);
        my $most_severe = ".";
        if (@$m_s_c) {
            $_ = lc for @$m_s_c;
		    $most_severe = (sort { $rank{$a} <=> $rank{$b} } @$m_s_c)[0];
        }
        push @add_info_field, "most_severe_consequence=".$most_severe;

        #Add new info field information
        push @info_field, @add_info_field;
        #print new and old information
        print join ";", @info_field;
        print "\t";
        #print everything after info field
        print join "\t", @VARIANTS[8..$#VARIANTS];
    }
}

sub findmax {
    my @in = @_;
    my $high = 0;
    foreach my $val (@in) {
        if ($val && $val ne '.') {
           if($val > $high) {
               $high = $val;
           }
       }
       else {

       }
    }
    return $high;
}

# Parse VCF meta info line (only FORMAT and INFO)
sub parse_metainfo {
    my $comment = shift;

    $comment =~ s/^##//;
    my( $type, $data ) = ( $comment =~ /^(.*?)=(.*)$/ );


    if( $type eq "FORMAT" or $type eq "INFO" or $type eq "SAMPLE" or $type eq "FILTER" ) {
	$data = remove_surrounding( $data, '<', '>' );
	my $pairs = keyval( $data, '=', ',' );
	return $type, $pairs;
    }

    return undef, undef;
}


# Parse VCF variant line
sub parse_variant {
    my( $var_str, $head, $meta ) = @_;
    my @var_data = split /\t/, $var_str;
    my %var;

    $var{ vcf_str } = $var_str;

    # First seven fields
    for ( 0..6 ) {
	$var{ $head->[$_] } = $var_data[$_];
    }

    # Eigth field, INFO
    $var{ INFO } = parse_info( $var_data[7] );

    # Parse VEP annotation field, if any
    if( $var{ INFO }->{ CSQ } ) {
	$var{ INFO }->{ CSQ } = parse_VEP_CSQ( $var{INFO}->{CSQ}, $meta->{INFO}->{CSQ} );
    }

    # Genotypes for each sample
    for ( 9 .. (@var_data-1) ) {
	$var{ GT } -> { $head->[$_] } = parse_genotype( $var_data[8], $var_data[$_] );
    }

    return \%var;
}


# Parse genotype field of VCF
sub parse_genotype {
    my( $format, $data ) = @_;

    my @format = split ':', $format;
    my @data   = split ':', $data;

    my %gt;
    @gt{@format} = @data;

    return \%gt;
}


# Parse info column of VCF file
sub parse_info {
    my $str = shift;
    my $info = keyval( $str, "=", ";" );

    return $info;
}


sub parse_VEP_CSQ {
    my( $CSQ_var, $CSQ_meta ) = @_;

    $CSQ_meta->{Description} =~ /Consequence annotations from Ensembl VEP\. Format: (.*?)$/;

    my @field_names = split /\|/, $1;

    my @transcripts = split /,/, $CSQ_var;

    my @data_transcripts;
    foreach my $transcript_CSQ ( @transcripts ) {
	my @values = split /\|/, $transcript_CSQ;

	my %data;
	for( 0 .. $#field_names ) {
	    if( $field_names[$_] eq "Consequence" ) {
		my @conseq_array = split '&', $values[$_];
		$data{ $field_names[$_] } = \@conseq_array;
	    }
	    else {
		$data{ $field_names[$_] } = ( $values[$_] or "" );
	    }
	}

	push( @data_transcripts, \%data )
    }
    return \@data_transcripts;
}

# Removes character(s) defined in arg2 if first in string, and arg3 if last in string.
sub remove_surrounding {
    my( $str, $before, $after ) = @_;
    $str =~ s/^$before//;
    $str =~ s/$after$//;
    return $str;
}


# Parse string with key value pairs. Return hash.
#  * Keys and values separated by 2nd argument.
#  * Pairs separated by 3rd argument
#  * Handles commas in values if surrounded by double quotes
sub keyval {
    my( $str, $keyval_sep, $pair_sep ) = @_;

    my @pair_str = split /$pair_sep/, $str;
    my %pairs;
    foreach( @pair_str ) {

	# If key-value separator exists, save the value for the key
	if( /$keyval_sep/ ) {
	    my( $key, $val ) = split /$keyval_sep/;
        $val = remove_surrounding( $val, '"', '"' );
	    $pairs{$key} = $val;
	}

	# Otherwise treat the whole string as a flag and set it to one (true).
	else {
	    $pairs{$_} = 1;
	}
    }
    return \%pairs;
}

sub excel_float {
    my $val = shift;

    return 0 if $val eq ".";

    $val =~ s/\./,/;
    return $val;
}


# THESE ARE ALL COUNTED AS "other":

#      4 Affects
#     14 association
#      2 _association
#     92 Conflicting_interpretations_of_pathogenicity
#      3 other
#      4 _other
#      6 protective
#      2 _protective
#     11 _risk_factor
#     37 risk_factor


sub CSQ {
	my ($csq) = shift;
    my @all_csq;
    foreach my $b (@$csq) {
        ## if canon pick consensus conseqeunce or if equal most severe ^
            my $tmp = $b->{Consequence};
            foreach my $conq (@$tmp) {
               push @all_csq, $conq;
            }
    
    }
    return \@all_csq; 

}
