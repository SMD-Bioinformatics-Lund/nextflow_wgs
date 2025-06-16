#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#use JSON::MaybeXS qw(encode_json decode_json);

# NOTE: This script is currently run as part of a container, not from the bin folder
# Keep here for reference, and eventually for extraction from the container

my $RSCRIPT = '/opt/plotting_POD.R'; # needs for a nice plot of the data
my $DISTANCE = 1000000; # distance allowed between connected segments

our %opt;
GetOptions( \%opt, 'ped=s', 'proband=s', 'snv=s', 'cnv=s', 'tmp=s', 'out=s', 'help' );

my %threshold = ('lower' => '0.54', 'upper' => '1.85'); 
my $AF_LIMIT = '0.05';

print_usage("ERROR: --snv, --ped and --proband ID required!\n") unless $opt{ped} and $opt{snv} and $opt{proband}; 

unless( -e $opt{snv} ) { print "ERROR: VCF file not found\n"; exit; }

unless( -e $opt{ped} ) { print "ERROR: PED file not found\n"; exit; }

if( $opt{cnv} eq 'NA' ){ $opt{cnv} = 0; }

unless( defined($opt{cnv}) ){ print "Running without cnv file!\n"; }

if( not defined($opt{tmp}) ){ $opt{tmp} = './tmp'; }else{ $opt{tmp} .= '/tmp'; }
unless( -e $opt{tmp} ) { system( "mkdir $opt{tmp}" ); }

unless( defined($opt{out}) ){ $opt{out} = './'; }
 
#print Dumper( \%opt ); exit;

# get PED info
my ( $PED_info, $found_p ) = fetch_PEDinfo( $opt{ped}, $opt{proband} );

unless( $found_p ){ print "ERROR: Proband not found in PED\n"; exit; }

if( $opt{cnv} ){

    print "Running with cnv file!\n";

    if( $opt{cnv} =~ /gz$/ ){
        open( BED, "zcat $opt{cnv} |" ) or die "gunzip $opt{cnv}: $!";
    }else{
        open( BED, "$opt{cnv}" ) or die "gunzip $opt{cnv}: $!";
    }

    open OUT, ">$opt{tmp}/dup.bed" or die "cannot write to file dup.bed!";
    while ( my $row = <BED> ){
        chomp( $row );
        next if( $row=~/^#/ );
        my ( $svtype, $end );
        my $print = 0;
        if( $row =~/SVTYPE=(.?DUP);/ ){ $svtype = $1; $print = 1; }
        if( $row =~/END=(\d+);/ ){ $end = $1; }
        my @row = split /\t/, $row;
        if( $print ){
            print OUT "$row[0]\t$row[1]\t$end\t$svtype".'_'."$row[0]".':'."$row[1]".'_'."$end\n";
        }
    }
    close BED;
    close OUT;
    system( "bedtools intersect -header -a $opt{snv} -b $opt{tmp}/dup.bed -wa > $opt{tmp}/intersected_vcf.vcf" );
    $opt{snv} = "$opt{tmp}/intersected_vcf.vcf";
}

if( $opt{snv} =~ /gz$/ ){
    open( IN, "zcat $opt{snv} |" ) or die "gunzip $opt{snv}: $!";
}else{
    open( IN, "$opt{snv}" ) or die "gunzip $opt{snv}: $!";
}


my ( %labels, %format, %ped, %summary, %data, %allele_ratio, @csq_labels, $counter );

while ( my $row = <IN> ){
    chomp( $row );
    next if( $row =~/^##/ );  #skip initial lines
    if( $row =~/^#CHROM/ ){
        my @row = split /\t/, $row;
        for my $i( 0..$#row ){ $labels{ $row[$i] } = $i; }
        for my $fam_id ( sort keys %$PED_info ){
            unless( $labels{ $fam_id } ){ print "ERROR: VCF file does not match $fam_id from PED file\n"; exit; }
        }

        if( defined $$PED_info{$opt{'proband'}}{'maternal_id'} ){ $ped{ 'Mother' } = $$PED_info{$opt{'proband'}}{'maternal_id'}; }
        if( defined $$PED_info{$opt{'proband'}}{'paternal_id'} ){ $ped{ 'Father' } = $$PED_info{$opt{'proband'}}{'paternal_id'}; }

    }else{
       $counter ++;
       my $pop_af = 0;
       my @row = split /\t/, $row;
       my @info = split /;/, $row[ $labels{'INFO'} ];
       for my $i(0..$#info){
           if( $info[$i] =~/^GNOMADAF=(\S+)$/ ){
               $pop_af = $1;
           }
       }
       #print "$row[$labels{'#CHROM'}]\t$row[$labels{'POS'}]\t$row[$labels{'REF'}]\t$row[$labels{'ALT'}]\t$pop_af\n";
        
       next if ( $pop_af <= $AF_LIMIT ); # if not common snp, move on
       
       my @format = split /:/, $row[ $labels{'FORMAT'} ];
       for my $i(0..$#format){ $format{ $format[$i] } = $i; }
       my @pb = split /:/, $row[ $labels{$opt{'proband'}} ];
       my ( $AD_0, $AD_1 ) = split /,/, $pb[ $format{'AD'}];
 
       if( ( $pb[ $format{'GT'} ] eq '0/1' ) && ((( $AD_0/$AD_1 ) >= $threshold{'upper'} ) || (( $AD_0/$AD_1 ) <= $threshold{'lower'} )) ){

#           print "$row[$labels{'#CHROM'}]\t$row[$labels{'POS'}]\t$row[$labels{'REF'}]\t$row[$labels{'ALT'}]\t$pb[ $format{'GT'} ]\t$pb[ $format{'AD'} ]\t".($AD_0/$AD_1)."\n";
           
           if(( $AD_0/$AD_1 ) >= $threshold{'upper'} ){ # barnet är av genotypen 0/0/1 (AAB)
               
               $summary{ $row[$labels{'#CHROM'}] }{'total'}++;
               
               my ( @m, @f );
               if( $ped{'Mother'} ){ @m = split /:/, $row[ $labels{ $ped{ 'Mother' } } ]; }
               if( $ped{'Father'} ){ @f = split /:/, $row[ $labels{ $ped{ 'Father' } } ]; }
               if( $m[ $format{'GT'} ] eq '1/1' ){ 
                   $summary{ $row[$labels{'#CHROM'}] }{'F-AAB'}++; 
                   $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = 'F-AAB';
                   $allele_ratio{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = ($AD_0/$AD_1);
               }    
               if( $f[ $format{'GT'} ] eq '1/1' ) { 
                   if( defined $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} and $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} eq 'F-AAB' ){ 
                       $data{$row[$labels{'#CHROM'}]} {$row[$labels{'POS'}] } = '0';
                       $summary{$row[$labels{'#CHROM'}]}{'F-AAB'}--; 
                   }else{ 
                       $summary{$row[$labels{'#CHROM'}]}{'M-AAB'}++; 
                       $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = 'M-AAB';
                       $allele_ratio{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = ($AD_0/$AD_1);
                   }
               }

                                 
           }else{ # barnet är av genotypen 0/1/1 (ABB)

               $summary{ $row[$labels{'#CHROM'}] }{'total'}++;

               my ( @m, @f );
               if( $ped{'Mother'} ){ @m = split /:/, $row[ $labels{ $ped{ 'Mother' } } ]; }
               if( $ped{'Father'} ){ @f = split /:/, $row[ $labels{ $ped{ 'Father' } } ]; }
               if( $m[ $format{'GT'} ] eq '0/0' ){ 
                   $summary{$row[$labels{'#CHROM'}]}{'F-ABB'}++;
                   $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = 'F-ABB';
                   $allele_ratio{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = ($AD_0/$AD_1);
               }
               if( $f[ $format{'GT'} ] eq '0/0' ){ 
                   if( defined $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} and $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} eq 'F-ABB' ){
                       $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = '0';
                       $summary{$row[$labels{'#CHROM'}]}{'F-ABB'}--;
                   }else{
                       $summary{$row[$labels{'#CHROM'}]}{'M-ABB'}++;
                       $data{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = 'M-ABB';
                       $allele_ratio{$row[$labels{'#CHROM'}]}{$row[$labels{'POS'}]} = ($AD_0/$AD_1);
                   }
               }
           }
       }
   }
}
close IN;

open OUT, ">$opt{tmp}/all_chromosomes.txt" or die "cannot print to file!";
open RES, ">$opt{out}".'/'.$opt{proband}."_POD_results.html" or die "cannot print to file!";
print RES "<table style=\"font-size:0.9em;text-align:center;width:600px;\"><tr style=\"background-color:lightgrey;\">";
for my $pod_value( 'Chrom', 'Skewed', 'POD-T', 'M-POD', 'P-POD', 'MP-score', 'M-AAB', 'P-AAB', 'M-ABB', 'P-ABB' ){ print RES "<th>$pod_value</th>"; }
print RES "</tr>\n"; 

my $seg_nr = 0;
my ( %segments, %seg_size, %origin );
for my $chr( 1..22,'X','Y' ){
    my $prev_parent = 'INIT';
    for my $pos ( sort {$a<=>$b} keys %{$data{$chr}} ){
        if( $data{$chr}{$pos} ne '0'){ 
            print OUT "$chr\t$pos\t$data{$chr}{$pos}\t$allele_ratio{$chr}{$pos}\n";
            my $parent = $data{ $chr }{ $pos };
            $parent =~s/-\S+//;
            if( $parent eq $prev_parent ){ 
                push( @{$segments{ $seg_nr }}, "$chr:$pos" );
                $seg_size{ $seg_nr } ++;
                $origin{ $seg_nr } = $parent;
            }else{
                $seg_nr ++;
                push( @{$segments{ $seg_nr }}, "$chr:$pos" );
                $seg_size{ $seg_nr } ++;
                $origin{ $seg_nr } = $parent;
            }
            $prev_parent = $parent;
        }
    }
    print RES "<tr style=\"background-color:lightgrey;\">";
    if( $summary{ $chr }{ total } ){ 
        print RES "<td>$chr</td><td>$summary{ $chr }{ total }</td>";
    }else{
        print RES "<td>$chr</td><td>0</td>";
    }
    my ( $maab, $faab, $mabb, $fabb ) = ( '0', '0', '0', '0' );
    if( $summary{$chr}{'M-AAB'} ){ $maab = $summary{$chr}{'M-AAB'}; }
    if( $summary{$chr}{'F-AAB'} ){ $faab = $summary{$chr}{'F-AAB'}; }
    if( $summary{$chr}{'M-ABB'} ){ $mabb = $summary{$chr}{'M-ABB'}; }
    if( $summary{$chr}{'F-ABB'} ){ $fabb = $summary{$chr}{'F-ABB'}; }

    print RES "<td>".( $maab + $faab + $mabb + $fabb )."</td>";
 
    print RES "<td>".( $maab + $mabb )."</td>";

    print RES "<td>".( $faab + $fabb )."</td>";
    
    if( ($maab + $mabb) > 0 and ($faab + $fabb) > 0 ){ 
        print RES "<td>".sprintf( "%.0f", 100 * abs((( $maab + $mabb ) - ( $faab + $fabb )))/(( $maab + $mabb ) + ( $faab + $fabb )))."</td>";
    }else{
        print RES "<td>0</td>";
    }
    for my $pod_value( 'M-AAB', 'F-AAB', 'M-ABB', 'F-ABB' ){
        if( $summary{ $chr }{ $pod_value } ){
            print RES "<td>$summary{ $chr }{ $pod_value }</td>";
        }else{
            print RES "<td>0</td>";
        }
    } 
    print RES "</tr>\n";
}
print RES "</table>\n";
close OUT;

# connect segments nearby, $DISTANCE decides inclusion
#######################################################################
my %connected_segments;
my %consumed_segments;
for my $nr( 1..($seg_nr - 1) ){
    for my $nr2( ($nr + 1)..$seg_nr ){
        my ($chr, $pos ) = split( /:/, $segments{$nr}[-1] );
        my ($chr2, $pos2 ) = split( /:/, $segments{$nr2}[0] );
        next if( $chr ne $chr2 );
        next if( $origin{$nr} ne $origin{$nr2} );
        if( $pos2 - $pos <= $DISTANCE ){
            #print "$nr\t$origin{$nr}\t$seg_size{$nr}\t$segments{$nr}[-1] <- $nr2\t$origin{$nr2}\t$seg_size{$nr2}\t$segments{$nr2}[0]\n";
            $connected_segments{$nr}{$nr2} = 1;
            $consumed_segments{$nr} = 1;
            for my $n( sort {$a<=>$b} keys %consumed_segments ){ if($connected_segments{$n}{$nr}){ $connected_segments{$n}{$nr2} = 1; $connected_segments{$nr}{$nr2} = 0; } }
        }
    }
}

my ( %connected_sizes, %connected_coords );
#for my $n(sort {$a<=>$b} keys %connected_segments){
for my $n( 1..$seg_nr ){
    $connected_sizes{$n} = $seg_size{$n};
    $connected_coords{$n}{'start'} = $segments{$n}[0];
    for my $n2(sort {$a<=>$b} keys %{$connected_segments{$n}}){
        if( $connected_segments{$n}{$n2} ){
            #print "$n -> $n2\t$origin{$n}\t$connected_segments{$n}{$n2}\n";
            $connected_sizes{$n} += $seg_size{$n2};
            $connected_coords{$n}{'end'} = $segments{$n2}[-1];
        }
    }
}

my %once;
#for my $n( 1..$seg_nr ){
print RES "<br><table style=\"font-size:0.9em;text-align:center;width:600px;\"><tr style=\"background-color:lightgrey;\">";
print RES "<th># Rank</th><th># Markers</th><th>Origin</th><th>Start</th><th>End</th></tr>";
$counter = 0;
for my $n( sort { $connected_sizes{$b}<=>$connected_sizes{$a} } keys %connected_sizes ){
    my $p_seg = "";
    for my $n2( 1..$seg_nr ){
        if( $connected_segments{$n}{$n2} ){ $p_seg = "$connected_sizes{$n}</td><td>$origin{$n}</td><td>$connected_coords{$n}{'start'}</td><td>$connected_coords{$n}{'end'}"; $once{$n} = 1; $once{$n2} = 1; }
    }
    if( $p_seg ne "" ){
        $counter ++; 
        print RES "<tr style=\"background-color:lightgrey;\">";
        print RES "<td>$counter</td><td>$p_seg</td></tr>\n";
    }
    if( !$once{$n} ){
        $counter ++; 
        print RES "<tr style=\"background-color:lightgrey;\">";
        print RES "<td>$counter</td><td>$seg_size{$n}</td><td>$origin{$n}</td><td>$segments{$n}[0]</td><td>$segments{$n}[-1]</td></tr>\n";
    }
    if( $counter >= 25 ){ last; }
}
print RES "</table>\n";
close RES;
#######################################################################

#print RES "<br><table style=\"font-size:0.9em;text-align:center;width:600px;\"><tr style=\"background-color:lightgrey;\">";
#print RES "<th># Rank</th><th># Markers</th><th>Origin</th><th>Start</th><th>End</th></tr>";
#$counter = 0;
#for my $seg_nr ( sort { $seg_size{$b}<=>$seg_size{$a} } keys %seg_size ){
#    next if( $seg_size{$seg_nr} <= 10 );
#    $counter ++;
#    $origin{$seg_nr} =~s/F/P/; 
#    print RES "<tr style=\"background-color:lightgrey;\">"; 
#    print RES "<td>$counter</td><td>$seg_size{$seg_nr}</td><td>$origin{$seg_nr}</td><td>$segments{$seg_nr}[0]</td><td>$segments{$seg_nr}[-1]</td></tr>\n";
#    if( $counter >= 25 ){ last; }
#}
#print RES "</table>\n";
#close RES;

system("Rscript --vanilla $RSCRIPT $opt{tmp}/all_chromosomes.txt $opt{out} $opt{proband} >/dev/null 2>&1");
system( "rm -r $opt{tmp}" );

sub fetch_PEDinfo {
    # sex: 1 = male, 2 = female
    # phenotype: 1 = unaffected, 2 = affected
    # relations: proband, mother, father and other
    my ( $ped_path, $pid ) = @_;
    my %sex = ( 'other' => 'U', '1' => 'M', '2' => 'F' );
    my %phenotype = ( '-9' => 'missing', '0' => 'missing', '1' => 'unaffected', '2' => 'affected' );
    my %type = ( 'P' => 'proband', 'M' => 'mother', 'F' => 'father', 'O' => 'other' );
    my %relations;
    my $proband_match = 0;
    open PED, "$ped_path" or die "can't open for reading";
    while( my $row = <PED> ){
        chomp($row);

        my @col = split /\t/, $row;
        $relations { $col[1] } -> {'fam_id'} = $col[0];
        $relations { $col[1] } -> {'paternal_id'} = $col[2];
        $relations { $col[1] } -> {'maternal_id'} = $col[3];
        $relations { $col[1] } -> {'sex'} = $sex{ $col[4] };
        $relations { $col[1] } -> {'phenotype'} = $phenotype{ $col[5] };
        if( $col[6] ) {
            $relations { $col[1] } -> {'type'} = $type{ $col[6] };
            #$relations { $col[1] } -> {'run'} = $col[7];
        } else {
        # determine type if its not given by PED
            if(( $col[2] eq '0' ) && ( $col[3] eq '0' )){
                if( $col[4] eq '1'){
                    $relations { $col[1] } -> {'type'} = 'father'; # funkar inte om där finns syskon till föräldrarna!
                }else{
                    $relations { $col[1] } -> {'type'} = 'mother'; # funkar inte om där finns syskon till föräldrarna!
                }
            }else{
                if( $col[1] eq $pid ){ # Matcha mot manuellt inslaget id
                    $relations { $col[1] } -> {'type'} = 'proband';
                    $proband_match = 1;
                }else{
                    $relations { $col[1] } -> {'type'} = 'other';
                }
            }
        }
        if( $col[7] ) {
            $relations { $col[1] } -> {'run'} = $col[7];
        }
    }
    close PED;
    return ( \%relations, $proband_match );
} # fetch_PEDinfo

sub print_usage {
    print "$_[0]\n\n" if $_[0];
    print "USAGE: create_pipeline_json_from_pedfile.pl --ped <PED FILE>\n\n";
    print "    --proband SAMPLE ID  ID referring to the proband in the PED file\n\n";
    print "    --ped     FILE       Path to file describing the family relationship\n\n";
    print "    --snv     FILE       Path to SNV VCF file \(required\)\n\n";
    print "    --cnv     FILE       Path to CNV VCF file \(optional\)\n\n";
    print "    --tmp     DIR        Where to put temporary files \(default is ./tmp\)\n\n";
    print "    --out     DIR        Path to output dir \(default is ./\)\n\n";
    print "    --help               Will print this message\n\n";
    print "PED file contains 6 columns:\n";
    print "  * Family ID\n";
    print "  * Individual ID\n";
    print "  * Paternal ID\n";
    print "  * Maternal ID\n";
    print "  * Sex\n";
    print "  * Phenotype\n";
    print "Example:\n";
    print "666-01\t666-01\t666-02\t666-03\t1\t2\n";
    print "666-01\t666-02\t0\t0\t1\t1\n";
    print "666-01\t666-03\t0\t0\t2\t1\n\n";
    exit(0);
} # print_usage
