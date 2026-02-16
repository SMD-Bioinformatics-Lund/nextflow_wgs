#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname(__FILE__);
use vcf2;
use Data::Dumper;
use List::Util qw( max min );
use Getopt::Long;
use GD::Graph;
use GD::Graph::lines;

## TODO:
# 1: What happens when tumor variants = 0? Plot empty?
# 2: Make a suggested genotype for contaminant
# 3: Clean up code, and make more subroutines

# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcf=s', 'case-id=s', 'normal', 'detect-level=s',
    'ADfield-name=s', 'high=s', 'binsize-cutoff=s' );
my $vcf            = vcf2->new( 'file' => $opt{vcf} );
my $id             = $opt{"case-id"};
my $check_normal   = $opt{normal};
my $high           = 0.3;
my $low            = 0;
my $ADfield        = "VD";
my $binsize_cutoff = 80;

if ( $opt{'detect-level'} ) {
    $low = $opt{'detect-level'};
}
if ( $opt{'ADfield-name'} ) {
    $ADfield = $opt{'ADfield-name'};
}
if ( $opt{'high'} ) {
    $high = $opt{'high'};
}
if ( $opt{'binsize-cutoff'} ) {
    $binsize_cutoff = $opt{'binsize-cutoff'};
}

my %dist_int;
my %dist_dec;
my $paired = 0;

## save all SNV variants present in gnomadgenomes > 5% and their allelefrequencies
while ( my $v = $vcf->next_var() ) {

    my $t_VAF;
    my $t_VD;
    my $t_DP;
    my $n_VAF;
    my $n_VD;
    my $n_DP;

    next if $v->{CHROM} =~ /X|Y/;
    next unless exists $v->{INFO}->{CSQ}->[0]->{gnomADg_AF};
    my $gnomad = $v->{INFO}->{CSQ}->[0]->{gnomADg_AF};
    next unless ( $gnomad =~ /\d+/ ) {
        $gnomad = $v->{INFO}->{CSQ}->[0]->{gnomADg_AF};
        $gnomad = find_max($gnomad);
    }
    next unless ( $gnomad >= 0.05 );

    ## only snvs
    next if length( $v->{REF} ) > 1 || length( $v->{ALT} ) > 1;
    my $var =
      $v->{CHROM} . "\t" . $v->{POS} . "\t" . $v->{REF} . "\t" . $v->{ALT};

    ## if paired returns normal sample-id
    $paired = paired( $id, $v );
    if ( $opt{normal} && $paired eq "0" ) {
        print "no normal sample in vcf\n";
    }

    for my $gt ( @{ $v->{GT} } ) {
        if ( $gt->{_sample_id} eq $id ) {
            $t_VD = $gt->{$ADfield};
            my @ads = split( /,/, $t_VD );
            my $min = min(@ads);
            $t_DP = $gt->{DP};
            if ( $ADfield eq 'VD' ) {
                $t_VAF = $gt->{VAF};
            }
            else {
                if ( $t_DP == 0 ) {
                    $t_VAF = 0;
                }
                else {
                    $t_VAF = $min / $t_DP;
                }
            }
        }
        else {
            $n_VD = $gt->{$ADfield};
            my @ads = split( /,/, $t_VD );
            my $min = min(@ads);
            $n_DP = $gt->{DP};
            if ( $ADfield eq 'VD' ) {
                $n_VAF = $gt->{VAF};
            }
            else {
                if ( $n_DP == 0 ) {
                    $n_VAF = 0;
                }
                else {
                    $n_VAF = $min / $n_DP;
                }
            }
        }
    }
    next if $t_DP <= 5;
    if ($paired) {
        if ($check_normal) {
            check_vaf( $n_VAF, 0, $var );
        }
        else {
            check_vaf( $n_VAF, $t_VAF, $var );
        }
    }
    else {
        check_vaf( $t_VAF, 0, $var );
    }

}

## if normal was asked for set id to other sample ##
if ( $opt{normal} && $paired ) {
    $id = $paired;
}

## get number of variants within each window set by high low and windowsize (0.005)
my ( $distri, $num_bins, $num_vars_bin ) = get_distibution( $low, $high );
my %distri = %$distri;

#print Dumper(%distri);
## if no bins, no variants, no contamination
if ( $num_bins == 0 ) {
    print "0.0\n";
    system("touch $id.png");
    exit;
}

## find heterozygous highpoint and plot distribution
my ( $af_at_highpoint, $bin_at_highpoint ) =
  find_heterozygous_peak( $num_bins, $num_vars_bin );

## if a highpoint is found, i.e. above nominal cutoff set by assay-hash ##
## this is likely the contamination score for the sample, try to find homozygous peak (should be a peak around 2x the AF-average of hetero peak)
## print the genotypes of the heterozygous loci and, if found, the homozygous loci too
if ($af_at_highpoint) {
    print $af_at_highpoint. "\n";
    ## try to find homozygous peak, starting from hetero highpoint
    my ( $bin_at_homo_highpoint, $homo_highpoint, $af_at_homo ) =
      find_homozygous_peak( $bin_at_highpoint, $af_at_highpoint );
    print_genotypes( $distri{$bin_at_highpoint}{VARS}, "hetero" );
    if ( $bin_at_homo_highpoint ne "cannot find homo" ) {
        print_genotypes( $distri{$bin_at_homo_highpoint}{VARS}, "homo" );
    }

}
else {
    print "0.0\n";
}

sub meanAF {
    my $afs = shift;
    my @afs = @$afs;

    my $sum = 0;
    foreach my $value (@afs) {
        $sum = $sum + $value;
    }
    my $mean = $sum / scalar(@afs);

    return $mean;
}

sub paired {
    my ( $id, $v ) = @_;

    my $n_id = 0;
    if ( scalar( @{ $v->{GT} } ) > 1 ) {
        for my $gt ( @{ $v->{GT} } ) {
            if ( $gt->{_sample_id} eq $id ) {

            }
            else {
                $n_id = $gt->{_sample_id};
            }

        }
    }
    return $n_id;
}

sub check_vaf {
    # Groups VAFs into bins, both on INTEGER and FLOAT level
    # Go from low to high, default 0.01 -> 0.30
    # If paired tumor sample, make sure it is not a tumor specific variant
    my $vaf      = shift;
    my $othervaf = shift;
    my $var      = shift;
    ## only one VAF, not paired
    if ( !$othervaf ) {
        if ( $vaf >= $low && $vaf <= $high ) {
            my $matchingbin = $vaf * 100;
            my @dec = split( '\.', $matchingbin );
            if ( scalar(@dec) > 1 ) {
                $matchingbin = $dec[0]++ if $dec[1] >= 0.5;
                $matchingbin = $dec[0]   if $dec[1] < 0.5;
            }
            else {
                $matchingbin = $dec[0];
            }

            $dist_int{$matchingbin}++;
            $dist_dec{$vaf}{COUNT}++;
            push @{ $dist_dec{$vaf}{VAR} }, $var;
        }
    }
    ## extra check that the variant is indeed germline, this is mostly fixed by gnomad_check earlier.
    else {
        if ( $vaf <= $low && $othervaf >= $low && $othervaf <= $high ) {

            my $matchingbin = $othervaf * 100;
            my @dec = split( '\.', $matchingbin );
            if ( scalar(@dec) > 1 ) {
                $matchingbin = $dec[0]++ if $dec[1] >= 0.5;
                $matchingbin = $dec[0]   if $dec[1] < 0.5;
            }
            else {
                $matchingbin = $dec[0];
            }

            $dist_int{$matchingbin}++;
            $dist_dec{$othervaf}{COUNT}++;
            push @{ $dist_dec{$othervaf}{VAR} }, $var;
        }
    }

}

sub find_homozygous_peak {
    my ( $first_peak, $af_at_hetero ) = @_;
    my $homo_highpoint        = 0;
    my $af_at_homo            = 0;
    my $bin_at_homo_highpoint = 0;
    for ( my $i = $first_peak + 1 ; $i <= $num_bins - 1 ; $i++ ) {
        if (   $distri{$i}{COUNT} < $distri{ $i + 1 }{COUNT}
            && $distri{$i}{COUNT} > $homo_highpoint )
        {
            $homo_highpoint        = $distri{ $i + 1 }->{COUNT};
            $af_at_homo            = $distri{ $i + 1 }->{MEAN};
            $bin_at_homo_highpoint = $i + 1;
        }
    }
    if (   $af_at_homo / $af_at_hetero <= 2.2
        && $af_at_homo / $af_at_hetero >= 1.8 )
    {
        return $bin_at_homo_highpoint, $homo_highpoint, $af_at_homo;
    }
    else {
        return "cannot find homo";
    }
}

sub print_genotypes {
    my ( $vars, $type ) = @_;
    my @vars      = @$vars;
    my $dist_file = "$id.genotypes.txt";
    open( GENOTYPES, '>>', $dist_file );
    foreach my $var (@vars) {
        my @tmp = split( "\t", $var );

        print GENOTYPES $tmp[0] . "\t" . $tmp[1] . "\t";
        if ( $type eq "hetero" ) {
            print GENOTYPES $tmp[2] . "/" . $tmp[3] . "\n";
        }
        elsif ( $type eq "homo" ) {
            print GENOTYPES $tmp[3] . "/" . $tmp[3] . "\n";
        }
    }
    close GENOTYPES;
}

sub get_distibution {
    my ( $low, $high ) = @_;
    my %distri;
    my %dist;
    my $inc          = 0.005;
    my $start        = 0;
    my $window       = 0.005;
    my $num_vars_bin = 0;
    my $num_bins     = 0;
    for ( my $i = $low ; $i <= $high ; $i = $i + $inc ) {
        my $sum = 0;
        my @afs;
        my @vars;
        foreach my $af ( sort { $a <=> $b } keys %dist_dec ) {
            if ( $af >= $i && $af <= $i + $window ) {
                $sum = $dist_dec{$af}{COUNT} + $sum;
                push @afs,  $af;
                push @vars, @{ $dist_dec{$af}{VAR} };
            }
        }
        next if $sum == 0;
        $num_bins++;
        my $j = $i + $window;
        $distri{$num_bins}{AFspan} = "$i-$j";
        $distri{$num_bins}->{COUNT} = $sum;
        my $mean = meanAF( \@afs );
        $distri{$num_bins}->{MEAN} = $mean;
        $distri{$num_bins}->{VARS} = \@vars;
        $num_vars_bin              = $sum + $num_vars_bin;

    }
    return \%distri, $num_bins, $num_vars_bin;
}

sub find_heterozygous_peak {
    my ( $num_bins, $num_vars_bin ) = @_;
    my $dist_file = "$id.dist.txt";
    open( DIST, '>', $dist_file );
    my $mean_bincount    = $num_vars_bin / $num_bins;
    my $highpoint        = 0;
    my $af_at_highpoint  = 0;
    my $bin_at_highpoint = 0;
    my @xdata;
    my @ydata;

    foreach my $af_c ( sort { $a <=> $b } keys %distri ) {
        push @ydata, $distri{$af_c}->{COUNT};
        my $xdata = sprintf( "%.3f", $distri{$af_c}->{MEAN} );
        push @xdata, $xdata;
        if ( $distri{$af_c}->{COUNT} / $mean_bincount >= 3.5 ) {
            if (   $distri{$af_c}->{COUNT} > $highpoint
                && $distri{$af_c}->{COUNT} > $binsize_cutoff )
            {
                $highpoint        = $distri{$af_c}->{COUNT};
                $af_at_highpoint  = $distri{$af_c}->{MEAN};
                $bin_at_highpoint = $af_c;
            }
        }
        print DIST $af_c . "\t"
          . $distri{$af_c}->{COUNT} . "\t"
          . $distri{$af_c}->{MEAN} . "\n";
    }
    close DIST;
    my @data;
    push @data, \@xdata, \@ydata;
    plot_distribution( \@data );
    return $af_at_highpoint, $bin_at_highpoint;
}

sub plot_distribution {
    my $data = shift;
    my @data = @$data;
    ## PLOT DISTRIBUTION ##
    #######################

    my $graph = GD::Graph::lines->new( 1500, 1000 );
    $graph->set_text_clr("#BFBF00");
    my $gd = $graph->plot( \@data ) or die $graph->error;
    open( IMG, ">$id.png" ) or die $!;
    binmode IMG;
    print IMG $gd->png;
    close IMG;
    #######################
}

sub find_max {
    my $af = shift;
    my @af = split( '&', $af );

    my $max = 0;
    for my $a (@af) {
        if ( $a > $max ) {
            $max = $a;
        }
    }
    return $max;
}
