#!/usr/bin/perl

use strict; 
use warnings; 
use IO::Zlib;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

## v 3.0 - oct 2015
## requires GPL.list with 3col files 

my $gzfile = shift @ARGV; 
my $gse=$1 if ($gzfile =~ m/(GSE\d+)/);

############################################
my $nsample_min = 8; 
my $nsample_max = 200;

my ( @lines,@values,$min,$max,$ave,$stdev,$logmin,$logmax,$logave,$logstdev,$nsample,$nprobes );
my $gpl=undef;              ## platform ID learnt from the file contents (not the file name!)
my $title="Probe_id,Symbol,Entrez_ID";
my $if_log2=undef;
my $gplok=0;                   ## whether file is ok to be processed 
my ( @idref,@pid,@tit,@geo,@gpls ); 

my %symbol;              ## maps probe id -> gene symbol 
my %entrez;              ## and probe id -> entrez id  
############################################

print STDOUT "Processing file $gzfile, dataset ID $gse\n"; 
open FILE,"gzip -dc $gzfile |" or die "ERROR: Failed to open GZIP archive file $gzfile.";
open GPL,"<","GPL.list" or die "ERROR: You must provide file GPL.list containing all of the acceptable platforms."; 

############################################
#       PARSING THE EXPRESSION MATRIX      #
############################################

while (<FILE>) { 
  if ( ! m/^!/ && ! m/^\"ID_REF\"/ && ! m/^\s+$/) {
    ## these are all of the numerical expression values, with probe IDs 
    s/\"//g; 
    chomp;
    push @lines,$_; 
    my @line = split(/\t+/);     
    shift @line; 
    #print "@line\n";
    push @values,@line; 
  } elsif (m/^\"ID_REF\"/) {
    ## ID_REF string with all of the GSM values #1 
    s/\"//g; 
    s/^ID_REF\t//g; 
    chomp;
    @idref=split(/\t/); 
  } elsif (m/^\!Sample_geo_accession/) {
    ## ID_REF string with all of the GSM values #2 - will compare with #1 for a sanity check
    s/\"//g; 
    s/^\!Sample_geo_accession\t//g; 
    chomp;
    @geo=split(/\t/); 
  } elsif (m/^\!Sample_platform_id/) {
    ## all of the platforms expressed as GPL\d*  
    s/\"//g; 
    s/^\!Sample_platform_id\t//g; 
    chomp;
    @pid=split(/\t/); 
    $gpl=$pid[0];
    foreach my $val (@pid) { 
      if ($val ne $gpl) { 
        print STDERR "ERROR: The expression matrix contains multiple platform IDs!\n";
        exit 1; 
      } 
    } 
  } elsif (m/^\!Sample_title/) {
    ## sample titles 
    s/^\!Sample_title\t//g; 
    s/(\"|\n)//g;
    s/\t+/XYXYXYX/g;
    s/\W+/_/g;
    s/_+/_/g;
    s/Na_ve/Naive/g; 
    s/XYXYXYX/\t/g;
    @tit=split(/\t/); 
  }
}

$nsample = $#pid+1;
$nprobes = $#lines+1;
printf STDOUT "Expression matrix parsed: platform ID %s, number of samples: %d, number of probes: %d\n",$gpl,$nsample,$nprobes; 

## Make appropriate title string, including both GEO ids and sample titles. 
for (my $i=0; $i<=$#tit; $i++) { 
  my $val=",".$tit[$i]."_".$geo[$i];
  $title=$title.$val; 
}
$title=$title."\n";

############################################
#         SOME BASIC SANITY CHECKS         #
############################################

while (<GPL>) {
  if (m/$gpl/) { 
    printf STDOUT "The platform used ($gpl) is on the list of accepted platforms! Proceeding...\n";
    $gplok=1;
    last;  
  } 
}

if (!$gplok) {  
  print STDERR "ERROR: The platform used ($gpl) is not in the accepted list.\n";
  exit 1; 
}

if ($nsample < $nsample_min) { 
  print STDERR "ERROR: The number of samples ($nsample) is less than the pre-defined minimum ($nsample_min)\n";
  exit 1; 
} elsif ($nsample > $nsample_max) {
  print STDERR "ERROR: The number of samples ($nsample) is more than the pre-defined maximum ($nsample_max)\n";
  exit 1; 
} 

############################################
#           MAKE ANNOTATION HASH           #
############################################

open ANN,"<","$gpl.3col" or die "$!"; 

while (<ANN>) {
  chomp;
  my @tt=split(/\t+/);
  
  $symbol{$tt[0]}=$tt[1];
  $entrez{$tt[0]}=$tt[2];
  
  if ($tt[1] =~ m/\/\/\// || $tt[2] =~ m/\/\/\//) { 
    my @sa=split(/\s*\/\/\/\s*/,$tt[1]);
    my @ea=split(/\s*\/\/\/\s*/,$tt[2]);
    my $mini=0;                        ## we associate the probes that have multiple genes in them with the one with the smallest Entrez ID. 
    for (my $i=0; $i<=$#ea; $i++) {
      $mini = $i if ($ea[$i] <= $ea[$mini]);
    }
    $symbol{$tt[0]}=$sa[$mini];
    $entrez{$tt[0]}=$ea[$mini];
  }
}


############################################
#     CLASSIFY THE EXPRESSION MATRIX       #
############################################

$min=$values[0];
$max=$values[0];
$ave=0;
$stdev=0; 

if ($#values != -1) { 
  for (my $i=0; $i<=$#values; $i++) {
    if ($values[$i]) {
      $max=$values[$i] if ($values[$i]>$max);
      $min=$values[$i] if ($values[$i]<$min);
      $ave+=$values[$i];
      $stdev+=$values[$i]**2;
    }
  }
  $stdev=sqrt(($stdev-($ave**2/($#values+1)))/($#values+1));           ## these are not column-wise, but overall stats 
  $ave/=($#values+1); 
} else { 
  print STDERR "ERROR: The expression matrix is empty!\n"; 
  exit 1; 
}

if ($max > 10 && $max < 25 && $ave > 4 && $ave < 10) {                                     ## these values should also be flexible 
  $if_log2=1; 
} elsif ($max > 5000 && $max < 1000000) {
  $if_log2=0;  
} else {  
  printf STDERR "ERROR: Metrics not within range. MIN %.3f, MAX %.3f, AVERAGE %.3f, STDEV %.3f\n",$min,$max,$ave,$stdev; 
  exit 1; 
}
my $type=($if_log2)?"LOG2-transformed":"RAW intensities";
printf STDOUT "Expression matrix metrics: MIN %.3f, MAX %.3f, AVERAGE %.3f, STDEV %.3f; matrix type is: %s\n",$min,$max,$ave,$stdev,$type; 


############################################
#  TRANSFORM AND PRINT THE FINAL CSV FILE  #
############################################

if ($gplok) { 
  my $csvname=$gse."_".$gpl."_preprocessed.csv"; 
  open CSV,">",$csvname or die "ERROR: Failed to open CSV file for writing."; 
  print CSV $title; 

  foreach my $line (@lines) { 
    my @vals=split(/\t+/,$line); 
    my $probe=shift @vals; 
    if ($if_log2) { 
      if (defined $probe && defined $symbol{$probe} && defined $entrez{$probe}) { 
        printf CSV ("%s,%s,%s",$probe,$symbol{$probe},$entrez{$probe});
        for (my $i=0; $i<=$#vals; $i++) {
          printf CSV ",%.2f",2**$vals[$i];
        } 
        print CSV "\n"; 
      } 
    } else { 
      if (defined $probe && defined $symbol{$probe} && defined $entrez{$probe}) { 
        printf CSV ("%s,%s,%s",$probe,$symbol{$probe},$entrez{$probe});
        for (my $i=0; $i<=$#vals; $i++) {
          printf CSV ",%.2f",$vals[$i];
        } 
        print CSV "\n";  
      } 
    }  
  }
  close CSV;
  print STDOUT "Pre-processed CSV file is successfully written!\n"; 
} else { 
  print STDERR "ERROR: Sanity check failed; no CSV file is written!\n"; 
  exit 1;
}  


close FILE;
close GPL;
close ANN;
