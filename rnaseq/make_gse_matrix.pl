#!/usr/bin/env perl 

use strict; 
use warnings; 
use Data::Dumper; 

my $gse = shift @ARGV; 
my $triple = shift @ARGV;   ## "triple" table of SRR-GSM-GSE
my $gene_tr = shift @ARGV;  ##  gene-to-transcript table 
my $gene_entrez = shift @ARGV; 
my $gsm_names = shift @ARGV; 

my $dir = "/pub37/alexp/tmp/sra/all_rn/abund/"; 

my $genes = {}; 
my $entrez = {}; 
my $gse_matrix = {};

open TRIPLE,"<",$triple or die "$!";
open GENE_TR,"<",$gene_tr or die "$!"; 
open ENTREZ,"<",$gene_entrez or die "$!"; 
open GSM_NAMES,"<",$gsm_names or die "$!"; 

while (<GENE_TR>) { 
  chomp; 
  my @t = split /\t+/; 
  $genes->{$t[1]} = $t[0]; 
}

while (<ENTREZ>) {
  chomp; 
  my @t = split /\t+/; 
  $entrez->{$t[0]} = "$t[1]\t$t[2]"; 
} 

while (<TRIPLE>) { 
  chomp; 
  my @t = split /\t+/; 
  push @{$gse_matrix->{$t[1]}->{SRR}},$t[0] if ($t[2] eq $gse); 
} 

foreach my $gsm (keys %{$gse_matrix}) { 
  my @srrs = @{$gse_matrix->{$gsm}->{SRR}}; 
  foreach my $srr (@srrs) { 
    my $srr_count_sum = 0; 

    my $abund = join '',$dir,$srr,".abundance.tsv";
    open ABUND,"<",$abund or die "$!"; 
    <ABUND>; ## skip header 
    while (<ABUND>) { 
      my $tr = (split /\t+/)[0]; 
      my $count = (split /\t+/)[3];
      my $gene = $genes->{$tr}; 
      if (defined $gse_matrix->{$gsm}->{exp}->{$gene}) { 
        $gse_matrix->{$gsm}->{exp}->{$gene} += $count; 
      } else {
        $gse_matrix->{$gsm}->{exp}->{$gene} = $count;
      }
      $srr_count_sum += $count;  
    } 

    close ABUND; 
    if (defined $gse_matrix->{$gsm}->{sum}) { 
      $gse_matrix->{$gsm}->{sum} += $srr_count_sum; 
    } else { 
      $gse_matrix->{$gsm}->{sum} = $srr_count_sum;
    } 
  }
  printf STDERR "Processed series %s, sample %s, count sum = %d\n",$gse,$gsm,$gse_matrix->{$gsm}->{sum}; 
}       

while (<GSM_NAMES>) { 
  chomp; 
  my @t = split /\t+/; 
  my $gsm = $t[1]; 
  my $name = $t[5]; 
  $gse_matrix->{$gsm}->{name} = $name if (defined $gse_matrix->{$gsm}); ## sometimes this will be defined twice or more, that's OK 
}  

# print Dumper $gse_matrix;

my @gsms = sort keys %{$gse_matrix}; 
my @unique_genes = do { my %seen; grep { !$seen{$_}++ } values %{$genes}};
my @genes = sort @unique_genes; 

my $header = "Ensembl\tName\tEntrez"; 

foreach my $gsm (@gsms) {
  my $title = join "-",$gse_matrix->{$gsm}->{name},$gsm;  
  $header = join "\t",$header,$title; 
} 

print "$header\n"; 

foreach my $gene (@genes) { 
  my $out = join "\t",$gene,$entrez->{$gene}; 
  foreach my $gsm (@gsms) { 
    $out = join "\t",$out,int($gse_matrix->{$gsm}->{exp}->{$gene}+0.5);
  } 
  print "$out\n"; 
}  

close GENE_TR;
close TRIPLE;
close ENTREZ; 
close GSM_NAMES;  
  
