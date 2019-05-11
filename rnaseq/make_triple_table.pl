#!/usr/bin/env perl 

use strict; 
use warnings; 

my $srr_gsm = shift @ARGV; 
my $gsm_gse = shift @ARGV; 

open SRR_GSM,"<",$srr_gsm or die "$!"; 
open GSM_GSE,"<",$gsm_gse or die "$!"; 

my %gsm2gse; 

while (<GSM_GSE>) { 
  chomp; 
  my @t = split /\t/; 
  $gsm2gse{$t[0]} = $t[1]; 
}

while (<SRR_GSM>) { 
  chomp; 
  my @t = split /\t/; 
  printf "%s\t%s\t%s\n",$t[0],$t[1],$gsm2gse{$t[1]}; 
} 

close SRR_GSM; 
close GSM_GSE; 
