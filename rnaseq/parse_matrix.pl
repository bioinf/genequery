#!/usr/bin/env perl 
#
#

use strict;
use warnings;

my $gz = shift @ARGV;

my $gse;
my @gsm;
my @lib;
my @org;
my @str;
my @ttl;

open GZ, "zcat $gz |" or die "gunzip $gz: $!";

while (<GZ>) {
  $gse = $1 if (m/!Series_geo_accession\s+"(GSE.*?)"/);
  @gsm = split /\t+/ if (m/"ID_REF"\t/);
  @lib = split /\t+/ if (m/!Sample_library_selection\t/);
  @org = split /\t+/ if (m/!Sample_organism_ch1\t/);
  @str = split /\t+/ if (m/!Sample_library_strategy\t/);
  @ttl = split /\t+/ if (m/!Sample_title\t/);
}

shift @gsm;
shift @lib;
shift @org;
shift @str;
shift @ttl;

open GZ, "zcat $gz |" or die "gunzip $gz: $!";

while (<GZ>) {
  $gse = $1 if (m/!Series_geo_accession\s+"(GSE.*?)"/);
  @gsm = split /\t+/ if (m/"ID_REF"\t/);
  @lib = split /\t+/ if (m/!Sample_library_selection\t/);
  @org = split /\t+/ if (m/!Sample_organism_ch1\t/);
  @str = split /\t+/ if (m/!Sample_library_strategy\t/);
  @ttl = split /\t+/ if (m/!Sample_title\t/);
}

shift @gsm;
shift @lib;
shift @org;
shift @str;
shift @ttl;

my $n_smp = scalar @gsm;

for (my $i=0; $i<$n_smp; $i++) {
  $gsm[$i] = "NONE" if (!defined $gsm[$i]);
  $lib[$i] = "NONE" if (!defined $lib[$i]);
  $org[$i] = "NONE" if (!defined $org[$i]);
  $str[$i] = "NONE" if (!defined $str[$i]);
  $ttl[$i] = "NONE" if (!defined $ttl[$i]);
  chomp $gsm[$i];
  chomp $lib[$i];
  chomp $org[$i];
  chomp $str[$i];
  chomp $ttl[$i];
  $gsm[$i] =~ s/"//g;
  $lib[$i] =~ s/"//g;
  $org[$i] =~ s/"//g;
  $str[$i] =~ s/"//g;
  $ttl[$i] =~ s/"//g;
  print "$gse\t$gsm[$i]\t$lib[$i]\t$org[$i]\t$str[$i]\t$ttl[$i]\n";
}
close GZ;

