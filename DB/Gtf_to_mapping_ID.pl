#!/usr/bin/perl
# gtfファイルからtranscriptIDとgeneIDのマッピングファイルを作成する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 2 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [gtf file]\n";
}
my %opt;
getopts("i:", \%opt);
my $FILE_GTF = $opt{i};


my %DATA;
open(IN, "< $FILE_GTF") or die "Can't open $FILE_GTF.\n";
while(<IN>){
  next if(/^#/); #ignore header
  chomp;
  my %attribs = ();
  my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");
  my @add_attributes = split(";", $attributes);

  foreach my $attr ( @add_attributes ) {
     next unless $attr =~ /^\s*(.+)\s"(.+)"$/;
     my $c_type  = $1;
     my $c_value = $2;
     $attribs{$c_type} = $c_value;
  }
	if(exists $attribs{"transcript_id"}){
		my $output = sprintf "%s_%s\t%s\n", $attribs{"gene_id"}, $attribs{"gene_name"}, $attribs{"transcript_id"};
		$DATA{$attribs{"transcript_id"}} = $output;
	}
}
close IN;

foreach my $key (keys %DATA){
	print $DATA{$key};
}

