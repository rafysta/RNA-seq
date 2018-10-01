#!/usr/bin/perl
# Jul 6, 2017 ver.1

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 6 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [original quality file] -m [sample mapping] -o [converted quality file] \n";
}
my %opt;
getopts("i:m:o:", \%opt);

my $FILE_in = $opt{i};
my $FILE_mapping = $opt{m};
my $FILE_out = $opt{o};

#---------------------------------------
# mapping
#---------------------------------------
my %num2name;
my $fh_map = IO::File->new($FILE_mapping) or die "cannot open $FILE_mapping: $!";
while($_ = $fh_map->getline()){
	s/\r?\n//;
	my ($num, $name) = split /\t/;
	$num2name{$num} = $name;
}
$fh_map->close();


#---------------------------------------
# load quality file
#---------------------------------------
my $fh_in = IO::File->new($FILE_in) or die "cannot open $FILE_in: $!";
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
$fh_out->print(join("\t", qw(sample Barcode State Signal IntegSignal Circularity Size Confidence DropIndex)) . "\n");
my %title_anno = &Array2hash(&SplitTitle($fh_in->getline()));
while($_ = $fh_in->getline()){
	s/\r?\n//;
	my @data = split /\t/;
	for(my $i=0; $i < @data; $i++){
		if($data[$i] eq ""){
			$data[$i] = "NA";
		}
	}
	my $sampleNum = $data[$title_anno{'Sample'}];
	unless(exists $num2name{$sampleNum}){
		next;
	}
	
	my @output;
	push @output, $num2name{$sampleNum};
	push @output, $data[$title_anno{'Barcode'}];
	push @output, $data[$title_anno{'State'}];
	push @output, $data[$title_anno{'Signal1'}];
	push @output, $data[$title_anno{'Integ Signal1'}];
	push @output, $data[$title_anno{'Circularity1'}];
	push @output, $data[$title_anno{'Size1'}];
	push @output, $data[$title_anno{'Confidence'}];
	push @output, $data[$title_anno{'Drop index'}];

	$fh_out->print(join("\t", @output) . "\n");
}
$fh_in->close();
$fh_out->close();

sub SplitTitle{
	my ($char) = @_;
	$char =~ s/\r?\n//;
	my @data = split /\t/, $char;
	return @data;
}

sub Array2hash{
	my @data = @_;
	my %hash;
	for(my $i = 0; $i < @data; $i++){
		$hash{$data[$i]} = $i;
	}	
	return %hash;
}




