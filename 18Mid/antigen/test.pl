#!/usr/bin/perl
use warnings;
use strict;
use tool;


my $dir="./modelData/referSeq/";
opendir(IN,$dir)or die"$!";
my @files;
while(my $file=readdir(IN)){
	if($file=~/^H/){
		push(@files,$file);
	}
}
closedir IN or die"$!";

foreach my$file(@files){
	my $infile=$dir.$file;
	my $seqRef=readSeq($infile);
	my @id=keys %{$seqRef};
	my $type=$file;
	$type=~s/\_.+$//;
	my $outfile=$dir.$type;
	open(OUT,">$outfile")or die"$!";
	for(my $i=0;$i<scalar @id;$i++){
		my $newName="standardSeq".($i+1);
		print OUT ">$newName\n",join("",@{$seqRef->{$id[$i]}}),"\n";
	}
	close OUT or die"$!";

	print "$file is over\n";
}

			
