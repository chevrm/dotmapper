#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

## Setup script path
my $script = abs_path($0);
my @sp = split(/\//, $script);
my $script_dir = join('/', @sp[0..$#sp-1]);

my ($ref, @tocompare) = @ARGV;

my $threads = 8; ## default = 8
my $evalue = '1e-3'; ## default = '1e-3'
my $idthresh = 0; ## default = 0
my $maxstep = 300; ## default = 300bp
my $nwin = 300;  # default = 300bp
my $avgwin = 3000;  ## default = 3000 windows
my $min_len = 3000; # default = 3kb

my $reflen = makefna($ref, 'ref.fna');
die "ABORT: $ref not above $min_len bp.\n" unless($reflen >= $min_len);
my @r = split(/\//, $ref);
my @rf = split(/\./, $r[-1]);
my $rpref = join('.', @rf[0..$#rf-1]);

foreach my $comp (@tocompare){
    my $complen = makefna($comp, 'comp.fna');
    if($complen < $min_len){
	print STDERR "SKIPPING: $comp not above $min_len bp.\n";
	next;
    }else{
	my @c = split(/\//, $comp);
	my @cf = split(/\./, $c[-1]);
	my $cpref = join('.', @cf[0..$#cf-1]);

	my $step = int( (($reflen+$complen)/2)/$avgwin );
	$step = 1 if ($step==0);
	$step = $maxstep if ($step > $maxstep);
	print "Ref=$ref\tComp=$comp\n";
	print "\tStepsize = $step bp\n";
	print "\tRef len: $reflen\tComp len: $complen\n";

	## Make tmp subject
	my $refwincount = 0;
	open my $ts, '>', 'tmp.s' or die $!;
	open my $win, '>', 'tmp.win' or die $!;
	my $ac = new Bio::SeqIO(-file=>'ref.fna', -format=>'fasta');
	while(my $seq = $ac->next_seq){
	    for (my $i=1;$i+$nwin<=$seq->length;$i+=$step){
		print $ts ">refcat_window$refwincount ".$seq->id."\n".$seq->subseq($i,$i+$nwin)."\n";
		print $win join("\t", "refcat_window$refwincount", $seq->id, $i, $i+$nwin)."\n";
		$refwincount+=1;
	    }
	}
	close $ts;

	## Make blast db
	system("makeblastdb -in tmp.s -out tmp.db -dbtype nucl > /dev/null");

	## Make tmp query
	open my $tq, '>', 'tmp.q' or die $!;
	my $compwincount = 0;
	my $bc = new Bio::SeqIO(-file=>'comp.fna', -format=>'fasta');
	while(my $seq = $bc->next_seq){
	    for (my $i=1;$i+$nwin<=$seq->length;$i+=$step){
		print $tq ">compcat_window$compwincount ".$seq->id."\n".$seq->subseq($i,$i+$nwin)."\n";
		print $win join("\t", "compcat_window$compwincount", $seq->id, $i, $i+$nwin)."\n";
		$compwincount+=1;
	    }
	}
	close $tq;
	close $win;

	## Do the blast!
	print "Performing blastn...";
	system("blastn -query tmp.q -db tmp.db -outfmt 6 -evalue $evalue -max_target_seqs 1000000 -num_threads $threads -out ba.bn");
	print "DONE!\n";

	## Parse the blast
	if(-z 'ba.bn'){
	    print STDERR "SKIPPING: blastn returned zero hits.\n";
	    next;
	}else{
	    my %bahit = ();
	    my $bh = 0;
	    open my $bl, '<', 'ba.bn' or die $!;
	    while(<$bl>){
		chomp;
		my ($query, $hit, $pctid, $alilen, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $_);
		if($pctid>=$idthresh){
		    $bahit{$query}{$hit} = $pctid;
		    $bh +=1;
		}
	    }
	    close $bl;

	    ## Write the python-ready csv
	    print "Generating matrix...";
	    open my $of, '>', 'r.csv' or die $!;
	    ## Header
	    foreach my $awn (0..$refwincount-1){
		print $of ",$awn";
	    }
	    print $of "\n";
	    ## Data rows
	    my $dp = 0;
	    my $hi = 0;
	    foreach my $bwn (0..$compwincount-1){
		my $bid = 'compcat_window' . $bwn;
		print $of $bwn;
		if(exists $bahit{$bid}){
		    foreach my $awn (0..$refwincount-1){
			my $aid = 'refcat_window' . $awn;
			if(exists $bahit{$bid}{$aid}){
			    print $of "," . $bahit{$bid}{$aid};
			    $hi = $bahit{$bid}{$aid} if($hi < $bahit{$bid}{$aid});
			    $dp += 1;
			}else{
			    print $of ",$idthresh";
			}
		    }
		}else{
		    foreach my $awn (0..$refwincount-1){
			print $of ",$idthresh";
		    }
		}
		print $of "\n";
	    }
	    close $of;
	    print "DONE\n";
	    
	    ## Make the dotplot heatmap
	    print "Generating heatmap...";
	    my $png = join(".", $rpref, $cpref, "png");
	    my ($aw, $bw) = ($refwincount-1,$compwincount-1);
	    system("python $script_dir/clustheat.py $rpref $cpref $png $aw $bw");
	    print "DONE\n\n";
	}
    }
}

## Cleanup!
my @torm = ('ba.bn', 'comp.fna', 'r.csv', 'ref.fna', 'tmp*');
foreach my $r (@torm){
    system("rm $r 2> /dev/null");
}


sub makefna{
    my ($f, $fna) = (shift, shift);
    my $len = 0;
    if($f =~ m/\.gb(k)?$/){
	open my $ofh, '>', $fna or die $!;
	my $gb = new Bio::SeqIO(-file=>$f, -format=>'genbank');
	my $seq=$gb->next_seq;
	print $ofh '>'.$seq->id."\n".$seq->seq."\n";
	close $ofh;
	$len = $seq->length;
    }elsif($f =~ m/\.f(n)?a(sta)?$/){
	open my $ofh, '>', $fna or die $!;
	my $fa = new Bio::SeqIO(-file=>$f, -format=>'fasta');
	my $seq=$fa->next_seq;
	print $ofh '>'.$seq->id."\n".$seq->seq."\n";
	close $ofh;
	$len = $seq->length;
    }
    return($len);
}

