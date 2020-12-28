#!/usr/bin/env perl 
#By C. Alex Buerkle, Univ. of Wyoming

## This program removes barcode sequences from fastq reads and places
## an identifier in the info line.  This now includes a function for
## correcting barcodes that are off by 1.

## cab 1oct19 -- took previous GBS 768 barcode script and modified it to run on paired end reads 
## cab 6may20 -- anchor sequences at primer sequence for amplicon, build hash of corrected mids
## cab 13may20 -- modified to work with primers for SARS-CoV-2
## cab 14May20 -- add locus name to output, many other changes over the summer
use warnings;

## to run on teton
##  module load swset/2018.05 gcc/7.3.0 perl-text-levenshtein-xs
## parse_count.pl mid_key.csv forward.fq reverse.fq machinename

use Text::Levenshtein::XS 'distance';

unless (scalar @ARGV > 3){
    die "I need four arguments to run: barcodefile fwdfile revfile MACHINENAME ";
}
my $barcodes = shift(@ARGV);
my $forwardfile = shift (@ARGV);
my $reversefile = shift (@ARGV);
my $devicename = shift (@ARGV);

unless($devicename){
    die "Please provide a machine name that appears after \@ on info lines of fastq file";
}

open (FORWARD, $forwardfile) or die "Could not open FORWARD file";
open (REVERSE, $reversefile) or die "Could not open REVERSE file";
open (MIDS, $barcodes) or die "Could not open MID file";

my %mids;
my %midsctr;
my %forwardbarcodeshash;
my %reversebarcodeshash;
my %forwardcorrections;
my %reversecorrections;
my @forwardbarcodearray;
my @reversebarcodearray;

# my %forwardprimers = ('GCGTTCTCCATTCTGGTTACT' => 'N1',
# 	       'TGACTTCCATGCCAATGCG' => 'N2',
# 	       'TGAGCGGCTGTCTCCAC' => 'RNAseP'); 

# my %reverseprimers = ('..CCCCAAAATCAGCGAAATG' => 'N1', # GACCCCAAAATCAGCGAAATG
# 	       '..GGAACTGATTACAAACATTGGC' => 'N2', # AAGGAACTGATTACAAACATTGGC
# 	       '..TTTGGACCTGCGAGCG' => 'RNAseP'); # GATTTGGACCTGCGAGCG
# my @locusnames = ('N1', 'N2', 'RNAseP');

my %forwardprimers = ('GTGCCAGCAGCCGCGGTAA' => '16S', # 515f GTG[CT]CAGC[AC]GCCGCGGTAA
		      'GTGCCAGCCGCCGCGGTAA' => '16S', # 16S are 19bp, ITS is 22bp
		      'GTGTCAGCAGCCGCGGTAA' => '16S',
		      'GTGTCAGCCGCCGCGGTAA' => '16S',
		      'CTTGGTCATTTAGAGGAAGTAA' => 'ITS'); ## coligos seems to have 2 Ts rather than 3 Ts in the primer region at CATTTA, made search more permissive on 26 August 2020, but am now using Levenshtein distance for primer overall
my %reverseprimers = ( 'GGACTACAAGGGTATCTAAT' => '16S', # 808r GGACTAC[ACT][ACG]GGGT[AT]TCTAAT
		       'GGACTACCAGGGTATCTAAT' => '16S', # all are 20bp
		       'GGACTACTAGGGTATCTAAT' => '16S',
		       'GGACTACACGGGTATCTAAT' => '16S',
		       'GGACTACCCGGGTATCTAAT' => '16S',
		       'GGACTACTCGGGTATCTAAT' => '16S',
		       'GGACTACAGGGGTATCTAAT' => '16S',
		       'GGACTACCGGGGTATCTAAT' => '16S',
		       'GGACTACTGGGGTATCTAAT' => '16S',
		       'GGACTACAAGGGTTTCTAAT' => '16S',
		       'GGACTACCAGGGTTTCTAAT' => '16S',
		       'GGACTACTAGGGTTTCTAAT' => '16S',
		       'GGACTACACGGGTTTCTAAT' => '16S', 
		       'GGACTACCCGGGTTTCTAAT' => '16S',
		       'GGACTACTCGGGTTTCTAAT' => '16S', 
		       'GGACTACAGGGGTTTCTAAT' => '16S',
		       'GGACTACCGGGGTTTCTAAT' => '16S', 
		       'GGACTACTGGGGTTTCTAAT' => '16S',
		       'GCTGCGTTCTTCATCGATGC' => 'ITS');
my %forwardprimerMisspellings = ();
my %reverseprimerMisspellings = ();

# in R:
# reverse16<-apply(expand.grid("'GGACTAC", c('A','C', 'T'), c('A','C','G'),
#                             'GGGT', c('A', 'T'), "TCTAAT'"), 1, paste0, collapse="")
# paste(reverse16, "=> '16S'")

my %locusnames = ('16S' => 0, 
		  'ITS' => 0);
my @locusnamesArray = sort keys %locusnames;

<MIDS>; ## get rid of top line in MIDS file
while(<MIDS>){
    chomp;
    @line = split ',', $_;
    $line[0] =~ tr/[a-z]/[A-Z]/; ## catch lower-case barcode input and make it uppercase
    $line[1] =~ tr/[a-z]/[A-Z]/; ## catch lower-case barcode input and make it uppercase

    $forwardbarcodeshash{$line[0]}++;
    $reversebarcodeshash{$line[1]}++;
    $mids{$line[0]}{$line[1]}{$line[2]} = $line[3]; ## added $line[2] which indicateds locus
    $midsctr{$line[0]}{$line[1]}{$line[2]} = 0; ## initialize counters to zero
    unless(exists $locusnames{$line[2]} ){
	die "$line[2] in mids input is not in hard-coded \%locusnames";
    }
}
close (MIDS) or die "Could not close MIDS\n";

@forwardbarcodearray = sort keys %forwardbarcodeshash;
@reversebarcodearray = sort keys %reversebarcodeshash;

## simplify file names by dropping leading folder name
$forwardfile =~ s/.*\/([\w.]+fq)$/$1/; 
$forwardfile =~ s/.*\/([\w.]+fa)$/$1/;
$forwardfile =~ s/.*\/([\w.]+fastq)$/$1/;
$reversefile =~ s/.*\/([\w.]+fq)$/$1/;  
$reversefile =~ s/.*\/([\w.]+fa)$/$1/;
$reversefile =~ s/.*\/([\w.]+fastq)$/$1/;

open (FWDSEQ, "> parsed_"."$forwardfile") or 
    die "Could not open FWDSEQ\n";
open (REVSEQ, "> parsed_"."$reversefile") or 
    die "Could not open REVSEQ\n";
open (FWDCRAP, "> truemiderrors_"."$forwardfile") or die "Could not open CRAP\n";
open (REVCRAP, "> truemiderrors_"."$reversefile") or die "Could not open CRAP\n";

open (FWDOTHER, "> phixOther_"."$forwardfile") or die "Could not open FWDOTHER\n";
open (REVOTHER, "> phixOther_"."$reversefile") or die "Could not open REVOTHER\n";

my $getit = 0;
my $seqcnt = 0;
my $fwdgoodmid = 0;
my $revgoodmid = 0;
my $fwdcorrectedmid = 0;
my $revcorrectedmid = 0;
my $fwdlocus = 0;
my $revlocus = 0;
my $fwdbarcodelength = 0;
my $revbarcodelength = 0;

my $goodmidpairctr = 0;
my $badmidpairctr = 0;
my $otherpairctr = 0;
my $correctedmidctr = 0;

my $goodfwdmidctr = 0;
my $goodrevmidctr = 0;
my $goodmidpairwithcorrectdirectionprimer = 0;
my $goodfwdmidwithcorrectdirectionprimer = 0;
my $goodrevmidwithcorrectdirectionprimer = 0;
my $goodfwdmidwithincorrectdirectionprimer = 0;
my $goodrevmidwithincorrectdirectionprimer = 0;
my $badfwdmidwithcorrectdirectionprimer = 0;
my $badrevmidwithcorrectdirectionprimer = 0;
my $badfwdmidwithincorrectdirectionprimer = 0; 
my $badrevmidwithincorrectdirectionprimer = 0;

my $adlen;
my $adrem = 0;
my @ad;
my $bclen = 10 ; # 10 bp barcode to begin (consider also 9 and 8)
my $qline;

my $forward = '';
my $reverse = '';
while (defined($forward = <FORWARD>) && defined($reverse = <REVERSE>)){
    if( !($forward =~ /^\@$devicename/) || 
	!($reverse =~ /^\@$devicename/) ){  ## input file is foobarred here, with data not in sets of 4 lines
	print "parse error -- $seqcnt\n$_\n";
	while(!($forward =~ /^\@$devicename/)){
	    $forward = <FORWARD>;
	    $reverse = <REVERSE>;
	}
	while(!($reverse =~ /^\@$devicename/)){
	    $forward = <FORWARD>;
	    $reverse = <REVERSE>;
	}
	print "Back on track --> $_\n";
    }

    $seqcnt++;
    if(!($seqcnt % 50000)){
	print "$seqcnt\n";
    }
    chomp($forward = <FORWARD>);
    chomp($reverse = <REVERSE>);

    $fwdlocus = checkforprimers(\$forward, \%forwardprimers, 
				\%forwardprimerMisspellings, (19, 22)); ## get 0, or string id for locus
    $revlocus = checkforprimers(\$reverse, \%reverseprimers, 
				\%reverseprimerMisspellings, (20)); ## get 0, or string id for locus

    my $fwdlocusrevprimer = checkforprimers(\$forward, \%reverseprimers, 
					    \%reverseprimerMisspellings, (20)); 
    my $revlocusfwdprimer = checkforprimers(\$reverse, \%forwardprimers, 
					    \%forwardprimerMisspellings, (19, 22)); 


    ($fwdgoodmid,  $fwdcorrectedmid, 
     $fwdbarcodelength, $forwardmid) = lookupmid(\$forward, 
						 \%forwardbarcodeshash, 
						 \@forwardbarcodearray,
						 \%forwardcorrections,
						 \$fwdlocus);
    ($revgoodmid, $revcorrectedmid, 
     $revbarcodelength, $reversemid) = lookupmid(\$reverse, 
						 \%reversebarcodeshash,
						 \@reversebarcodearray,
						 \%reversecorrections, 
						 \$revlocus);

# $forward and $reverse are edited in place by lookupmid
    $correctedmidctr = $correctedmidctr + $fwdcorrectedmid + $revcorrectedmid;
#     print "$fwdgoodmid $forwardmid ($fwdbarcodelength), $revgoodmid $reversemid ($revbarcodelength)\n";

    if(exists $mids{$forwardmid} && 
       exists $mids{$forwardmid}{$reversemid} &&
       exists $mids{$forwardmid}{$reversemid}{$fwdlocus} && 
       $fwdgoodmid && $revgoodmid && $fwdlocus && $revlocus && ($fwdlocus eq $revlocus)){
	$goodmidpairctr++;
	$midsctr{$forwardmid}{$reversemid}{$fwdlocus}++; # fwdlocus and revlocus are the same
	
	print FWDSEQ "@"."$fwdlocus $forwardmid $reversemid $mids{$forwardmid}{$reversemid}{$fwdlocus}\n"; ##  removed $seqcnt
	print REVSEQ "@"."$revlocus $forwardmid $reversemid $mids{$forwardmid}{$reversemid}{$revlocus}\n";
	
	print FWDSEQ "$forward\n";
	print REVSEQ "$reverse\n";
    }
    elsif($revlocus eq 0 && $fwdlocus eq 0){ ## both ends lack primer
	print FWDOTHER "$forward\n";
	print REVOTHER "$reverse\n";
	$otherpairctr++;
    }
    else{ ## everything else: bad mid on one or both ends, one end lacks primer
	print FWDCRAP "$forward\n";	
	print REVCRAP "$reverse\n";
	$badmidpairctr++;
    }    

    chomp($forward = <FORWARD>);  ### should be info line that separates seq and qual, often +
    chomp($reverse = <REVERSE>);  ### should be info line that separates seq and qual, often +

    if(exists $mids{$forwardmid} && 
       exists $mids{$forwardmid}{$reversemid} &&
       exists $mids{$forwardmid}{$reversemid}{$fwdlocus} && 
       $fwdgoodmid && $revgoodmid && $fwdlocus && $revlocus && ($fwdlocus eq $revlocus) ){
	print FWDSEQ "$forward\n";
	print REVSEQ "$reverse\n";
    }
    chomp($forward = <FORWARD>);  ### quality line
    chomp($reverse = <REVERSE>);  ### 

    if(exists $mids{$forwardmid} && 
       exists $mids{$forwardmid}{$reversemid} &&
       exists $mids{$forwardmid}{$reversemid}{$fwdlocus} && 
       $fwdgoodmid && $revgoodmid && $fwdlocus && $revlocus && ($fwdlocus eq $revlocus)){
	makeQline(\$forward, $fwdbarcodelength); ## edit qline in place (pass by reference)
	makeQline(\$reverse, $revbarcodelength);	
	print FWDSEQ "$forward\n";
	print REVSEQ "$reverse\n";	
    }
    ### more counters for diagnostics and the parse report
    $goodfwdmidctr += $fwdgoodmid;
    $goodrevmidctr += $revgoodmid;
    $goodmidpairwithcorrectdirectionprimer += $fwdgoodmid * $revgoodmid * string2logical($fwdlocus);
    $goodfwdmidwithcorrectdirectionprimer += $fwdgoodmid * string2logical($fwdlocus);
    $goodrevmidwithcorrectdirectionprimer += $revgoodmid * string2logical($revlocus);
    $goodfwdmidwithincorrectdirectionprimer += $fwdgoodmid * string2logical($fwdlocusrevprimer);
    $goodrevmidwithincorrectdirectionprimer  += $revgoodmid * string2logical($revlocusfwdprimer);
    $badfwdmidwithcorrectdirectionprimer += !$fwdgoodmid * string2logical($fwdlocus);
    $badrevmidwithcorrectdirectionprimer += !$revgoodmid * string2logical($revlocus);
    $badfwdmidwithincorrectdirectionprimer += !$fwdgoodmid * string2logical($fwdlocusrevprimer);
    $badrevmidwithincorrectdirectionprimer += !$revgoodmid * string2logical($revlocusfwdprimer);

}

open(REPORT, "> parsereport_$forwardfile") or die "Failed to open report file";
print REPORT "Good barcode pair count: $goodmidpairctr (corrected $correctedmidctr individual mids)\n";
print REPORT "Pairs that lacked primer sequence: $otherpairctr\n";
print REPORT "Bad barcode pair count: $badmidpairctr\n";

print REPORT "goodfwdmidctr: $goodfwdmidctr\n";
print REPORT "goodrevmidctr: $goodrevmidctr\n";
print REPORT "goodmidpairwithcorrectdirectionprimer: $goodmidpairwithcorrectdirectionprimer\n";
print REPORT "goodfwdmidwithcorrectdirectionprimer: $goodfwdmidwithcorrectdirectionprimer\n";
print REPORT "goodrevmidwithcorrectdirectionprimer: $goodrevmidwithcorrectdirectionprimer\n";
print REPORT "goodfwdmidwithincorrectdirectionprimer:  $goodfwdmidwithincorrectdirectionprimer\n";
print REPORT "goodrevmidwithincorrectdirectionprimer: $goodrevmidwithincorrectdirectionprimer\n";
print REPORT "badfwdmidwithcorrectdirectionprimer:   $badfwdmidwithcorrectdirectionprimer\n";
print REPORT "badrevmidwithcorrectdirectionprimer:   $badrevmidwithcorrectdirectionprimer\n";
print REPORT "badfwdmidwithincorrectdirectionprimer: $badfwdmidwithincorrectdirectionprimer\n";
print REPORT "badrevmidwithincorrectdirectionprimer: $badrevmidwithincorrectdirectionprimer\n";
print REPORT "Number misspellings of forward primers: ", scalar keys %forwardprimerMisspellings, "\n";
print REPORT "Number misspellings of reverse primers: ", scalar keys %reverseprimerMisspellings, "\n";

my $outcount = '';
my @tmploci = ();
foreach my $fwd (sort keys %mids){
    foreach my $rev (sort keys %{$mids{$fwd}}){
	@tmploci = sort keys %{$mids{$fwd}{$rev}}; 
	## use hash above to get a locus name, so we can use below to
	## retrieve sample name
	print REPORT "$fwd,$rev,$mids{$fwd}{$rev}{$tmploci[0]},"; 
	foreach $i (0..($#locusnamesArray-1)){
	    if(exists $midsctr{$fwd}{$rev}{$locusnamesArray[$i]}){
		$outcount = $midsctr{$fwd}{$rev}{$locusnamesArray[$i]};
	    }
	    else{
		$outcount = 'NA';
	    }
	    print REPORT "$outcount,";
	}
	if(exists $midsctr{$fwd}{$rev}{$locusnamesArray[$#locusnamesArray]}){
	    $outcount = $midsctr{$fwd}{$rev}{$locusnamesArray[$#locusnamesArray]};
	}
	else{
	    $outcount = 'NA';
	}
	print REPORT "$outcount\n";
    }
}
close(REPORT) or die "Could not close REPORT\n";
close(FWDSEQ) or die "Could not close FWDSEQ\n";
close(REVSEQ) or die "Could not close REVSEQ\n";
close(FWDCRAP) or die "Could not close FWDCRAP\n";
close(REVCRAP) or die "Could not close FWDCRAP\n";
close(FWDOTHER) or die "Could not close FWDOTHER\n";
close(REVOTHER) or die "Could not close REVOTHER\n";

##--------- sub routines -----------------##

sub string2logical{
    if($_[0] eq 0){
	return(0);
    }
    else{
	return(1);
    }
}

sub correctmid{
    my $testseq = $_[0];
    my $arrayref = $_[1];

    my $corrected;
    my $min = 100;
    my $mindex = -1;
    for my $i (0 .. $#{$arrayref}){
	$dist = distance($testseq, ${$arrayref}[$i]);
	if($min > $dist){
	    $min = $dist;
	    $mindex = $i;
	}
    }
    return(${$arrayref}[$mindex], $min);
}

sub checkforprimers{
    my $target = $_[0];
    my $primerhashref = $_[1];
    my $altprimerhashref = $_[2]; ## hash to store edits/misspellings of expected primer sequence
    my $primerlengtharray = $_[3];
    my $queryseq = '';
    my $primer = '';

    ## test for exact match to a known primer sequence in one of two
    ## hashes the canonical primer sequences in primerhashref, and the
    ## misspellings in altprimerhasref.  If neither produces a hit,
    ## search for possible edits in the second for loop
    foreach my $midlength (10, 9, 8){
	foreach my $primerlength ( @{$primerlengtharray} ){
	    $queryseq = substr($$target, $midlength, $primerlength);
	    if( exists( ${$primerhashref}{$queryseq} )){
		return( ${$primerhashref}{$queryseq});
	    }
	    elsif(exists( ${$altprimerhashref}{$queryseq} )){
		return(${$altprimerhashref}{$queryseq});
	    }
	}
    }	
    foreach $primer (keys %{$primerhashref}){
	## implement distance measure, rather than exact match that is used above.
	foreach $midlength (10, 9, 8){
	    $queryseq = substr($$target, $midlength, length $primer);
	    if( distance($queryseq, $primer) < 5){ 
		## edit distance of 4 or smaller
		## store the match in another hash, so it can be used in lookup above
		${$altprimerhashref}{$queryseq} = ${$primerhashref}{$primer};
		return(${$primerhashref}{$primer});
	    }
	}
    }    
    return(0);  
}
    # previously:
    # foreach my $primer (keys %{$primerhashref}){
    # 	if($$target =~ /\w{8,10}$primer/){
    # 	    return(${$primerhashref}{$primer}); ## returns locus name, based on extact match to primer
    # 	}
    # }


sub lookupmid{
    my $targetref = $_[0];
    my $hashref = $_[1];
    my $arrayref = $_[2];
    my $correctionshashref = $_[3];
    my $locusprimerref = $_[4];
    my $goodmid = 0;
    my $correctedmid = 0;
    my $mid = '';

    $line10 = substr($$targetref, 0, $bclen); ## substr EXPR,OFFSET,LENGTH
    $line9 = substr($$targetref, 0, ($bclen - 1)); # 9 bp barcode
    $line8 = substr($$targetref, 0, ($bclen - 2)); # 8 bp barcode

    if(${$hashref}{$line10} || ${$correctionshashref}{$line10}){ ##  barcode is in $targetref, remove barcode
	$whichL = 10;
	$goodmid = 1;
	$$targetref =~ s/^$line10//;
	if(${$hashref}{$line10}){
	    $mid = $line10;
	}
	else{
	    $mid = ${$correctionshashref}{$line10};
	    $correctedmid = 1;
	}
    }	
    elsif(${$hashref}{$line9} || ${$correctionshashref}{$line9}){ ##  barcode is in $targetref, remove barcode 
	$whichL = 9;
	$goodmid = 1;
	$$targetref =~ s/^$line9//;
	if(${$hashref}{$line9}){
	    $mid = $line9;
	}
	else{
	    $mid = ${$correctionshashref}{$line9};
	    $correctedmid = 1;
	}
    }
    elsif(${$hashref}{$line8} || ${$correctionshashref}{$line8}){ ##  barcode is in $targetref, remove barcode 
	$whichL = 8;
	$goodmid = 1;
	$$targetref =~ s/^$line8//;
	if(${$hashref}{$line8}){
	    $mid = $line8;
	}
	else{
	    $mid = ${$correctionshashref}{$line8};
	    $correctedmid = 1;
	}
    }
    else { ## potential mid error, we already know barcode is there
	if($$locusprimerref){
	    ($line10b, $n10) = correctmid($line10, $arrayref);
	    $minN = $n10;
	    $whichL = 10;
	    if ($minN > 1){
		($line9b, $n9) = correctmid($line9, $arrayref);
		if ($n9 < $minN){
		    $minN = $n9;
		    $whichL = 9;
		}
		if ($minN > 1){
		    ($line8b, $n8) = correctmid($line8, $arrayref);
		    if ($n8 < $minN){
			$minN = $n8;
			$whichL = 8;
		    }
		}
	    }
	    if ($minN > 1){ ## can't correct
		$goodmid = 0;
	    }
	    else { ## has been corrected
		if ($whichL == 10){
		    ${$correctionshashref}{$line10} = $line10b;
		    ## corrections should only need to happen only once
		    $mid = $line10b;
		    $$targetref =~ s/^$line10//;
		}
		elsif ($whichL == 9){
		    ${$correctionshashref}{$line9} = $line9b;
		    $mid = $line9b;
		    $$targetref =~ s/^$line9//;
		}
		elsif ($whichL == 8){
		    ${$correctionshashref}{$line8} = $line8b;
		    $mid = $line8b;
		    $$targetref =~ s/$line8//;
		}
		$goodmid = 1;
		$correctedmid = 1;
	    }
	}
	else{
	    $goodmid = 0;
	}
    }
    return($goodmid, $correctedmid, $whichL, $mid);
}

sub makeQline{
    my $qlineref = $_[0];
    my $barcodelength = $_[1];
    
    if ($barcodelength == 10){
	$seqlength = (length ${$qlineref}) - $bclen ;
	${$qlineref} = substr(${$qlineref}, $bclen, $seqlength);
    }
    elsif ($barcodelength == 9){
	$seqlength = (length ${$qlineref}) - $bclen  + 1;
	${$qlineref} = substr(${$qlineref}, ($bclen - 1), $seqlength);
    }
    elsif ($barcodelength == 8){
	$seqlength = (length ${$qlineref}) - $bclen + 2;
	${$qlineref} = substr(${$qlineref}, ($bclen - 2), $seqlength);
    }
}
