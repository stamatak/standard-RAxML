#!/usr/bin/env perl

print $ARGV[0]." ".$ARGV[1]." ".$ARGV[2]." ".$#ARGV."\n";

if($#ARGV < 0 || $#ARGV > 3)
  {
    print "ERROR: usage: \"perl bsBranchLengths.pl alignmentFileName treeFileName numberOfReplicates\"\n";
    exit;
  }


$alignment = $ARGV[0];
$tree      = $ARGV[1];
$reps      = $ARGV[2];


system("./raxmlHPC -f j -m GTRCAT -b 12345 -# ".$reps." -s ".$alignment." -n REPLICATES");


for($i = 0; $i < $reps; $i++)
  {
    system("./raxmlHPC -f e -m GTRGAMMA -t ".$tree." -s ".$alignment.".BS".$i." -n REP".$i);
  }

$concat = "";

for($i = 0; $i < $reps; $i++)
  {
    $concat = $concat." RAxML_result.REP".$i." ";
  }

system("cat ".$concat." > bsTrees");




