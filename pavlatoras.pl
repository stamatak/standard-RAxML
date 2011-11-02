#!/usr/bin/perl
use strict;
use warnings;

# those are the ref alignment and ref tree

my $referenceData = "44.phy";
my $referenceTree = "44.tree";

# this is the subalignment, i.e., placements on a branch

my $subAlignment = "10.phy";

# estimate model parames on reference tree to obtain a binary model parameter file

system("./raxmlHPC-SSE3 -f e -m GTRGAMMA -t ".$referenceTree." -s ".$referenceData." -n MODEL\n");

# do a ML search on the subalignment using the binary model data file, that is, without optimizing 
# ML model params at all

system("./raxmlHPC-SSE3 -m GTRGAMMA -R RAxML_binaryModelParameters.MODEL -p 12345 -s ".$subAlignment." -n ML_TREE\n");

my $fileNames = "cat RAxML_bestTree.ML_TREE ";
my $range = 1000000;

# now generate 10 or more random trees  on the subalignment of reads

for(my $i = 0; $i < 10; $i++)
{
    my $r = rand($range); 
	
    system("./raxmlHPC-SSE3 -y -d -m GTRCAT -p ".$r." -s ".$subAlignment." -n R".$i."\n");
    
    $fileNames = $fileNames." RAxML_randomTree.R".$i." ";
}

# concatenate all tree files, i.e., the ML tree and the 10 random trees

system($fileNames. " > allTrees");

# compute per site log likelihoods on all those trees, the per-site log likes can then be used with CONSEL for a significance 
# test 

system("./raxmlHPC-SSE3 -m GTRGAMMA -R RAxML_binaryModelParameters.MODEL -s ".$subAlignment." -z allTrees -f g -n allLikes\n");
