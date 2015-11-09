#!/usr/bin/perl

# Script to add CONECT records to pdb files containing coarse grain particles
#
# Usage: ./connect_pdb_atom.pl [input pdbfile] > [output pdbfile]
#
# The sscript was built for a POPC bilayer with a coarse grain protein
#
# The coarse grain protein follows MARTINI coarse grain force field
#
# Modified by Xiaofeng Yu, from 
# Author: Vlad Cojocaru (vlad.cojocaru@eml-r.villa-bosch.de)


use strict;
use warnings;

my @aminoacids=("GLY", "ALA", "PRO", "VAL", 
                "LEU", "ILE", "CYS", "MET",
	        "SER", "THR", "ASN", "GLN", 
                "PHE", "HIS", "TRP", "TYR", 
	        "LYS", "ARG", "ASP", "GLU"); 

my $pdbfile=$ARGV[0];

my @data=&store_data();
&write_pdb(@data);
&write_connect(@data);


sub store_data {
# subroutine to read in the pdb file; should work for most pdb files
# can be used in other scripts to modify pdb files if required

my $prev_res=0;
my $curr_res=0;
my $prev_rid=0;
my $atom=0;

my %residues;
my %columns;

open PDBFILE, "$pdbfile" or die "Cannot open file $pdbfile for reading: $!";

while ( my $line=<PDBFILE> ) {
 if ( $line =~ /^ATOM/ || $line =~ /^HETATM/ ) {
  chomp $line;
  $atom++;
  $columns{"atom"}=$atom;
  $columns{"ind"}=substr $line, 6, 5; $columns{"ind"} =~ s/^\s+//; $columns{"ind"} =~ s/\s+$//;
  $columns{"name"}=substr $line, 11, 5; $columns{"name"} =~ s/^\s+//; $columns{"name"} =~ s/\s+$//;
  $columns{"rname"}=substr $line, 16, 5; $columns{"rname"} =~ s/^\s+//; $columns{"rname"} =~ s/\s+$//;
  $columns{"chain"}=substr $line, 21, 1; $columns{"chain"} =~ s/^\s+//; $columns{"chain"} =~ s/\s+$//;
  $columns{"rid"}=substr $line, 22, 4; $columns{"rid"} =~ s/^\s+//; $columns{"rid"} =~ s/\s+$//;
  $columns{"x"}=substr $line, 26, 12; $columns{"x"} =~ s/^\s+//; $columns{"x"} =~ s/\s+$//;
  $columns{"y"}=substr $line, 38, 8; $columns{"y"} =~ s/^\s+//; $columns{"y"} =~ s/\s+$//;
  $columns{"z"}=substr $line, 46, 8; $columns{"z"} =~ s/^\s+//; $columns{"z"} =~ s/\s+$//;
  $columns{"occ"}=substr $line, 54, 6; $columns{"occ"} =~ s/^\s+//; $columns{"occ"} =~ s/\s+$//;
  $columns{"b"}=substr $line, 60, 6; $columns{"b"} =~ s/^\s+//; $columns{"b"} =~ s/\s+$//;

  if ( $columns{"rid"} != $prev_rid) { $curr_res++; } 

  foreach my $column (keys %columns) {
   push @{ $residues{$curr_res}{$column} }, $columns{$column} unless ( $column eq "rname" || $column eq "chain" );
   $residues{$curr_res}{$column} = $columns{$column} if ( $column eq "rname" && $curr_res != $prev_res );
   $residues{$curr_res}{$column} = $columns{$column} if ( $column eq "chain" && $curr_res != $prev_res );  
  }
  
  $prev_rid=$columns{"rid"};
  $prev_res=$curr_res;
 }
}
close PDBFILE or die "Cannot close file $pdbfile: $!";

return (\%residues, \%columns);

}

sub write_pdb {
# subroutine that writes a new pdb file with the data read in 
# the "store_data" subroutine; also useful if modification of the 
# original pdb file without introducing CONECT records is required;
# it needs adaptation depending on pdb file and modification desired

my %residues = %{$_[0]};
my %columns = %{$_[1]};

foreach my $residue (sort {$a <=> $b} keys %residues) {
 my $chain=$residues{$residue}{"chain"}; 
 my $rname=$residues{$residue}{"rname"};
 my $frname;
 if ($rname =~ /[A-Za-z1-9][A-Za-z1-9][A-Za-z1-9][A-Za-z1-9]/) {
  $frname="$rname";
 } elsif ($rname =~ /[A-Za-z1-9][A-Za-z1-9][A-Za-z1-9]/) {
  $frname="$rname ";
 } elsif ($rname =~ /[A-Za-z1-9][A-Za-z1-9]/) {
  $frname="$rname  ";
 } elsif ($rname =~ /[A-Za-z1-9]/) {
  $frname="$rname   ";
 }
 my @residue_atoms = @{ $residues{$residue}{"atom"} };
 for (my $i=0; $i <= $#residue_atoms; $i++) {
  my $name=$residues{$residue}{"name"}[$i];
  my $fname;
  if ($name =~ /[A-Za-z1-9][A-Za-z1-9][A-Za-z1-9][A-Za-z1-9]/) {
   $fname="$name";
  } elsif ($name =~ /[A-Za-z1-9][A-Za-z1-9][A-Za-z1-9]/) { 
   $fname=" $name";
  } elsif ($name =~ /[A-Za-z1-9][A-Za-z1-9]/) { 
   $fname=" $name ";
  } elsif ($name =~ /[A-Za-z1-9]/) { 
   $fname=" $name  ";
  }

  printf "ATOM%7s%5s%5s%1s%4s%12s%8s%8s%6s%6s\n",
         $residues{$residue}{"atom"}[$i], 
	 $fname, 
	 $frname, 
	 $chain, 
	 $residues{$residue}{"rid"}[$i],
	 $residues{$residue}{"x"}[$i],
	 $residues{$residue}{"y"}[$i],
	 $residues{$residue}{"z"}[$i],
	 $residues{$residue}{"occ"}[$i], 
	 $residues{$residue}{"b"}[$i];
 }
}

}


sub write_connect {
# subroutine to write the correct CONECT records;
# this was built for POPC lipids and Sansom's coarse grain proteins
# it needs adaptation if used with other force fields or other lipids
 
my %residues = %{$_[0]};
my %columns = %{$_[1]};

my $prev_res=0;
foreach my $curr_res (sort {$a <=> $b} keys %residues) {
  
 my $curr_rname = $residues{$curr_res}{"rname"};
 my $prev_rname = $residues{$prev_res}{"rname"} unless ( $prev_res == 0 );

 my @curr_res_atoms = @{ $residues{$curr_res}{"atom"} };
 my @prev_res_atoms = @{ $residues{$prev_res}{"atom"} } unless ( $prev_res == 0 );

 my (%curr_data, %prev_data);
 
 for (my $i=0; $i <= $#curr_res_atoms; $i++) {
  my $curr_name = $residues{$curr_res}{"name"}[$i];
  my $curr_atom = $residues{$curr_res}{"atom"}[$i];
  my $curr_ind = $residues{$curr_res}{"ind"}[$i];
  push @{ $curr_data{$curr_name} }, $curr_atom;
  push @{ $curr_data{$curr_name} }, $curr_ind;
 }

 if ( $prev_res != 0 ) {
  for (my $i=0; $i<= $#prev_res_atoms; $i++) {
   my $prev_name = $residues{$prev_res}{"name"}[$i];
   my $prev_atom = $residues{$prev_res}{"atom"}[$i];
   my $prev_ind = $residues{$prev_res}{"ind"}[$i];
   push @{ $prev_data{$prev_name} }, $prev_atom;
   push @{ $prev_data{$prev_name} }, $prev_ind;
  } 
 }
 
 if ( $curr_rname =~ "POP.*" ) {

  printf "CONECT%5s%5s\n", $curr_data{"NC3"}[0], $curr_data{"PO4"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"PO4"}[0], $curr_data{"GL1"}[0];
  printf "CONECT%5s%5s%5s\n", $curr_data{"GL1"}[0], $curr_data{"GL2"}[0], $curr_data{"C1A"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"C1A"}[0], $curr_data{"C2A"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"C2A"}[0], $curr_data{"C3A"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"C3A"}[0], $curr_data{"C4A"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"GL2"}[0], $curr_data{"C1B"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"C1B"}[0], $curr_data{"C2B"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"C2B"}[0], $curr_data{"D3B"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"D3B"}[0], $curr_data{"C4B"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"C4B"}[0], $curr_data{"C5B"}[0];

 }
 
 my @check_prev_rname = grep { /$prev_rname/ } @aminoacids unless ( $prev_res == 0 );
 my @check_curr_rname = grep { /$curr_rname/ } @aminoacids;
 if ( @check_prev_rname && @check_curr_rname ) {
  printf "CONECT%5s%5s\n", $prev_data{"BB"}[0], $curr_data{"BB"}[0];  
 }  
 
 if ( $curr_rname =~ /ASN|ASP|MET|LEU|ILE|GLU|GLN|CYS|SER|PRO|THR|VAL/ ) {
  printf "CONECT%5s%5s\n", $curr_data{"BB"}[0], $curr_data{"SC1"}[0];  
 }
 if ( $curr_rname =~ /PHE/ ) {
  printf "CONECT%5s%5s\n", $curr_data{"BB"}[0], $curr_data{"SC1"}[0];  
  printf "CONECT%5s%5s\n", $curr_data{"SC1"}[0], $curr_data{"SC2"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"SC2"}[0], $curr_data{"SC3"}[0];
 }
 if ( $curr_rname =~ /HIS/ ) {
  printf "CONECT%5s%5s\n", $curr_data{"BB"}[0], $curr_data{"SC1"}[0];  
  printf "CONECT%5s%5s\n", $curr_data{"SC1"}[0], $curr_data{"SC2"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"SC2"}[0], $curr_data{"SC3"}[0];
 }
 if ( $curr_rname =~ /TRP/ ) {
  printf "CONECT%5s%5s\n", $curr_data{"BB"}[0], $curr_data{"SC1"}[0];  
  printf "CONECT%5s%5s\n", $curr_data{"SC1"}[0], $curr_data{"SC2"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"SC2"}[0], $curr_data{"SC3"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"SC2"}[0], $curr_data{"SC4"}[0];
 }
 if ( $curr_rname =~ /TYR/ ) {
  printf "CONECT%5s%5s\n", $curr_data{"BB"}[0], $curr_data{"SC1"}[0];  
  printf "CONECT%5s%5s\n", $curr_data{"SC1"}[0], $curr_data{"SC2"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"SC2"}[0], $curr_data{"SC3"}[0];
 }
 if ( $curr_rname =~ /LYS/ ) {
  printf "CONECT%5s%5s\n", $curr_data{"BB"}[0], $curr_data{"SC1"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"SC1"}[0], $curr_data{"SC2"}[0];  
 }
 if ( $curr_rname =~ /ARG/ ) {
  printf "CONECT%5s%5s\n", $curr_data{"BB"}[0], $curr_data{"SC1"}[0];
  printf "CONECT%5s%5s\n", $curr_data{"SC1"}[0], $curr_data{"SC2"}[0];  
 }

 $prev_res=$curr_res;
 
}

}
