### &&&&&&&&&&&&& SUBROUTINE put array in a phylip sequential ####################################################################
sub array_to_phylip_seq {
my (@phylip) =(@_);
my $SEQ='';
`rm conversion.phy`;
open(OUT,">>conversion.phy");
($tu,$sequencetest) = split(/,/,$phylip[0]);
$length = length ($sequencetest);
$n_taxa=(@phylip);
#print "$n_taxa,$length,$datatype\n";
       my $a=1; 
       foreach $phylip (@phylip) {           
            my ($taxa,$seq) = split (/,/,$phylip);
            $taxa_complete = $taxa;
	    if (length($taxa) > 10) {($taxa) = substr($taxa,0,7)} 	
              
            print OUT "$taxa_complete	$taxa$a\n";             
            if (length($taxa_complete) > 10) {push (@phylip_array," $taxa$a   $seq\n")}
		else {push (@phylip_array," $taxa   $seq\n")}
            $a++;
    }
$header = " $n_taxa  $length\n";
unshift (@phylip_array,$header);
return (@phylip_array);
}

### &&&&&&&&&&&&& SUBROUTINE put array in a phylip  ####################################################################
sub array_to_phylip {
my (@phylip) =(@_);
my $SEQ='';
`rm conversion.phy`;
open(OUT,">>conversion.phy");
($tu,$sequencetest) = split(/,/,$phylip[0]);
$length = length ($sequencetest);
$n_taxa=(@phylip);
#print "$n_taxa,$length,$datatype\n";
 
    for ($pos=0;$pos<=$length;$pos+=50) {
    #   print "$pos\n";
       my $a=1; 
       foreach $phylip (@phylip) {           
            my ($taxa,$seq) = split (/,/,$phylip);
            $taxa_complete = $taxa;
            if (length($taxa) > 10) {($taxa) = substr($taxa,0,7)}  
            my $sub = substr($seq,$pos,50);
            if ($pos ==0){
		if (length($taxa_complete) > 10) {push (@phylip_array," $taxa$a   $sub\n")}
		else {push (@phylip_array," $taxa   $sub\n")}                
                print OUT "$taxa_complete	$taxa$a\n"; 
                $a++;
                }else{push (@phylip_array,"             $sub\n")}
            }
    push (@phylip_array,"\n");
    }
$header = " $n_taxa  $length\n";
unshift (@phylip_array,$header);
return (@phylip_array);
}
### &&&&&&&&&&&&& SUBROUTINE put phylip file in array  ####################################################################
sub phylip_to_array {
my ($infile) =(@_);
my $SEQ='';
my $seq='';
my @array=();
my @taxa = ();
my $FLAG= 0;
open(PHY, "$infile")
or die "sub_phy_to_array can not open infile\n"; 



my$a=0;
while (<PHY>) {
      
      chomp;
      if ($a == 0) {
          ($taxa_n,$character_n) = $_ =~ /\s+(\d+)\s+(\d+)/;
          $a++;
          }else{
               ($all) = $_ =~ /(.+)/;
               if ($all)  { push (@all,$all)}
               }
      }

for ($i=0;$i<$taxa_n;$i++) {
    my $a=0;  
    for ($j=0;($i+$j) < @all; $j=$j+$taxa_n) {
        if ($a ==0) {
           ($tax,$seq) = $all[$i+$j] =~ /(\S+)\s+(\S+.+)/; 
           }else{
                ($seq) = $all[$i+$j];
                }
        $a++;
        $SEQ .= $seq;
        }
    $SEQ =~ s/\s//g;
    push ( @array,join(",",$tax,$SEQ));
    $SEQ ='';
    $tax='';
    }

return(@array);
}
 


### &&&&&&&&&&&&& SUBROUTINE put clustal file in array ######################################################################################################
sub clustal_to_array {
my ($infile) =(@_);
my $SEQ='';
my $seq='';
my @array=();
my @taxa = ();
my $FLAG= 0;
open(CLS, "$infile")
or die "sub_clustal_to_array can not open infile\n"; 

###    fetch the list of taxa #######
my $a=0;
while (<CLS>) {
      
      chomp;
      if ($a > 1) {
          ($tax) = $_ =~ /(\S+)\s*\S+/;
          if ($tax) {
          foreach $taxa (@taxa){
                  if ($tax eq $taxa) {$FLAG=1}
          }
      if ($FLAG == 0){
         push (@taxa,$tax);
         }else{$FLAG = 0}
         }
      }
      $a++;
}
close (CLS);

#print @taxa;
###############    fetch seq  


 
foreach $taxa (@taxa) {
    open(CLS, "$infile");
    while (<CLS>) {   
       chomp;
            my ($taxa2,$seq) = $_ =~ /^(\S+)\s+(\S+)/;
#print "$taxa2 , $taxa\n";
                  if ($taxa eq $taxa2) {
                     $SEQ .= $seq;
#print $SEQ;
                     }
                  }
  push (@array,join(",",$taxa,$SEQ));
  #print $SEQ;
  $SEQ = '';
  }
return(@array);
}
 
### &&&&&&&&&&&&& SUBROUTINE put array in A nexus  ####################################################################
sub array_to_nexus {
my (@nexus) =(@_);
my $SEQ='';

($tu,$sequencetest) = split(/,/,$nexus[0]);
$length = length ($sequencetest);
($datatype) = $sequencetest =~ /(^[ATCG-]+$)/;
#print "--$datatype\n";
if ($datatype) {$datatyp = "DNA"}else{$datatyp='protein'}
$n_taxa=(@nexus);
#print "$n_taxa,$length,$datatype\n";

$header = "#NEXUS\nBEGIN DATA;\nDimensions ntax=$n_taxa nchar=$length;\nformat missing=?\ninterleave datatype=$datatyp gap= -;\n\nmatrix\n";
$end = ";\nend;";
 
    for ($pos=0;$pos<=$length;$pos+=50) {
    #   print "$pos\n";
       foreach $nexus (@nexus) {
            my ($taxa,$seq) = split (/,/,$nexus);
            my $sub = substr($seq,$pos,50);
            push (@nexus_array,"$taxa\t$sub\n")
            }
    push (@nexus_array,"\n");
    }

unshift (@nexus_array,$header);
push (@nexus_array,$end); 
return (@nexus_array);
}
### &&&&&&&&&&&&& SUBROUTINE put nexus file in array with name and seq ####################################################################
sub nexus_to_array {
my ($infile) =(@_);
my $SEQ='';
my @array=();
my @taxa=();
my $FLAG=0;
open(NXS, "$infile")
or die "sub_fasta_to_array can not open infile\n"; 


###    fetch the list of taxa #######
my $a=0;
while (<NXS>) {
      chomp;
      my ($n_taxa) = $_ =~ /\S+\s+ntax=(\d+).+/;
      if ($n_taxa) {$num_taxa = $n_taxa}
      if (($FLAG == 1) && ($a < $num_taxa)) {
         $a++;
         my ($taxa) = $_ =~ /(\S+)\s+.+/;
         push (@taxa,$taxa); 
      }
      my ($flag) = $_ =~ /(matrix)/;
      if ($flag eq 'matrix') {$FLAG =1}
  }
close (NXS);
  
###############    fetch seq   ############
       
foreach $tax (@taxa) { 
             open(NXS, "$infile");
             while (<NXS>) {
                  chomp;
                  my ($taxa,$seq) = $_ =~ /(\S+)\s+(\S+)/;
                  if ($tax eq $taxa) {
                     $SEQ .= $seq;
                     }
                  }
push (@array,join(",",$tax,$SEQ));
$SEQ = '';
}
return(@array);
}
### &&&&&&&&&&&&& SUBROUTINE put fasta file in array with name and seq   ####################################################################

sub fasta_to_array {
my ($infile) =(@_);
my $SEQ='';
my @fasta=();

open(FASTA, "$infile")
or die "sub_fasta_to_array can not open infile\n"; 

while (<FASTA>) {
      chomp;
      my ($line) = $_ =~ /^>(.+$)/;
      if ($line) {
         if ($SEQ) {
            push (@fasta, join(",",$name_old,$SEQ));
            $SEQ ='';
            my $name_old='';
         } 
            $name_old = $line;
      }else{
      my ($line_seq) = $_ =~ /(^\S+\**)/;
      if ($line_seq) {
         $SEQ .= $line_seq;
         $line_seq = '';
      }}
}
push (@fasta, join(",",$name_old,$SEQ));
return (@fasta);

}
######

### &&&&&&&&&&&&& SUBROUTINE from array to fasta ==========>   print join ("",@fasta_format);  ######################################################################################################

sub array_to_fasta {
my @fasta_format=();
my (@fasta) =@_;



foreach $fasta (@fasta) {
        my ($name,$seq)= split (/,/,$fasta);
        ($name1) = $name =~ />(\S+)/;
        if ($name1) { $name = $name1}
        $name = ">$name";  
        for($i=0; $i <= length($seq);$i+=50) {
           $sub = substr($seq,$i,50);
           $back .= $sub;
           $back .= "\n";
        }
        push (@fasta_format, join("\n",$name,$back));
        $back='';
}
return (@fasta_format);

}
######

### &&&&&&&&&&&&& SUBROUTINE parse_fasta_array_length

sub parse_flat_length {
my @length=();
my (@array) =@_;

foreach $seq (@array) {
        my ($name,$pureseq)= split (/,/,$seq);
        $length = length ($pureseq);
        push (@length, join(",",$name,$length));
        }

return (@length);
}

######

### &&&&&&&&&&&&& SUBROUTINE parse_gff ==================>>> it gives back to you an array with an element for every gene  in this element there is in order contig,start,end,orientation,name,path,coverage,identity,n of exons,start_exon1,end_exon1,start_exon2,end_exon2 etc..

sub parse_gff {
my ($infile) =(@_);
my @gff=();
my $exons_count =0;

open(GFF, "$infile")
or die "sub parse_gff can not open infile\n"; 

while (<GFF>) {
      chomp;
      my($ctg,$spec,$gene_start,$gene_end,$gene_ori,$gene_name,$path,$coverage,$identity)= $_ =~ /^(\S+)\s\S+\s(\S+)\s(\d+)\s(\d+)\s\S*\s(.)\s.\sID=(\S+)\..+(path\d+).+coverage=(\d*.\d);identity=(\d*.\d)/;     
      if ($spec) {
         push(@gene_des, join (",",$ctg,$gene_start,$gene_end,$gene_ori,$gene_name,$path,$coverage,$identity));
      }
      my($exon_start,$exon_end)= $_ =~ /^\S+\s\S+\sexon\s(\d+)\s(\d+).+/;  
      if ($exon_start) {
         push (@exons, join(";",$exon_start,"$exon_end;"));
         $exon_count++;
      }
      my ($hashtag) = $_ =~ /(^###$)/;
      if ($hashtag) {
         push (@gff, join("",@gene_des,",$exon_count,",@exons));
         @exons =();
         @gene_des=();
         $exon_count =0;
      }
}
         push (@gff, join("",@gene_des,",$exon_count,",@exons));
return (@gff);
}

### &&&&&&&&&&&&& SUBROUTINE TRANSPOsE A MATRIX (require 2 arg: infile, delimitator  ########################################################################################################################################
sub transpose_file {
#my $infile =(@_);
my (@array) =(@_);
my $infile=$array[0];
my $del=$array[1]; #delimitator if empty = tab

    open(FILE, "$infile")
    or die "sub transpone_file can not open infile\n";
    fresse : while (<FILE>) {
       chomp;
       my $line=$_; 
#print "$line\n";
       if ($del){print "CAZZO";
            (@num_of_column) = $line =~ /$del/g;
            $num_of_column = (@num_of_column);
          }else{
                (@num_of_column) = $line =~ /\t/g;
                $num_of_column = (@num_of_column);
               }
       last fresse;
       }
    $num_of_column++;   
   #print $num_of_column;
    for ($i=1;$i<=$num_of_column;$i++){
                      if ($del){`cut -d "$del" -f $i $infile > column_temp`}else{`cut -f $i $infile > column_temp`}
                      open (ORCO,"column_temp");
                      while  (<ORCO>) {
                         chomp;
                           ($lett) = $_ =~ /(.+)/;
#print "$lett\n";
                           $LETT.=$lett;
                           $LETT.="\t"; 
                          }
 #                     `rm column_temp`;
                      push (@transponed,$LETT);
                      $LETT='';
                      }                        
`rm column_temp`;
return (@transponed);
}

1
