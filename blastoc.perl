# This Perl scripts processes BLAS/LAPACK files, which are written in
# Fortran and then generates corresponding headers for C.  The
# BLAS/LAPACK files have a regular format which can be exploted to
# extract the relevant information.

# The file name.
$file = $ARGV[0];

# Open the file.
open(FORTRAN, "<$file") or die("Could not open file");

# Get the first line.  This should be the func/sub call.
$line = <FORTRAN>;
$line =~ s/^\s*//;     # Remove white space at begining of line.
$line =~ s/\( /\(/;    # Replace "( " and " ) with "(" and ")".
$line =~ s/ \)/\)/;
$line =~ s/, /,/g;     # Replace ", " with ",".

# Split the line along white space.
@tokens = split(/\s+/, $line);

# Set the return type.
$return_type="void";
if($line =~ /FUNCTION/){
    $return_type=$tokens[0];
}

# Get the function name and args.
$names  = @tokens[-1];              
@vars = split(/[,\(\)]/, $names);
$func = shift(@vars);

# Hash to keep track of the type of variable.
%typehash = ();

# We will change when we have pointers.
$pntrchar = "";
%pntrhash = ();

# Parse the file until we get to PURPOSE .
while($line = <FORTRAN>){
    $line   =~ s/^\s*//;   # Remove some gunk from the line.
    $line =~ s/, /,/g;     # Replace ", " with ",".
    $line   =~ s/PRECISION//;
    $line   =~ s/\([A-Z,\* ]*\)//g;
    # print $line, "\n";
    $fchar  = substr($line, 0, 1);  # Get the first character.
    if($fchar eq '*'){
	if ($line =~ /Array/){      # If we are to "Array" make pointers.
	    $pntrchar = "\*";
	}
    }
    else{
	# Take the first token and set it to the type.
	# The rest of the tokens are the variables of that type.
	@id   = split(/[\s,]+/, $line);
	$type = shift @id;
	foreach(@id){
	    $pntrhash{ $_ } = $pntrchar;
	    $typehash{ $_ } = $type;
	    # print $typehash{ $_ }, "\n";
	}
    }
    if ($line =~ /Purpose/){last;}  # Break of of loop if we're done.
}

# Create our Fortran function header for C.
# \L changes all to lower case. \U changes all to upper case.
$header = "\L$return_type $func\_\(";
foreach(@vars){
    # $header .= "$typehash{ $_ }$pntrhash{ $_ } $_, ";
    $header .= "\L$typehash{ $_ }\* \U$_, ";
}
$header =~ s/,\ $/\)/;   # Swap comma for ).
$header =~ s/integer/int/g;
$header =~ s/character/char/g;
# $header =~ tr/A-Z/a-z/;  # Change all to lower case.

# Check what we've got.
print $header, "\n";

# Make corresponding C wrapper.
$wrapper = "$return_type $func\(";
foreach(@vars){
    $wrapper .= "$typehash{ $_ }$pntrhash{ $_ } $_, ";
}
$wrapper =~ s/,\ $/\)/; # Swap comma for ).
$wrapper =~ tr/A-Z/a-z/;  # Change all to lower case.
$wrapper =~ s/integer/int/g;
$wrapper =~ s/character/char/g;

print $wrapper, "\n";

# Make wrapper body.
$call = "\L\treturn $func\_\(";
foreach(@vars){
    $ref = "&";
    if($pntrhash{ $_ } eq "\*"){ $ref = "" }
    $call .= "$ref$_, ";
}
$call =~ s/,\ $/\)/; # Swap comma for ).
$call =~ tr/A-Z/a-z/;  # Change all to lower case.

print $call, "\n";

# Useful to check the script.
# while ( my ($key, $value) = each(%typehash) ) {
#    print "$key => $value\n";
#}

# Close the file.
close(FORTRAN);
